import os
import warnings
from multiprocessing import Pool
from pathlib import Path
from typing import Union, List

import numpy as np
import scipy as sc
from Bio import BiopythonWarning
from numba import njit

from pymembrane.parameters.thylakoid.protein_default import prepare_thylakoid_protein_default
from pymembrane.structure.atomic_protein import DynamicProteinAtomic
from pymembrane.structure.cg_particle import CGProtein
from pymembrane.structure.membrane import calculate_atomic_clashes_sparse, calculate_atomic_clashes
from pymembrane.util.physical_constants import precision

warnings.simplefilter('ignore', BiopythonWarning)


def get_r_range(
        protein_a: Union[DynamicProteinAtomic, CGProtein],
        protein_b: Union[DynamicProteinAtomic, CGProtein],
        step: float = 5.0,
        cushion: float = 5.0
) -> np.ndarray:
    """
    This gets the range of radii that can be used in the parameterization. The
    lower bound is dynamically adjusted based on the step size and the max value
    between the inner radii of either protein. The upper bound is also dynamically
    adjusted based on the step size and the sum of the outer radii of both proteins.
    The bounds are then used to generate a range of radii with the specified step.

    Parameters
    ----------
    1. protein_a : DynamicProteinAtomic or CGProtein
                   The first protein.

    2. protein_b : DynamicProteinAtomic or CGProtein
                   The second protein.

    3. step : float (Units: Angstroms)
              The step size to use in the range.

    4. cushion : float (Units: Angstroms)
                 The safety margin to use in the range to increase both bounds by.

    Returns
    -------
    1. r_range : np.ndarray
                 The range of radii to use in the parameterization.
    """

    min_r = np.floor(-cushion + max(protein_a.inner_radius, protein_b.inner_radius))
    max_r = np.ceil(cushion + (protein_a.radius + protein_b.radius))
    adjusted_max_r = max_r + (step - 0.01)

    return np.arange(min_r, adjusted_max_r, step)


def protein_clashes(
        protein_a: DynamicProteinAtomic,
        protein_b: DynamicProteinAtomic,
        phi_b: float,
        overlap_limit=None,
        sparse = True
):
    """
    This function is used to calculate the number of clashes between two
    proteins with set degrees of freedom with the inputs processed.

    Parameters
    ----------
    1. protein_a : DynamicProteinAtomic
                    The first protein in the pair.

    2. protein_b : DynamicProteinAtomic
                    The second protein in the pair.

    3. r : float
          The distance to set the proteins apart.

    4. phi_b : float (Units: Radians)
               The angle to set the second protein.

    5. overlap_limit : int
                       Maximum number of atoms in overlapped chains

    Returns
    -------
    1. clashes_list : list
                      [radius, phi_a, phi_b, clash_count]
    """
    # Rotate Protein_b correctly
    protein_b.rotate_by_axis(phi_b, 'z')

    if sparse:
        clash = calculate_atomic_clashes_sparse(protein_a, protein_b, overlap_limit)
    else:
        clash = calculate_atomic_clashes(protein_a, protein_b)

    # Undo Protein_b Rotation
    protein_b.rotate_by_axis(-phi_b, 'z')

    return clash


def _parallelize_protein_default_clashes(
        protein_a: str,
        protein_b: str,
        species: str,
        r: float,
        phi_a: float,
        list_phi_b: np.array,
        scratch: str,
) -> None:
    """
    This is a helper function to parallelize protein clashes over r. This means
    that each core will work on a specific r over the ranges of angles.

    Parameters
    ----------
    1. protein_a : str
                   The first protein in the pair.

    2. protein_b : str
                  The second protein in the pair.

    3. species : str
                 the species associated with both proteins

    3. r : float
          The distance to set the proteins apart.

    4. phi_range : np.ndarray
                   The range of angles to set the proteins.

    5. scratch : str
                 The path to the scratch directory.

    Returns
    -------
    """
    print(f'Running: {r} {phi_a}')
    protein_a_atomic = prepare_thylakoid_protein_default(species, protein_a)
    protein_a_atomic.rotate_by_axis(phi_a, 'z')
    protein_b_atomic = prepare_thylakoid_protein_default(species, protein_b)
    protein_b_atomic.translate(np.array([r, 0, 0]))
    results = [protein_clashes(protein_a_atomic, protein_b_atomic, phi_b, overlap_limit=10**4) for phi_b in list_phi_b]
    if not scratch is None:
        np.save(Path(scratch, f'{r}_{np.rad2deg(phi_a)}.npy'), [list_phi_b, results])


def parameterize_default(
        protein_a: str,
        protein_b: str,
        species: str,
        phi_steps: float = 5,
        r_steps: float = 5,
        r_range: Union[List[float], np.ndarray] = None,
        cores: int = None,
        path: str = None,
        cleanup: bool = True
) -> np.ndarray:
    """
    Parameterize a pair of proteins across a range of distances and angles for
    use with the force-field.

    Parameters
    ----------
    1. protein_a : str
                    The first protein in the pair.

    2. protein_b : str
                    The second protein in the pair.

    3. species : str
                 Species type to use for both proteins

    3. phi_steps : float [unit: degrees]
                   The number of steps to use in the range of angles.

    4. r_steps : float [unit: Angstroms]
                 The step size to use in the range of radii.

    5. r_range : List[float] or np.ndarray [unit: Angstroms]
                 The range of radii to use in the parameterization.

    6. cores : int
               The number of cores to parallelize over.

    7. path : str
              The path to save the parameterization to.

    8. cleanup : bool
                 Whether to remove the scratch directory after run.

    Returns
    -------
    1. clash_list : np.ndarray
                    A list of [r, phi_a, phi_b, clashes] for the ranges of
                    each degree of freedom.
    """
    if r_range is None:
        protein_a_atomic = prepare_thylakoid_protein_default(species, protein_a)
        protein_b_atomic = prepare_thylakoid_protein_default(species, protein_b)
        r_range = get_r_range(protein_a_atomic, protein_b_atomic, r_steps)
        print(f'r_range: {r_range}')
    phi_range = np.arange(0, 2 * np.pi, np.deg2rad(phi_steps))

    if path is None:
        path = Path.cwd()
    scratch = Path(path, f'{protein_a}_{protein_b}_parameterization')
    os.makedirs(scratch, exist_ok=True)

    tasks = [(protein_a, protein_b, species, r, phi_a, phi_range, scratch) for r in r_range for phi_a in phi_range]

    with Pool(processes=cores) as pool:
        for result in pool.starmap(_parallelize_protein_default_clashes, tasks):
            pass

    clash_list = np.array([np.load(f) for f in scratch.glob('*.npy')])

    if cleanup:
        for f in scratch.glob('*.npy'):
            f.unlink()
        if not any(scratch.iterdir()):
            scratch.rmdir()
        else:
            warnings.warn(f'{scratch} could not be deleted as it is non-empty.')

    phi_range = np.rad2deg(phi_range)
    final_data = np.array([r_range, phi_range, phi_range, clash_list], dtype=object)

    np.save(Path(path, f'clash_1d_{protein_a}_{protein_b}.npy'), final_data, allow_pickle=True)
    return final_data


# ========================= SIGMA R PARAMETERIZATION =========================


@njit
def _compute_sigma12(r_1: float, r_2: float):
    """
    The sigma in the LJ potential represents the diameter of an atom as it is
    when V(r) = 0 where r = σ. However, when describing interactions between
    two different atoms, one can compute the effective interaction mimicking
    the form of the LJ potential.

    Here, we will use the simple Lorentz-Berthelot rule. In our function,
    r_1 and r_2 are the RADII of the two atoms.

    Parameters
    ----------
    1. r_1 : float
             The radius of the first atom.
    2. r_2 : float
             The radius of the second atom.
    """
    return r_1 + r_2


@njit
def _compute_epsilon12(e_1: float, e_2: float):
    """
    See above doc string. Here, e_1 and e_2 are the LJ depths of the two atoms.
    We will use the simple Lorentz-Berthelot rule.

    Parameters
    ----------
    1. e_1 : float
             The LJ depth of the first atom.
    2. e_2 : float
             The LJ depth of the second atom.

    Returns
    -------
    1. epsilon : float
                 The effective epsilon between the two atoms.
    """
    return np.power(e_1 * e_2, 0.5)


@njit
def _total_wca_potential(
        coords_a: np.ndarray,
        coords_b: np.ndarray,
        radii_a: np.ndarray,
        radii_b: np.ndarray,
        epsilons_a: np.ndarray,
        epsilons_b: np.ndarray,
        cutoff: float
):
    potential = 0.0

    for coord_b, radius_b, epsilon_b in zip(coords_b, radii_b, epsilons_b):

        sigmas = _compute_sigma12(radii_a, radius_b)
        epsilons = _compute_epsilon12(epsilons_a, epsilon_b)

        distances = np.sqrt(np.sum((coords_a - coord_b) ** 2, axis=1))

        potential_matrix = epsilons * (4 * (np.power(sigmas / distances, 12) - np.power(sigmas / distances, 6)) + 1)
        potential_matrix = np.where(distances < np.power(2, 1/6) * sigmas, potential_matrix, 0)
        potential += np.sum(potential_matrix)

        if potential > cutoff:
            return cutoff

    return potential


def wca_potential(
        protein_a: Union[DynamicProteinAtomic, CGProtein],
        protein_b: Union[DynamicProteinAtomic, CGProtein],
        r: float,
        offset: float,
        cutoff: float = 1e9
):

    if np.sum(np.abs(protein_a.location)) > precision:
        protein_a.translate(-protein_a.location)

    protein_b.translate(np.array([r, 0, 0])-protein_b.location)

    coords_a = protein_a.R2_atomic_xyz
    coords_b = protein_b.R2_atomic_xyz

    sigma_a = protein_a.list_wca_sigmas
    sigma_b = protein_b.list_wca_sigmas

    epsilons_a = protein_a.list_wca_epsilons
    epsilons_b = protein_b.list_wca_epsilons

    return _total_wca_potential(coords_a, coords_b,
                                sigma_a, sigma_b,
                                epsilons_a, epsilons_b,
                                cutoff) - offset


def __root_potential(r, cg_lhcii_1, cg_lhcii_2, offset, cutoff):
    return wca_potential(cg_lhcii_1, cg_lhcii_2, r, offset, cutoff)


def determine_sigma_r(
        protein_a: Union[DynamicProteinAtomic, CGProtein],
        protein_b: Union[DynamicProteinAtomic, CGProtein],
        offset: float,
        root_guess: int = 80,
        cutoff: float = 1e9
):
    root_search = sc.optimize.root_scalar(__root_potential,
                                          args=(protein_a, protein_b, offset, cutoff),
                                          bracket=(30, 200), x0=root_guess, xtol=0.01, rtol=0.01)

    return root_search.root


def scan_sigmar_over_phi_b(
        protein_a: Union[DynamicProteinAtomic, CGProtein],
        protein_b: Union[DynamicProteinAtomic, CGProtein],
        phi_a: float,
        list_phi_b: Union[np.ndarray, List[float]],
        offset: float,
        root_guess: int,
        cutoff: float
):
    protein_a.rotate_by_axis(phi_a, 'z')

    list_sigma_r = []
    for phi_b in list_phi_b:
        protein_b.rotate_by_axis(phi_b, 'z')
        sigma_r = determine_sigma_r(protein_a, protein_b, offset, root_guess, cutoff)
        protein_b.rotate_by_axis(-phi_b, 'z')
        list_sigma_r.append(sigma_r)

    return list_sigma_r


def save_scan_sigmar_overphib(path_save,
                              atomic_protein_a, atomic_protein_b,
                              phi_a, list_phib,
                              offset, root_guess, cutoff):
    list_sigma_r = scan_sigmar_over_phi_b(atomic_protein_a, atomic_protein_b,
                                          phi_a, list_phib,
                                          offset, root_guess, cutoff)
    np.save(Path(path_save, f'sigmar_phia{phi_a}.npy'), list_sigma_r)


def construct_sigmar_scan(atomic_protein_a, atomic_protein_b,
                          list_phia, list_phib, offset,
                          path=None, root_guess=80, cutoff=10**9,
                          n_cores=None, cleanup=True):
    if path is None:
        path = Path.cwd()
    scratch = Path(path, f'{atomic_protein_a.name}_{atomic_protein_b.name}_parameterization')
    os.makedirs(scratch, exist_ok=True)
    print(f'Directory made:{scratch}')

    tasks = [(scratch, atomic_protein_a, atomic_protein_b,
              phi_a, list_phib, offset, root_guess, cutoff)
             for phi_a in list_phia]

    with Pool(processes=n_cores) as pool:
        pool.starmap(save_scan_sigmar_overphib, tasks)

    S2_sigma_r = np.array([np.load(f) for f in scratch.glob('*.npy')])

    if cleanup:
        for f in scratch.glob('*.npy'):
            f.unlink()
        if not any(scratch.iterdir()):
            scratch.rmdir()
        else:
            warnings.warn(f'{scratch} could not be deleted as it is non-empty.')

    np.savez(Path(path, f'sigmar_{atomic_protein_a.name}_{atomic_protein_b.name}.npz'),
             list_phia=list_phia, list_phib=list_phib, S2_sigma_r=S2_sigma_r)
    return S2_sigma_r
