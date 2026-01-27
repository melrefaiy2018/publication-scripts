import numpy as np
from Bio.PDB.Residue import Residue as BioResidue

from pymembrane.structure.atomic_pigment import ResidueAtomic
from pymembrane.util.physical_constants import dipole_CC, tresp_CC
from Bio.PDB.Residue import Residue as BioResidue
from pymembrane.structure.atomic_pigment import ResidueAtomic


def scale_charges(list_partial_charges, list_atom_positions, vacuum_mag):
    """
    Calculates the scaling factor for the transition charges of a single pigment. The equation used to calculate the
    transition dipole moment from the atomic partial charges is found in text in this paper:
    https://doi.org/10.1007/s11120-007-9248-z two paragraphs above equation 5.

    Parameters
    ----------
    1. list_partial_charges : list(float)
                            List of the partial charges that correspond to each
                            atom in the pigment.

    2. list_atom_positions : list(tuple(float))
                            List containing the x, y, z position of each atom in
                            the pigment.

    3. vacuum_mag : float [unit: e]
                    Magnitude of the dipole vector in vacuum.

    Returns
    -------
    1. scaling_factor : float
                        Scaling factor that adjusts all of the transition charges
                        so the dipole moment of the pigment matches the experimental
                        dipole moment.
    """
    dipole_mag = np.linalg.norm(sum([list_partial_charges[atom]
                                     * list_atom_positions[atom] for atom in
                                        range(len(list_atom_positions))]))
    return dipole_mag / vacuum_mag

def calculate_tresp_coupling(pigA, pigB, f_val, dielectric, coupling_const=tresp_CC):
    """
    Calculates the coupling between two pigments using TrEsp.

    PARAMETERS
    ----------
    1. pigA : object(PigmentAtomic)

    2. pigB : object(PigmentAtomic)

    3. f_val : float
               Screening/local field correction factor. Renger introduces this constant in this paper,
               https://doi.org/10.1021/ja9072222, under equation 9: "a factor f was introduced to take into account
               the electronic polarizability of the environment in an implicit way as discussed above."

    4. dielectric : int
                    Relative dielectric constant of the protein environment.

    RETURNS
    --------
    1. pigment_coupling_cm : float
                            Coupling between pigments A and B (cm^-1).

    RAISES
    ------
    ValueError : If required pigment data is missing, invalid, or atoms are too close
    KeyError : If required dictionary keys are not present
    """
    # Helper function to get pigment name for error messages
    def get_pig_name(pig):
        return getattr(pig, 'name', str(pig))
    
    pigA_name = get_pig_name(pigA)
    pigB_name = get_pig_name(pigB)
    
    # Validate inputs - check for None first
    if f_val is None:
        raise ValueError("f_val cannot be None")
    if f_val <= 0:
        raise ValueError(f"f_val must be positive, got {f_val}")
    
    if dielectric is None:
        raise ValueError("dielectric constant cannot be None. "
                        "Please set self.dielectric_tresp attribute before calling this function.")
    if dielectric <= 0:
        raise ValueError(f"dielectric constant must be positive, got {dielectric}")
    
    # Validate pigment A data
    try:
        if 'tresp_atoms' not in pigA.dict_data:
            raise KeyError(f"Pigment '{pigA_name}' missing 'tresp_atoms' in dict_data")
        if 'tresp_pc' not in pigA.dict_data:
            raise KeyError(f"Pigment '{pigA_name}' missing 'tresp_pc' in dict_data")
        if 'vacuum_mag' not in pigA.dict_data:
            raise KeyError(f"Pigment '{pigA_name}' missing 'vacuum_mag' in dict_data")
            
        n_pigA_atoms = len(pigA.dict_data['tresp_atoms'])
        if n_pigA_atoms == 0:
            raise ValueError(f"Pigment '{pigA_name}' has no TrEsp atoms defined")
        if len(pigA.dict_data['tresp_pc']) != n_pigA_atoms:
            raise ValueError(f"Pigment '{pigA_name}' has {n_pigA_atoms} atoms but "
                           f"{len(pigA.dict_data['tresp_pc'])} charges")
        if pigA.dict_data['vacuum_mag'] == 0:
            raise ValueError(f"Pigment '{pigA_name}' has vacuum_mag = 0")
            
    except (AttributeError, TypeError) as e:
        raise ValueError(f"Pigment A ('{pigA_name}') has invalid or missing dict_data: {e}") from e
    
    # Validate pigment B data
    try:
        if 'tresp_atoms' not in pigB.dict_data:
            raise KeyError(f"Pigment '{pigB_name}' missing 'tresp_atoms' in dict_data")
        if 'tresp_pc' not in pigB.dict_data:
            raise KeyError(f"Pigment '{pigB_name}' missing 'tresp_pc' in dict_data")
        if 'vacuum_mag' not in pigB.dict_data:
            raise KeyError(f"Pigment '{pigB_name}' missing 'vacuum_mag' in dict_data")
            
        n_pigB_atoms = len(pigB.dict_data['tresp_atoms'])
        if n_pigB_atoms == 0:
            raise ValueError(f"Pigment '{pigB_name}' has no TrEsp atoms defined")
        if len(pigB.dict_data['tresp_pc']) != n_pigB_atoms:
            raise ValueError(f"Pigment '{pigB_name}' has {n_pigB_atoms} atoms but "
                           f"{len(pigB.dict_data['tresp_pc'])} charges")
        if pigB.dict_data['vacuum_mag'] == 0:
            raise ValueError(f"Pigment '{pigB_name}' has vacuum_mag = 0")
            
    except (AttributeError, TypeError) as e:
        raise ValueError(f"Pigment B ('{pigB_name}') has invalid or missing dict_data: {e}") from e
    
    # Get atom positions with error handling
    try:
        R2_pigA_pos = np.array([pigA.get_atom_xyz(pigA.dict_data['tresp_atoms'][atom]) for
                        atom in range(n_pigA_atoms)])
    except KeyError as e:
        raise KeyError(f"Pigment '{pigA_name}' missing atom in structure: {e}") from e
    except Exception as e:
        raise ValueError(f"Error getting atom positions for '{pigA_name}': {e}") from e
        
    try:
        R2_pigB_pos = np.array([pigB.get_atom_xyz(pigB.dict_data['tresp_atoms'][atom]) for
                        atom in range(n_pigB_atoms)])
    except KeyError as e:
        raise KeyError(f"Pigment '{pigB_name}' missing atom in structure: {e}") from e
    except Exception as e:
        raise ValueError(f"Error getting atom positions for '{pigB_name}': {e}") from e
    
    # Check for NaN or infinite coordinates
    if not np.all(np.isfinite(R2_pigA_pos)):
        raise ValueError(f"Pigment '{pigA_name}' has invalid (NaN/Inf) coordinates")
    if not np.all(np.isfinite(R2_pigB_pos)):
        raise ValueError(f"Pigment '{pigB_name}' has invalid (NaN/Inf) coordinates")
    
    # Scale charges
    try:
        pigA_scale = scale_charges(pigA.dict_data['tresp_pc'], R2_pigA_pos,
                                    pigA.dict_data['vacuum_mag'])
        if pigA_scale == 0 or not np.isfinite(pigA_scale):
            raise ValueError(f"Invalid charge scaling for '{pigA_name}': {pigA_scale}")
    except Exception as e:
        raise ValueError(f"Error scaling charges for '{pigA_name}': {e}") from e
        
    try:
        pigB_scale = scale_charges(pigB.dict_data['tresp_pc'], R2_pigB_pos,
                                    pigB.dict_data['vacuum_mag'])
        if pigB_scale == 0 or not np.isfinite(pigB_scale):
            raise ValueError(f"Invalid charge scaling for '{pigB_name}': {pigB_scale}")
    except Exception as e:
        raise ValueError(f"Error scaling charges for '{pigB_name}': {e}") from e
    
    # Calculate scaled charges
    Q2_pigA = np.array([[pigA.dict_data['tresp_pc'][atom_a] / pigA_scale for
                            atom_a in range(n_pigA_atoms)]])
    Q2_pigB = np.array([[pigB.dict_data['tresp_pc'][atom_b] / pigB_scale for
                            atom_b in range(n_pigB_atoms)]])

    Q2_pigA_tile = np.tile(Q2_pigA.T, (1, n_pigB_atoms))
    Q2_pigB_tile = np.tile(Q2_pigB, (n_pigA_atoms, 1))

    R3_pigA = np.tile(R2_pigA_pos, (1, 1, n_pigB_atoms)).reshape(n_pigA_atoms,
                                                                n_pigB_atoms, 3)
    R3_pigB = np.tile(R2_pigB_pos, (n_pigA_atoms, 1, 1))

    R2_ab_norm = np.linalg.norm((R3_pigA-R3_pigB), axis=2)
    
    # Check for atoms that are too close (distance < 0.1 Angstrom)
    min_distance = np.min(R2_ab_norm)
    if min_distance < 0.1:
        raise ValueError(f"Atoms between '{pigA_name}' and '{pigB_name}' are too close "
                        f"(min distance: {min_distance:.4f} Å). This may indicate "
                        f"overlapping structures or invalid coordinates.")
    
    # Check for zero distances (would cause division by zero)
    if np.any(R2_ab_norm == 0):
        raise ValueError(f"Zero distance found between atoms of '{pigA_name}' and '{pigB_name}'. "
                        f"Pigments may be overlapping.")

    # Calculate coupling
    V_ab = coupling_const * f_val/dielectric * \
           np.sum((Q2_pigA_tile * Q2_pigB_tile)/R2_ab_norm)
    
    # Check for invalid result
    if not np.isfinite(V_ab):
        raise ValueError(f"Calculated coupling between '{pigA_name}' and '{pigB_name}' "
                        f"is invalid (NaN/Inf). Check input parameters and coordinates.")
    
    return V_ab


def calculate_dipole_coupling(pigA, pigB, dielectric, coupling_const=dipole_CC):
    """
    Calculates the dipole coupling between two pigments.

    PARAMETERS
    -----------
    1. pigA : object(PigmentAtomic)

    2. pigB : object(PigmentAtomic)

    3. coupling_constant : float
                            Coulomb constant ((Nm^2)/C^2).

    4. dielectric : int
                    Relative dielectric constant of the protein environment.

    RETURNS
    --------
    1. J_AB_cm : float
                    Coupling between pigments A and B (cm^-1).
    """
    pigA_dipole = pigA.get_dipole_dir('qy') * pigA.dict_data['dipole_mag']
    pigB_dipole = pigB.get_dipole_dir('qy') * pigB.dict_data['dipole_mag']
    distance_AB = (pigA.location - pigB.location)
    J_AB_cm = coupling_const/dielectric * ((np.dot(pigA_dipole, pigB_dipole) / (
            np.linalg.norm(distance_AB) ** 3)) - (3 * (np.dot(pigA_dipole, distance_AB)
            * (np.dot(pigB_dipole, distance_AB))) / (np.linalg.norm(distance_AB) ** 5)))
    return J_AB_cm


def name_background_object(background_atomic):
    if isinstance(background_atomic, ResidueAtomic):
        return background_atomic.name.strip()

    elif isinstance(background_atomic, BioResidue):
        return f"{background_atomic.full_id[2]}_{background_atomic.resname}" \
                f"_{background_atomic.id[1]}"
    
    # Handle tuple case (background_atomic might be a tuple of (name, object))
    elif isinstance(background_atomic, tuple) and len(background_atomic) == 2:
        return str(background_atomic[0])  # Use the name part of the tuple
    
    # Handle case where it might be a simple value or other object with string representation
    elif hasattr(background_atomic, '__str__'):
        return str(background_atomic)

    else:
        raise ValueError(f'Background type {type(background_atomic)} is not supported!')


def pigment_site_energy_shift_by_component(pigA, background_atomic, dielectric_eff):
    """
    Calculates the contribution of the protein on the site energy for the pigment.

    Parameters
    ----------
    pigA
    background_atomic

    Returns
    -------

    """
    # pigA number atoms, delta charge, and position
    # ----------------------------------------------
    n_pigA_atoms = len(pigA.dict_data['q_0011_atom'])
    q_pigA_delta = np.array([pigA.dict_data['q_11'][atom_a]
                            for atom_a in range(n_pigA_atoms)]) \
                            - np.array([pigA.dict_data['q_00'][atom_a] for atom_a in
                            range(n_pigA_atoms)])

    # For later transpose actions, needs to be a 2D array
    Q2_pigA_delta = np.reshape(q_pigA_delta, newshape=[1, n_pigA_atoms])

    R2_pigA_pos = []
    for atom in range(n_pigA_atoms):
        try:
            pos = pigA.get_atom_xyz(pigA.dict_data['q_0011_atom'][atom])
            R2_pigA_pos.append(pos)
        except KeyError as e:
            raise KeyError(
                f"Missing key '{e.args[0]}' in pigment '{getattr(pigA, 'name', pigA)}' at atom index {atom}. "
                f"Check pigment data: {pigA.dict_data.get('q_0011_atom', {})}"
            ) from e
    R2_pigA_pos = np.array(R2_pigA_pos)

    # Handle tuple case first
    if isinstance(background_atomic, tuple) and len(background_atomic) == 2:
        # Extract the actual object from the tuple
        background_atomic = background_atomic[1]

    # Background number atom, charge, and position (multiple cases)
    # -------------------------------------------------------------
    if isinstance(background_atomic, ResidueAtomic):
        # Case: Atomic Pigment Instance
        n_protein_atoms = len(background_atomic.dict_data['q_0011_atom'])

        Q2_protein_00 = np.reshape(np.array([background_atomic.dict_data['q_00'][atom_b] for atom_b in
                                             range(n_protein_atoms)]),
                                   newshape=[1, n_protein_atoms])
        R2_protein_pos = np.array([background_atomic.get_atom_xyz(background_atomic.dict_data[
                                                                      'q_0011_atom'][atom]) for
                                   atom in range(n_protein_atoms)])

    elif isinstance(background_atomic, BioResidue):
        # Case: Individual Protein Residue
        n_protein_atoms = len([atom for atom in background_atomic.get_atoms() if
                               atom.pqr_charge is not None])
        if n_protein_atoms == 0:
            raise ValueError(f"The atoms count of {background_atomic.resname} " +
                             f"molecule can't be equal zero")

        Q2_protein_00 = np.reshape(np.array([atom.pqr_charge for atom in
                                             background_atomic.get_atoms() if
                                             atom.pqr_charge is not None]),
                                   newshape=[1, n_protein_atoms])

        R2_protein_pos = np.array([atom.coord for atom in
                                   background_atomic.get_atoms() if
                                   atom.pqr_charge is not None])

    else:
        raise ValueError('Background type is not supported!')
    # Tile charge values: [n_pigA, n_protein]
    Q2_pigA_tile = np.tile(Q2_pigA_delta.T, (1, n_protein_atoms))  # pigA delta q
    Q2_background_tile = np.tile(Q2_protein_00, (n_pigA_atoms, 1))  # pigB q_00

    # tile pigA and protein_residues distance: [n_pigA, n_protein, 3]
    R3_pigA = np.tile(R2_pigA_pos, (1, 1, n_protein_atoms)).reshape(n_pigA_atoms,
                                                                    n_protein_atoms, 3)
    R3_pigB = np.tile(R2_protein_pos, (n_pigA_atoms, 1, 1))

    # Calculate distance between atoms pigA and protein atoms: [n_pigA, n_protein]
    R2_ab_norm = np.linalg.norm((R3_pigA - R3_pigB), axis=2)  # R2 vector

    # delta E
    return (tresp_CC / dielectric_eff) * np.sum(Q2_pigA_tile
                                                * Q2_background_tile
                                                / R2_ab_norm)


def detailed_residue_contribution(pigA, background_atomic, dielectric_eff, 
                                 significance_threshold=0.1, max_dominant_atoms=10):
    """
    Calculates detailed residue contribution with atom-level breakdown and identifies dominant atoms.
    
    Parameters
    ----------
    pigA : PigmentAtomic
        The pigment for which to calculate contributions
    background_atomic : ResidueAtomic or BioResidue
        The background residue/molecule contributing to the site energy
    dielectric_eff : float
        Effective dielectric constant
    significance_threshold : float, optional
        Threshold for considering an atom contribution significant (default: 0.1)
    max_dominant_atoms : int, optional
        Maximum number of dominant atoms to identify (default: 10)
        
    Returns
    -------
    dict
        Detailed contribution breakdown with format:
        {
            'total_contribution': float,
            'atom_breakdown': {atom_name: contribution_value, ...},
            'dominant_atoms': {atom_name: contribution_value, ...},
            'num_atoms': int,
            'num_significant_atoms': int
        }
    """
    
    # Handle special case for vacuum energy (CLA molecules)
    residue_name = name_background_object(background_atomic)
    if 'CLA' in residue_name:
        # For CLA molecules, return vacuum energy format
        vacuum_energy = pigment_site_energy_shift_by_component(pigA, background_atomic, dielectric_eff)
        return {
            'vacuum': {
                'total_contribution': float(vacuum_energy),
                'atom_breakdown': {'vacuum_energy': float(vacuum_energy)}
            }
        }
    
    # pigA number atoms, delta charge, and position
    n_pigA_atoms = len(pigA.dict_data['q_0011_atom'])
    q_pigA_delta = np.array([pigA.dict_data['q_11'][atom_a]
                            for atom_a in range(n_pigA_atoms)]) \
                            - np.array([pigA.dict_data['q_00'][atom_a] for atom_a in
                            range(n_pigA_atoms)])

    R2_pigA_pos = np.array([pigA.get_atom_xyz(pigA.dict_data['q_0011_atom'][atom])
                            for atom in range(n_pigA_atoms)])

    # Handle tuple case first
    if isinstance(background_atomic, tuple) and len(background_atomic) == 2:
        # Extract the actual object from the tuple
        background_atomic = background_atomic[1]
    
    # Background atom information
    if isinstance(background_atomic, ResidueAtomic):
        # Case: Atomic Pigment Instance
        n_protein_atoms = len(background_atomic.dict_data['q_0011_atom'])
        protein_charges = np.array([background_atomic.dict_data['q_00'][atom_b] 
                                   for atom_b in range(n_protein_atoms)])
        protein_positions = np.array([background_atomic.get_atom_xyz(
                                     background_atomic.dict_data['q_0011_atom'][atom]) 
                                     for atom in range(n_protein_atoms)])
        atom_names = [background_atomic.dict_data['q_0011_atom'][atom] 
                     for atom in range(n_protein_atoms)]
        
    elif isinstance(background_atomic, BioResidue):
        # Case: Individual Protein Residue
        charged_atoms = [atom for atom in background_atomic.get_atoms() 
                        if atom.pqr_charge is not None]
        n_protein_atoms = len(charged_atoms)
        
        if n_protein_atoms == 0:
            raise ValueError(f"The atoms count of {background_atomic.resname} "
                           f"molecule can't be equal zero")
        
        protein_charges = np.array([atom.pqr_charge for atom in charged_atoms])
        protein_positions = np.array([atom.coord for atom in charged_atoms])
        atom_names = [atom.name.strip() for atom in charged_atoms]
        
    else:
        raise ValueError('Background type is not supported!')
    
    # Calculate per-atom contributions
    atom_contributions = {}
    
    for j, atom_name in enumerate(atom_names):
        # Calculate contribution from this specific atom to all pigA atoms
        atom_contribution = 0.0
        
        for i in range(n_pigA_atoms):
            # Distance between pigA atom i and protein atom j
            distance = np.linalg.norm(R2_pigA_pos[i] - protein_positions[j])
            
            # Contribution from this atom pair
            contribution = (tresp_CC / dielectric_eff) * \
                          (q_pigA_delta[i] * protein_charges[j] / distance)
            atom_contribution += contribution
        
        atom_contributions[atom_name] = np.float64(atom_contribution)
    
    # Calculate total contribution
    total_contribution = sum(atom_contributions.values())
    
    # Identify significant atoms (above threshold)
    significant_atoms = {name: contrib for name, contrib in atom_contributions.items() 
                        if abs(contrib) >= significance_threshold}
    
    # Identify dominant atoms (largest absolute contributions)
    sorted_atoms = sorted(atom_contributions.items(), 
                         key=lambda x: abs(x[1]), reverse=True)
    dominant_atoms = dict(sorted_atoms[:max_dominant_atoms])
    
    # Filter dominant atoms to only include significant ones
    dominant_atoms = {name: contrib for name, contrib in dominant_atoms.items() 
                     if abs(contrib) >= significance_threshold}
    
    return {
        'total_contribution': np.float64(total_contribution),
        'atom_breakdown': atom_contributions,
        'dominant_atoms': dominant_atoms,
        'num_atoms': n_protein_atoms,
        'num_significant_atoms': len(significant_atoms)
    }


def calculate_cdc_detailed_contributions(pigA, background_list, dielectric_eff, 
                                       significance_threshold=0.1, max_dominant_atoms=10):
    """
    Calculate detailed CDC contributions for all background residues/molecules.
    
    Parameters
    ----------
    pigA : PigmentAtomic
        The pigment for which to calculate contributions
    background_list : list
        List of background residues/molecules (ResidueAtomic or BioResidue objects)
    dielectric_eff : float
        Effective dielectric constant
    significance_threshold : float, optional
        Threshold for considering an atom contribution significant (default: 0.1)
    max_dominant_atoms : int, optional
        Maximum number of dominant atoms to identify per residue (default: 10)
        
    Returns
    -------
    dict
        Dictionary with detailed contributions for each residue:
        {
            'residue_id': {
                'total_contribution': float,
                'atom_breakdown': {atom_name: contribution, ...},
                'dominant_atoms': {atom_name: contribution, ...},
                'num_atoms': int,
                'num_significant_atoms': int
            },
            ...
        }
    """
    
    detailed_contributions = {}
    
    for background_atomic in background_list:
        residue_id = "unknown"  # Initialize with default value
        try:
            # Generate unique residue identifier
            residue_id = name_background_object(background_atomic)
            
            # Calculate detailed contribution for this residue
            contribution_data = detailed_residue_contribution(
                pigA, background_atomic, dielectric_eff, 
                significance_threshold, max_dominant_atoms
            )
            
            detailed_contributions[residue_id] = contribution_data
            
        except Exception as e:
            print(f"Warning: Could not calculate contribution for {residue_id}: {e}")
            print(f"Background object type: {type(background_atomic)}")
            continue
    
    return detailed_contributions


def calculate_detailed_contributions_from_cdc_results(wscp_atomic, dict_total_contribution_site_energy, 
                                                     pigment_name=None, 
                                                     significance_threshold=0.1, max_dominant_atoms=10):
    """
    Calculate detailed CDC contributions using the results from calculate_cdc_site_shift.
    
    This is a convenience function that works directly with the output of 
    ElectrostaticProteinAtomic.calculate_cdc_site_shift().
    
    Parameters
    ----------
    wscp_atomic : ElectrostaticProteinAtomic
        The protein object containing pigments and residues
    dict_total_contribution_site_energy : dict
        Output from calculate_cdc_site_shift()[1] - dictionary mapping pigment names 
        to their background contribution dictionaries
    pigment_name : str, optional
        Name of specific pigment to analyze. If None, analyzes the first pigment.
    significance_threshold : float, optional
        Threshold for considering an atom contribution significant (default: 0.1)
    max_dominant_atoms : int, optional
        Maximum number of dominant atoms to identify per residue (default: 10)
        
    Returns
    -------
    dict
        Detailed contribution breakdown for the specified pigment
    """
    
    # Get the pigment object
    if pigment_name is None:
        pigment_obj = wscp_atomic.list_pigments[0]
    else:
        pigment_obj = None
        for pig in wscp_atomic.list_pigments:
            if pigment_name in pig.name:
                pigment_obj = pig
                break
        if pigment_obj is None:
            raise ValueError(f"Pigment {pigment_name} not found")
    
    # Build the background list (same as what calculate_cdc_site_shift uses)
    background_list = []
    
    # Add other pigments as background objects
    for pig in wscp_atomic.list_pigments:
        if pig != pigment_obj:
            background_list.append(pig)
    
    # Add residues that are not pigments as background objects
    from pymembrane.structure.atomic_protein import pigment_name_from_residue
    for res in wscp_atomic.get_residues():
        if pigment_name_from_residue(res) not in wscp_atomic.dict_pigments.keys():
            background_list.append(res)
    
    # Calculate detailed contributions
    detailed_results = calculate_cdc_detailed_contributions(
        pigA=pigment_obj,
        background_list=background_list,
        dielectric_eff=wscp_atomic.dielectric_cdc,
        significance_threshold=significance_threshold,
        max_dominant_atoms=max_dominant_atoms
    )
    
    return detailed_results

