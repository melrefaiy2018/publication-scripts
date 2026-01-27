import copy
import csv
import os
from pathlib import Path
from typing import Literal, List, Optional, Type, Union,Dict
from warnings import warn

import numpy as np
import pandas as pd
from Bio import PDB
from matplotlib import pyplot as plt
from matplotlib.pyplot import Axes
from mpl_toolkits.mplot3d import art3d
from scipy.spatial import ConvexHull
from shapely import MultiPoint, Polygon
from shapely import concave_hull as shapely_concave_hull
from shapely import convex_hull as shapely_convex_hull

from pymembrane.exciton.hamiltonian import Hamiltonian
from pymembrane.structure.electrostatics import (calculate_tresp_coupling,
                                                 calculate_dipole_coupling,
                                                 name_background_object,
                                                 pigment_site_energy_shift_by_component)
from pymembrane.util.fake_objects_for_testing import FakePDB
from pymembrane.util.linear_algebra import construct_rotation_for_flatten, construct_rotation_matrix
from pymembrane.util.load_hamiltonian import load_hamiltonian, match_pigment_names, update_pigment_name
from pymembrane.util.physical_constants import dict_covalent_radii, dict_vdw_radii, dict_wca_sigma, dict_wca_epsilon


def load_protein(name, prot_type, parameter_mod, file_pdb, class_type=None):
    """Load a protein from a PDB file.

    Args:
        name: The name of the protein.
        prot_type: The type of the protein.
        parameter_mod: The path to the parameter files.
        file_pdb: The name of the PDB file.
        class_type: The class type of the protein. Options are 'StaticProteinAtomic',
            'DynamicProteinAtomic', 'PigmentProteinAtomic', or 'ElectrostaticProteinAtomic'.
            Defaults to 'DynamicProteinAtomic' if None.

    Returns:
        An instance of the specified protein class (StaticProteinAtomic, DynamicProteinAtomic,
        PigmentProteinAtomic, or ElectrostaticProteinAtomic).
    """
    path = Path(parameter_mod)
    path_pdb = str(path.joinpath(file_pdb))
    pdb_atomic = PDB.PDBParser().get_structure(name, path_pdb)
    if class_type == 'StaticProteinAtomic':
        return StaticProteinAtomic(pdb_atomic, name=name, prot_type=prot_type)
    elif class_type is None or class_type == 'DynamicProteinAtomic':
        return DynamicProteinAtomic(pdb_atomic, name=name, prot_type=prot_type)
    elif class_type == 'PigmentProteinAtomic':
        return PigmentProteinAtomic(pdb_atomic, name=name, prot_type=prot_type)
    elif class_type == 'ElectrostaticProteinAtomic':
        return ElectrostaticProteinAtomic(pdb_atomic, name=name, prot_type=prot_type)
    else:
        #TODO: Replace with correct exception
        print(f'WARNING: load_protein() from atomic_protein.py class_type={class_type} is not allowed.')



def pigment_name_from_residue(res):
    """Return the name of a pigment from a residue.

    Args:
        res: Bio.PDB.Residue.Residue object.

    Returns:
        The name of the pigment as a string in the format:
        '{structure_name}_{chain_id}_{resname}_{res_id}'.
    """
    return f'{res.get_parent().get_parent().full_id[0]}_{res.get_parent().id}_{res.resname}_{res.id[1]}'


class StaticProteinAtomic:
    __slots__ = ['_pdb', '__name', '__type', 'dict_data', 'dict_pigments',
                 'H2_hamiltonian', '_atom_coords', '_atom_masses', '_mass',
                 '_center_of_mass', '_R2_residue_center', '__list_vdw_radii',
                 '__list_covalent_radii', '__list_wca_sigmas', '__list_wca_epsilons',
                 '__hull', '__hull2d', '__hull2d_radius', '__concave_hull',
                 '__dict_hull2d_by_chain', '__inner_radius', '_dict_index_atom_by_chain',
                 '_dict_index_res_by_chain',
                 '_fval_tresp', '_dielectric_tresp', '_dielectric_dipole']

    def __init__(
            self,
            pdb: Union[str, PDB.Structure.Structure, FakePDB],
            name: str,
            prot_type: str = None
    ) -> None:
        """Initialize a StaticProteinAtomic object.

        This is the base class for the set of atomic protein classes used to represent
        proteins. The base of the atomic representation is the BioPython PDB object.

        Roles:
            1. PDB Data Access: All calls that access data FROM the pdb object should use
               this class without calling to PDB.
            2. Protein Identity: This class is responsible for the identity of the atomic
               protein.
            3. Protein Shape: This class will be responsible for defining the protein shape
               (convex hull, etc.).

        Args:
            pdb: Either a PDB.Structure.Structure object instantiated by the Structure class
                of the Biopython package, a path (str) pointing to the PDB file, or a
                FakePDB object for testing.
            name: The name of the ProteinAtomic object.
            prot_type: Identifies the 'type' of protein for membrane simulations. If None,
                defaults to the value of name.

        Raises:
            ValueError: If pdb is not a valid type (PDB.Structure.Structure, str, or FakePDB).
        """
        if isinstance(pdb, str):
            self._pdb = PDB.PDBParser(QUIET=True).get_structure(str(name), pdb)
        elif isinstance(pdb, PDB.Structure.Structure):
            self._pdb = pdb
        # FakePDB is used for testing.
        elif isinstance(pdb, FakePDB):
            self._pdb = pdb
        else:
            raise ValueError('Incorrect type of PDB passed. Please input a PDB'
                             ' object or a path to a PDB file.')

        if prot_type is None:
            prot_type = name

        # Define Name and type of Protein
        # -------------------------------
        self.__name = str(name)
        self.__type = prot_type

        # Define General Data Storage Structure
        # -------------------------------------
        # Note: This data cannot depend on the position/orientation of particle
        self.dict_data = {'MONTECARLO_TRANSLATE': True}

        # Define General Storage Structure for Pigments
        # ---------------------------------------------
        self.dict_pigments = {}

        # Define Excitonic Hamiltonian Storage
        # ------------------------------------
        self.H2_hamiltonian = None

        # Define Center-of-Mass
        # ---------------------
        atom_coords = np.array([atom.coord for atom in self._pdb.get_atoms()])
        self._atom_masses = np.array([atom.mass for atom in self._pdb.get_atoms()])
        self._mass = np.sum(self._atom_masses)

        self._center_of_mass = np.sum(atom_coords * self._atom_masses[:, None], axis=0) / self.mass

        self._R2_residue_center = np.array(
            [np.sum([atom.coord * atom.mass for atom in res.get_atoms()], axis=0)
             / np.sum([atom.mass for atom in res.get_atoms()])
             for res in self.get_residues()])

        self._dict_index_atom_by_chain = {chain.id : [index for (index, atom) in enumerate(self.get_atoms())
                                                      if atom.get_parent().get_parent().id == chain.id]
                                          for chain in self._pdb.get_chains()}

        self._dict_index_res_by_chain = {chain.id : [index for (index, res) in enumerate(self.get_residues())
                                                     if res.get_parent().id == chain.id]
                                         for chain in self._pdb.get_chains()}

        # Define Atomic Properties
        # ------------------------
        self.__list_vdw_radii = None
        self.__list_covalent_radii = None
        self.__list_wca_sigmas = None
        self.__list_wca_epsilons = None

        # Define Shape Properties
        # ------------------------
        self._reset_pdb()

        # Define Coupling Parameters
        # --------------------------
        self._fval_tresp = None
        self._dielectric_tresp = None
        self._dielectric_dipole = None

    def __eq__(self, other) -> bool:
        if isinstance(other, self.__class__):
            return self.name == other.name
        else:
            return False

    def __str__(self) -> str:
        return f'{type(self).__name__} {self.name}'

    def __repr__(self) -> str:
        return f'''{type(self).__name__} {self.name}
                   Mass: {self.mass}
                   Location: {self.center_of_mass}'''

    def __hash__(self):
        return hash(self.type)

    def _reset_pdb(self):
        # Some of these depend on the position and orientation of
        # the pdb object and need to be accesible to classes that
        # modify the pdb structure.
        self.__hull = None
        self.__hull2d = None
        self.__dict_hull2d_by_chain = None
        self.__hull2d_radius = None
        self.__concave_hull = None
        self.__inner_radius = None

    @property
    def mass(self) -> float:
        """Mass of the Protein_atomic object.

        Returns:
            Mass in atomic mass units (a.m.u.).
        """
        return self._mass

    @property
    def name(self) -> str:
        """Name of the protein.

        Returns:
            The protein name.
        """
        return self.__name

    @name.setter
    def name(self, new_name: str) -> None:
        """Set the name of the protein.

        Args:
            new_name: The new protein name.
        """
        self.__name = new_name

    @property
    def type(self) -> str:
        """Type of the protein.

        Returns:
            The protein type identifier.
        """
        return self.__type

    @type.setter
    def type(self, new_type: str) -> None:
        """Set the type of the protein.

        Args:
            new_type: The new protein type identifier.
        """
        self.__type = new_type

    @property
    def center_of_mass(self) -> np.ndarray:
        """Get the location of the protein (as determined by its center of mass).

        Returns:
            Center of mass of the protein in Angstroms.
        """
        return self._center_of_mass

    @property
    def location(self) -> np.ndarray:
        """Get the location of the protein (as determined by its center of mass).

        Returns:
            Center of mass of the protein in Angstroms.
        """
        return self._center_of_mass

    @property
    def inertia_tensor(self) -> np.ndarray:
        """Constructs the inertia tensor matrix for the protein.
        
        The inertia tensor summarizes all moments of inertia of an object with one quantity.
        For a rigid object of N point masses m_k, the moment of inertia tensor is defined as:
        
            I_ij = Σ_(k=1)^N [m_k * (||r_k||² * δ_ij - x_i^k * x_j^k)]
        
        where the tensor has the form:
            | I_11  I_12  I_13 |
            | I_21  I_22  I_23 |
            | I_31  I_32  I_33 |

        See: https://en.wikipedia.org/wiki/Moment_of_inertia#Inertia_tensor

        Returns:
            Inertia tensor matrix in units of a.m.u. * Angstrom^2.
        """
        # self._atom_coords = np.array([atom.coord for atom in self._pdb.get_atoms()])
        atom_pos = self._atom_coords - self.center_of_mass
        atom_masses = self._atom_masses

        # Construct Inertial Matrix
        # -------------------------
        xx = np.sum(atom_masses * (atom_pos[:, 1] ** 2 + atom_pos[:, 2] ** 2))
        yy = np.sum(atom_masses * (atom_pos[:, 0] ** 2 + atom_pos[:, 2] ** 2))
        zz = np.sum(atom_masses * (atom_pos[:, 0] ** 2 + atom_pos[:, 1] ** 2))
        xy = -np.sum(atom_masses * atom_pos[:, 0] * atom_pos[:, 1])
        xz = -np.sum(atom_masses * atom_pos[:, 0] * atom_pos[:, 2])
        yz = -np.sum(atom_masses * atom_pos[:, 1] * atom_pos[:, 2])

        return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]])

    @property
    def fval_tresp(self):
        warn('To properly use fval_tresp, please use PigmentProteinAtomic instead.')
        return None

    @property
    def dielectric_tresp(self):
        warn('To properly use dielectric_tresp, please use PigmentProteinAtomic instead.')
        return None

    @property
    def dielectric_dipole(self):
        warn('To properly use dielectric_dipole, please use PigmentProteinAtomic instead.')
        return None

    def get_models(self) -> List[PDB.Model.Model]:
        """Get a list of models from the PDB.

        Returns:
            List of PDB.Model objects from the PDB.
        """
        return self._pdb.get_models()

    def get_chains(self) -> List[PDB.Chain.Chain]:
        """Get a list of chains from the PDB.

        Returns:
            List of PDB.Chain objects from the PDB.
        """
        return self._pdb.get_chains()

    def get_residues(self) -> List[PDB.Residue.Residue]:
        """Get a list of residues from the PDB.

        Returns:
            List of PDB.Residue objects from the PDB.
        """
        return self._pdb.get_residues()

    def get_atoms(self) -> List[PDB.Atom.Atom]:
        """Get a list of atoms from the PDB.

        Returns:
            List of PDB.Atom objects from the PDB.
        """
        return self._pdb.get_atoms()

    def construct_subcomplex(
            self,
            list_chains_keep: List[str],
            name: str,
            class_type: Type['StaticProteinAtomic'],
            model_index: Optional[int] = None,
            center_atomic: bool = True
    ) -> 'StaticProteinAtomic':
        """Construct a new Protein_atomic object containing only specified chains.
        
        This does not alter the original ProteinAtomic object.

        Args:
            list_chains_keep: Chain names from the PDB to keep in the subcomplex.
            name: Name of the new subcomplex.
            class_type: Type of subcomplex (StaticProteinAtomic, DynamicProteinAtomic,
                PigmentProteinAtomic, or ElectrostaticProteinAtomic).
            model_index: Model index from the Biopython.PDB object. For structures from
                x-ray crystallography, this can be ignored.
            center_atomic: Whether to center ProteinAtomic on the origin (True) or to
                keep the same position (False).

        Returns:
            Instance of the specified ProteinAtomic class type.

        Raises:
            ValueError: If model_index is not specified when multiple models exist,
                if an invalid model_index is provided, if specified chains are not
                in the PDB, or if an unsupported class_type is provided.
        """
        new_pdb = copy.deepcopy(self._pdb)

        if len(new_pdb.child_list) != 1 and model_index is None:
            raise ValueError(f'Must specify a model index.')

        if len(new_pdb.child_list) == 1 and model_index is not None and model_index > 0:
            raise ValueError(f'Your PDB only has one model. An index of {model_index} is invalid.')

        model = (new_pdb.child_list[model_index]
                 if model_index is not None
                 else new_pdb.child_list[0])

        available_chain_ids = [chain.id for chain in model.child_list]
        for chain in list_chains_keep:
            if chain not in available_chain_ids:
                raise ValueError(f'Chains {list_chains_keep} not in PDB.')

        for chain_id in available_chain_ids:
            if chain_id not in list_chains_keep:
                model.detach_child(chain_id)

        if (class_type == DynamicProteinAtomic
                or class_type == PigmentProteinAtomic
                or class_type == ElectrostaticProteinAtomic):
            return class_type(new_pdb, name, center_atomic)
        elif class_type == StaticProteinAtomic:
            return class_type(new_pdb, name)
        else:
            raise ValueError(f'Class type {class_type} not supported.')

    @property
    def R2_atomic_xyz(self) -> np.ndarray:
        """Get coordinates for all atoms in the protein.

        Returns:
            Array of coordinates for all atoms in Angstroms.
        """
        return np.array([atom.coord for atom in self._pdb.get_atoms()])

    @property
    def R2_residue_center(self) -> np.ndarray:
        """Get the center of masses of each residue in the protein.

        Returns:
            Array of center of masses of each residue in Angstroms.
        """
        return self._R2_residue_center

    @property
    def list_vdw_radii(self) -> np.ndarray:
        """Get a list of Van der Waals radii for each atom in the protein.

        Returns:
            Van der Waals radii for each atom in Angstroms.
        """
        if self.__list_vdw_radii is None:
            self.__list_vdw_radii = np.array([dict_vdw_radii[atom.element]
                                              for atom in self._pdb.get_atoms()])
        return np.array(self.__list_vdw_radii)

    @property
    def list_covalent_radii(self) -> np.ndarray:
        """Get a list of covalent radii for each atom in the protein.

        Returns:
            Covalent radii for each atom in Angstroms.
        """
        if self.__list_covalent_radii is None:
            radii = [dict_covalent_radii[atom.element]
                     for atom in self._pdb.get_atoms()]
            self.__list_covalent_radii = np.array(radii)
        return self.__list_covalent_radii

    @property
    def list_wca_sigmas(self) -> np.ndarray:
        """Get a list of the Weeks-Chandler-Anderson (WCA) potential radii for each atom.

        Returns:
            WCA radii for each atom in Angstroms.
        """
        if self.__list_wca_sigmas is None:
            self.__list_wca_sigmas = np.array([dict_wca_sigma[atom.element]
                                               if atom.element in dict_wca_sigma.keys()
                                               else 0
                                               for atom in self._pdb.get_atoms()
                                               ])
        return np.array(self.__list_wca_sigmas)

    @property
    def list_wca_epsilons(self) -> np.ndarray:
        """Get a list of the Weeks-Chandler-Anderson (WCA) potential epsilon for each atom.

        Returns:
            WCA epsilon for each atom in cm^-1.
        """
        if self.__list_wca_epsilons is None:
            self.__list_wca_epsilons = np.array([dict_wca_epsilon[atom.element]
                                                 if atom.element in dict_wca_epsilon.keys()
                                                 else 0
                                                 for atom in self._pdb.get_atoms()
                                                 ]
                                                )
        return np.array(self.__list_wca_epsilons)

    @property
    def hull3d(self) -> ConvexHull:
        """Get the 3D convex hull of the protein using SciPy.

        Returns:
            3D convex hull of the protein.
        """
        if self.__hull is None:
            self.__hull = ConvexHull(self.R2_residue_center)
        return self.__hull

    @property
    def hull2d(self) -> Polygon:
        """Get the 2D convex hull of the protein using shapely.

        Returns:
            2D convex hull of the protein as a shapely.Polygon.
        """
        if self.__hull2d is None:
            self.__hull2d = shapely_convex_hull(MultiPoint(self._R2_residue_center))
        return self.__hull2d

    @property
    def dict_hull2d_by_chain(self) -> dict:
        """Get a dictionary of the 2D convex hull of each chain in the protein.

        Returns:
            Dictionary mapping chain IDs to their 2D concave hulls (shapely.Polygon).
        """
        if self.__dict_hull2d_by_chain is None:
            self.__dict_hull2d_by_chain = {}
            for chain in self._pdb.get_chains():
                coords = np.array([atom.get_coord() for atom in chain.get_atoms()])
                self.__dict_hull2d_by_chain[chain.id] = shapely_concave_hull(MultiPoint(coords), ratio=0.20)
        return self.__dict_hull2d_by_chain

    @property
    def convex_area(self) -> float:
        """Get the area of the 2D convex hull of the protein.

        Returns:
            Area of the 2D convex hull in Angstroms^2.
        """
        return self.hull2d.area

    @property
    def hull2d_vertices(self) -> np.ndarray:
        """Get the vertices of the 2D convex hull of the protein.

        Returns:
            Vertices of the 2D convex hull.
        """
        return np.array(self.hull2d.exterior.coords)

    @property
    def radius(self) -> float:
        """Get the radius of the largest circle that inscribes the 2D convex hull.
        
        The circle is centered at the protein's center of mass.

        Returns:
            Radius of the circle in Angstroms.
        """

        if self.__hull2d_radius is None:
            diffs = self.hull2d_vertices - self.center_of_mass
            norms = np.linalg.norm(diffs, axis=1)
            self.__hull2d_radius = np.max(norms)

        return self.__hull2d_radius

    @property
    def concave_hull2d(self) -> Polygon:
        """Get the 2D concave hull of the protein.

        Returns:
            2D concave hull as a shapely.Polygon.
        """
        if self.__concave_hull is None:
            self.__concave_hull = shapely_concave_hull(MultiPoint(self.R2_residue_center), ratio=0.35)

        return self.__concave_hull

    @property
    def area(self) -> float:
        """Get the area of the 2D concave hull of the protein.

        Returns:
            Area of the 2D concave hull in Angstroms^2.
        """
        return self.concave_hull2d.area

    @property
    def concave_hull2d_vertices(self) -> np.ndarray:
        """Get the vertices of the 2D concave hull of the protein.

        Returns:
            Vertices of the 2D concave hull.
        """
        return np.array(self.concave_hull2d.exterior.coords)

    @property
    def inner_radius(self) -> float:
        """Get the radius of the largest circle that fits inside the concave hull.
        
        The circle is centered at the protein's center of mass.

        Returns:
            The inner radius in Angstroms.
        """
        if self.__inner_radius is None:
            diffs = np.linalg.norm(self.concave_hull2d_vertices - self.center_of_mass, axis=1)
            self.__inner_radius = np.min(diffs)
        return self.__inner_radius

    def visualize_protein_hull3d(self, ax: Axes, color: str = 'g') -> None:
        """Visualize the 3D convex hull of the protein using matplotlib.

        Args:
            ax: The matplotlib Axes to plot the 3D convex hull on.
            color: The color of the 3D convex hull (matplotlib color code). Defaults to 'g'.
        """
        list_vtx = []
        R2_residue_center = self.R2_residue_center
        for list_v in self.hull3d.simplices:
            vtx = np.array([R2_residue_center[list_v[0]],
                            R2_residue_center[list_v[1]],
                            R2_residue_center[list_v[2]]])
            list_vtx.append(vtx)
        tri = art3d.Line3DCollection(list_vtx)
        tri.set_color(color)
        tri.set_edgecolor(color)
        ax.add_collection3d(tri)
        ax.scatter([R2_residue_center[vtx][0] for vtx in self.hull3d.vertices],
                    [R2_residue_center[vtx][1] for vtx in self.hull3d.vertices],
                    [R2_residue_center[vtx][2] for vtx in self.hull3d.vertices],
                   color
                   )


    def visualize_protein_hull2d(
            self,
            ax: Axes =None,
            hull_type: Literal['convex', 'concave'] = 'convex',
            color: str = 'g',
            plot_label: str = None
    ) -> Axes:
        """Visualize the 2D hull of the protein.

        Args:
            ax: The matplotlib Axes to plot on. If None, a new figure and axes are created.
            hull_type: Type of hull to plot ('convex' or 'concave'). Defaults to 'convex'.
            color: The color of the hull (matplotlib color code). Defaults to 'g'.
            plot_label: Label for the plot legend.

        Returns:
            The matplotlib Axes object with the plot.

        Raises:
            ValueError: If hull_type is not 'convex' or 'concave'.
        """

        if ax is None:
            fig, ax = plt.subplots()
        ax.set_aspect('equal')

        if hull_type == 'convex':
            vertices = np.vstack([self.hull2d_vertices, self.hull2d_vertices[0]])
        elif hull_type == 'concave':
            vertices = np.vstack([self.concave_hull2d_vertices, self.concave_hull2d_vertices[0]])
        else:
            raise ValueError("hull_type must be 'convex' or 'concave'")

        ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=1)
        ax.fill(vertices[:, 0], vertices[:, 1], color=color, alpha=0.65, label=plot_label)

        return ax

    def visualize_protein_hull2d_by_chain(
            self,
            dict_chain_color: dict,
            ax: Axes =None,
            plot_label: str = None,
    ) -> Axes:
        """Visualize the 2D hull of each chain in the protein with different colors.

        Args:
            dict_chain_color: Dictionary mapping chain IDs to colors.
            ax: The matplotlib Axes to plot on. If None, a new figure and axes are created.
            plot_label: Label for the plot legend.

        Returns:
            The matplotlib Axes object with the plot.
        """

        if ax is None:
            fig, ax = plt.subplots()
        ax.set_aspect('equal')

        for (chain, hull) in self.dict_hull2d_by_chain.items():
            if chain not in dict_chain_color.keys():
                continue

            color = dict_chain_color[chain]

            hull2d_vertices = hull.exterior.coords
            vertices = np.vstack([hull2d_vertices, hull2d_vertices[0]])

            ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=1)
            ax.fill(vertices[:, 0], vertices[:, 1], color=color, alpha=0.65, label=plot_label)

        return ax


class DynamicProteinAtomic(StaticProteinAtomic):
    __slots__ = []

    def __init__(
            self,
            pdb: Union[str, PDB.Structure.Structure, FakePDB],
            name: str,
            prot_type: str = None,
            center_atomic: bool = True
    ) -> None:
        """Initialize a DynamicProteinAtomic object.

        This is an extension of StaticProteinAtomic that allows modifications to be made
        on the PDB object through wrapper methods.

        Roles:
            1. PDB Modifications: Instead of editing the PDB directly, you can use wrapper
               methods to modify the PDB object. This enforces the principle that the PDB
               object should be inaccessible under most use cases.

        Args:
            pdb: Either a PDB.Structure.Structure object instantiated by the Structure class
                of the Biopython package, a path (str) pointing to the PDB file, or a
                FakePDB object for testing.
            name: The name of the ProteinAtomic object.
            prot_type: Identifies the 'type' of protein for membrane simulations. If None,
                defaults to the value of name.
            center_atomic: Whether to set the origin at the center of mass.
        """
        super().__init__(pdb, name, prot_type)

        # Center Atomic Protein
        # ---------------------
        if center_atomic:
            self._transform_pdb(rot=np.eye(3), tran=-self.center_of_mass)

    def _transform_pdb(
            self,
            rot: np.ndarray = None,
            tran: np.ndarray = None
    ) -> None:
        """Perform rotation and translation on the PDB object and update properties.

        Args:
            rot: Right multiplying 3D rotation matrix in radians. If None, defaults to
                identity matrix.
            tran: Translation vector [x,y,z] in Angstroms. If None, defaults to zero vector.
        """
        if rot is None:
            rot = np.eye(3)
        if tran is None:
            tran = np.zeros(3)

        self._pdb.transform(rot=rot, tran=tran)
        self._center_of_mass = np.dot(self.center_of_mass, rot) + tran
        self._R2_residue_center = np.dot(self._R2_residue_center, rot) + tran

        self._reset_pdb()


    def translate(self, translate: np.ndarray):
        """Translate the location of the protein by the specified vector.

        Args:
            translate: Translation vector in Angstroms.
        """
        self._transform_pdb(tran=translate)

    def rotate_by_axis(self, angle, axis='z'):
        """Rotate the protein by a specified angle around a given axis.

        Args:
            angle: The angle to rotate the protein in radians.
            axis: The axis of rotation ('x', 'y', or 'z'). Defaults to 'z'.
        """
        center_of_mass = self.center_of_mass
        self.translate(-center_of_mass)
        self._transform_pdb(rot=construct_rotation_matrix(angle, axis))
        self.translate(center_of_mass)

    def rotate_by_matrix(self, matrix: np.ndarray):
        """Rotate the orientation of the PDB using a rotation matrix.

        Args:
            matrix: Left multiplying rotation matrix.
        """
        center_of_mass = self.center_of_mass
        self._transform_pdb(tran=-center_of_mass)
        self._transform_pdb(rot=matrix)
        self._transform_pdb(tran=center_of_mass)

    def _flatten(self, Z1_unit=None):
        """Orient the protein such that a special moment of inertia points in the z-direction.

        This constructs one 'special' moment of inertia that points in the direction of
        the membrane normal and orients the protein accordingly.

        Args:
            Z1_unit: Unit vector in the z-direction.
        """
        center_of_mass = self.center_of_mass
        self._pdb.transform(np.eye(3), -center_of_mass)
        R2_rot = construct_rotation_for_flatten(self.inertia_tensor, Z1_unit)
        self._transform_pdb(rot=np.transpose(R2_rot), tran=np.array([0, 0, 0]))
        self._transform_pdb(rot=np.eye(3), tran=center_of_mass)

    def _flip_up_down(self):
        """Flip the protein upside down by rotating 180 degrees around the z-axis."""
        com_save = self.center_of_mass
        self._transform_pdb(tran=-com_save)
        self._transform_pdb(rot=np.array([[1, 0, 0], [0, 1, 0], [0, 0, -1]]))
        self._transform_pdb(tran=com_save)


class PigmentProteinAtomic(DynamicProteinAtomic):

    def __init__(
            self,
            pdb: Union[str, PDB.Structure.Structure, FakePDB],
            name: str,
            prot_type: str = None,
            center_atomic: bool = True,
            dict_pigment_param: dict = None,
            dict_pigment_by_chain_param: dict = None
    ) -> None:
        """Initialize a PigmentProteinAtomic object.

        This is an extension of DynamicProteinAtomic that allows the presence of pigments
        within the protein, along with associated analyses regarding those pigments.

        Roles:
            1. Pigment Management: Handle pigments within the protein using associated
               methods and classes from atomic_pigment.py, including data analysis.
            2. Spectroscopy Calculations: Perform spectroscopy calculations using the
               Hamiltonian, which includes site energies and electronic coupling.

        Args:
            pdb: Either a PDB.Structure.Structure object instantiated by the Structure class
                of the Biopython package, a path (str) pointing to the PDB file, or a
                FakePDB object for testing.
            name: The name of the ProteinAtomic object.
            prot_type: Identifies the 'type' of protein for membrane simulations. If None,
                defaults to the value of name.
            center_atomic: Whether to set the origin at the center of mass.
            dict_pigment_param: Dictionary mapping residue names to pigment parameters:
                {resname: (ClassPigment, lineshape, dict_coupling_data, dict_q0011_charges,
                disorder, dipole_mag)}. Adds data to every pigment matching the resname.
            dict_pigment_by_chain_param: Dictionary mapping (chain, resname) tuples to
                pigment parameters: {(chain, resname): (ClassPigment, lineshape,
                dict_coupling_data, dict_q0011_charges, disorder, dipole_mag)}.
                Adds data to pigments of specific chains in the PDB.
        """
        super().__init__(pdb, name, prot_type, center_atomic)

        # Prepare Pigments by Dictionary
        # ------------------------------
        # First we perform all the of the pigment definitions by chain, then
        # we 'fill in' anything missing by finding any residue not yet in
        # self.dict_pigments with a resname in dict_pigment and constructing
        # that object.
        if dict_pigment_by_chain_param is not None:
            for ((chain, resname), param) in dict_pigment_by_chain_param.items():
                self.prepare_pigments(resname, *param)

        if dict_pigment_param is not None:
            for (resname, param) in dict_pigment_param.items():
                self.prepare_pigments(resname, None, *param)

        self._list_pigments = list(self.dict_pigments.values())
        self._central_frequency = None

    @property
    def list_pigments(self):
        return self._list_pigments

    @list_pigments.setter
    def list_pigments(self, new_list_pigments):
        self._list_pigments = new_list_pigments
        self.H2_hamiltonian = None

    @property
    def central_frequency(self):
        if self._central_frequency is not None:
            return self._central_frequency
        else:
            return 0

    @property
    def fval_tresp(self):
        return self._fval_tresp

    @fval_tresp.setter
    def fval_tresp(self, new_fval_tresp):
        self._fval_tresp = new_fval_tresp

    @property
    def dielectric_tresp(self):
        return self._dielectric_tresp

    @dielectric_tresp.setter
    def dielectric_tresp(self, new_dielectric_tresp):
        self._dielectric_tresp = new_dielectric_tresp

    @property
    def dielectric_dipole(self):
        return self._dielectric_dipole

    @dielectric_dipole.setter
    def dielectric_dipole(self, new_dielectric_dipole):
        self._dielectric_dipole = new_dielectric_dipole


    def prepare_pigments(
            self,
            resname: str,
            ClassPigments,
            list_chain: list[str] = None,
            lineshape=None,
            dict_coupling_data=None,
            dict_q0011_charges=None,
            disorder=None,
            **kwargs
    ) -> None:
        """Find and create pigment objects for residues with the given resname.

        Args:
            resname: Three letter residue name from the PDB.
            ClassPigments: Class from atomic_pigment.py that will prepare the pigment object.
            list_chain: List of chain IDs to search in. If None, searches all chains.
            lineshape: Lineshape function for the pigment.
            dict_coupling_data: Dictionary with atom names as keys and q_0011 charges as values.
            dict_q0011_charges: Dictionary of q_0011 charge data for the pigment.
            disorder: Disorder parameter (sigma_e) for the pigment.
            **kwargs: Additional keyword arguments passed to the ClassPigments constructor.
        """
        for res in self.get_residues():
            if (list_chain is not None and res.get_parent().id in list_chain) or (list_chain is None):
                pigment_name = pigment_name_from_residue(res)
                if res.get_resname() == resname and pigment_name not in self.dict_pigments.keys():
                    # Add Pigment to Dictionary
                    self.dict_pigments[pigment_name] = ClassPigments(res, lineshape, **kwargs)

                    # Update Coupling Data
                    if dict_coupling_data is not None:
                        self.dict_pigments[pigment_name].get_coupling_data(dict_coupling_data)

                    # Update Charge Data
                    if dict_q0011_charges is not None:
                        self.dict_pigments[pigment_name].get_q_0011_charge(dict_q0011_charges)

                    if disorder is not None:
                        self.dict_pigments[pigment_name].dict_data['sigma_e'] = disorder

        # If you are adding pigments, then you need a new Hamiltonian
        self.list_pigments = list(self.dict_pigments.values())

    def load_hamiltonian(self, dict_or_path, map_chains_from_saved=None):
        """Load the Hamiltonian representation for a system.

        Loads from a file or a dictionary of file paths, and sets up structure and couplings
        based on the provided data. Handles various configurations of Hamiltonian data,
        allowing flexibility in how structure and coupling information is defined and
        mapped to pigments.

        Args:
            dict_or_path: Either a dictionary mapping chain identifiers to file paths for
                Hamiltonian data, or a single file path to a Hamiltonian dataset.
            map_chains_from_saved: Optional mapping function to adjust chain indices during
                Hamiltonian loading.
        """
        # Remove central frequency
        self.set_central_frequency(0)

        if type(dict_or_path) is not dict:
            list_name_saved, H2_ham_0 = load_hamiltonian(dict_or_path, map_chains_from_saved)
            list_name_site = [pig.name for pig in self.list_pigments]
            list_name = match_pigment_names(list_name_saved, list_name_site)
            if list_name is not None:
                list_sigma = [self.dict_pigments[name].dict_data['sigma_e']
                              if 'sigma_e' in self.dict_pigments[name].dict_data.keys() else 0
                              for name in list_name]
                self.H2_hamiltonian = Hamiltonian(H2_ham_0, list_name, list_sigma)
                self.set_central_frequency()
            else:
                print(f'No Hamiltonian added to chain {self.name}: pigment names and structure names do not match.')

        else:
            list_site_names = [pig.name for pig in self.list_pigments]
            list_sigma = [self.dict_pigments[name].dict_data['sigma_e']
                          if 'sigma_e' in self.dict_pigments[name].dict_data.keys() else 0
                          for name in list_site_names]

            # Initialize Zero Array
            # ---------------------
            self.H2_hamiltonian = Hamiltonian(np.zeros([self.n_pigments, self.n_pigments]),
                                              list_site_names, list_sigma)

            # Add diagonal sub-blocks
            # -----------------------
            # The tricky part here is figuring out how to generalize this matching process
            # so that we only need to store one copy of the LHCII hamiltonian even though it
            # has many different chain indices in different structures.
            dict_pigname_by_chain = {}
            for (tuple_chain, path_to_ham) in dict_or_path.items():
                list_name_site = [pig.name for pig in self.list_pigments if pig.residue.get_parent().id in tuple_chain]
                dict_pigname_by_chain[tuple_chain] = list_name_site
                list_name_saved, h2_ham = load_hamiltonian(path_to_ham, map_chains_from_saved)
                # The names in the hamiltonian csv should not depend on the atomic protein name we assign when preparing
                # the atomic protein, However, when we load the hamiltonian we are expecting the row/column name should
                # exactly match the pigment name, that is why I am using 'update_pigment_name' function, that preserve
                # the chain_name,pigment_type and pigment_index, and just update the name prefix with the atomic protein name.
                list_name_saved_updated = update_pigment_name(list_name_saved,self.name)
                # print(f'{list_name_saved_updated=}')
                if list_name_saved_updated is not None:
                    self.H2_hamiltonian.update_subset(h2_ham,
                                                      list_name_saved_updated,
                                                      list_name_saved_updated)

            # Add off-diagonal blocks: TrEsp or Dipole-Dipole Coupling
            # --------------------------------------------------------
            # In this section, we will default to using TrESP unless the
            # the data is not available, in which case we will use
            # dipole-dipole coupling.
            list_chains = list(dict_or_path.keys())
            for (index_a, chain_a) in enumerate(list_chains):
                for (index_b, chain_b) in enumerate(list_chains[index_a:]):
                    index_b = index_b + index_a
                    if index_b > index_a:
                        list_pig_a = [self.dict_pigments[name] for name in dict_pigname_by_chain[chain_a]]
                        list_pig_b = [self.dict_pigments[name] for name in dict_pigname_by_chain[chain_b]]
                        for pig_a in list_pig_a:
                            for pig_b in list_pig_b:
                                if ((pig_a.dict_data['tresp_pc'] is not None)
                                        and (pig_b.dict_data['tresp_pc'] is not None)):
                                    v_ab = calculate_tresp_coupling(pig_a, pig_b, self.fval_tresp,
                                                                    self.dielectric_tresp)
                                else:
                                    v_ab = calculate_dipole_coupling(pig_a, pig_b, self.dielectric_dipole)
                                self.H2_hamiltonian.update_coupling(v_ab, pig_a.name, pig_b.name)
            self.set_central_frequency()

    def set_central_frequency(self, central_frequency=None):
        if (central_frequency is None) and (self.H2_hamiltonian is not None):
            central_frequency = np.max(np.diag(self.H2_hamiltonian.H2_ham_0))

        self._central_frequency = central_frequency
        if self.H2_hamiltonian is not None:
            self.H2_hamiltonian.set_central_frequency(central_frequency)

    def build_hamiltonian_domains(self, calc_method=None, hamiltonian_filter=15, N_ens=1000,
                                  domain_cutoff=0.1, coupling_cutoff=20, precomputed_domains=None):
        """Build domains for the loaded Hamiltonian in this AtomicProtein instance.

        Args:
            calc_method: The calculation method ('Energy' or 'Coupling'). Ignored if
                precomputed_domains is provided.
            hamiltonian_filter: Filter for Hamiltonian matrix (Energy method).
            N_ens: Number of disorder realizations (Energy method).
            domain_cutoff: Cutoff for domain formation (Energy method).
            coupling_cutoff: Cutoff for electronic coupling in cm^-1 (Coupling method).
            precomputed_domains: Preloaded domain list for direct use. If provided,
                calc_method is ignored.

        Returns:
            The list of domains (list_pigments_by_domain).

        Raises:
            ValueError: If Hamiltonian has not been loaded, if calc_method is not specified
                when precomputed_domains is not provided, or if an invalid calc_method is given.
        """
        # Use precomputed_domains directly if provided
        if precomputed_domains is not None:
            print("Loading precomputed domain data...")
            self.list_pigments_by_domain = precomputed_domains
            print(f'precomputed_domains= {precomputed_domains}')
            return self.list_pigments_by_domain

        # Ensure Hamiltonian is loaded if precomputed_domains is not provided
        if self.H2_hamiltonian is None:
            raise ValueError("Hamiltonian has not been loaded.")

        # Check that calc_method is provided if precomputed_domains is not
        if calc_method is None:
            raise ValueError("`calc_method` must be specified unless `precomputed_domains` is provided.")

        # Build domains based on the calculation method
        if calc_method.lower() == "energy":
            print('Building domains based on Energy ...')
            self.list_pigments_by_domain = self.H2_hamiltonian.build_energetic_domains(
                domain_cutoff=domain_cutoff,
                hamiltonian_filter=hamiltonian_filter,
                N_ens=N_ens
            )
            return self.list_pigments_by_domain

        elif calc_method.lower() == "coupling":
            print('Building domains based on coupling ...')
            self.list_pigments_by_domain = self.H2_hamiltonian.build_coupling_domains(
                coupling_cutoff=coupling_cutoff
            )
            return self.list_pigments_by_domain

        else:
            raise ValueError("Invalid calculation method. Choose 'Energy' or 'Coupling'.")

    def __build_dict_domain_by_name(self, list_pigments_by_domain):
        """Build a dictionary mapping pigment names to their domain indices.

        Args:
            list_pigments_by_domain: A list of lists where each inner list contains
                pigment names for a specific domain.

        Returns:
            Dictionary mapping pigment names to their respective domain index.
        """
        return {name: index for index, domain in enumerate(list_pigments_by_domain) for name in domain}

    def define_domains(self, list_pigments_by_domain):
        """Assign domain indices to pigments based on the provided list.

        Args:
            list_pigments_by_domain: A list of pigments grouped by domain.
        """
        dict_domain_by_name = self.__build_dict_domain_by_name(list_pigments_by_domain)
        for name, domain_index in dict_domain_by_name.items():
            self.dict_pigments[name].domain_index = domain_index

    def visualize_pigment_location(self, ax, color='k'):
        """Visualize pigment locations in 3D.

        Args:
            ax: Matplotlib 3D Axes object for visualization.
            color: Color for the pigment markers. Defaults to 'k' (black).
        """
        R3_pigment_positions = np.zeros([self.n_pigments, 3], dtype=np.float32)
        for (index, pigment) in enumerate(self.dict_pigments.values()):
            R3_pigment_positions[index, :] = pigment.location

        ax.scatter(R3_pigment_positions[:, 0],
                   R3_pigment_positions[:, 1],
                   R3_pigment_positions[:, 2],
                   color=color, marker='o')

    @property
    def n_pigments(self):
        """Count the number of pigments in dict_pigments.

        Returns:
            Number of pigments.
        """
        return len(self.dict_pigments)

    def delete_pigments(self, list_remove_name):
        """Delete specified pigments from the protein.

        Args:
            list_remove_name: List of pigment names to remove.
        """
        for name in list_remove_name:
            if name in self.dict_pigments.keys():
                list_names = [pig.name for pig in self.list_pigments]
                index_name = list_names.index(name)
                self.dict_pigments.pop(name)
                self.list_pigments.pop(index_name)
                if self.H2_hamiltonian is not None:
                    self.H2_hamiltonian.delete_pigments(list_remove_name)
            else:
                # Do we want to raise an error?
                print(f'ERROR: {name} is not in AtomicProteion.dict_pigments.keys() and could not be removed.')

class ElectrostaticProteinAtomic(PigmentProteinAtomic):
    """Extension of PigmentProteinAtomic for managing electrostatic charges on atoms.

    This class manages the presence of electrostatic charges on atoms associated with
    protein residues and provides CDC (Charge Density Coupling) capabilities.
    """

    def __init__(self, _pdb, name, prot_type=None, center_atomic=True,
                 dict_pigment_param=None, dict_pigment_by_chain_param=None,
                 path_extended_pdb=None):
        """Initialize an ElectrostaticProteinAtomic object.

        Args:
            _pdb: Either a PDB.Structure.Structure object, a path pointing to the PDB file,
                or a FakePDB object for testing.
            name: The name of the ProteinAtomic object.
            prot_type: Identifies the 'type' of protein for membrane simulations. If None,
                defaults to the value of name.
            center_atomic: Whether to set the origin at the center of mass.
            dict_pigment_param: Dictionary mapping residue names to pigment parameters.
            dict_pigment_by_chain_param: Dictionary mapping (chain, resname) tuples to
                pigment parameters.
            path_extended_pdb: Path to the extended PDB file containing atomic partial charges
                in columns 82-90.
        """
        super().__init__(_pdb, name, prot_type, center_atomic,
                         dict_pigment_param, dict_pigment_by_chain_param)

        if path_extended_pdb is not None:
            self.set_atomic_charges(path_extended_pdb)

        self.dielectric_cdc = None

    def set_atomic_charges(self, path_extended_pdb):
        """Add atomic charges to each atom in the protein from an extended PDB file.

        Uses the pqr_charge property of the Bio.PDB.Atom module.

        Args:
            path_extended_pdb: Path to the extended PDB file that contains atomic partial
                charges in columns 82-90.

        Raises:
            ValueError: If path_extended_pdb is not a string path to a PDB file.
        """
        # freeze columns to 82-90
        # gromacs removes chain names
        # set flag at beginning of function to test if user actually passed in a regular PDB with no charges
        # append charges to a list and check if all the values of a list are none to warn user
        if type(path_extended_pdb) != type('str'):
            raise ValueError(
                'This function requires an extended PDB file with charges in columns 82-90, you cannot '
                'build the atomic protein object using a PDB object, it must be a PDB path')
        else:
            list_atomic_ids = [
                f'{atom.get_full_id()[2]}_{atom.get_parent().resname}_{atom.get_parent().id[1]}_' \
                f'{atom.get_id()}' for atom in self._pdb.get_atoms()]
            list_atomic_atoms = [atom for atom in self._pdb.get_atoms()]
            dict_atomic_ids = dict(zip(list_atomic_ids, list_atomic_atoms))
            for pdb_line in open(path_extended_pdb).readlines():
                if pdb_line[0:6].strip() in ['ATOM', 'HETATM']:
                    atom_id = f'{pdb_line[21].strip()}_{pdb_line[17:20].strip()}_' \
                              f'{int(pdb_line[22:26].strip())}' \
                              f'_{pdb_line[12:17].strip()}'
                    if pdb_line[80:90].strip() == '':
                        charge = None
                    elif pdb_line[80:90].strip() == 'None':
                        charge = None
                    else:
                        charge = float(pdb_line[80:90].strip())
                    dict_atomic_ids[atom_id].set_charge(charge)

    def calculate_cdc_site_shift(self,
                                 verbose=False):
        """Calculate the site energy shift for each pigment due to CDC.

        Optionally shows the contribution of individual residues on the total shift for
        each site. All data is stored in the pigment dictionaries and also returned.

        Args:
            verbose: Whether to print detailed information about contributions. Defaults
                to False.

        Returns:
            Tuple containing:
                - dict_total_site_energy: Dictionary mapping pigment names to total shift
                  in site energy (cm^-1).
                - dict_total_contribution_site_energy: Dictionary mapping pigment names to
                  dictionaries of residue contributions.
        """
        dict_total_site_energy = {}
        dict_total_contribution_site_energy = {}
        if verbose: print('---------------\n Calculating CDC \n---------------')

        for pigA in self.dict_pigments.values():
            if verbose: print(f'{pigA.name}')

            dict_background_contribution = {}
            delta_E_sum = 0
            list_background_objects = [pigB for pigB in
                                       self.dict_pigments.values()
                                       if not pigA == pigB]

            list_background_objects.extend([res for res in
                                            self.get_residues()
                                            if pigment_name_from_residue(res)
                                            not in self.dict_pigments.keys()])

            # Calculate contributions individually
            for background_object in list_background_objects:
                name = name_background_object(background_object)
                delta_E = pigment_site_energy_shift_by_component(pigA,
                                                                 background_object,
                                                                 self.dielectric_cdc)
                delta_E_sum += delta_E
                dict_background_contribution[name] = delta_E
                if verbose: print(f"Interaction between {pigA.name} and {name}: {delta_E}")

            # Store data in pigment
            pigA.dict_data['deltaE_shift'] = delta_E_sum
            pigA.dict_data['deltaE_shift_components'] = dict_background_contribution

            # Construct Contribution Dictionary by pigment
            dict_total_contribution_site_energy[pigA.name.strip()] = \
                dict_background_contribution
            dict_total_site_energy[pigA.name.strip()] = delta_E_sum

        return dict_total_site_energy, dict_total_contribution_site_energy

    def calculate_total_site_energy(self, dict_site_energy_shift, E_0a, E_0b):
        """Calculate the total site energy for each pigment.

        Args:
            dict_site_energy_shift: Dictionary containing the site energy shifts for each
                pigment in cm^-1.
            E_0a: The base energy for pigment type 'CLA' in cm^-1.
            E_0b: The base energy for pigment type 'CHL' in cm^-1.

        Returns:
            Dictionary mapping each pigment to its corresponding total site energy in cm^-1.

        Raises:
            ValueError: If the pigment type is not supported (not 'CLA' or 'CHL').
        """
        E_total = []
        for key, value in dict_site_energy_shift.items():
            if 'CLA' in key:
                E_total.append(value + E_0a)
            elif 'CHL' in key:
                E_total.append(value + E_0b)
            else:
                raise ValueError(f'Pigment type ({key}) not supported.')

        return dict(zip(list(dict_site_energy_shift.keys()), E_total))



    def construct_hamiltonian(self , e0_a , e0_b , site_energy_calc='cdc' , coupling_calc='tresp' , dir_save=None) :
        """Construct a Hamiltonian matrix for a given set of pigments.

        Args:
            e0_a: Base energy level for site A in cm^-1.
            e0_b: Base energy level for site B in cm^-1.
            site_energy_calc: Method for site energy calculation. Currently, only 'cdc'
                (charge density coupling) is supported.
            coupling_calc: Method for coupling calculation. Options are 'dipole' or 'tresp'
                (transition electrostatic potential).
            dir_save: Directory path to save the Hamiltonian matrix as a CSV file. If None,
                the matrix is not saved.

        Returns:
            Hamiltonian matrix as a nested dictionary with site energies on the diagonal
            and coupling values off-diagonal.

        Raises:
            ValueError: If an unsupported site energy calculation or coupling calculation
                method is provided.
        """
        if not (coupling_calc == 'dipole' or coupling_calc == 'tresp') :
            raise ValueError(f"Coupling calculation {coupling_calc} not recognized (expected: 'dipole' or 'tresp').")

        if site_energy_calc == 'cdc' :
            dict_shifts = self.calculate_cdc_site_shift()[0]
            dict_site_energies = self.calculate_total_site_energy(dict_shifts , e0_a , e0_b)
        else :
            raise ValueError(f'Site energy calculation {site_energy_calc} currently not supported')

        # Initialize Hamiltonian as a nested dictionary
        H2_hamiltonian = { pigA.name : { pigB.name : 0.0 for pigB in self.dict_pigments.values() } for pigA in
                           self.dict_pigments.values() }

        # Populate the Hamiltonian matrix
        for pigA in self.dict_pigments.values() :
            for pigB in self.dict_pigments.values() :
                if pigA.name == pigB.name :
                    H2_hamiltonian[pigA.name][pigB.name] = dict_site_energies[pigA.name]
                else :
                    if coupling_calc == 'tresp' :
                        H2_hamiltonian[pigA.name][pigB.name] = calculate_tresp_coupling(pigA ,  pigB ,
                                                                                                self.fval_tresp ,
                                                                                                self.dielectric_tresp)
                    elif coupling_calc == 'dipole' :
                        H2_hamiltonian[pigA.name][pigB.name] = calculate_dipole_coupling(pigA , pigB ,
                                                                                                self.dielectric_dipole)
        return H2_hamiltonian

    @staticmethod
    def _save_hamiltonian_to_csv(H2_hamiltonian , dir_save) :
        """Save the Hamiltonian matrix as a CSV file.

        Args:
            H2_hamiltonian: The Hamiltonian matrix as a nested dictionary.
            dir_save: Directory path to save the Hamiltonian matrix as a CSV file.

        Raises:
            ValueError: If dir_save is not provided.
        """
        if not dir_save :
            raise ValueError("The 'dir_save' parameter must be provided to save the file.")

        # Ensure the directory exists, if not, create it
        if not os.path.exists(os.path.dirname(dir_save)) :
            os.makedirs(os.path.dirname(dir_save))

        # Write the Hamiltonian matrix to the CSV file
        with open(dir_save , 'w' , newline='') as csvfile :
            writer = csv.writer(csvfile)
            # Write header (pigment names)
            header = [''] + list(H2_hamiltonian.keys())
            writer.writerow(header)
            # Write rows (pigment name and corresponding values)
            for pigA , row in H2_hamiltonian.items() :
                writer.writerow([pigA] + [row[pigB] for pigB in H2_hamiltonian.keys()])

    def __classify_contribution_by_type(self, dict_data: Dict[str, float], list_pigments: List[str], list_cofactors: List[str]) -> Dict[str, str]:
        """Classify each contribution as 'Pigments', 'Cofactors', 'Water', or 'Residues'.

        Args:
            dict_data: Dictionary of residue contributions.
            list_pigments: List of substrings to identify pigments.
            list_cofactors: List of substrings to identify cofactors.

        Returns:
            Dictionary where keys are the same as in dict_data and values are the
            classified types.
        """
        return {
            key: ('Pigments' if any(pigment in key for pigment in list_pigments)
                    else 'Cofactors' if any(cofactor in key for cofactor in list_cofactors)
                    else 'Water' if 'HOH' in key
                    else 'Residues')
                    for key in dict_data
        }

    def __calculate_sum_of_contributions(self, contributions: Dict[str, float], classified_contribution: Dict[str, str]) -> Dict[str, float]:
        """Sum the contributions for each category and calculate the total.

        Args:
            contributions: Dictionary of residue contributions.
            classified_contribution: Dictionary classifying the contributions by type.

        Returns:
            Dictionary with summed contributions for each category ('Residues', 'Pigments',
            'Cofactors', 'Water') and the total sum.
        """
        sums = {category: sum(value for res, value in contributions.items()
                              if classified_contribution[res] == category)
                for category in ['Residues', 'Pigments', 'Cofactors', 'Water']}  #['protein', 'cofactor', 'water', 'pigment']
        sums['total'] = sum(sums.values())  # Add total sum
        return sums

    def __plot_aggregated_contributions(self, pigment: str, sum_contributions: Dict[str, float], save_dir: str) -> None:
        """Save aggregated contributions for each category as a plot.

        Args:
            pigment: The pigment being plotted.
            sum_contributions: Dictionary of summed contributions.
            save_dir: Directory to save the plot.
        """
        # Ensure the save directory exists
        os.makedirs(save_dir, exist_ok=True)
        
        # Create the bar plot
        categories, values = zip(*sum_contributions.items())  # Extract keys and values
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Improved publication-quality colors
        plt.figure(figsize=(10, 6))
        plt.bar(categories, values, color=colors)
        plt.ylabel('wave number (cm-1)')
        plt.title(f'Site Energy Decomposition of Pigment: {pigment}')
        plt.tight_layout()

        # Save the plot as an image file
        plot_path = os.path.join(save_dir, f"{pigment}_contributions.png")
        plt.savefig(plot_path)
        plt.close()  # Close the plot to avoid displaying it

    def plot_agregated_energy_contribution(self, list_pigments: Optional[List[str]] = None, list_cofactors: Optional[List[str]] = None, save_dir: str = "plots") -> None:
        """Classify contributions for each pigment and save aggregated contribution plots.

        If no lists for pigments or cofactors are provided, default lists are used.

        Args:
            list_pigments: List of substrings to identify pigments. Defaults to
                ['CLA', 'CHL', 'PHO'].
            list_cofactors: List of substrings to identify cofactors. Defaults to
                ['LHG', 'LUT', 'NEX', 'XAT', 'LMG'].
            save_dir: Directory to save the plots. Defaults to "plots".
        """
        list_pigments = list_pigments or ['CLA', 'CHL', 'PHO']
        list_cofactors = list_cofactors or ['LHG', 'LUT', 'NEX', 'XAT', 'LMG']

        for pigment, data in self.dict_pigments.items():
            residue_contribution = data.dict_data['deltaE_shift_components']
            classified_contribution = self.__classify_contribution_by_type(residue_contribution, list_pigments, list_cofactors)

            # Store classified contributions back into the object if needed
            self.dict_pigments[pigment].dict_data['classified_contribution_by_type'] = classified_contribution

            # Calculate sum and save plot contributions
            sum_contributions = self.__calculate_sum_of_contributions(residue_contribution, classified_contribution)
            self.__plot_aggregated_contributions(pigment, sum_contributions, save_dir)