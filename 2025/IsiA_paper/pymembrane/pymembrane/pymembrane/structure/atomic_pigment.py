import copy
import itertools
import math
from abc import ABC, abstractmethod
from itertools import product

import numpy as np

from pymembrane.util.physical_constants import dict_covalent_radii
from pymembrane.util.procrustes import procrustes

default_colors = {'C': 'k',
                  'O': 'r',
                  'N': 'b'}


def nth(iterable, n, default=None):
    """Returns the nth item or a default value.
    
    Args:
        iterable: An iterable to extract from.
        n: The index of the item to extract.
        default: The default value to return if the iterable is exhausted.
        
    Returns:
        The nth item from the iterable, or the default value if exhausted.
    """
    return next(itertools.islice(iterable, n, None), default)


def dist_molecular(atoms_small, atoms_large):
    """Calculates the total distance between two sets of points.
    
    This function accepts two sets of points (without labels) and defines the total
    distance between these vectors by finding the minimum distance from each point
    in the smaller set to any point in the larger set.
    
    Args:
        atoms_small: Array of coordinates for the smaller set of atoms.
        atoms_large: Array of coordinates for the larger set of atoms.
        
    Returns:
        float: The total distance between the two sets of atoms.
    """
    dist = 0
    for p1 in atoms_small:
        dist += np.min(np.linalg.norm(atoms_large - p1, axis=1))
    return dist

class ResidueAtomic(ABC):
    """Abstract base class for atomic residue objects.
    
    This class provides a template to construct atomic residue classes, which contain
    the atomic level information of residues such as position, mass, and bonds.
    
    Attributes:
        residue: A residue object from the Biopython Residue class.
        list_atom_names: List of atom names in the residue.
        bonds: Set of bonds between atoms within the residue.
        dict_data: Dictionary to store additional data.
    """
    def __init__(self, residue):
        """Initializes the ResidueAtomic object.

        Args:
            residue: A residue object from the Biopython Residue class.
        """
        self.residue = residue
        self.list_atom_names = [atom.name for atom in residue.get_atoms()]
        self.bonds = self.define_bonds()
        self.dict_data = {}

    @property
    def n_atoms(self):
        """Returns the number of atoms in the residue.
        
        Returns:
            int: The number of atoms.
        """
        return len([atom for atom in self.residue.get_atoms()])

    @property
    @abstractmethod
    def location(self):
        """Returns the location of the residue.
        
        This is an abstract method that must be implemented by subclasses.
        
        Returns:
            array: The X,Y,Z coordinates of the residue location.
        """
        pass

    @property
    def mass(self):
        """Calculates the mass of the residue.
        
        The mass is calculated by summing the mass of each atom in the residue.

        Returns:
            float: The total mass of the residue.
        """
        mass = 0
        for atom in self.residue.get_atoms():
            mass += atom.mass

        return mass

    def define_bonds(self, scale_factor=1.05):
        """Defines the bonds between atoms within the residue.
        
        Bonds are defined based on the covalent radii of the atoms. Two atoms are
        considered bonded if their distance is less than the sum of their covalent
        radii multiplied by the scale factor.

        Args:
            scale_factor: A slight scaling factor to ensure that all bonds are
                accounted for. Default is 1.05.

        Returns:
            set: A set containing tuples that define the bonds by the names of
                both atoms bonded.
        """
        bonds = set([])
        list_atom_pairs = [(a1, a2) for a1, a2 in product(self.residue.get_atoms(), self.residue.get_atoms())
                           if a1 != a2]
        for atom1, atom2 in list_atom_pairs:
            r1 = dict_covalent_radii[atom1.element]
            r2 = dict_covalent_radii[atom2.element]
            distance = np.linalg.norm(atom1.coord - atom2.coord)
            if distance <= scale_factor * (r1 + r2):
                name = tuple(sorted([atom1.id, atom2.id]))
                if not name in bonds:
                        bonds.add(name)
        return bonds

    def visualize(self, ax, set_color=None):
        """Constructs a visualization of the residue.

        Args:
            ax: The matplotlib axes object for the visualization.
            set_color: Optional color to override default atom colors.
        """
        self._visualize_atoms(ax, set_color)
        self._visualize_bonds(ax)

    def _visualize_atoms(self, ax, set_color):
        """Constructs a visualization of the atoms in a residue.

        Args:
            ax: The matplotlib axes object for the visualization.
            set_color: Optional color to override default atom colors.
        """
        for atom in self.residue.get_atoms():
            if atom.element in default_colors.keys():
                color = default_colors[atom.element]
                if set_color:
                    color = set_color
            else:
                color = 'purple'

            ax.scatter([atom.coord[0]], [atom.coord[1]], [atom.coord[2]],
                       marker='o',
                       s=30,
                       color=color)

    def _visualize_bonds(self, ax):
        """Constructs a visualization of the bonds within a residue.

        Args:
            ax: The matplotlib axes object for the visualization.
        """
        for (atom1_id, atom2_id) in self.bonds:
            try:
                atom1 = self.residue[atom1_id]
                atom2 = self.residue[atom2_id]
                ax.plot([atom1.coord[0], atom2.coord[0]],
                        [atom1.coord[1], atom2.coord[1]],
                        [atom1.coord[2], atom2.coord[2]],
                        linewidth=3, color='k', alpha=0.5)
            except:
                pass

    def find_overlap_transform(self, residue, element):
        """Finds the optimal transformation to overlap two residues.
        
        This function accepts two residue entries and an element symbol to be considered.
        It extracts all the corresponding elements in both residues and attempts all possible
        orderings to match them. For each possible ordering it calculates the procrustes
        mapping. The final selection is based on the minimum residue mapping across all
        atoms (ignoring labels).

        Args:
            residue: The target residue to align to.
            element: The element symbol to use for alignment (e.g., 'N', 'C').

        Returns:
            tuple: A tuple containing:
                - z_norm (float): The residual sum of squared errors, normalized according
                    to a measure of the scale of X, ((X - X.mean(0))**2).sum().
                - Z (array): The matrix of transformed residue coordinates.
                - tform (dict): A dict specifying the rotation, translation and scaling
                    that maps self -> residue.
        """
        # Determine which Chlorophyll has more atoms
        # ------------------------------------------
        ncla_to = residue.n_atoms
        ncla_from = self.n_atoms
        if ncla_to <= ncla_from:
            from_large = True
        else:
            from_large = False

        # Select Atoms by Element Type
        # ----------------------------
        atoms_to = [atom.coord for atom in residue.residue.get_atoms() if (atom.element == element)]
        atoms_to = [[x, y, z] for (x, y, z) in atoms_to]
        # print(atoms1)

        atoms_from = [atom.coord for atom in self.residue.get_atoms() if (atom.element == element)]
        atoms_from = [[x, y, z] for (x, y, z) in atoms_from]
        # print(atoms2)

        atoms_to_full = [atom.coord for atom in residue.residue.get_atoms()]
        atoms_from_full = [atom.coord for atom in self.residue.get_atoms()]
        if not (len(atoms_from) == len(atoms_to)):
            print("Error: Two pigments have a different number of " + element + " atoms.")
            pass
        else:
            print("Number of permutational mappings: " + str(math.factorial(len(atoms_from))))

        # Calculate Mapping for All Permutations of "atoms2"
        # --------------------------------------------------
        atom_perms = itertools.permutations(atoms_from)
        z_perms = []
        z_Nperms = []
        count = 0
        for atoms_test in atom_perms:
            z_Nnorm, _, z_trans = procrustes(atoms_to, atoms_test, scaling=False, reflection=False)
            atoms_from_trans = np.dot(np.array(atoms_from_full), np.array(z_trans['rotation'])) + np.array(
                z_trans['translation'])
            if from_large:
                z_norm = dist_molecular(np.array(atoms_to_full), np.array(atoms_from_trans))
            else:
                z_norm = dist_molecular(np.array(atoms_from_trans), np.array(atoms_to_full))
            z_perms.append(z_norm)
            z_Nperms.append(z_Nnorm)
            count += 1

        atom_perms = itertools.permutations(atoms_from)
        i_min = z_perms.index(min(z_perms))
        perm_min = nth(atom_perms, i_min)
        # print("Index of Best Permutation: " + str(i_min))
        # print(perm_min)
        z_norm, Z, tform = procrustes(atoms_to, perm_min, scaling=False, reflection=False)

        return z_norm, Z, tform




class PigmentAtomic(ResidueAtomic):
    """Atomic pigment class.
    
    This class inherits from the ResidueAtomic class to create an atomic pigment object.
    It contains important information specific to pigments such as the dipole moment.
    
    Attributes:
        residue: A residue object from the Biopython Residue class.
    """
    def __init__(self, residue):
        """Initializes the PigmentAtomic object.

        Args:
            residue: A residue object from the Biopython Residue class.
        """
        super().__init__(residue)

    def __str__(self):
        return f'{type(self).__name__} {self.name}'


    @property
    @abstractmethod
    def name(self):
        """Returns the name of the residue.
        
        This is an abstract method that must be implemented by subclasses.

        Returns:
            str: The name of the residue.
        """
        pass

    @abstractmethod
    def get_dipole_dir(self):
        """Calculates the direction of dipole vector of the pigment.
        
        This is an abstract method that must be implemented by subclasses.

        Returns:
            array: A unit vector describing the direction of the dipole moment.
        """
        pass

    @abstractmethod
    def visualize_dipole(self, ax):
        """Visualizes the dipole vector of the pigment.
        
        This is an abstract method that must be implemented by subclasses.

        Args:
            ax: The matplotlib axes object for the visualization.
        """
        pass

    def get_coupling_data(self, dict_coupling):
        """Adds coupling data necessary for TrEsp coupling calculations.
        
        This function adds the data from coupling_data.py to the pigment's data_dict.

        Args:
            dict_coupling: Dictionary of coupling data containing:
                - tresp_atoms: Atom names for TrEsp calculation.
                - tresp_pc: Point charges for TrEsp.
                - vacuum_mag: Vacuum dipole magnitude.
                - dipole_mag: Dipole magnitude.
        """
        self.dict_data['tresp_atoms'] = dict_coupling['tresp_atoms']
        self.dict_data['tresp_pc'] = dict_coupling['tresp_pc']
        self.dict_data['vacuum_mag'] = dict_coupling['vacuum_mag']
        self.dict_data['dipole_mag'] = dict_coupling['dipole_mag']

    def get_q_0011_charge(self, dict_q0011_charges):
        """Adds q00 and q11 charge data to the pigment's data_dict.

        Args:
            dict_q0011_charges: Dictionary of charge data containing:
                - atom: List of atom names for charge assignment.
                - q_00: Ground state charges.
                - q_11: Excited state charges.
        """
        self.dict_data['q_0011_atom'] = dict_q0011_charges['atom']
        self.dict_data['q_00'] = dict_q0011_charges['q_00']
        self.dict_data['q_11'] = dict_q0011_charges['q_11']


class DefaultPigment(PigmentAtomic):
    """Base class for managing pigments with lineshape functionality.
    
    This class inherits from PigmentAtomic class to create the base class for managing
    pigments. The core of this class is to implement methods to manage the pigment lineshape
    that describes the electronic coupling to the bath of vibrational degrees of freedom.
    
    Attributes:
        _lineshape: The lineshape object describing vibrational coupling.
        domain: The domain to which this pigment belongs.
    """
    def __init__(self, residue, lineshape):
        """Initializes the DefaultPigment object.
        
        Args:
            residue: A residue object from the Biopython Residue class.
            lineshape: The lineshape object for this pigment.
        """
        super().__init__(residue)
        self._lineshape = lineshape
        if self._lineshape is None:
            print(f'WARNING: Lineshape of {self} is None. No spectral calculations can be performed until this updated.')
        self.domain = None
        self.domain_index = 0  # Default domain index

    def update_lineshape(self, lineshape):
        """Updates the lineshape of the pigment.
        
        Args:
            lineshape: The new lineshape object.
        """
        self._lineshape = lineshape

    @property
    def explicit_vib(self):
        """Returns whether explicit vibrational modes are included.
        
        Returns:
            bool: True if explicit vibrational modes are included.
        """
        return self._lineshape.explicit_vib

    @property
    def t_axis(self):
        """Returns the time axis for lineshape calculations.
        
        Returns:
            array: The time axis array.
        """
        return self._lineshape.t_axis

    def g_t(self, temp):
        """Returns the lineshape function g(t) at a given temperature.
        
        Args:
            temp: Temperature in Kelvin.
            
        Returns:
            array: The lineshape function values.
        """
        return self._lineshape.g_t(temp)

    def vibronic_absorption(self, temp):
        """Returns the vibronic absorption spectrum at a given temperature.
        
        Args:
            temp: Temperature in Kelvin.
            
        Returns:
            array: The vibronic absorption spectrum.
        """
        return self._lineshape.vibronic_absorption(temp)

    def vibronic_fluorescence(self, temp):
        """Returns the vibronic fluorescence spectrum at a given temperature.
        
        Args:
            temp: Temperature in Kelvin.
            
        Returns:
            array: The vibronic fluorescence spectrum.
        """
        return self._lineshape.vibronic_fluorescence(temp)

    @property
    def e_lambda(self):
        """Returns the reorganization energy.
        
        Returns:
            float: The reorganization energy.
        """
        return self._lineshape.e_lambda

    @property
    def get_dipole_dir(self, **kwargs):
        """Returns the dipole direction.
        
        Args:
            **kwargs: Additional keyword arguments.
            
        Returns:
            array: The dipole direction vector.
        """
        return self.get_dipole_dir(**kwargs)

    @property
    def dipole(self):
        """Returns the full dipole moment vector.
        
        Returns:
            array: The dipole moment vector (magnitude * direction).
        """
        return self.dict_data['dipole_mag']*self.get_dipole_dir()

    @property
    def _explicit_vib_s(self):
        """Returns the Huang-Rhys factors for explicit vibrational modes.
        
        Returns:
            array: The Huang-Rhys factors.
        """
        return self._lineshape._explicit_vib_s

    @property
    def _explicit_vib_w(self):
        """Returns the frequencies for explicit vibrational modes.
        
        Returns:
            array: The vibrational frequencies.
        """
        return self._lineshape._explicit_vib_w

    @property
    def fc_00(self):
        """Returns the Franck-Condon factor for 0-0 transition.
        
        Returns:
            float: The Franck-Condon factor.
        """
        if self.explicit_vib:
            return np.prod([np.exp(-svib/2) for svib in self._explicit_vib_s])
        else:
            return 1


class PorphyrinAtomic(DefaultPigment):
    """Atomic porphyrin ring pigment class.
    
    This class inherits from the DefaultPigment class to create an atomic porphyrin ring
    object. It is intended to act as a base class for the large set of pigment types which
    derive from porphyrin rings (e.g., chlorophylls, pheophytins).
    
    Attributes:
        qy_atoms: Tuple of atoms defining the Qy dipole vector direction.
        qx_atoms: Tuple of atoms defining the Qx dipole vector direction.
    """

    def __init__(self, residue, lineshape=None, qy_atoms=None, qx_atoms=None):
        """Initializes the PorphyrinAtomic object.

        Args:
            residue: A residue object from the Biopython Residue class.
            lineshape: The lineshape object for this pigment. Defaults to None.
            qy_atoms: Tuple of atom names (str) that define the direction of the Qy
                dipole vector. Defaults to None, which auto-detects NB-ND or N1B-N1D.
            qx_atoms: Tuple of atom names (str) that define the direction of the Qx
                dipole vector. Defaults to None, which auto-detects NA-NC or N1A-N1C.
        """
        super().__init__(residue, lineshape)
        # REMEMBER: atomic nomenclature is reversed in the pdb files
        if qy_atoms is None:
            if 'NB' in self.residue.child_dict.keys():
                self.qy_atoms = (self.residue['NB'], self.residue['ND'])
            elif 'N1B' in self.residue.child_dict.keys():
                self.qy_atoms = (self.residue['N1B'], self.residue['N1D'])
            else:
                print(f'{self.name}: Qy dipole not defined')
        else:
            self.qy_atoms = (self.residue[qy_atoms[0]], self.residue[qy_atoms[1]])
        if qx_atoms is None:
            if 'NA' in self.residue.child_dict.keys():
                self.qx_atoms = (self.residue['NA'], self.residue['NC'])
            elif 'N1A' in self.residue.child_dict.keys():
                self.qx_atoms = (self.residue['N1A'], self.residue['N1C'])
            else:
                print(f'{self.name}: Qx dipole not defined')
        else:
            self.qx_atoms = (self.residue[qx_atoms[0]], self.residue[qx_atoms[1]])

    def __repr__(self):
        return f'{type(self).__name__} {self.name}: ' \
               f'qy dipole direction {self.get_dipole_dir()}: ' \
               f'location {self.location}'

    @property
    def name(self):
        """Returns the unique name identifier for this pigment.
        
        Returns:
            str: A name string constructed from PDB ID, chain ID, residue name,
                and residue number.
        """
        return f'{self.residue.get_parent().get_parent().full_id[0]}_{self.residue.get_parent().id}_{self.residue.resname}_{self.residue.id[1]}'

    def get_dipole_dir(self, tdm='qy'):
        """Returns a unit vector in the direction of the dipole moment.

        Args:
            tdm: The axis that defines the transition dipole moment. Can be
                either 'qy' or 'qx'. Defaults to 'qy'.

        Returns:
            array: A unit vector describing the X,Y,Z coordinates of the transition
                dipole moment direction.
                
        Raises:
            ValueError: If tdm is not 'qy' or 'qx'.
        """
        if tdm == 'qy':
            v = self.qy_atoms[0].coord - self.qy_atoms[1].coord
            return v / np.linalg.norm(v)
        elif tdm == 'qx':
            v = self.qx_atoms[0].coord - self.qx_atoms[1].coord
            return v / np.linalg.norm(v)
        else:
            raise ValueError('Unknown dipole for chlorin-based pigment. Allowed values: qy, qx')

    def get_coupling_data(self, dict_coupling):
        """Adds coupling data with automatic nitrogen atom name correction.
        
        This method intercepts the dictionaries early and corrects the most common naming
        error by mapping the nitrogen atom names to the correct version for this instance
        of the pigment (NA/NB/NC/ND vs N1A/N1B/N1C/N1D).

        Args:
            dict_coupling: Dictionary of coupling data.
        """
        new_dict_coupling = None

        if (dict_coupling['tresp_atoms'] is not None) and ('N1A' in self.list_atom_names) and ('NA' in dict_coupling['tresp_atoms']):
            new_dict_coupling = copy.deepcopy(dict_coupling)
            for (atom_name_tresp, new_atom_name) in zip(['NA', 'NB', 'NC', 'ND'],
                                                  ['N1A', 'N1B', 'N1C', 'N1D']):
                index_atom = new_dict_coupling['tresp_atoms'].index(atom_name_tresp)
                new_dict_coupling['tresp_atoms'][index_atom] = new_atom_name
        elif (dict_coupling['tresp_atoms'] is not None) and ('NA' in self.list_atom_names) and ('N1A' in dict_coupling['tresp_atoms']):
            new_dict_coupling = copy.deepcopy(dict_coupling)
            for (atom_name_tresp, new_atom_name) in zip(['N1A', 'N1B', 'N1C', 'N1D'],
                                                  ['NA', 'NB', 'NC', 'ND']):
                index_atom = new_dict_coupling['tresp_atoms'].index(atom_name_tresp)
                new_dict_coupling['tresp_atoms'][index_atom] = new_atom_name
        else:
            new_dict_coupling = dict_coupling

        super().get_coupling_data(new_dict_coupling)

    @property
    def location(self):
        """Returns the location of the porphyrin ring.
        
        The location is calculated as the average position of the four nitrogen atoms
        in the porphyrin ring.

        Returns:
            array: An array containing the X,Y,Z coordinates of the porphyrin center.
        """
        return (self.residue['NA'].coord + self.residue['NB'].coord +
                self.residue['NC'].coord + self.residue['ND'].coord) / 4

    def visualize_dipole(self, ax, tdm='qy', color='k', scale_factor=1.):
        """Constructs a visualization of the dipole vector.

        Args:
            ax: The matplotlib axes object for the visualization.
            tdm: The transition dipole moment to visualize. Can be 'qy' or 'qx'.
                Defaults to 'qy'.
            color: The color of the dipole moment on the plot. Defaults to 'k' (black).
            scale_factor: Scaling factor for the dipole vector length. Defaults to 1.
        """
        vector = self.get_dipole_dir(tdm)

        x0 = [-0.5 * vector[0] + self.location[0]]
        y0 = [-0.5 * vector[1] + self.location[1]]
        z0 = [-0.5 * vector[2] + self.location[2]]
        dx = [vector[0] * scale_factor]
        dy = [vector[1] * scale_factor]
        dz = [vector[2] * scale_factor]
        ax.quiver(x0, y0, z0, dx, dy, dz, color=color, lw=1)

    def get_atom_xyz(self, atom):
        """Returns the x,y,z coordinates of a specified atom.

        Args:
            atom: Atom name (str) that you want to get the coordinates of.

        Returns:
            array: An array containing the x,y,z position of the atom.
        """
        return self.residue[atom].coord



class ChlorophyllAtomic(PorphyrinAtomic):
    """Atomic chlorophyll pigment class.
    
    This class inherits from the PorphyrinAtomic class to create an atomic chlorophyll
    object. Chlorophylls contain a central magnesium atom.
    """
    def __init__(self, residue, lineshape=None, qy_atoms=None, qx_atoms=None):
        """Initializes the ChlorophyllAtomic object.

        Args:
            residue: A residue object from the Biopython Residue class.
            lineshape: The lineshape object for this pigment. Defaults to None.
            qy_atoms: Tuple of atom names (str) that define the direction of the Qy
                dipole vector. Defaults to None.
            qx_atoms: Tuple of atom names (str) that define the direction of the Qx
                dipole vector. Defaults to None.
        """
        super().__init__(residue, lineshape, qy_atoms, qx_atoms)

    @property
    def location(self):
        """Returns the location of the chlorophyll.
        
        The location is defined as the coordinates of the central magnesium atom.

        Returns:
            array: An array containing the X,Y,Z coordinates of the central
                magnesium atom.
        """
        return self.residue['MG'].coord


class PheophytinAtomic(PorphyrinAtomic):
    """Atomic pheophytin pigment class.
    
    This class inherits from the PorphyrinAtomic class to create an atomic pheophytin
    object. Pheophytins are like chlorophylls but without the central magnesium atom.
    """
    def __init__(self, residue, lineshape=None, qy_atoms=None, qx_atoms=None):
        """Initializes the PheophytinAtomic object.

        Args:
            residue: A residue object from the Biopython Residue class.
            lineshape: The lineshape object for this pigment. Defaults to None.
            qy_atoms: Tuple of atom names (str) that define the direction of the Qy
                dipole vector. Defaults to None.
            qx_atoms: Tuple of atom names (str) that define the direction of the Qx
                dipole vector. Defaults to None.
        """
        super().__init__(residue, lineshape, qy_atoms, qx_atoms)
