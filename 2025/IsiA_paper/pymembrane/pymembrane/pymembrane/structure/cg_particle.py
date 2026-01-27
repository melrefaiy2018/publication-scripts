from typing import Literal, Union

import matplotlib.pyplot as plt
import numba
import numpy as np
from mpl_toolkits.mplot3d import art3d
from numpy.typing import ArrayLike, NDArray
from shapely import MultiPoint
from shapely import convex_hull as shapely_convex_hull

from pymembrane.structure.atomic_pigment import (ChlorophyllAtomic, DefaultPigment, PheophytinAtomic, PigmentAtomic,
                                                 PorphyrinAtomic)
from pymembrane.structure.atomic_pigment import default_colors
from pymembrane.structure.atomic_protein import (DynamicProteinAtomic, PigmentProteinAtomic, StaticProteinAtomic)
from pymembrane.util.linear_algebra import construct_rotation_matrix

AtomicProteinObject = Union[
    StaticProteinAtomic,
    DynamicProteinAtomic,
    PigmentProteinAtomic
]

AtomicPigmentObject = Union[
    PigmentAtomic,
    DefaultPigment,
    PorphyrinAtomic,
    ChlorophyllAtomic,
    PheophytinAtomic
]


@numba.njit()
def array_map_position(
        R2_pos: NDArray[float],
        R2_orient: NDArray[float],
        atomic_center_of_mass: NDArray[float],
        atomic_offset: NDArray[float],
        cg_location: NDArray[float]
) -> NDArray[float]:
    """
    Maps an (N,3) array of atom positions to a new array of coarse-grained positions.

    Parameters
    ----------
    1.) R2_pos: np.array(floats)
                An Nx3 array of each atom's position in the protein.

    2.) R2_orient: np.array(floats)
                   A matrix that defines the orientation of the protein.

    3.) atomic_center_of_mass: np.array(floats)
                               The point of the atomic protein object's center of mass.

    4.) atomic_offset: np.array(floats)
                       The offset between the positions of the atomic protein and the
                       coarse-grained protein.

    5.) cg_location: np.array(floats)
                     The location of the coarse-grained protein object's center of mass.

    Returns
    -------
    1.) list_map_pos: np.array(floats)
                      An array of positions for each atom's position in the
                      coarse-grained protein object.
    """
    list_map_pos = np.zeros_like(R2_pos)
    for (index, pos) in enumerate(R2_pos):
        list_map_pos[index,:] = ((pos - (atomic_center_of_mass+atomic_offset)) @
                                 R2_orient + cg_location)
    return list_map_pos


class CGProtein:
    """
    This is a class that constructs a coarse-grained protein object. This class will not
    contain all the atomic information about the protein, but instead will point back
    to the ProteinAtomic object if necessary.
    """
    __slots__ = [
        'atomic',
        '__name',
        '__hull',
        '_Z1_norm',
        '__location',
        '__original_location',
        'R2_orient',
        '_atomic_offset',
        'dict_pigments',
        '_restraint_potential',
        'list_pigments',
        '_fval_tresp',
        '_dielectric_tresp',
        '_dielectric_dipole'
    ]
    def __init__(
            self,
            atomic_protein: AtomicProteinObject,
            name: str,
            location: NDArray[float] = None,
            R2_orient: NDArray[float] = None,
            atomic_offset: NDArray[float] = None,
    ) -> None:
        """
        This function instantiates the CGProtein object and initializes its
        attributes.

        Parameters
        ----------
        1.) atomic_protein : object
                             An atomic protein object created from the ProteinAtomic
                             class.

        2.) name: str
                  The name of the coarse-grained protein object

        3.) location: np.array(floats), optional
                      An array that defines the location of the protein using X,Y,Z
                      coordinates.

        3.) R2_orient : array(floats)
                        A matrix that defines the orientation of the protein.

        5.) atomic_offset: np.array(floats), optional
                           The offset between the positions of the atomic protein and
                           the coarse-grained protein.

        Returns
        -------
        None
        """
        self.atomic = atomic_protein
        self.__name = name
        self.__hull = None
        self._Z1_norm = None
        self._restraint_potential = None

        # Define Coupling Parameters
        # --------------------------
        self._fval_tresp = None
        self._dielectric_tresp = None
        self._dielectric_dipole = None

        if location is None:
            self.__location = np.zeros(3, dtype=np.float64)
        else:
            self.__location = np.array(location, dtype=np.float64)

        # Initialize the original location
        self.__original_location = self.__location.copy()

        if R2_orient is None:
            self.R2_orient = np.eye(3)
        else:
            self.R2_orient = R2_orient

        if atomic_offset is None:
            self._atomic_offset = np.zeros(3, dtype=np.float64)
        else:
            self._atomic_offset = atomic_offset

        self.dict_pigments = {f'{self.name}_{key}': CGPigment(value, self)
                              for (key, value) in self.atomic.dict_pigments.items()}
        self.list_pigments = list(self.dict_pigments.values())

    def __repr__(self):
        return f'CGParticle({type(self.atomic).__name__}) {self.name}: Loc({self.__location}), ' \
               f'Orientation({self.R2_orient}), AtomicOffset({self._atomic_offset})\n'\
               f'* AtomicProperties: {self.atomic.__repr__()}'

    def __str__(self):
        return f'CGParticle({type(self.atomic).__name__}, {self.atomic.name}) {self.name}: Loc({self.__location}),' \
               f'Orientation({self.R2_orient}), AtomicOffset({self._atomic_offset})\n'

    @property
    def name(self):
        """
        Access the name of the CGProtein object.
        """
        return self.__name

    @property
    def center_of_mass(self):
        """
        Access the center of mass of the CGProtein object.
        """
        return self.atomic.center_of_mass

    @property
    def type(self):
        """
        Access the type of protein.
        """
        return self.atomic.type

    @property
    def Z1_norm(self):
        """
        Access the Z1_norm of the protein. If not given set to the normalized k-vector.
        """
        if self._Z1_norm is None:
            return np.array([0, 0, 1])
        else:
            return self._Z1_norm

    def update_original_location(self):
        """
        Save the current location as the original location.
        """
        self.__original_location = self.__location.copy()

    @property
    def original_location(self):
        """
        Access the original location.
        """
        return self.__original_location

    def set_original_location(self, new_location: ArrayLike) -> None:
        """
        Explicitly override the original location.

        This method should only be used when you intentionally want to reset the
        original location.

        Parameters
        ----------
        1.) new_location: ArrayLike
                          The new original location to set.

        Returns
        -------
        None
        """
        if (not isinstance(new_location, (list, tuple, np.ndarray)) or
                len(new_location) != 3):
            raise ValueError("Original location must be a 3-element array-like "
                             "structure.")

        self.__original_location = np.array(new_location, dtype=np.float64)

    def translate(
            self,
            vector: NDArray[float],
            update_original: bool = False
    ) -> None:
        """
        Translates the location of a coarse-grained protein  object by the
        coordinates specified in the input vector. If specified, the input vector can
        become the original location.

        Parameters
        ----------
        1.) vector: np.array(floats)
                    A vector that defines the translation.

        2.) update_original: bool, optional
                             If True, update the original location to the new location.

        Returns
        -------
        None
        """
        self.__location += np.array(vector, dtype=np.float64)
        if update_original:
            self.update_original_location()

    def rotate_by_axis(self, angle: float, axis: str = 'z') -> None:
        """
        This function rotates a coarse-grained protein object by the angle
        specified in the input along the specified axis.

        Parameters
        ---------
        1.) angle : float
                    The angle of rotation. Angle should be in radians.

        2.) axis: str, optional
                  The axis along which the rotation will take place. If not specified,
                  rotation will take place about the z-axis.

        Returns
        -------
        None
        """
        self.rotate_by_matrix(construct_rotation_matrix(angle, axis))

    def rotate_by_matrix(self, R2_rotate: NDArray[float]) -> None:
        """
        This function transforms the orientation of the coarse-gained protein
        using the input rotation matrix.

        Parameters
        ----------
        1.) R2_rotate : array(floats)
                        A rotation matrix defining the transformation.

        Returns
        -------
        None
        """
        self.R2_orient = self.R2_orient @ R2_rotate

    @property
    def dict_data(self):
        """
        Access the user-specified dictionary of data pertaining to the atomic protein
        object that is used to initialize the CGProtein.
        """
        return self.atomic.dict_data

    @property
    def H2_hamiltonian(self):
        return self.atomic.H2_hamiltonian

    def visualize_pigment_location(
            self,
            ax: plt.Axes,
            color: str ='k'
    ) -> None:
        """
        Creates a 2D scatter plot of the pigment locations for a CGProtein object. This
        will only work if the atomic object is a PigmentProteinAtomic & has a nonempty
        dict_pigments.

        Parameters
        ----------
        1.) ax: plt.Axes
                The axes along which the pigment will appear.

        2.) color: str, optional
                   The color of each point on the scatter plot.

        Returns
        -------
        1.) dict_data : dict
                        A dictionary of data pertaining to the protein.
        """
        return self.atomic.dict_data

    @property
    def H2_hamiltonian(self):
        return self.atomic.H2_hamiltonian

    def visualize_pigment_location(self, ax, color='k'):
        R3_pigment_positions = np.zeros([self.atomic.n_pigments,3], dtype=np.float32)
        for (index, pigment) in enumerate(self.dict_pigments.values()):
            R3_pigment_positions[index,:] = pigment.location

        ax.scatter(R3_pigment_positions[:, 0],
                   R3_pigment_positions[:, 1],
                   R3_pigment_positions[:, 2],
                   color=color, marker = 'o')

    def _map_position(self, pos: NDArray[float]) -> NDArray[float]:
        """
        This function will map the position of the pigment from the atomic
        coordinates to the coarse-grained coordinates.

        Parameters
        ----------
        1.) pos : array
                  An array containing x,y,z coordinates that define the atomic
                  position of the pigment.
                  # Note: if using the array version it should be [N_atoms, 3]

        Returns
        -------
        1.) cg_pos : array
                     An array containing x,y,z coordinates that represent the atomic
                     position of the pigment in the coarse-grained coordinate system.
        """
        if np.shape(np.shape(pos))[0] == 2:
            return array_map_position(pos, self.R2_orient, self.atomic.center_of_mass,
                                      self._atomic_offset, self.__location)
        elif np.shape(np.shape(pos))[0] == 1:
            return (pos - (self.atomic.center_of_mass+self._atomic_offset))@self.R2_orient + self.__location
        else:
            # TODO: This needs to be a raise error call.
            print('ERROR!')

    def get_models(self):
        """
        Access the models of the atomic object's PDB file.
        """
        return self.atomic.get_models()

    def get_chains(self):
        """
        Access the polypeptide chains from the atomic object's PDB file.
        """
        return self.atomic.get_chains()

    def get_residues(self):
        """
        Access the residues from the atomic object's PDB file.
        """
        return self.atomic.get_residues()

    def get_atoms(self):
        """
        Access the atoms from the atomic object's PDB file.
        """
        return self.atomic.get_atoms()

    @property
    def R2_residue_center(self):
        """
        Access the array of each residue's center of mass.
        """
        return self._map_position(self.atomic.R2_residue_center)

    @property
    def R2_atomic_xyz(self):
        """
        Access the array of each atom's xyz coordinates.
        """
        return self._map_position(np.array([atom.coord for atom in
                                            self.atomic.get_atoms()]))

    @property
    def list_covalent_radii(self):
        """
        Access the list of each atom's covalent radii.
        """
        return self.atomic.list_covalent_radii

    @property
    def list_vdw_radii(self):
        """
        Access the list of each atom's Van der Waals radii.
        """
        return self.atomic.list_vdw_radii

    @property
    def list_wca_epsilons(self):
        """
        Access the list of each atom's Weeks-Chandler-Andersen (WCA) potential well
        depths (epsilon).
        https://webff-documentation.readthedocs.io/en/latest/Reference/NonBond-Weeks-Chandler-Anderson.html
        """
        return self.atomic.list_wca_epsilons

    @property
    def list_wca_sigmas(self):
        """
        Access the list of each atom's Weeks-Chandler-Andersen (WCA) zero-potential
        radii (sigmas).
        https://webff-documentation.readthedocs.io/en/latest/Reference/NonBond-Weeks-Chandler-Anderson.html
        """
        return self.atomic.list_wca_sigmas

    @property
    def location(self):
        """
        Access the location of the CGProtein object's center of mass.
        """
        return self.__location

    @property
    def hull3d(self):
        """
        Access the 3D convex hull of the atomic object.
        """
        return self.atomic.hull3d

    @property
    def hull3d_vertices(self):
        """
        Access the vertices of the 3D convex hull.
        """
        return self.R2_residue_center[self.hull3d.vertices]

    @property
    def hull2d(self):
        """
        Access the 2D convex hull of the atomic object.
        """
        return shapely_convex_hull(MultiPoint(self.R2_residue_center))

    @property
    def hull2d_vertices(self):
        """
        Access the vertices of the 2D convex hull.
        """
        return self._map_position(self.atomic.hull2d_vertices)

    @property
    def radius(self):
        """
        Access the radius of the circle centered at the CGProtein object center of mass
        which will inscribe the entire 2D convex hull.
        """
        return self.atomic.radius

    @property
    def concave_hull2d_ratio(self):
        """
        Get the alpha parameter for the concave hull.

        Alpha determines the maximum radius allowed for the circumcircles of the
        Delaunay triangles made from the set of points. In simple terms, this determines
        the concavity of the hull.

        https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-018-1305-z
        https://gwlucastrig.github.io/TinfourDocs/DelaunayIntro/index.html
        """
        return self.atomic.concave_hull2d_ratio

    @property
    def concave_hull2d(self):
        """
        Access the 2D concave hull of the atomic object with concave ratio set by
        concave_hull2d_alpha.
        """
        return self.atomic.concave_hull2d

    @property
    def concave_hull2d_vertices(self):
        """
        Access the vertices of the 2D concave hull with concave ratio set by
        concave_hull2d_alpha.
        """
        return self._map_position(self.atomic.concave_hull2d_vertices)

    @property
    def inner_radius(self):
        """
        Access the radius of the smallest possible circle that can fit entirely inside the concave hull.
        """
        return self.atomic.inner_radius

    def visualize_protein_hull3d(
            self,
            ax: plt.Axes = None,
            color: str ='g'
    ) -> None:
        """
        Creates a plot of the 3D convex hull of the atomic object as a triangulated surface.

        Parameters
        ----------
        1.) ax: plt.Axes, optional
                The axes along which the 3D protein hull will appear.

        2.) color: str, optional
                   The color of each triangle in the figure.

        Returns
        -------
        None
        """
        #TODO: consistent default values of ax
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

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

    def visualize_protein_hull2d(
            self,
            ax: plt.Axes = None,
            hull_type: Literal['convex', 'concave'] = 'convex',
            color: str = 'g',
            plot_label: str = None
    ) -> plt.Axes:
        """
        Creates a plot of the 2D hull (convex or concave) of the atomic
        object.

        Parameters
        ----------
        1.) ax: plt.Axes, optional
                The axes along which the 2D hull will appear.

        2.) hull_type: Literal['convex', 'concave'], optional
                       Determines whether the convex or concave hull will be plotted.

        3.) color: str, optional
                   The color of the polygon that represents the hull.

        4.) plot_label: str, optional
                        Label for the plot.


        """
        if ax is None:
            fig, ax = plt.subplots()
        ax.set_aspect('equal')

        if hull_type == 'convex':
            vertices = np.vstack([
                self.hull2d_vertices,
                self.hull2d_vertices[0]
            ])
        elif hull_type == 'concave':
            vertices = np.vstack([
                self.concave_hull2d_vertices,
                self.concave_hull2d_vertices[0]
            ])
        else:
            raise ValueError("hull_type must be 'convex' or 'concave'")

        ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=1)
        ax.fill(
                vertices[:, 0],
                vertices[:, 1],
                color=color,
                alpha=0.65,
                label=plot_label
        )

        return ax

    def visualize_orientation(
            self,
            ax: plt.Axes,
            axis: str ,
            color: str = 'g',
            scale_factor: float = 100.,
            line_width: float = 5
    ) -> None:
        """
        Creates a plot of a 3D arrow representing the orientation
        of the CGProtein object with respect to a specified axis.

        1.) ax: plt.Axes
                 The axes of a Matplotlib plot in 3D

        2.) axis: str, optional
                   The axis about which the orientation will be visualized.

        3.) color: str, optional
                    The color of the arrow on the plot.

        4.) scale_factor: float, optional
                          The factor by which the length of the arrow will be scaled on
                          the plot.

        5.) line_width: float, optional
                        The width of the arrow on the plot.

        Returns
        -------
        None
        """
        index_axis = 'xyz'.index(axis)
        unit_vector = np.zeros(3)
        unit_vector[index_axis] = 1
        vector = scale_factor * (self.R2_orient @ unit_vector)
        ax.quiver(self.__location[0], self.__location[1], self.__location[2],
                  vector[0], vector[1], vector[2],
                  color=color,
                  lw=line_width)

    @property
    def disk_radius(self):
        if self._disk_radius is not None:
            return self._disk_radius
        else:
            return self.radius


class CGPigment:
    """
    This is a class that constructs a course-grained pigment object. This class will
    not contain all of the atomic information about the pigment, but instead will
    point back to the PigmentAtomic object if necessary.
    """

    def __init__(
            self,
            atomic_pigment: AtomicPigmentObject,
            cg_protein: CGProtein
    ) -> None:
        """
        Instantiates the CGPigment object and initializes its attributes.

        Parameters
        ----------
        1.) atomic_pigment: AtomicPigmentObject
                            An atomic pigment object created in one of the classes that
                            inherit from the ResidueAtomic class.

        2.) cg_protein: CGProtein
                        A coarse-grained protein object created in the CGProtein class.

        Returns
        -------
        None
        """
        self.atomic = atomic_pigment
        self.cg_protein = cg_protein

    def __repr__(self):
        return f'CGPigment({type(self.atomic).__name__}, {self.name}): ' \
               f'qy dipole direction {self.get_dipole_dir()}: ' \
               f'location {self.location}'

    def __str__(self):
        return f'CGPigment({type(self.atomic).__name__}, {self.name})'

    @property
    def name(self):
        """
        Access the name of the atomic protein and the atomic pigment objects.
        """
        return f'{self.cg_protein.name}_{self.atomic.name}'
    def get_dipole_dir(self, *args, **kwargs) -> NDArray[float]:
        """
        Returns the coordinates that represent the direction of the
        transition dipole moment of the pigment.

        Returns
        -------
        1.) tdm: np.array(floats)
                    A vector describing the X,Y,Z coordinates of the transition dipole
                    moment.
        """
        return self.cg_protein.R2_orient @ (self.atomic.get_dipole_dir(*args, **kwargs))

    @property
    def location(self):
        """
        This function returns the location of the pigment (defined by the position of the center of mass)

        Parameters
        ----------
        None

        Returns
        -------
        1.) location : array(floats)
                       An array containing the X,Y,Z coordinates that define the
                       location of the pigment.
        """
        return self._map_position(self.atomic.location)

    @property
    def R2_atomic_xyz(self):
        """
        This function returns all atomic coordinates of the pigment mapped to the
        coarse-grained coordinate system.

        Parameters
        ----------
        None

        Returns
        -------
        1.) R2_atomic_xyz : array(floats)
                            An array of shape [N_atoms, 3] containing the X,Y,Z
                            coordinates of all atoms in the pigment, mapped to the
                            coarse-grained coordinate system.
        """
        return self._map_position(np.array([atom.coord for atom in self.atomic.residue.get_atoms()]))

    @property
    def mass(self):
        """
        This function returns the mass of the PigmentAtomic object

        Parameters
        ----------
        None

        Returns
        -------
        1.) mass : float
                   The mass of the pigment.
        """
        return self.atomic.mass

    @property
    def dict_data(self):
        """
        This function returns the dictionary of user-specified data pertaining to the
        atomic pigment.

        Parameters
        ----------
        None

        Returns
        -------
        1.) dict_data : dict
                        A dictionary of data pertaining to the pigment object.
        """
        return self.atomic.dict_data

    def visualize(self, ax):
        """
        This function constructs a visualization of the pigment.

        Parameters
        ----------
        1.) ax : object
                 The axes along which the visualization will appear.

        Returns
        -------
        None
        """
        self._visualize_atoms(ax)
        self._visualize_bonds(ax)

    def _map_position(self, pos):
        """
        This function will map the position of the pigment from the atomic
        coordinates to the coarse-grained coordinates.

        Parameters
        ----------
        1.) pos : array
                  An array containing x,y,z coordinates that define the atomic
                  position of the pigment.
                  # Note: if using the array version it should be [N_atoms, 3]

        Returns
        -------
        1.) cg_pos : array
                     An array containing x,y,z coordinates that represent the atomic
                     position of the pigment in the coarse-grained coordinate system.
        """
        if np.shape(np.shape(pos))[0] == 2:
            return array_map_position(pos, self.cg_protein.R2_orient, 
                                     self.cg_protein.atomic.center_of_mass,
                                     self.cg_protein._atomic_offset, 
                                     self.cg_protein.location)
        elif np.shape(np.shape(pos))[0] == 1:
            return self.cg_protein.R2_orient @ (pos - (self.cg_protein.atomic.center_of_mass+self.cg_protein._atomic_offset)) + self.cg_protein.location
        else:
            # TODO: This needs to be a raise error call.
            print('ERROR!')


    def _visualize_atoms(self, ax):
        """
        This function constructs a visualization of the atoms in a pigment.

        Parameters
        ----------
        1.) ax : object
                 The axes along which the visualization will appear.

        Returns
        -------
        None
        """
        for atom in self.atomic.residue.get_atoms():
            if atom.element in default_colors.keys():
                color = default_colors[atom.element]
            else:
                color = 'purple'
            coord = self._map_position(atom.coord)
            ax.scatter([coord[0]], [coord[1]], [coord[2]],
                       marker='o',
                       s=30,
                       color=color)

    def _visualize_bonds(self, ax: plt.Axes) -> None:
        """
        Constructs a visualization of the bonds within a pigment.

        Parameters
        ----------
        1.) ax: plt.Axes
                The axes along which the visualization of bonds will appear.

        Returns
        -------
        None
        """
        for (atom1_id, atom2_id) in self.atomic.bonds:
            try:
                atom1 = self._map_position(self.atomic.residue[atom1_id].coord)
                atom2 = self._map_position(self.atomic.residue[atom2_id].coord)
                ax.plot([atom1[0], atom2[0]],
                        [atom1[1], atom2[1]],
                        [atom1[2], atom2[2]],
                        linewidth=3, color='k', alpha=0.5)
            except:
                pass
    def visualize_dipole(
            self,
            ax: plt.Axes,
            color: str ='k',
            scale_factor: float = 1.,
            *args,
            **kwargs
    ) -> None:
        """
        Constructs a visualization of the dipole vector of a pigment.

        Parameters
        ----------
        1.) ax: plt.Axes
                The axes along which the visualization will appear.

        2.) color: str, optional
                   The color of the dipole moment on the plot.

        3.) scale_factor: float, optional
                          The factor to which the dipole vector will be scaled on the
                          plot.

        Returns
        -------
        None
        """
        vector = self.get_dipole_dir(*args, **kwargs)
        location = self.location
        x0 = [-0.5*vector[0] + location[0]]
        y0 = [-0.5*vector[1] + location[1]]
        z0 = [-0.5*vector[2] + location[2]]
        dx = [vector[0] * scale_factor]
        dy = [vector[1] * scale_factor]
        dz = [vector[2] * scale_factor]

        ax.quiver(x0, y0, z0, dx, dy, dz, color=color)

    def get_atom_xyz(self, atom):
        """
        This function returns a list giving the x,y,z coordinates for each atom given in
        the atom_list input.

        Parameters
        ----------
        1.) atom : str
                    Atom that you want to get the coordinates of.

        Returns
        -------
        1.) atom_position : tuple(float)
                            A list containing the x,y,z position the atom.
        """
        return self._map_position(self.atomic.residue[atom].coord)

    def get_coupling_data(self, dictionary_name):
        """
        This function will add the data that is necessary to perform a tresp coupling
        calculation to the pigments data_dict. The data comes from coupling_data.py

        Parameters
        ----------
        1.) dictionary_name : str
                              dictionary you want to retrieve data from

        Returns
        -------
        None

        """
        self.atomic.get_coupling_data(dictionary_name)


class CGProteinFlat(CGProtein):
    """
    This class inherits from CGProtein and flattens a coarse-grained protein into 
    2 dimensions. 
    """
    def __init__(
            self,
            atomic_protein: AtomicProteinObject,
            name: str,
            location: NDArray = None,
            angle: float = None,
            atomic_offset: NDArray = None,
    ) -> None:
        """
        Instantiates the CGProteinFlat object and initializes its attributes.

        Parameters
        ----------
        1.) atomic_protein: AtomicProteinObject
                            An atomic protein object created from atomic_protein.

        2.) name: str
                  The name of the CGProteinFlat object.

        3.) location: array(floats), optional
                      An array specifying the center of mass of the CGProteinFlat object.

        4.) angle: float, optional
                   The angle by which to rotate the CGProteinFlat object about the
                   z-axis in radians.

        5.) atomic_offset: array(floats), optional
                            An array specifying the offset between the atomic protein's
                            position and the coarse-grained protein's position.

        Returns
        -------
        None

        """
        #TODO: set np.eye(3) to np.eye(3, dtype=float)
        super().__init__(
            atomic_protein,
            name,
            location,
            np.eye(3),
            atomic_offset,
        )

        if angle is not None:
            self.__angle = angle
            self.rotate(angle)
        else:
            self.__angle = 0

    def rotate(self, angle: float) -> None:
        """
        Rotates the CGProteinFlat object by a given angle using a rotation matrix
        constructed from a quaternion. Below is a derivation of the formula used to get
        the rotation matrix.

        (Note that this function uses a minus sign for the imaginary components because
        rotate_by_matrix multiplies from the right side)

        https://www.songho.ca/opengl/gl_quaternion.html#:~:text=To%20rotate%20the%20quaternion%20p,back%20of%20the%20previous%20multiplication

        Parameters
        ----------
        1.) angle: float
                   The angle by which to rotate about the z-axis. This should be in
                   radians.

        Returns
        -------
        None
        """
        self.__angle += angle

        # Construct Matrix for rotation around z-axis
        # -------------------------------------------
        axis = np.array([0, 0, 1])
        a = np.cos(angle / 2)
        b, c, d = -axis * np.sin(angle / 2)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
        super().rotate_by_matrix(np.array([
            [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
            [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
            [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
        ]))


    def rotate_by_axis(self, angle: float, axis: str = 'z') -> None:
        """
        Implemented to avoid the use of rotate_by_axis on the flat object.
        """
        #TODO: this should be a proper raised error
        print('ERROR: rotate_by_axis() is not supported for CGProteinFlat '
              'objects.')

    def rotate_by_matrix(self, R2_rotate: NDArray) -> None:
        """
        Implemented to avoid the use of rotate_by_matrix on the flat object.
        """
        print('ERROR: rotate_by_matrix() is not supported for CGProteinFlat objects.')

    @property
    def angle(self):
        """
        Access the angle at which the CGProteinFlat object is rotated about the z-axis
        in radians.
        """
        return self.__angle

    @angle.setter
    def angle(self, new_angle: float) -> None:
        """
        Sets a new angle for the flattened object to be rotated about the z-axis.

        Parameters
        ----------
        1.) new_angle: float
                       The new angle to rotate the CGProteinFlat object to in radians.

        Returns
        -------
        None

        """
        while new_angle<0:
            new_angle += 2*np.pi

        while new_angle>2*np.pi:
            new_angle -= 2*np.pi

        self.__angle = new_angle




