import copy
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import ConnectionPatch
from scipy.sparse import coo_array
from shapely import concave_hull
from shapely.geometry import MultiPoint
from functools import cached_property
from tqdm import tqdm
from pymembrane.exciton.exciton_aggregate import ExcitonicAggregate, ExcitonDomain
from pymembrane.exciton.lineshape import PigmentLineshape_Kubo
from pymembrane.exciton.linear_spectra import LinearSpectra
from pymembrane.structure.membrane import ProteinAggregate
from pymembrane.exciton.forster import generalized_forster_thermal


def eigenprop(k_ij, p0, t_axis):
    """Calculate population dynamics using eigenvalue decomposition of the rate matrix.

    This function accepts an initial population vector and a time axis that defines all the time 
    points of interest. The calculation proceeds by calculating the eigenvectors of the rate matrix, 
    which is a memory intensive calculation but makes the calculation scale favorably with time 
    points (since each time point does not depend on the previous time points).

    Note:
        Definition of k_ij:
        - The element k_ij[i][j] defines the rate i<--j in units of 1/fs
        - The columns of k_ij should sum to zero.

    Args:
        k_ij (numpy.ndarray): Rate matrix where element [i][j] defines rate i<--j. Units: 1/fs.
        p0 (numpy.ndarray): Initial population vector (Nx1 column vector).
        t_axis (numpy.ndarray or list): Time points of interest. Units: fs.

    Returns:
        numpy.ndarray: Population matrix p_t with shape (len(t_axis), N) where N is the number of states.
    """
    # Calculate Dynamics
    # ------------------
    t_axis = np.array(t_axis)
    e_kin, u_kin = np.linalg.eig(k_ij)  # u_kin are column-wise eigenvectors
    u_inv = np.linalg.inv(u_kin)
    p_t = np.zeros([len(t_axis), len(k_ij[0])])

    for (index_t, t) in enumerate(t_axis):
        p_t[index_t, :] = np.dot(u_kin,
                                 np.dot(np.diagflat(np.exp(t * e_kin)),
                                        np.dot(u_inv, p0)))

    return p_t


def numericalprop(k_ij, p0, t_axis):
    '''
    This function will numerically integrate the population dynamics in a series of finite
    time steps (dt) for a total of Nt steps according to the rate matrix (k_ij) and the
    initial population (p0).

    INPUTS: p0, dt, Nt
    OUTPUT: t_axis, p_t
    '''
    pi = p0
    t_old = t_axis[0]
    p_t = [pi]

    for t_i in t_axis[1:]:
        pi = pi + k_ij@pi * (t_i - t_old)
        t_old = t_i
        p_t.append(pi)
    
    return p_t


def runge_kutta_prop(K2_rate, P1_0, t_axis):
    """Propagate population dynamics using 4th order Runge-Kutta method.

    Args:
        K2_rate (numpy.ndarray): Rate matrix. Units: 1/fs.
        P1_0 (numpy.ndarray): Initial population vector.
        t_axis (numpy.ndarray or list): Time points at which to evaluate populations. Units: fs.

    Returns:
        list: List of population vectors at each time point in t_axis.
    """
    P1_i = P1_0
    P2_t_pop = [P1_0]
    t_old = t_axis[0]

    for t_i in t_axis[1:]:
        P1_i = runge_kutta_step(K2_rate, P1_i, t_i-t_old)
        P2_t_pop.append(P1_i)
        t_old = t_i

    return P2_t_pop


def runge_kutta_step(K2_rate, P1_pop, tau):
    """Perform a single 4th order Runge-Kutta integration step.

    Performs a single Runge-Kutta step from the current time to a time tau forward using
    the classical RK4 method with coefficients [0.0, 0.5, 0.5, 1.0].

    Args:
        K2_rate (numpy.ndarray): Rate matrix. Units: 1/fs.
        P1_pop (numpy.ndarray): Current population vector.
        tau (float): Timestep of the calculation. Units: fs.

    Returns:
        numpy.ndarray: Updated population vector after time tau.
    """
    # Calculation constants
    # ---------------------
    k = [[] for i in range(4)]
    c_rk = [0.0, 0.5, 0.5, 1.0]

    for i in range(4):
        # Update system values: phi_tmp
        if i == 0:
            P1_tmp = copy.deepcopy(P1_pop)
        else:
            P1_tmp = P1_pop + c_rk[i] * k[i - 1] * tau

        # Calculate system derivatives
        k[i] = K2_rate@P1_tmp

    # Actual Integration Step
    P1_pop = P1_pop + tau * (k[0] + 2.0 * k[1] + 2.0 * k[2] + k[3]) / 6.0

    return P1_pop


class RateDomain:
    """A domain representing pigments that are not electronically coupled.

    This class represents a domain of pigments connected through energy transfer rates
    rather than electronic coupling. It provides a simplified interface compatible with
    ExcitonDomain for rate matrix calculations.

    Attributes:
        t_axis (numpy.ndarray): Time axis from PigmentLineshape_Kubo class.
        domain_type (str): The type identifier of the domain.
        connected_to (list): List of [indices, rates] for outgoing connections.
        connected_from (list): List of [indices, rates] for incoming connections.
    """
    t_axis = PigmentLineshape_Kubo.t_axis

    def __init__(self, domain_type, connectivity, location=None):
        """Initialize a RateDomain.

        Args:
            domain_type (str): The type identifier of the domain (e.g., 'FL', 'ISC', 'RP1').
            connectivity (list): A nested list structure [[to_indices, to_rates], [from_indices, from_rates]]
                               specifying connections to and from this domain.
            location (numpy.ndarray, optional): The 3D spatial location of the domain. Defaults to None.
        """
        self.domain_type = domain_type
        self.connected_to = connectivity[0]
        self.connected_from = connectivity[1]
        self.__location = location
        self._explicit_vib = False

    def __repr__(self):
        """Return string representation of the RateDomain.

        Returns:
            str: String representation in the format 'Rate Domain (domain_type)'.
        """
        return f'Rate Domain ({self.name})'

    def __len__(self):
        """Return the number of states in the domain.

        Returns:
            int: Always returns 1 for RateDomain.
        """
        return 1

    @property
    def location(self):
        """Get the spatial location of the domain.

        Returns:
            numpy.ndarray or None: The 3D spatial location of the domain.
        """
        return self.__location

    @cached_property
    def name(self):
        """Get the name of the domain.

        Returns:
            str: The domain type as the name.
        """
        return f'{self.domain_type}'

    @property
    def list_names(self):
        """Get a list of domain names.

        Returns:
            list: A single-element list containing the domain type.
        """
        return [self.domain_type]

    def thermal_pop(self, temp):
        """Calculate thermal population.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            numpy.ndarray: Array containing [1] (no thermal distribution for rate domains).
        """
        return np.array([1])

    def thermal_mu2(self, temp):
        """Calculate thermal average of transition dipole moment squared.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            float: Always returns 0 for rate domains.
        """
        return 0

    def _raw_thermal_pop(self, temp):
        """Calculate raw thermal population.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            list: List containing [1].
        """
        return [1]

    def G2_gt_site(self, temp):
        """Calculate site lineshape function.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            numpy.ndarray: Zero array with shape of t_axis.
        """
        return np.zeros_like(self.t_axis)

    @property
    def e_lambda_exciton(self):
        """Get exciton reorganization energies.

        Returns:
            None: Rate domains have no reorganization energy.
        """
        return None

    def calculate_lifetime_by_temp(self, temp):
        """Calculate lifetime by temperature.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            list: List containing [0].
        """
        return [0]

    def abs_and_fl_exc_by_temp(self, temp):
        """Calculate absorption and fluorescence spectra by temperature.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            tuple: (t_axis, zero array, zero array) representing time axis and zero spectra.
        """
        return (self.t_axis, np.zeros_like(self.t_axis), np.zeros_like(self.t_axis))

    def calculate_abs_fl_exc(self, temp):
        """Calculate absorption and fluorescence spectra.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            tuple: (t_axis, zero array, zero array) representing time axis and zero spectra.
        """
        return (self.t_axis, np.zeros_like(self.t_axis), np.zeros_like(self.t_axis))

    def calculate_abs_vib_wdipole(self, temp):
        """Calculate absorption spectrum with vibrational structure weighted by dipole.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            numpy.ndarray: Zero array with shape of t_axis.
        """
        return np.zeros_like(self.t_axis)

    def calculate_fl_vib_wdipole(self, temp):
        """Calculate fluorescence spectrum with vibrational structure weighted by dipole.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            numpy.ndarray: Zero array with shape of t_axis.
        """
        return np.zeros_like(self.t_axis)

    def G2_gt_exciton(self, temp):
        """Calculate exciton lineshape function.

        Args:
            temp (float): Temperature. Units: K.

        Returns:
            numpy.ndarray: Zero array with shape of t_axis.
        """
        return np.zeros_like(self.t_axis)

    @property
    def n_exc(self):
        """Get the number of exciton states.

        Returns:
            int: Always returns 1 for RateDomain.
        """
        return 1

    @property
    def U2_exc_site(self):
        """Get transformation matrix from exciton to site basis.

        Returns:
            numpy.ndarray: Identity matrix of shape (1, 1).
        """
        return np.array([[1]])

    @property
    def U2_site_exc(self):
        """Get transformation matrix from site to exciton basis.

        Returns:
            numpy.ndarray: The stored transformation matrix.
        """
        return self._U2_site_exc

    @property
    def dipole_exciton(self):
        """Get exciton transition dipole moment.

        Returns:
            numpy.ndarray: Zero vector of shape (3,).
        """
        return np.zeros(3)

    @property
    def dipole_site(self):
        """Get site transition dipole moment.

        Returns:
            numpy.ndarray: Zero vector of shape (3,).
        """
        return np.zeros(3)

    @property
    def H2_domain(self):
        """Get domain Hamiltonian.

        Returns:
            numpy.ndarray: Zero matrix of shape (1, 1).
        """
        return np.array([[0]])

    @property
    def is_rc(self):
        """Check if this domain is a reaction center.

        Returns:
            bool: Always returns False for RateDomain.
        """
        return False


class RateMatrix(LinearSpectra):
    """Rate matrix for energy and charge transfer dynamics in photosynthetic systems.

    This class constructs and manages a rate matrix for excitation energy transfer (EET),
    charge separation, fluorescence, and intersystem crossing. It extends LinearSpectra to
    incorporate spectroscopic properties with kinetic modeling.

    Attributes:
        structure (ProteinAggregate): The protein aggregate structure.
        temp (float): Temperature in Kelvin.
        dict_protein_to_index (dict): Maps protein names to domain indices.
    """

    def __init__(self,
                 structure: ProteinAggregate,
                 temp,
                 r_neighbor=70,
                 k_ij=None,
                 ):
        super().__init__(structure)

        if k_ij is None:
            self.__protein_rate_matrix = self.construct_rate_matrix(temp, r_neighbor)
        else:
            self.__protein_rate_matrix = k_ij

        self.__rate_matrix = self.protein_rate_matrix
        self.__temp = temp
        self.__list_index_fl = []
        self.__list_index_rp2 = []
        self.__list_index_isc = []

    def construct_rate_matrix(self, temp, neighbor_dist=70, use_progress=True):
        list_don = []
        list_acc = []
        list_rates = []

        # Create progress bar for donor domains if available
        has_tqdm = True  # tqdm is imported at module level
        domain_iter = tqdm(enumerate(self.list_domain), total=len(self.list_domain),
                          desc="Computing rates (serial)", unit="domain") if (has_tqdm and use_progress) else enumerate(self.list_domain)

        for (index_don, domain_donor) in domain_iter:
            list_acceptor_index = self.construct_neighbor_list(index_don, neighbor_dist)
            for index_acc in list_acceptor_index:
                list_don.append(index_don)
                list_acc.append(index_acc)
                domain_acc = self.list_domain[index_acc]
                list_rates.append(generalized_forster_thermal(domain_acc,
                                                              domain_donor,
                                                              temp,
                                                              self.H2_hamiltonian
                                                              )
                                  )
        return coo_array((list_rates, (list_acc, list_don)),
                         shape=(self.n_dom, self.n_dom))

    def construct_neighbor_list(self, index_donor, r_cutoff):
        list_loc = np.array([dom.location for dom in self.list_domain])
        list_dist = np.linalg.norm(list_loc - self.list_domain[index_donor].location, axis=1)
        list_close = list(np.where(list_dist<r_cutoff)[0])
        list_close.remove(index_donor)
        return list_close

    @property
    def temp(self):
        """Get the temperature.

        Returns:
            float: Temperature. Units: K.
        """
        return self.__temp

    @property
    def protein_rate_matrix(self):
        """Get the protein rate matrix.

        Returns:
            scipy.sparse.coo_array: The protein rate matrix.
        """
        return self.__protein_rate_matrix

    def time_evolve(self, P1_0, N_t, dt):
        """Time evolve the population using numerical propagation.

        Args:
            P1_0 (numpy.ndarray): Initial population vector.
            N_t (int): Number of time steps.
            dt (float): Time step size. Units: fs.

        Returns:
            tuple: (t_axis, P1_t) where t_axis is the time array and P1_t is the list of populations.
        """
        t_axis = np.arange(0, N_t*dt, dt)
        P1_t = numericalprop(self.rate_matrix, P1_0, t_axis)
        return t_axis, P1_t

    def _add_domain(self, new_domain):
        """Add a new domain to the rate matrix.

        Args:
            new_domain (RateDomain): The domain to add with connectivity information.
        """
        self.list_domain.append(new_domain)
        index_new_domain = len(self.list_domain)-1
        row = list(self.__rate_matrix.row)
        col = list(self.__rate_matrix.col)
        data = list(self.__rate_matrix.data)
        for (index_to, rate) in zip(new_domain.connected_to[0],
                                    new_domain.connected_to[1]):
            row.append(index_to)
            col.append(index_new_domain)
            data.append(rate)

        for (index_from, rate) in zip(new_domain.connected_from[0],
                                    new_domain.connected_from[1]):
            row.append(index_new_domain)
            col.append(index_from)
            data.append(rate)

        self.__rate_matrix = coo_array((data, (row, col)),
                                       shape=(self.n_dom,self.n_dom))

    def _correct_diagonal(self):
        """Correct diagonal elements to ensure column sums equal zero.

        This function checks that the sum of rates out of all states (as controlled by the 
        negative diagonal entry) are equal to 0. If they are not, then the diagonal term is 
        updated to ensure the sum of columns is 0.
        """
        row = list(self.__rate_matrix.row)
        col = list(self.__rate_matrix.col)
        data = list(self.__rate_matrix.data)
        diag_correction = -np.sum(self.__rate_matrix, 0)
        for index_col in np.arange(self.__rate_matrix.shape[1]):
            row.append(index_col)
            col.append(index_col)
            data.append(diag_correction[index_col])

        self.__rate_matrix = coo_array((data, (row, col)),
                                       shape=(self.n_dom,self.n_dom))

    def add_isc(self, k_isc):
        """Add intersystem crossing (ISC) pathway to all exciton domains.

        Creates a sink domain for intersystem crossing that is connected to all ExcitonDomain 
        objects with the specified rate constant.

        Args:
            k_isc (float): Intersystem crossing rate constant. Units: 1/fs.
        """
        list_rates_from = []
        list_connections_from = []
        for (index, domain) in enumerate(self.list_domain):
            if type(domain) == ExcitonDomain:
                list_rates_from.append(k_isc)
                list_connections_from.append(index)

        isc_domain = RateDomain(domain_type='ISC',
                                connectivity=[[[], []],
                                              [list_connections_from, list_rates_from]])
        self.__list_index_isc.append(len(self.list_domain))
        self._add_domain(isc_domain)

    def add_cs_rp2(self, k_cs, k_rc, k_irr):
        '''
        The key trick here is that we need to know how to identify the RC proteins
        and the specific domain that is responsible for charge separation.
        '''
        for (index_rc, domain) in enumerate(self.list_domain):
            # Check if the protein is the RC
            if domain.is_rc:
                index_rp1 = len(self.list_domain)
                rp1_domain = RateDomain(domain_type='RP1',
                                        connectivity=[[[index_rc], [k_rc]],
                                                      [[index_rc], [k_cs]]])
                self._add_domain(rp1_domain)
                rp2_domain = RateDomain(domain_type='RP2',
                                        connectivity=[[[], []],
                                                      [[index_rp1], [k_irr]]])
                self.__list_index_rp2.append(len(self.list_domain))
                self._add_domain(rp2_domain)

    def add_fluorescence(self, k_fl):
        """Add fluorescence decay pathway to all exciton domains.

        Creates a fluorescence sink domain connected to all ExcitonDomain objects. The rates 
        are weighted by the thermal average of transition dipole moment squared to account for 
        temperature-dependent emission.

        Args:
            k_fl (float): Base fluorescence rate constant. Units: 1/fs.
        """
        list_therm_mu2_from = []
        list_index_fl_from = []
        for (index, domain) in enumerate(self.list_domain):
            if type(domain) == ExcitonDomain:
                list_therm_mu2_from.append(domain.thermal_mu2(self.temp))
                list_index_fl_from.append(index)

        list_therm_rate_from = k_fl*np.array(list_therm_mu2_from)/np.mean(list_therm_mu2_from)
        fl_domain = RateDomain(domain_type='FL',
                               connectivity=[[[], []],
                                             [list_index_fl_from, list_therm_rate_from]])
        self.__list_index_fl.append(len(self.list_domain))
        self._add_domain(fl_domain)

    def add_rp1_to_ct(self, ct_pigment_name, k_ct, k_rp1):
        """Add an RP1 domain coupled to a charge transfer (CT) state with reversible rates.

        Creates a radical pair domain that is reversibly connected to a specific CT pigment 
        state. This is used to model charge separation from CT states.

        Args:
            ct_pigment_name (str): Name of the CT pigment to identify its domain.
            k_ct (float): Forward rate from CT to RP1 (charge separation). Units: 1/fs.
            k_rp1 (float): Backward rate from RP1 to CT (charge recombination). Units: 1/fs.

        Raises:
            ValueError: If the CT pigment name is not found in any domain.
        """
        # Find the CT domain
        ct_index = None
        for i, domain in enumerate(self.list_domain):
            if hasattr(domain, 'list_pigments'):
                pigment_names = [p.name for p in domain.list_pigments]
                if ct_pigment_name in pigment_names:
                    ct_index = i
                    break
        
        if ct_index is None:
            raise ValueError(f"CT pigment {ct_pigment_name} not found in any domain.")
        
        # Create and add the RP1 domain        
        rp1_domain = RateDomain(domain_type='CTRP1',
                                connectivity=[[[ct_index], [k_rp1]],
                                              [[ct_index], [k_ct]]])

        self.__list_index_ct_rp1.append(len(self.list_domain))
        self._add_domain(rp1_domain)

    def calculate_thermal_rate(self, list_acceptor_index, list_donor_index):
        """Calculate thermally averaged rate from donor set to acceptor set.

        This method calculates the rate of transport from a set of donor domains to a set of 
        acceptor domains assuming that excitation is initially thermalized within the set of 
        donor domains.

        Args:
            list_acceptor_index (list): List of acceptor domain indices.
            list_donor_index (list): List of donor domain indices.

        Returns:
            float: Effective thermal rate from donors to acceptors. Units: 1/fs.

        Note:
            TODO: Is there a better way to slice a sparse array without casting it to dense?
        """
        thermal_pop = np.array([np.sum(self.list_domain[index_donor]._raw_thermal_pop(self.temp)) for index_donor in list_donor_index])
        thermal_pop = thermal_pop / np.sum(thermal_pop)  # normalized P_d
        rate_eff = np.sum([P_d * np.sum([self.rate_matrix.tocsr()[index_acceptor, index_donor]
                                         for index_acceptor in list_acceptor_index])
                           for (P_d, index_donor) in zip(thermal_pop, list_donor_index)])
        return rate_eff

    @property
    def list_index_fl(self):
        """Get indices of fluorescence domains.

        Returns:
            list: List of indices for fluorescence sink domains.
        """
        return self.__list_index_fl

    @property
    def list_index_isc(self):
        """Get indices of intersystem crossing domains.

        Returns:
            list: List of indices for ISC sink domains.
        """
        return self.__list_index_isc

    @property
    def list_index_rp2(self):
        """Get indices of radical pair 2 (RP2) domains.

        Returns:
            list: List of indices for RP2 sink domains.
        """
        return self.__list_index_rp2

    @property
    def list_index_ct_rp1(self):
        """Get indices of CT-coupled RP1 domains.

        Returns:
            list: List of indices for CT-RP1 domains.
        """
        return self.__list_index_ct_rp1

    @property
    def list_index_rate_domains(self):
        """Get indices of all rate sink domains.

        Returns:
            list: Combined list of indices for all rate domains (RP2, ISC, and FL).
        """
        index_rate_domains = self.__list_index_rp2 + self.__list_index_isc + self.__list_index_fl
        return index_rate_domains
        
    @property
    def rate_matrix(self):
        """Get the full rate matrix including all added domains.

        Returns:
            scipy.sparse.coo_array: The complete rate matrix.
        """
        return self.__rate_matrix

    @property
    def H2_hamiltonian(self):
        """Get the Hamiltonian from the structure.

        Returns:
            numpy.ndarray: The Hamiltonian matrix.
        """
        return self.structure.H2_hamiltonian

    def plot_protein_connectivity(
        self,
        dict_chains_to_domains: Dict[str, List[int]],
        dict_chain_to_protein: Dict[str, str],
        dict_protein_to_legend: Dict[str, str],
        threshold: float = 100,
        ax: plt.Axes = None
    ):
        """Plot a protein connectivity map showing EET timescales between proteins.

        Creates a 2D visualization of protein complexes with arrows indicating energy transfer
        pathways and their characteristic timescales. Proteins are shown as filled convex hulls
        with labeled connections between them.

        Args:
            dict_chains_to_domains (Dict[str, List[int]]): Chain-to-domain map that groups domains
                                                           by the chain they're located in.
            dict_chain_to_protein (Dict[str, str]): Chain-to-protein map that labels each PDB chain
                                                    with a protein name.
            dict_protein_to_legend (Dict[str, str]): Protein-to-color map specifying the color for
                                                     each protein in the plot.
            threshold (float, optional): The maximum transfer time to display in the map.
                                        Units: ps. Defaults to 100.
            ax (plt.Axes, optional): Matplotlib axes object to plot on. If None, creates new figure.
                                    Defaults to None.

        Returns:
            plt.Axes: Matplotlib axes object with the connectivity map.
        """

        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        chain_centers = {}

        for chain, _ in dict_chains_to_domains.items():
            coords = []

            # Check if structure has get_chains method (ProteinAggregate)
            if hasattr(self.structure, 'get_chains'):
                for pdb_chain in self.structure.get_chains():
                    if pdb_chain.id == chain:
                        coords.extend([atom.coord for atom in pdb_chain.get_atoms()])

            # If structure is a Membrane with list_proteins, use transformed CG coordinates
            elif hasattr(self.structure, 'list_proteins'):
                # For membranes, we need to match by protein name + chain to avoid confusion
                # between chains with same ID in different proteins (e.g., LHCII in C2S2 vs standalone LHCII)

                # Parse the chain identifier to see if it includes protein name
                # Expected format could be: "proteinname_chain" or just "chain"
                if '_' in chain:
                    # Format: "proteinname_chain" - match specific protein
                    target_protein_name, target_chain_id = chain.rsplit('_', 1)

                    for protein in self.structure.list_proteins:
                        if protein.name == target_protein_name:
                            if hasattr(protein, 'atomic') and hasattr(protein.atomic, 'get_chains'):
                                for pdb_chain in protein.atomic.get_chains():
                                    if pdb_chain.id == target_chain_id:
                                        # Get indices of atoms in this specific chain
                                        chain_atom_indices = []
                                        all_atoms = list(protein.atomic.get_atoms())
                                        for idx, atom in enumerate(all_atoms):
                                            if atom.get_parent().get_parent().id == target_chain_id:
                                                chain_atom_indices.append(idx)

                                        # Get only the transformed coordinates for this chain
                                        if chain_atom_indices:
                                            coords.extend(protein.R2_atomic_xyz[chain_atom_indices])
                                        break
                            break
                else:
                    # Legacy format: just "chain" - find first matching chain across all proteins
                    for protein in self.structure.list_proteins:
                        if hasattr(protein, 'atomic') and hasattr(protein.atomic, 'get_chains'):
                            for pdb_chain in protein.atomic.get_chains():
                                if pdb_chain.id == chain:
                                    # Get indices of atoms in this specific chain
                                    chain_atom_indices = []
                                    all_atoms = list(protein.atomic.get_atoms())
                                    for idx, atom in enumerate(all_atoms):
                                        if atom.get_parent().get_parent().id == chain:
                                            chain_atom_indices.append(idx)

                                    # Get only the transformed coordinates for this chain
                                    if chain_atom_indices:
                                        coords.extend(protein.R2_atomic_xyz[chain_atom_indices])
                                    break

            # Fallback: use domain pigment coordinates (which are already transformed for CG)
            if not coords:
                for index in dict_chains_to_domains[chain]:
                    pigments = self.list_domain[index].list_pigments
                    for pigment in pigments:
                        # Handle both CGPigment (has .atomic.residue) and atomic pigments (has .residue)
                        if hasattr(pigment, 'atomic') and hasattr(pigment.atomic, 'residue'):
                            # For CGPigment, use R2_atomic_xyz which are transformed coordinates
                            coords.extend(pigment.R2_atomic_xyz)
                        elif hasattr(pigment, 'residue'):
                            coords.extend([atom.coord for atom in pigment.residue.get_atoms()])
                        elif hasattr(pigment, 'R2_atomic_xyz'):
                            # If all else fails, use the mapped atomic coordinates directly
                            coords.extend(pigment.R2_atomic_xyz)

            coords = np.array(coords)
            hull = concave_hull(MultiPoint(coords), ratio=0.20)
            x, y = hull.exterior.xy

            # Get protein name for this chain
            protein_name = dict_chain_to_protein.get(chain, 'unknown')

            # Extract chain ID for labeling (last part after final underscore)
            if '_' in chain:
                chain_id = chain.rsplit('_', 1)[1]  # Get the actual chain letter
                display_label = chain_id  # Just show the chain letter
            else:
                chain_id = chain
                display_label = protein_name.upper().replace('_', ' ')

            # Try to get color by exact match first, then try matching by protein type
            color = dict_protein_to_legend.get(protein_name, None)

            # If no direct match, try to extract base protein type (e.g., 'lhcii' from 'lhcii_0')
            if color is None:
                if '_' in protein_name:
                    base_protein_type = protein_name.split('_')[0]
                    # Try uppercase version first
                    color = dict_protein_to_legend.get(base_protein_type.upper(), None)
                    # Try lowercase if uppercase didn't work
                    if color is None:
                        color = dict_protein_to_legend.get(base_protein_type, None)

            # Default to grey if still no match
            if color is None:
                print(f"Warning: No color found for protein '{protein_name}'")
                if '_' in protein_name:
                    print(f"  Base type: '{protein_name.split('_')[0]}'")
                print(f"  Available keys: {list(dict_protein_to_legend.keys())}")
                color = '#999999'

            ax.plot(x, y, linewidth=1, color=color)
            ax.fill(x, y, alpha=0.5, color=color)

            center = hull.centroid
            chain_centers[chain] = (center.x, center.y)
            ax.text(center.x, center.y, display_label, ha='center',
                    va='center', fontweight='bold', color='black', fontsize=10)

        chains = list(dict_chains_to_domains.keys())
        for i, chain1 in enumerate(chains):
            for j, chain2 in enumerate(chains):
                if i < j:
                    rate_forward = self.calculate_thermal_rate(dict_chains_to_domains[chain2],
                                                                    dict_chains_to_domains[chain1])
                    rate_reverse = self.calculate_thermal_rate(dict_chains_to_domains[chain1],
                                                                    dict_chains_to_domains[chain2])

                    # Convert rates from fs^-1 to ps^-1
                    rate_forward_ps = rate_forward * 1000
                    rate_reverse_ps = rate_reverse * 1000

                    time_forward = 1 / rate_forward_ps if rate_forward_ps > 0 else float('inf')
                    time_reverse = 1 / rate_reverse_ps if rate_reverse_ps > 0 else float('inf')

                    if time_forward < threshold or time_reverse < threshold:
                        x1, y1 = chain_centers[chain1]
                        x2, y2 = chain_centers[chain2]

                        angle = np.arctan2(y2 - y1, x2 - x1)
                        offset = 4

                        if time_forward < threshold:
                            dx = offset * np.sin(angle)
                            dy = -offset * np.cos(angle)
                            arrow = ConnectionPatch((x1 + dx, y1 + dy), (x2 + dx, y2 + dy), 'data', 'data',
                                                    arrowstyle='->', connectionstyle='arc3,rad=0.1',
                                                    color='#000000', alpha=0.7, linewidth=1.5)
                            ax.add_artist(arrow)

                            mid_x = (x1 + x2) / 2 + 1.5 * dx
                            mid_y = (y1 + y2) / 2 + 1.5 * dy
                            label_angle = np.degrees(angle)
                            if 90 < label_angle < 270:
                                label_angle += 180
                            ax.text(mid_x, mid_y, f'{time_forward:.0f} ps', ha='center', va='center',
                                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1),
                                    rotation=label_angle, fontsize=8, color='black')

                        if time_reverse < threshold:
                            dx = -offset * np.sin(angle)
                            dy = offset * np.cos(angle)
                            arrow = ConnectionPatch((x2 + dx, y2 + dy), (x1 + dx, y1 + dy), 'data', 'data',
                                                    arrowstyle='->', connectionstyle='arc3,rad=-0.1',
                                                    color='#000000', alpha=0.7, linewidth=1.5)
                            ax.add_artist(arrow)

                            mid_x = (x1 + x2) / 2 + 1.5 * dx
                            mid_y = (y1 + y2) / 2 + 1.5 * dy
                            label_angle = np.degrees(angle + np.pi)
                            if 90 < label_angle < 270:
                                label_angle += 180
                            ax.text(mid_x, mid_y, f'{time_reverse:.0f} ps', ha='center', va='center',
                                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1),
                                    rotation=label_angle, fontsize=8, color='black')

        return ax
    def plot_protein_rates(
            self,
            inter_protein_cutoff=40,
            intra_timescale_range=(0.1, 20),
            unit='ps',
            cmap='viridis_r',
            hull_alpha=0.3,
            show_intra_rates=True,
            show_inter_rates=True,
            inter_arrow_style=True,
            figsize=(14, 10),
            ax=None
    ):
        """Plot rates with clear separation of intra- and inter-protein transfers.

        Shows each protein as a colored hull with internal domain rates displayed
        inside, and inter-protein transfers shown as labeled arrows between proteins.

        Parameters
        ----------
        inter_protein_cutoff : float, optional
            Maximum timescale (in unit) to show for inter-protein transfers.
            Default is 40 ps.
        intra_timescale_range : tuple, optional
            (min, max) timescale range for coloring intra-protein rates.
            Default is (0.1, 20) ps.
        unit : str, optional
            Time unit. Options: 'ps' (default), 'fs', 'ns'.
        cmap : str, optional
            Colormap for intra-protein rate lines. Default is 'viridis_r'.
        hull_alpha : float, optional
            Transparency of protein hull fill. Default is 0.3.
        show_intra_rates : bool, optional
            Whether to show intra-protein rate lines. Default is True.
        show_inter_rates : bool, optional
            Whether to show inter-protein transfer arrows. Default is True.
        inter_arrow_style : bool, optional
            If True, show inter-protein rates as arrows with labels.
            If False, show as colored lines. Default is True.
        figsize : tuple, optional
            Figure size. Default is (14, 10).
        ax : plt.Axes, optional
            Matplotlib axes. If None, creates new figure.

        Returns
        -------
        fig : matplotlib.figure.Figure
            The figure object.
        ax : plt.Axes
            The axes object.
        """
        from matplotlib.collections import LineCollection
        from matplotlib.colors import LogNorm
        from matplotlib.patches import FancyArrowPatch

        # Create figure
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.get_figure()

        # Get domain locations
        locations = []
        for dom in self.list_domain:
            loc = dom.location
            if loc is not None:
                locations.append([float(loc[0]), float(loc[1])])
            else:
                locations.append([0.0, 0.0])
        locations = np.array(locations)

        # Group domains by protein
        protein_domains = {}  # protein_name -> list of domain indices
        domain_to_protein = {}  # domain_index -> protein_name

        for i, dom in enumerate(self.list_domain):
            if hasattr(dom, 'list_pigments') and len(dom.list_pigments) > 0:
                pig_name = dom.list_pigments[0].name
                parts = pig_name.split('_')
                if len(parts) >= 2:
                    protein_name = f"{parts[0]}_{parts[1]}"
                else:
                    protein_name = parts[0]
            else:
                protein_name = f'protein_{i}'

            if protein_name not in protein_domains:
                protein_domains[protein_name] = []
            protein_domains[protein_name].append(i)
            domain_to_protein[i] = protein_name

        # Get rate matrix and convert to timescales
        rate_dense = self.rate_matrix.toarray()
        with np.errstate(divide='ignore', invalid='ignore'):
            if unit == 'ps':
                timescales = np.where(rate_dense > 0, 1.0 / (rate_dense * 1000), np.inf)
            elif unit == 'fs':
                timescales = np.where(rate_dense > 0, 1.0 / rate_dense, np.inf)
            elif unit == 'ns':
                timescales = np.where(rate_dense > 0, 1.0 / (rate_dense * 1e6), np.inf)
            else:
                raise ValueError(f"Unknown unit '{unit}'. Use 'ps', 'fs', or 'ns'.")

        # Setup colormap for proteins
        protein_names = list(protein_domains.keys())
        n_proteins = len(protein_names)
        protein_cmap = plt.cm.get_cmap('tab20', max(n_proteins, 1))
        protein_colors = {name: protein_cmap(i) for i, name in enumerate(protein_names)}

        # Draw protein hulls and collect centroids
        protein_centroids = {}

        for prot_idx, (protein_name, domain_indices) in enumerate(protein_domains.items()):
            # Collect coordinates for hull
            coords = []
            for dom_idx in domain_indices:
                dom = self.list_domain[dom_idx]
                if hasattr(dom, 'list_pigments'):
                    for pig in dom.list_pigments:
                        if hasattr(pig, 'residue'):
                            coords.extend([atom.coord[:2] for atom in pig.residue.get_atoms()])
                        elif hasattr(pig, 'location'):
                            coords.append(pig.location[:2])

            if len(coords) < 3:
                coords = [locations[i] for i in domain_indices]

            if len(coords) >= 3:
                coords = np.array(coords)
                try:
                    hull = concave_hull(MultiPoint(coords), ratio=0.3)
                    x, y = hull.exterior.xy
                    color = protein_colors[protein_name]
                    ax.fill(x, y, alpha=hull_alpha, color=color, zorder=1)
                    ax.plot(x, y, linewidth=2, color=color, zorder=2)

                    cx, cy = hull.centroid.x, hull.centroid.y
                    protein_centroids[protein_name] = (cx, cy)

                    # Add protein label
                    ax.text(cx, cy, protein_name, ha='center', va='center',
                           fontsize=10, fontweight='bold', color='black',
                           bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=2),
                           zorder=15)
                except Exception:
                    # Fallback: use mean of domain locations
                    cx = np.mean([locations[i][0] for i in domain_indices])
                    cy = np.mean([locations[i][1] for i in domain_indices])
                    protein_centroids[protein_name] = (cx, cy)
            else:
                cx = np.mean([locations[i][0] for i in domain_indices])
                cy = np.mean([locations[i][1] for i in domain_indices])
                protein_centroids[protein_name] = (cx, cy)

        # Initialize lists for counting
        intra_segments = []
        inter_protein_rates = {}

        # Draw intra-protein rates
        if show_intra_rates:
            intra_timescales_list = []

            n_dom = len(self.list_domain)
            for i in range(n_dom):
                for j in range(i + 1, n_dom):
                    # Check if same protein
                    if domain_to_protein[i] == domain_to_protein[j]:
                        t_ij = timescales[i, j]
                        t_ji = timescales[j, i]
                        t_min = min(t_ij, t_ji)

                        if t_min < np.inf and intra_timescale_range[0] <= t_min <= intra_timescale_range[1]:
                            x1, y1 = locations[i]
                            x2, y2 = locations[j]
                            intra_segments.append([(x1, y1), (x2, y2)])
                            intra_timescales_list.append(t_min)

            if intra_segments:
                intra_timescales_arr = np.array(intra_timescales_list)
                norm = LogNorm(vmin=intra_timescale_range[0], vmax=intra_timescale_range[1])

                # Linewidths: faster = thicker
                lw_min, lw_max = 0.5, 2.5
                log_t = np.log10(intra_timescales_arr)
                log_range = np.log10(intra_timescale_range[1]) - np.log10(intra_timescale_range[0])
                log_t_norm = (np.log10(intra_timescale_range[1]) - log_t) / log_range
                linewidths = lw_min + (lw_max - lw_min) * log_t_norm

                lc = LineCollection(intra_segments, cmap=cmap, norm=norm, alpha=0.7,
                                   linewidths=linewidths, zorder=3)
                lc.set_array(intra_timescales_arr)
                ax.add_collection(lc)

                # Add colorbar for intra-protein rates
                cbar = fig.colorbar(lc, ax=ax, label=f'Intra-protein timescale ({unit})',
                                   shrink=0.6, pad=0.02)

            # Plot domain centers
            for protein_name, domain_indices in protein_domains.items():
                dom_locs = np.array([locations[i] for i in domain_indices])
                color = protein_colors[protein_name]
                ax.scatter(dom_locs[:, 0], dom_locs[:, 1], c=[color], s=40,
                          edgecolors='black', linewidths=0.5, zorder=5)

        # Draw inter-protein rates
        if show_inter_rates:
            # Calculate thermal rates between protein pairs
            protein_names_list = list(protein_domains.keys())
            for i, prot1 in enumerate(protein_names_list):
                for j, prot2 in enumerate(protein_names_list):
                    if i < j:  # Only compute once per pair
                        # Get domain indices for each protein
                        domains_prot1 = protein_domains[prot1]
                        domains_prot2 = protein_domains[prot2]

                        # Calculate thermal rates in both directions
                        rate_forward = self.calculate_thermal_rate(domains_prot2, domains_prot1)
                        rate_reverse = self.calculate_thermal_rate(domains_prot1, domains_prot2)

                        # Convert rates from fs^-1 to timescales in requested unit
                        if unit == 'ps':
                            rate_forward_ps = rate_forward * 1000
                            rate_reverse_ps = rate_reverse * 1000
                        elif unit == 'fs':
                            rate_forward_ps = rate_forward
                            rate_reverse_ps = rate_reverse
                        elif unit == 'ns':
                            rate_forward_ps = rate_forward * 1e6
                            rate_reverse_ps = rate_reverse * 1e6

                        time_forward = 1 / rate_forward_ps if rate_forward_ps > 0 else float('inf')
                        time_reverse = 1 / rate_reverse_ps if rate_reverse_ps > 0 else float('inf')

                        # Store rates if below cutoff
                        if time_forward < inter_protein_cutoff:
                            inter_protein_rates[(prot1, prot2)] = time_forward
                        if time_reverse < inter_protein_cutoff:
                            inter_protein_rates[(prot2, prot1)] = time_reverse

            # Draw inter-protein arrows
            for (prot_from, prot_to), timescale in inter_protein_rates.items():
                if prot_from in protein_centroids and prot_to in protein_centroids:
                    x1, y1 = protein_centroids[prot_from]
                    x2, y2 = protein_centroids[prot_to]

                    # Offset arrows slightly to avoid overlap for bidirectional
                    angle = np.arctan2(y2 - y1, x2 - x1)
                    offset = 5
                    dx = offset * np.sin(angle)
                    dy = -offset * np.cos(angle)

                    if inter_arrow_style:
                        # Draw as arrow with label
                        arrow = FancyArrowPatch(
                            (x1 + dx, y1 + dy), (x2 + dx, y2 + dy),
                            arrowstyle='-|>', mutation_scale=15,
                            color='red', linewidth=2, alpha=0.8, zorder=12
                        )
                        ax.add_patch(arrow)

                        # Add timescale label
                        mid_x = (x1 + x2) / 2 + 2 * dx
                        mid_y = (y1 + y2) / 2 + 2 * dy
                        ax.text(mid_x, mid_y, f'{timescale:.1f} {unit}',
                               fontsize=9, fontweight='bold', color='red',
                               bbox=dict(facecolor='white', edgecolor='red',
                                        alpha=0.9, pad=2, boxstyle='round'),
                               ha='center', va='center', zorder=13)
                    else:
                        # Draw as colored line
                        ax.plot([x1 + dx, x2 + dx], [y1 + dy, y2 + dy],
                               'r-', linewidth=2, alpha=0.8, zorder=12)

        ax.set_aspect('equal')
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')

        # Title with summary
        n_intra = len(intra_segments) if show_intra_rates else 0
        n_inter = len(inter_protein_rates) if show_inter_rates else 0
        ax.set_title(f'Protein Rate Map\n{n_proteins} proteins, {n_intra} intra-rates, {n_inter} inter-protein transfers (< {inter_protein_cutoff} {unit})')

        # Auto-set limits
        margin = 30
        ax.set_xlim(locations[:, 0].min() - margin, locations[:, 0].max() + margin)
        ax.set_ylim(locations[:, 1].min() - margin, locations[:, 1].max() + margin)

        plt.tight_layout()
        return fig, ax

# TODO: Implement Runge-Kutta 4th order
# TODO: Build classes that contains rates to simulate observables
