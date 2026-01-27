import copy

import networkx as nx
import numpy as np
import pandas as pd


class Hamiltonian:
    """Manages the electronic Hamiltonian indexed by pigment site names.
    
    This class encapsulates the electronic Hamiltonian matrix for a system of pigments,
    providing methods to manipulate, analyze, and extract information from the Hamiltonian.
    The core data structure is a numpy array with a bidirectional mapping between pigment
    site names and matrix indices.
    
    Attributes:
        list_site_names (list): List of pigment site names in the Hamiltonian.
        dict_index_by_label (dict): Mapping from site names to matrix indices.
        H2_ham_0 (np.ndarray): Base Hamiltonian matrix without disorder.
        h2_ham (np.ndarray): Current Hamiltonian matrix including any added disorder.
        central_frequency (float): Central frequency reference for the Hamiltonian.
    
    Note:
        TODO: This class should own the 'central_frequency' idea for better encapsulation.
    """

    def __init__(self, h2_ham, list_site_names, list_sigma=None):
        """Initialize the Hamiltonian object.

        Args:
            h2_ham (array-like): Initial Hamiltonian matrix as a 2D array.
            list_site_names (list): List of site names corresponding to matrix indices.
            list_sigma (list, optional): List of disorder standard deviations for each site.
                Defaults to None.
        """
        self.__h2_ham_0 = np.array(h2_ham)
        self.__h2_disord = np.zeros_like(h2_ham)
        self.__dict_index_by_label = {name: index for (index, name) in enumerate(list_site_names)}
        self.__list_site_names = list_site_names
        self._list_sigma = list_sigma
        # Central frequency shift is tracked separately so the raw Hamiltonian stays intact
        self._e_0 = 0

        #  Disjoint Set Union Properties for finding domains
        self._domains_parnet = np.arange(len(list_site_names))
        self._domains_size = np.ones(len(list_site_names), dtype=int)

    def __array__(self):
        """Return the Hamiltonian matrix as a numpy array.
        
        Returns:
            np.ndarray: The current Hamiltonian matrix including disorder.
        """
        return self.h2_ham

    def _make_domain(self, index):
        """Initialize a new domain for disjoint set union.
        
        Args:
            index (int): Index of the site to make into a new domain.
        """
        self._domains_parent[index] = index
        self._domains_size[index] = 1

    def _join_domains(self, index_a, index_b):
        """Join two domains in the disjoint set union structure.
        
        Uses union by size heuristic to maintain balanced tree structure.
        
        Args:
            index_a (int): Index of a site in the first domain.
            index_b (int): Index of a site in the second domain.
        """
        a = self._find_domain(index_a)
        b = self._find_domain(index_b)

        if a != b:
            if self._domains_size[a] < self._domains_size[b]:
                self._domains_parent[a] = b
                self._domains_size[b] += self._domains_size[a]
            else:
                self._domains_parent[b] = a
                self._domains_size[a] += self._domains_size[b]

    def _find_domain(self, index):
        """Find the root domain of a given site with path compression.
        
        Args:
            index (int): Index of the site to find the domain for.
            
        Returns:
            int: Index of the root domain.
        """
        if self._domains_parent[index] != index:
            self._domains_parent[index] = self._find_domain(self._domains_parent[index])
        return self._domains_parent[index]

    def get_domains(self):
        """Get all connected domains in the Hamiltonian.
        
        Returns:
            list: List of lists, where each inner list contains the indices of sites
                in a connected domain.
        """
        domains = {}
        for i in range(len(self._domains_parent)):
            root = self._find_domain(i)
            if root not in domains:
                domains[root] = []
            domains[root].append(i)
        return list(domains.values())

    @property
    def H2_ham_0(self):
        """Get the base Hamiltonian matrix without disorder.
        
        Returns:
            np.ndarray: The base Hamiltonian matrix.
        """
        return self.__h2_ham_0

    def update_subset(self, h2_subset, list_row_name, list_col_name):
        """Update a subset of the Hamiltonian matrix.
        
        Args:
            h2_subset (np.ndarray): New values for the subset.
            list_row_name (list): List of row site names for the subset.
            list_col_name (list): List of column site names for the subset.
        """
        self.__h2_ham_0[np.ix_([self.__dict_index_by_label[name] for name in list_row_name],
                               [self.__dict_index_by_label[name] for name in list_col_name])] = h2_subset

    def update_coupling(self, v_ab, name_a, name_b):
        """Update the coupling between two sites.
        
        Updates both the (a,b) and (b,a) elements to maintain hermiticity.
        
        Args:
            v_ab (float): New coupling value.
            name_a (str): Name of the first site.
            name_b (str): Name of the second site.
        """
        self.__h2_ham_0[self.__dict_index_by_label[name_a],
                        self.__dict_index_by_label[name_b]] = v_ab
        self.__h2_ham_0[self.__dict_index_by_label[name_b],
                        self.__dict_index_by_label[name_a]] = v_ab

    def subset_by_index(self, list_row_index, list_col_index):
        """Extract a subset of the Hamiltonian by matrix indices.
        
        Args:
            list_row_index (list): List of row indices to extract.
            list_col_index (list): List of column indices to extract.
            
        Returns:
            np.ndarray: Subset of the Hamiltonian matrix.
        """
        return self.h2_ham[np.ix_(list_row_index, list_col_index)]

    def subset_by_name(self, list_row_name, list_col_name):
        """Extract a subset of the Hamiltonian by site names.
        
        Args:
            list_row_name (list): List of row site names to extract.
            list_col_name (list): List of column site names to extract.
            
        Returns:
            np.ndarray: Subset of the Hamiltonian matrix.
        """
        return self.h2_ham[np.ix_([self.__dict_index_by_label[name] for name in list_row_name],
                                  [self.__dict_index_by_label[name] for name in list_col_name])]

    @property
    def list_site_names(self):
        """Get the list of site names in the Hamiltonian.
        
        Returns:
            list: List of pigment site names.
        """
        return self.__list_site_names

    @property
    def dict_index_by_label(self):
        """Get the mapping from site names to matrix indices.
        
        Returns:
            dict: Dictionary mapping site names to their corresponding indices.
        """
        return self.__dict_index_by_label

    def add_disorder(self, seed=None):
        """Add diagonal disorder to the Hamiltonian.
        
        Generates random disorder based on the standard deviations in list_sigma.
        
        Args:
            seed (int, optional): Random seed for reproducibility. Defaults to None.
        """
        rng = np.random.default_rng(seed=seed)
        self.__h2_disord = np.diag(rng.normal(0, self._list_sigma))

    @property
    def h2_ham(self):
        """Get the current Hamiltonian matrix including disorder.
        
        Returns:
            np.ndarray: Base Hamiltonian shifted by the central frequency and
                combined with any diagonal disorder.
        """
        identity = np.eye(len(self.__list_site_names))
        return self.__h2_ham_0 - identity * self._e_0 + self.__h2_disord

    def set_central_frequency(self, e_0):
        """Set a central frequency reference for the Hamiltonian.
        
        Stores a frequency offset that is subtracted from the diagonal when
        accessing the working Hamiltonian via `h2_ham`.
        
        Args:
            e_0 (float): Central frequency value to set.
        """
        # Store the central frequency without mutating the base Hamiltonian.
        # The shift is applied dynamically in the `h2_ham` property.
        self._e_0 = e_0

    @property
    def central_frequency(self):
        """Get the current central frequency reference.
        
        Returns:
            float: The central frequency value.
        """
        return self._e_0

    def delete_pigments(self, list_remove_name):
        """Remove pigments from the Hamiltonian.
        
        Removes the specified pigments and updates all internal structures accordingly,
        including the Hamiltonian matrices, site name lists, and index mappings.
        
        Args:
            list_remove_name (list): List of pigment site names to remove.
        """
        for name in list_remove_name:
            if name in self.__list_site_names:
                # Determine the index mapping
                index_old = self.__list_site_names.index(name)
                list_index_keep = list(np.arange(len(self.__list_site_names)))
                list_index_keep.pop(index_old)

                # Update the Hamiltonian Properties
                self.__list_site_names.pop(index_old)
                self.__dict_index_by_label = {name: index for (index, name) in enumerate(self.list_site_names)}
                self.__h2_ham_0 = self.__h2_ham_0[np.ix_(list_index_keep, list_index_keep)]
                self.__h2_disord = self.__h2_disord[np.ix_(list_index_keep, list_index_keep)]
                self._list_sigma.pop(index_old)
            else:
                print(f'ERROR: {name} is not in H2_hamiltonian.list_site_names and could not be removed.')

    def save(self, path, format='npz'):
        """Save the Hamiltonian to a file.
        
        Saves the base Hamiltonian matrix and site names to either a compressed numpy file
        or a CSV file format.
        
        Args:
            path (str): File path where the Hamiltonian should be saved.
            format (str, optional): Output format. Either 'npz' for numpy compressed format
                or 'csv' for CSV format. Defaults to 'npz'.
        
        Raises:
            ValueError: If format is not 'npz' or 'csv'.
        """
        if format == 'npz':
            np.savez_compressed(path, list_site_names=self.list_site_names, h2_ham_0=self.__h2_ham_0)
        elif format == 'csv':
            # Create a DataFrame with site names as both index and columns
            df = pd.DataFrame(
                self.__h2_ham_0,
                index=self.list_site_names,
                columns=self.list_site_names
            )
            df.to_csv(path)
        else:
            raise ValueError(f"Unsupported format '{format}'. Use 'npz' or 'csv'.")

    def diagonalize_hamiltonian(self, hamiltonian_0=None):
        """Diagonalize the provided Hamiltonian or the current internal matrix.
        
        Args:
            hamiltonian_0 (np.ndarray, optional): Matrix to diagonalize. If None,
                the object's current Hamiltonian (including disorder) is used.
        
        Returns:
            tuple: Eigenvalues and eigenvectors of the selected Hamiltonian.
        """
        matrix = self.h2_ham if hamiltonian_0 is None else np.array(hamiltonian_0)
        eigenvalues, eigenvectors = np.linalg.eig(matrix)
        return eigenvalues, eigenvectors

    def calculate_excitonic_overlap(self, hamiltonian_0):
        _, u_hat = self.diagonalize_hamiltonian(hamiltonian_0)
        num_sites_row, num_sites_col = u_hat.shape
        s_matrix = np.zeros((num_sites_row, num_sites_col))
        for mu in range(num_sites_row):
            numerator = (u_hat[mu, :] ** 2 * u_hat[:, :] ** 2)
            denominator = np.sum(u_hat ** 4, axis=0)
            s_matrix[mu, :] = np.sum(numerator / denominator, axis=1)
        return s_matrix

    def build_energetic_domains(self, domain_cutoff: float, hamiltonian_filter: float, N_ens: int = 100):
        """Build energetic domains based on excitonic overlap analysis.
        
        Identifies groups of pigments that are strongly coupled energetically by 
        analyzing excitonic overlap across multiple disorder realizations.
        
        Args:
            domain_cutoff (float): Threshold for excitonic overlap to consider sites
                as belonging to the same domain.
            hamiltonian_filter (float): Minimum absolute coupling value to retain in
                the Hamiltonian. Smaller couplings are set to zero.
            N_ens (int, optional): Number of disorder ensemble realizations to average
                over. Defaults to 100.
                
        Returns:
            list: List of lists, where each inner list contains the site names of
                pigments belonging to a strongly coupled energetic domain.
                
        Raises:
            ValueError: If matrix shapes are inconsistent during computation.
        """
        # Create a deep copy of the Hamiltonian
        hamiltonian = copy.deepcopy(self)

        # Ensure h2_ham is converted to NumPy array
        hamiltonian_0 = np.array(hamiltonian.h2_ham)
        hamiltonian_0[np.abs(hamiltonian_0) < hamiltonian_filter] = 0  # Apply filter
        # print(hamiltonian_0)
        # Initialize count matrix
        count_matrix = np.zeros_like(hamiltonian._Hamiltonian__h2_ham_0)

        # Add disorder and calculate overlap
        for index in np.arange(N_ens):
            hamiltonian.add_disorder(seed=index)
            # Use hamiltonian.h2_ham which includes both the filter and disorder
            exciton_overlap_matrix = hamiltonian.calculate_excitonic_overlap(hamiltonian.h2_ham)

            # Validate shape of overlap matrix
            if exciton_overlap_matrix.shape != count_matrix.shape:
                raise ValueError("Mismatch in matrix shapes while updating count_matrix")

            count_matrix += (exciton_overlap_matrix > domain_cutoff).astype(int)

        # Create DataFrame and validate list_site_names
        if len(self.list_site_names) != count_matrix.shape[0]:
            raise ValueError("Length of list_site_names does not match the shape of count_matrix")

        df_overlap = pd.DataFrame(count_matrix, index=self.list_site_names, columns=self.list_site_names)

        # Derive strongly coupled domains
        list_pigments_by_domain = self.find_strongly_coupled_domains(threshold=(N_ens / 2), dataframe=df_overlap)

        return list_pigments_by_domain

    def build_coupling_domains(self, coupling_cutoff):
        """Build coupling domains based on direct coupling strength.
        
        Identifies groups of pigments that are strongly coupled based on the
        magnitude of their coupling matrix elements.
        
        Args:
            coupling_cutoff (float): Minimum absolute coupling value to consider
                two sites as strongly coupled.
                
        Returns:
            list: List of lists, where each inner list contains the site names of
                pigments belonging to a strongly coupled domain.
                
        Raises:
            AttributeError: If the Hamiltonian matrix or site names are not defined.
            ValueError: If there is a mismatch between matrix size and site names.
        """
        if not hasattr(self, "h2_ham") or not hasattr(self, "list_site_names"):
            raise AttributeError("Hamiltonian matrix or site names not properly defined.")

        hamiltonian = copy.deepcopy(self)
        hamiltonian_0 = hamiltonian.h2_ham

        if hamiltonian_0.shape[0] != len(self.list_site_names):
            raise ValueError("Mismatch between Hamiltonian matrix size and site names length.")

        df = pd.DataFrame(
            hamiltonian_0,
            index=self.list_site_names,
            columns=self.list_site_names
        )
        return self.find_strongly_coupled_domains(dataframe=df, threshold=coupling_cutoff)

    @staticmethod
    def find_strongly_coupled_domains(threshold=20, dataframe=None):
        """Find strongly coupled domains using graph connectivity analysis.
        
        Uses graph theory to identify connected components where sites are connected
        if their coupling or overlap exceeds the threshold.
        
        Args:
            threshold (float, optional): Minimum absolute value for coupling/overlap
                to consider two sites as connected. Defaults to 20.
            dataframe (pd.DataFrame, optional): DataFrame containing coupling or overlap
                values with site names as row and column indices. Must be provided.
                
        Returns:
            list: List of lists, where each inner list contains the site names of
                pigments in a connected component (strongly coupled domain).
                
        Raises:
            ValueError: If dataframe is None or threshold is invalid.
        """
        if dataframe is None:
            raise ValueError("Dataframe must be provided.")

        if not isinstance(threshold, (int, float)) or threshold < 0:
            raise ValueError("Threshold must be a non-negative number.")

        pigment_names = dataframe.columns.tolist()
        mask = (dataframe.abs() > threshold)
        coupling_mask = mask.to_numpy()

        g = nx.Graph()
        g.add_nodes_from(pigment_names)

        for i in range(len(pigment_names)):
            for j in range(i + 1, len(pigment_names)):
                if coupling_mask[i, j]:
                    g.add_edge(pigment_names[i], pigment_names[j])

        return [list(component) for component in nx.connected_components(g)]

    def _subset_hamiltonian(self, list_pigment_names):
        """Subset the Hamiltonian to include only specified pigments.
        
        Updates the internal Hamiltonian matrix, disorder matrix, site names list,
        and index mapping to include only the specified pigments.
        
        Args:
            list_pigment_names (list): List of pigment site names to keep in the
                Hamiltonian.
        """
        # Determine the indices of the pigments to keep
        indices_to_keep = [self.__dict_index_by_label[name] for name in list_pigment_names]

        # Update the internal Hamiltonian matrices and site names
        self.__h2_ham_0 = self.__h2_ham_0[np.ix_(indices_to_keep, indices_to_keep)]
        self.__h2_disord = self.__h2_disord[np.ix_(indices_to_keep, indices_to_keep)]
        self.__list_site_names = [self.__list_site_names[i] for i in indices_to_keep]
        self.__dict_index_by_label = {name: idx for idx, name in enumerate(self.__list_site_names)}

        # Update the disorder list if it exists
        if self._list_sigma is not None:
            self._list_sigma = [self._list_sigma[i] for i in indices_to_keep]
