import copy
import os
from itertools import cycle

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from cycler import cycler
from numpy import linalg as LA
from pyvis.network import Network
from pymembrane.util.physical_constants import c, hbar
from tqdm import tqdm

def apply_default_plot_settings():
    """
    Function to apply default plot settings for all future figures.
    """
    rcParams = {
        'font.size': 16,
        'axes.labelsize': 20,
        'axes.titlesize': 22,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'figure.figsize': (12, 6),
        'figure.dpi': 300,
        'lines.linewidth': 4,
        'grid.color': 'gray',
        'grid.linestyle': '--',
        'grid.linewidth': 0.3,
        'grid.alpha': 0.5,
        'axes.linewidth': 1.5,
        'axes.edgecolor': 'black',
    }

    # Apply the settings to Matplotlib
    plt.rcParams.update(rcParams)
apply_default_plot_settings()


def get_sliced_data(total_array_data, lambda_axis, num_values):
    """
    Returns a sliced version of the data arrays abs_array and lambda_axis_a, centered around the maximum value.

    Parameters:
    total_array_data (ndarray): Input data array for the y-axis.
    lambda_axis (ndarray): Input data array for the x-axis.
    num_values (int): Number of values to include before and after the maximum.

    Returns:
    sliced_lambda_axis (ndarray): Sliced array for the x-axis.
    sliced_total_array (ndarray): Sliced array for the y-axis.
    """

    # Find the index of the maximum value
    max_index = np.argmax(total_array_data)

    # Get a range of values before and after the maximum
    start_index = max_index - num_values
    end_index = max_index + num_values
    sliced_total_array = total_array_data[start_index:end_index]
    sliced_lambda_axis = lambda_axis[start_index:end_index]

    return sliced_lambda_axis, sliced_total_array


# build a function that will create a directory if it does not exist:
def create_directory(directory_path):
    """
    Creates a directory if it does not exist.

    Parameters:
    - directory_path (str): The path to the directory to be created.

    Returns:
    - None
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    return None


def convert_to_nm(energy):
    """Convert energy (1/cm) to wavelength in nm."""
    return (2 * np.pi * hbar) * c / energy

def convert_nm_to_wavenumber(wavelength_nm):
    """Convert wavelength in nm to wavenumber (1/cm)."""
    return 10 ** 7 / wavelength_nm



def get_dict_diagonal_values_from_df(dataframe):
    """
    Function to extract the index and diagonal values from a DataFrame and return them as a dictionary.

    Args:
    dataframe (pandas.DataFrame): Input DataFrame.

    Returns:
    dict: Dictionary containing the index and corresponding diagonal values.
    """
    # Convert DataFrame to NumPy array
    array = dataframe.to_numpy()

    # Get diagonal values
    diagonal_values = np.diag(array)

    # Create a dictionary for index and diagonal values
    dict_diagonal = dict(zip(dataframe.index, diagonal_values))

    return dict_diagonal
    
def get_dict_diagonal_values_from_dict(dict_hamiltonian):
    """
    Function to extract diagonal values from a nested dictionary where the diagonal is defined as
    the value whose key in the inner dictionary matches the key of the outer dictionary.

    Args:
    big_dict (dict): Input nested dictionary.

    Returns:
    dict: Dictionary containing the outer keys and corresponding diagonal values.
    """
    dict_diagonal = {}

    # Iterate through the big_dict to find diagonal values
    for outer_key, inner_dict in dict_hamiltonian.items():
        # Check if the outer key is in the inner dictionary
        if outer_key in inner_dict:
            # The diagonal value is where the key of the outer dictionary matches the key in the inner dictionary
            dict_diagonal[outer_key] = inner_dict[outer_key]

    return dict_diagonal

def sort_dict_site_energy(dict_site_energy):
    """Return site energies sorted into a new dictionary by increasing energy."""
    sorted_keys = sorted(dict_site_energy, key=lambda k: dict_site_energy[k])
    return {key: dict_site_energy[key] for key in sorted_keys}

def sort_diagonal_values_df(df):
    """Return a dict of sorted diagonal values from a square DataFrame."""
    diagonal_values = pd.Series(np.diagonal(df), index=df.index)
    sorted_diagonal = diagonal_values.sort_values(ascending=True)
    return sorted_diagonal.to_dict()

def build_SitEnergy_df(data_dict):
    """Format a site energy dict and return a pivoted DataFrame."""
    check_key_format = list(data_dict.keys())[0]
    if len(check_key_format.split('_')) == 4:
        data = [
            {
                'Chain': key.split('_')[1],
                'Pigment': f"{key.split('_')[2]}_{key.split('_')[3]}",
                'site_shift': value,
            }
            for key, value in data_dict.items()
        ]
    elif len(check_key_format.split('_')) == 3:
        data = [
            {
                'Chain': key.split('_')[0],
                'Pigment': f"{key.split('_')[1]}_{key.split('_')[2]}",
                'site_shift': value,
            }
            for key, value in data_dict.items()
        ]
    else:
        raise ValueError('The key format is not recognized!')

    df = pd.DataFrame(data)
    df_pivoted = df.pivot(index='Pigment', columns='Chain', values='site_shift').reset_index()
    df_pivoted.columns.name = None
    df_pivoted.columns = ['Pigment'] + [f'Chain_{chain}' for chain in df_pivoted.columns[1:]]

    return df_pivoted




def build_domains_data(protein_atomic, list_pigments_by_domain, show_plot=False):
    """Compute pigment positions and domain metadata; optionally show a scatter plot."""
    dict_pigment_data = {}
    # Iterate over each pigment:
    for pigment, value in protein_atomic.dict_pigments.items():
        # Extract the xyz vector of the pigment
        xyz_vector = value.location[:2]
        # Add the pigment name and xyz vector to the dictionary
        dict_pigment_data[pigment] = xyz_vector

    # visualize pigment distributions:
    # ================================
    # Extract x, y coordinates for each pigment
    x = [v[0] for v in dict_pigment_data.values()]
    y = [v[1] for v in dict_pigment_data.values()]

    if show_plot:
        # Plot the dots
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111)  # Remove projection argument for 2D plot
        ax.scatter(x, y, c='b', label='Pigments')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('Pigment Data Visualization')
        # ax.legend()
        # plt.savefig(f'{md_scratch}Dynamics_{temp}K/Pigment_Data_Visualization.png')
        plt.show()

    # Dictionary to store the angles and xy coordinates of each domain
    dict_domain_info = {}
    # Iterate over each domain and its associated pigments
    for idx, domain_pigments in enumerate(list_pigments_by_domain):
        # Calculate the center of mass of pigments in the domain
        domain_center = np.mean([dict_pigment_data[f'{pigment}'] for pigment in domain_pigments], axis=0)
        # Calculate the vector from the origin to the domain center
        domain_vector = domain_center

        # Calculate the angle between the fixed vector and the domain vector using arctan2
        angle = np.arctan2(domain_vector[1], domain_vector[0])  # Angle in radians from -pi to pi
        angle_degrees = np.degrees(angle) % 360  # Convert to degrees and ensure the range from 0 to 360

        # get the chain name of the domain:
        chain_name = domain_pigments[0].split('_')[0]

        # Initialize the domain entry
        dict_domain_info[idx] = {
            'angle_degrees': angle_degrees,
            'xy_coordinates': domain_center.tolist(),  # Convert numpy array to list for better JSON compatibility
            'chain': chain_name,
            'pigments': domain_pigments
        }

    return dict_domain_info, dict_pigment_data

def format_pigment_label(pigment):
    """Shorten pigment label for plotting."""
    parts = pigment.split('_')
    if len(parts) >= 4:
        return f"{parts[1]}_{parts[3]}"
    else:
        return pigment

def visualize_pigments_network(list_pigments_by_domain, dict_pigment_info, save_path=None, view='2D', selected_chains=None):
    """
    Visualizes pigment networks either in 2D or 3D.

    Parameters:
    - list_pigments_by_domain: List of lists, where each sublist represents pigments in a specific domain.
    - dict_pigment_info: Dictionary containing positional information for pigments.
    - save_path: Optional path to save the visualization.
    - view: '2D' for 2D visualization, '3D' for 3D visualization.
    - selected_chains: Optional list of chains to filter the visualization by specific chains.
    """

    G = nx.Graph()
    color_cycle = cycle([
        'forestgreen', 'skyblue', 'gold', 'salmon', 'orchid', 'crimson', 'chocolate', 'slategrey',
        'teal', 'coral', 'indigo', 'olive', 'tomato', 'turquoise', 'maroon', 'peru',
        'navy', 'limegreen', 'magenta', 'sienna', 'steelblue', 'plum', 'darkkhaki', 'seagreen',
        'orangered', 'mediumvioletred', 'cadetblue', 'saddlebrown', 'darkslateblue', 'lightseagreen',
        'firebrick', 'mediumorchid', 'darkcyan'
    ])
    node_colors = []

    for pigments in list_pigments_by_domain:
        color = next(color_cycle)
        for pigment in pigments:
            chain = pigment.split('_')[0]
            if selected_chains is None or chain in selected_chains:
                G.add_node(pigment, color=color)
                node_colors.append(color)

        pigments_in_graph = [p for p in pigments if p in G]
        for i in range(len(pigments_in_graph)):
            for j in range(i + 1, len(pigments_in_graph)):
                G.add_edge(pigments_in_graph[i], pigments_in_graph[j])

    pos = {pigment: dict_pigment_info[f'{pigment}'] for pigment in G.nodes}
    labels = {pigment: format_pigment_label(pigment) for pigment in G.nodes}

    if view == '2D':
        hex_colors = [mcolors.rgb2hex(color) for color in node_colors]
        plt.figure(figsize=(12, 12))

        # Improved node style: black border, larger size
        nx.draw_networkx_nodes(G, pos, node_color=hex_colors, edgecolors='black', linewidths=1.5, node_size=3500)
        nx.draw_networkx_edges(G, pos, edge_color="#BBBBBB", width=2, alpha=0.5)
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=14, font_weight='bold', font_color='white')

        plt.title("Pigments Network Visualization")
        plt.axis('off')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Visualization saved as {save_path}")
        else:
            plt.show()

    elif view == '3D':
        pos_3d = {key: np.append(np.array(value), np.random.normal(scale=3)) for key, value in pos.items()}

        edge_x, edge_y, edge_z = [], [], []
        for edge in G.edges():
            x0, y0, z0 = pos_3d[edge[0]]
            x1, y1, z1 = pos_3d[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_z.extend([z0, z1, None])

        edge_trace = go.Scatter3d(x=edge_x, y=edge_y, z=edge_z, mode='lines',
                                  line=dict(color='#BBBBBB', width=2))

        node_x, node_y, node_z, node_color, node_text = [], [], [], [], []
        for node in G.nodes():
            x, y, z = pos_3d[node]
            node_x.append(x)
            node_y.append(y)
            node_z.append(z)
            node_color.append(mcolors.rgb2hex(G.nodes[node]['color']))
            node_text.append(format_pigment_label(node))

        node_trace = go.Scatter3d(
            x=node_x, y=node_y, z=node_z,
            mode='markers+text',
            marker=dict(size=15, color=node_color),
            text=node_text,
            textposition="top center"
        )

        fig = go.Figure(data=[edge_trace, node_trace])
        fig.update_layout(title="Pigments Network Visualization (3D)", showlegend=False)

        if save_path:
            fig.write_html(save_path)
            print(f"Visualization saved as {save_path}")
        else:
            fig.show()



def visualize_pyvis_pigments_network(list_pigments_by_domain, dict_pigment_info, save_path=None):
    """Visualize pigment domains with PyVis and optionally persist to HTML."""
    net = Network(height="800px", width="1000px", notebook=True)

    # Generate a unique color for each domain
    colors = plt.cm.get_cmap('hsv', len(list_pigments_by_domain))

    # Add nodes and edges with colors and coordinates
    for idx, pigments in enumerate(list_pigments_by_domain):
        color = colors(idx)
        hex_color = mcolors.rgb2hex(color)
        for pigment in pigments:
            x_coord = dict_pigment_info[f'{pigment}'][0]
            y_coord = dict_pigment_info[f'{pigment}'][1]
            net.add_node(pigment, label=pigment, title=pigment, color=hex_color, x=x_coord, y=y_coord)
        for i in range(len(pigments)):
            for j in range(i + 1, len(pigments)):
                net.add_edge(pigments[i], pigments[j])

    # Save or show the network
    if save_path:
        net.save_graph(save_path)
        print(f"Visualization saved as {save_path}")
    else:
        net.show_buttons(filter_=['physics'])
        net.show(f"{save_path}pyvis_pigments_network.html")

def prepare_domain_energy_data(raw_data):
    """
    Prepares domain energy data from raw text input.

    Parameters:
        raw_data (str): Raw input data as a string.

    Returns:
        dict: Dictionary with domain names as keys and energy values as lists.
    """
    dict_data = {}
    for index, line in enumerate(raw_data.split('\n')):
        line = line.strip()
        if line:
            parts = line.split(':')
            pigment_names = [name.strip(" '[]") for name in parts[1].split(',')]
            values = [float(v.strip()) for v in parts[2].strip('[]').split(',')]
            dict_data[f'Domain_{index}'] = dict(zip(pigment_names, values))
    return dict_data

def plot_energy_distribution(ax, energy_values, domain_names, lowest_energy):
    """
    Plots energy distribution on the given axis.

    Parameters:
        ax (matplotlib.axes.Axes): Axis to plot on.
        energy_values (list): List of energy values.
        domain_names (list): List of domain names.
        lowest_energy (float): Lowest energy value for reference.
    """
    ax.boxplot(energy_values, widths=0.7, medianprops=dict(color='red', linewidth=2), boxprops=dict(linewidth=2))
    ax.plot([1, len(domain_names)], [lowest_energy + 200, lowest_energy + 200], 'k--', lw=2)
    ax.plot([1, len(domain_names)], [lowest_energy, lowest_energy], 'k--', lw=2)
    ax.fill_between([1, len(domain_names)], lowest_energy, lowest_energy + 200, color='gray', alpha=0.3)
    ax.set_xticklabels(domain_names, rotation=90)
    ax.set_xlabel('Domain')
    ax.set_ylabel('Energy')
    ax.grid(True)

def plot_domain_energy(
    linear_protein=None,
    chains=None,
    multiple_figures=False,
    num_subplots=10,
    fig_size=(12, 8),
    save_path='domain_energy_visualization'
):
    """
    Plot energy distribution across domains with different options.

    Parameters:
        linear_protein (object): Object containing domain energy data.
        chains (list, optional): List of chain identifiers to filter by chains.
        multiple_figures (bool, optional): Flag to plot multiple subplots.
        num_subplots (int, optional): Number of subplots if multiple_figures is True.
        fig_size (tuple, optional): Size of the plot figure.
        save_path (str, optional): Path to save the plot.

    Returns:
        None
    """
    if linear_protein is None:
        raise ValueError("linear_protein cannot be None")

    domain_names = []
    energy_values = []
    lowest_energy = float('inf')

    # Collect all energy values and calculate the lowest energy
    for index in linear_protein.list_domain:
        domain_energy = index._list_e_exc
        min_domain_energy = np.min(domain_energy)
        if min_domain_energy < lowest_energy:
            lowest_energy = min_domain_energy

    # Filter domains based on chains if specified and collect domain names and energy values
    for index in linear_protein.list_domain:
        # Extract the list of pigment names (e.g., ['CP24_4_CLA_602', 'CP24_4_CLA_603'])
        list_names = index.list_names
        
        # Extract the index part from each pigment name (e.g., '602' from 'CP24_4_CLA_602')
        extracted_indices = [f"{name.split('_')[-3]}_{name.split('_')[-1]}" for name in list_names]
        
        # Join the indices with underscores (e.g., '602_603')
        domain_name = '_'.join(extracted_indices)
        
        # Access the domain energy or any other attribute as needed
        domain_energy = index._list_e_exc

        # Print or use the generated domain name and energy as needed
        print(f"Domain Name: {domain_name}, Domain Energy: {domain_energy}")

        if chains and not any(chain in domain_name for chain in chains):
            continue

        energy_values.append(domain_energy)
        domain_names.append(domain_name)

    # domain_names = transform_domain_names(domain_names)

    print(f'Lowest energy: {lowest_energy}')

    if multiple_figures:
        fig, axs = plt.subplots(num_subplots, 1, figsize=(15, 6 * num_subplots), sharex=True)
        for i, ax in enumerate(axs):
            start_index = i * (len(domain_names) // num_subplots)
            end_index = (i + 1) * (len(domain_names) // num_subplots) if i < num_subplots - 1 else len(domain_names)
            plot_energy_distribution(ax, energy_values[start_index:end_index], domain_names[start_index:end_index], lowest_energy)
    else:
        fig, ax = plt.subplots(figsize=fig_size)
        plot_energy_distribution(ax, energy_values, domain_names, lowest_energy)

    plt.title('Energy of Domains')
    plt.tight_layout()

    if not os.path.exists(save_path):
        os.makedirs(save_path)

    if chains is not None:
        chains_str = "_".join(chains)
        save_path = f'{save_path}/domain_energy_{chains_str}.png'
    else:
        save_path = f'{save_path}/domain_energy.png'
    plt.tight_layout()
    plt.savefig(save_path)
    plt.show()

def extract_matching_data(csv_file_path, desired_names):
    """
    Extracts data from a CSV matrix that matches the given row and column names.

    Parameters:
        csv_file_path (str): Path to the CSV file containing the matrix.
        desired_names (list): List of row and column names to be extracted.

    Returns:
        pd.DataFrame: DataFrame containing the extracted data.
    """
    # Read the CSV file into a Pandas DataFrame
    df = pd.read_csv(
        csv_file_path, index_col=0
    )  # Assuming row names are in the first column

    # Filter the DataFrame to keep only the desired row and column names
    filtered_df = df.loc[desired_names, desired_names]

    return filtered_df


def thermal_rate_time_scale(list_accpetor_domains, list_donor_domains, temperature):
    """
    Estimate thermal rate time scale between donor and acceptor domains.

    Expects `generalized_forster_exciton` and `atomic_c2s2m2` to be available in scope.
    """
    K_Dons = []
    for donor_index, domain in enumerate(list_donor_domains):
        K_alpha_summed = np.sum(np.concatenate([generalized_forster_exciton(acceptor_domain, domain, temperature, H2_hamiltonian=atomic_c2s2m2.H2_hamiltonian)
                                                for acceptor_domain in list_accpetor_domains],axis=0), axis=0)
        K_Dons.extend(K_alpha_summed)
    raw_thermal_pop = [domain._raw_thermal_pop(temperature) for domain in list_donor_domains] 
    raw_thermal_pop_flat = np.concatenate(raw_thermal_pop)
    P_d = raw_thermal_pop_flat/np.sum(raw_thermal_pop_flat) #normalized P_d
    K_eff = np.sum(P_d[don_ind] * K_Dons[don_ind] for don_ind in range(len(P_d)))
    K_eff = (1/ K_eff)/1000  # convert timescale to ps
    return K_eff

def prepare_subchain_exciton_calculation(chain_name, list_pigments_by_domain, hamiltonian_path):
    """
    Prepare the subchain exciton calculation by extracting the relevant data from the Hamiltonian.

    Args:
        chain_name (str): The name of the chain.
        list_pigments_by_domain (list): A list of pigments grouped by domain.
        hamiltonian_path (str): The file path to the Hamiltonian data.

    Returns:
        tuple: A tuple containing the list of sub-pigment domains and the extracted sub-chain Hamiltonian.
    """
    df_H_data = pd.read_csv(hamiltonian_path, index_col=0)
    # list_sub_pigment_domains = [pigments for pigments in list_pigments_by_domain if all(pigment.startswith(f'{chain_name}') for pigment in pigments)] # 'N_CLA_602'
    list_sub_pigment_domains = [pigments for pigments in list_pigments_by_domain if
                                all(pigment.split('_')[1] == chain_name for pigment in pigments)]  # 'LHCII_N_CLA_602'

    # prepare the hamiltonian:
    # list_sub_chain = [pigment for pigment in df_H_data.columns.to_list() if f'{chain_name}_' in pigment]
    list_sub_chain = [pigment for pigment in df_H_data.columns.to_list() if f'{chain_name}' in pigment.split('_')[1]]
    H_sub_chain = extract_matching_data(hamiltonian_path, list_sub_chain)
    return list_sub_pigment_domains, H_sub_chain





def build_domains_mo(protein_atomic, calc_method, path_hamiltonian, domain_cutoff=0.1, hamiltonian_filter=15, N_ens=1000, coupling_cutoff=20, saving_dir=f'{os.getcwd()}/'):
    """
    Build domains based on the specified calculation method.

    Parameters:
    - calc_method (str): Either 'Energy' or 'Coupling'.
    - domain_cutoff (float): Cutoff value for excitonic overlap based domain grouping.
    - hamiltonian_filter (int): Off-diagonal filter for the Hamiltonian before overlap calculation.
    - N_ens (int): Number of disorder realizations to average over.
    - coupling_cutoff (int): Coupling cutoff used for the 'Coupling' method.

    Returns:
    - list_pigments_by_domain: List of pigment groups per domain.
    - dict_domain_by_name: Mapping from pigment name to domain index.
    """
    df_H = pd.read_csv(path_hamiltonian, index_col=0)
    list_columns_name = df_H.columns.tolist()

    if calc_method == 'Energy':
        print('#### Build Domains based on Energy Overlap ####')
        list_pigments_by_domain = _build_energetic_domains(
            protein_atomic,
            N_ens,
            domain_cutoff,
            hamiltonian_filter,
            list_columns_name,
        )
        print(f'{list_pigments_by_domain=}')
        print(f'Number of domains: {len(list_pigments_by_domain)}')
        with open(f'{saving_dir}/list_pigments_by_domain_{domain_cutoff}.txt', 'w') as f:
            f.write(f'{list_pigments_by_domain}=')
    elif calc_method == 'Coupling':
        print('#### Build Domains based on their Electronic coupling based on Renger definition ####')
        list_pigments_by_domain = find_strongly_coupled_domains(
            path_hamiltonian=path_hamiltonian,
            threshold=coupling_cutoff,
        )
        print(f'{list_pigments_by_domain=}')
        print(f'Number of domains: {len(list_pigments_by_domain)}')
        with open(f'{saving_dir}/list_pigments_by_domain_{domain_cutoff}.txt', 'w') as f:
            f.write(f'{list_pigments_by_domain}=')
    else:
        raise ValueError('The provided method is not defined!')

    dict_domain_by_name = {}
    for index_domain, list_names in enumerate(list_pigments_by_domain):
        for name in list_names:
            dict_domain_by_name[name] = index_domain

    return list_pigments_by_domain, dict_domain_by_name


def _build_energetic_domains(protein_atomic, N_ens, domain_cutoff, hamiltonian_filter, list_columns_name):
    """
    Build energetic domains from disorder realizations and overlap cutoffs.

    Args:
        N_ens (int): Number of disorder realizations.
        domain_cutoff (float): Cutoff value for excitonic overlap matrix.
        hamiltonian_filter (float): Filter value for off-diagonal elements of the Hamiltonian.

    Returns:
        list: List of strongly coupled pigments by domain.
    """
    hamiltonian = copy.deepcopy(protein_atomic.H2_hamiltonian)
    hamiltonian_0 = hamiltonian.H2_ham
    hamiltonian_0[np.abs(hamiltonian_0) < hamiltonian_filter] = 0

    count_matrix = np.zeros_like(hamiltonian_0)

    for index in tqdm(np.arange(N_ens)):
        protein_atomic.H2_hamiltonian.add_disorder(seed=index)
        excitonic_overlap_matrix = __calculate_excitonic_overlap(protein_atomic.H2_hamiltonian.H2_ham)
        count_matrix += (excitonic_overlap_matrix > domain_cutoff).astype(int)

    df_overlap = pd.DataFrame(count_matrix, index=list_columns_name, columns=list_columns_name)
    return find_strongly_coupled_domains(dataframe=df_overlap, threshold=N_ens / 2)


def __diagonalize_hamiltonian(H2_ham):
    """
    Diagonalizes the given Hamiltonian matrix.

    Parameters:
    - H2_ham: The Hamiltonian matrix to be diagonalized.

    Returns:
    - U_hat: The unitary matrix containing the eigenvectors of the Hamiltonian matrix.
    """
    eigenvalues, eigenvectors = LA.eig(H2_ham)
    U_hat = eigenvectors
    return U_hat

def __calculate_excitonic_overlap(U):
    """
    Calculates the excitonic overlap matrix.

        Parameters:
        - U: numpy.ndarray
            The input matrix U.

        Returns:
        - S_matrix: numpy.ndarray
            The calculated excitonic overlap matrix.
        """
    U_hat = __diagonalize_hamiltonian(U)
    num_sites_row, num_sites_col = U_hat.shape
    S_matrix = np.zeros((num_sites_row, num_sites_col))
    for mu in range(num_sites_row):
        numerator = (U_hat[mu, :] ** 2 * U_hat[:, :] ** 2)
        denominator = np.sum(U_hat ** 4, axis=0)
        S_matrix[mu, :] = np.sum(numerator / denominator, axis=1)
    return S_matrix


def find_strongly_coupled_domains(threshold=20, path_hamiltonian=None, dataframe=None, output_name=None):
    """
    Finds strongly coupled domains in a Hamiltonian matrix.

    Parameters:
    - path_hamiltonian (str): The file path to the Hamiltonian matrix in CSV format.
    - threshold (float): The threshold value for coupling strength. Default is 20.

    Returns:
    - formatted_output (list): A list of lists representing the strongly coupled domains.

    The function loads the Hamiltonian matrix from the given CSV file, extracts the coupling values
    greater than the threshold or less than the negative threshold, builds a graph based on the pigment
    names, finds the connected nodes in the graph, and returns the strongly coupled domains as a list
    of lists. It also saves the formatted output as a text file named 'strongly_coupled_domains.txt'
    in the same directory as the Hamiltonian matrix file.
    """
    # Load data from CSV file
    if path_hamiltonian is None:
        if dataframe is not None:
            data = dataframe
        else:
            raise ValueError('Either path_hamiltonian or dataframe should be provided!')
    else:
        data = pd.read_csv(path_hamiltonian, index_col=0)

    pigment_names = data.columns.tolist()

    # Extract coupling values greater than threshold or less than -threshold
    mask = (data > threshold) | (data < -threshold)
    coupling_mask = mask.to_numpy()

    # Build the graph
    G = nx.Graph()
    G.add_nodes_from(pigment_names)

    for i in range(len(pigment_names)):
        for j in range(i + 1, len(pigment_names)):
            if coupling_mask[i, j]:
                G.add_edge(pigment_names[i], pigment_names[j])

    # Find connected nodes
    connected_nodes = list(nx.connected_components(G))

    # Format the output as a list of lists:
    formatted_output = []
    for component in connected_nodes:
        formatted_output.append(list(component))

    # save the output as a txt file:
    # ==============================
    # Extract the directory path
    dir_path = os.path.dirname(path_hamiltonian) if path_hamiltonian is not None else os.getcwd()
    if output_name is None:
        file_name = f'{dir_path}/strongly_coupled_domains_cutoff_{threshold}.txt'
    else:
        file_name = f'{dir_path}/{output_name}.txt'

    # # Open the file in write mode ("w")
    # with open(file_name, "w") as file:
    #     # Use repr() to convert the entire list to a string preserving the format
    #     file.write(repr(formatted_output))

    return formatted_output




# Step 13: Initialize the population for each chain (normalized)
def simulate_initial_population(nonlinear_system, condition):
    """
    Simulates and normalizes the initial population vector for different excitation conditions.
    
    Parameters:
    nonlinear_system : object
        The system containing the domains and pigment information.
    condition : str
        The excitation condition, can be 'all', 'CLA', or 'CHL'.
    
    Returns:
    np.ndarray
        The normalized initial population vector based on the excitation condition.
    """
    # Initialize the population array with zeros
    initial_population = np.zeros(len(nonlinear_system.list_domain), dtype=np.float64)

    # Loop over each domain and populate based on condition
    for i, domain in enumerate(nonlinear_system.list_domain):
        if condition == 'all':
            # Excite all pigments
            initial_population[i] = len(domain.list_names)
        elif condition == 'CLA':
            # Excite only CLA pigments
            initial_population[i] = sum(1 for name in domain.list_names if 'CLA' in name)
        elif condition == 'CHL':
            # Excite only CHL pigments
            initial_population[i] = sum(1 for name in domain.list_names if 'CHL' in name)
        else:
            raise ValueError("Invalid condition. Use 'all', 'CLA', or 'CHL'.")
    
    # Normalize the population
    total_population = np.sum(initial_population)
    if total_population > 0:
        initial_population /= total_population  # Normalize only if total population is non-zero
    
    return initial_population

def simulate_population(nonlinear_system, initial_population, N_t=100000, dt=100):
    """
    Simulates the time evolution of the population for the given initial population.
    
    Parameters:
    nonlinear_system : object
        The system containing domains and rate matrices.
    initial_population : np.ndarray
        The initial population vector.
    N_t : int
        Number of time steps (default: 100000).
    dt : int
        Time step size (default: 100).
        
    Returns:
    t_axis : np.ndarray
        The time axis array.
    domains_dict : dict
        Dictionary with keys as domain indices and values as population arrays for each domain.
    """
    # Step 14: Time axis and propagation of populations over time
    t_axis = np.arange(0, N_t * dt, dt)
    t_axis, list_pt_mixed = nonlinear_system.time_evolve(P1_0=initial_population, N_t=N_t, dt=dt)

    # Convert list_pt_mixed to a NumPy array for easier manipulation
    list_pt_mixed = np.array(list_pt_mixed)

    # Create a dictionary to store the populations for each domain
    domains_dict = {}

    # Sum the population for each domain at each time step
    for i in range(list_pt_mixed.shape[1]):  # list_pt_mixed.shape[1] gives the number of domains (columns)
        domains_dict[f'domain_{i}'] = list_pt_mixed[:, i]

    return t_axis, domains_dict

def plot_population_over_time(t_axis, domains_dict, nonlinear_system, condition_name, xlim=50, fig_x=10, fig_y=8):
    """
    Plots the populations over time for each domain.
    
    Parameters:
    t_axis : np.ndarray
        The time axis array.
    domains_dict : dict
        Dictionary with keys as domain indices and values as population arrays for each domain.
    nonlinear_system : object
        The system containing the list of domains.
    condition_name : str
        The name of the condition (e.g., 'all', 'CLA', 'CHL') for labeling the plot.
    xlim : int
        The x-axis limit for the plot (default: 50).
        
    Returns:
    None
    """
    plt.figure(figsize=(fig_x, fig_y))

    # Plot each domain dynamically
    for i in range(len(domains_dict)):  # len(domains_dict) gives the number of domains
        pigment_index = [name.split('_')[-1] for name in nonlinear_system.list_domain[i].list_names]  # Extract pigment index

        plt.plot(t_axis / 1000, domains_dict[f'domain_{i}'], 
                 label=pigment_index, linewidth=2)

    plt.xlim(0, xlim)
    plt.xlabel('Time (ps)')
    plt.ylabel('Population')
    plt.title(f'{condition_name}: Population of Excited States Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(f'{condition_name}_dynamic.png', dpi=600)
    plt.show()



def simulate_population_by_type(nonlinear_system, initial_population, N_t=100000, dt=100):
    """
    Simulates the time evolution of the population and sums the populations for CLA and CHL pigments.
    
    Parameters:
    nonlinear_system : object
        The system containing domains and rate matrices.
    initial_population : np.ndarray
        The initial population vector.
    N_t : int
        Number of time steps (default: 100000).
    dt : int
        Time step size (default: 100).
        
    Returns:
    t_axis : np.ndarray
        The time axis array.
    CLA_population : np.ndarray
        Summed population over time for all CLA pigments.
    CHL_population : np.ndarray
        Summed population over time for all CHL pigments.
    """
    # Step 14: Time axis and propagation of populations over time
    t_axis = np.arange(0, N_t * dt, dt)
    t_axis, list_pt_mixed = nonlinear_system.time_evolve(P1_0=initial_population, N_t=N_t, dt=dt)

    # Convert list_pt_mixed to a NumPy array for easier manipulation
    list_pt_mixed = np.array(list_pt_mixed)

    # Initialize arrays to store the summed populations for CLA and CHL
    CLA_population = np.zeros(list_pt_mixed.shape[0])
    CHL_population = np.zeros(list_pt_mixed.shape[0])

    # Loop over domains to sum populations for CLA and CHL
    for i, domain in enumerate(nonlinear_system.list_domain):
        domain_names = domain.list_names  # Get the names for this domain

        if any('CLA' in name for name in domain_names):  # Check if it's a CLA domain
            CLA_population += list_pt_mixed[:, i]  # Sum over time for CLA
        elif any('CHL' in name for name in domain_names):  # Check if it's a CHL domain
            CHL_population += list_pt_mixed[:, i]  # Sum over time for CHL

    return t_axis, CLA_population, CHL_population


def plot_population_by_type(t_axis, CLA_population, CHL_population, condition_name, xlim=50):
    """
    Plots the summed populations over time for CLA and CHL pigments.
    
    Parameters:
    t_axis : np.ndarray
        The time axis array.
    CLA_population : np.ndarray
        Summed population over time for all CLA pigments.
    CHL_population : np.ndarray
        Summed population over time for all CHL pigments.
    condition_name : str
        The name of the condition (e.g., 'all', 'CLA', 'CHL') for labeling the plot.
    xlim : int
        The x-axis limit for the plot (default: 50).
        
    Returns:
    None
    """
    plt.figure(figsize=(10, 6))

    # Plot CLA population over time
    plt.plot(t_axis / 1000, CLA_population, label="CLA", linewidth=2)

    # Plot CHL population over time
    plt.plot(t_axis / 1000, CHL_population, label="CHL", linewidth=2)

    plt.xlim(0, xlim)
    plt.xlabel('Time (ps)')
    plt.ylabel('Population')
    plt.title(f'{condition_name}: Summed Population of CLA and CHL Over Time')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(f'{condition_name}_CLA_CHL_dynamic.png', dpi=600)
    plt.show()

def plot_population_fits(t_axis, donor_population, acceptor_population, protein_name, q2_hardcode=None, save_path=None, xlim=100):
    """
    Plots the fitted population dynamics for donor and acceptor states based on a rate equation.

    Parameters
    ----------
    t_axis : np.array(float)
        The time points of the analyzed population dynamics.
    donor_population : np.array(float)
        The population of the donor state over time.
    acceptor_population : np.array(float)
        The population of the acceptor state over time.
    q2_hardcode : float or None, optional
        The thermalized population of the acceptor state. Set to the average of the 
        final 100 data points of acceptor_population if None.
    save_path : str, optional
        File path to save the plot. If None, the plot will not be saved.

    Returns
    -------
    k_rate : float
        The calculated transfer rate.
    t : float
        The decay time (tau) in ps.
    """
    # Get the mean population for the acceptor to use for fitting
    pop_mean_sum = np.mean(acceptor_population[-100:])  # Using last 100 points for mean

    # Calculate the rate of transport
    k_rate = _find_rate_of_transport(t_axis, acceptor_population, q2_hardcode)
    
    # Calculate decay time tau
    t = ((1 / k_rate) / 1000)  # Adjusting time scaling
    
    # Fit the population curves
    pop_accp_fit = _rate_fit(t_axis, k_rate, pop_mean_sum)
    don_pop_fit = 1 - pop_accp_fit

    # Set up the plot
    plt.figure(figsize=(8, 6))

    # Plot fitted acceptor and donor population
    plt.plot(t_axis / 1000, pop_accp_fit, color="black", linestyle="--", linewidth=3, label="Acc. Pop. Fit")
    plt.plot(t_axis / 1000, acceptor_population, color="black", linewidth=3, label="Acc. Pop")
    plt.plot(t_axis / 1000, don_pop_fit, color="green", linestyle="--", linewidth=3, label="Donor Fit")
    plt.plot(t_axis / 1000, donor_population, color="green", linewidth=3, label="Donor")

    # Set the title with LaTeX for tau
    plt.title(f'{protein_name}: $\\tau = {round(t)}$ ps', fontsize=16, weight='bold')
    plt.xlabel("Time (ps)", fontsize=14)
    plt.ylabel("Population", fontsize=14)

    # Set axis limits and ticks
    plt.xlim(-5, xlim)
    # plt.xticks([0, 50, 100], fontsize=12)
    plt.ylim(-0.05, 1.05)
    plt.yticks([0.0, 0.5, 1.0], fontsize=12)

    # Customize the legend
    plt.legend(frameon=False, fontsize=12, loc="upper right")

    # Add a grid for better readability
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)

    # Ensure layout is tight to avoid cutting off labels
    plt.tight_layout()

    # Save the figure if save_path is provided
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')

    # Display the plot
    plt.show()

    # Clear the figure to prevent overlap in future plots
    plt.clf()

    return k_rate, t


def _find_rate_of_transport(t_axis, population_2, q2_hardcode=None):
    """
    Finds the rate of transport, k, for the equation P1(0) * (1 - np.exp(-2*k * t_axis))
    by analyzing P2(t) in a system with two populations (P1 and P2).

    Parameters
    ----------
    t_axis : np.array(float)
        The time points of the analyzed population dynamics.
    population_2 : np.array(float)
        The population of the acceptor state at all times.
    q2_hardcode : float or None
        The thermalized population of the acceptor state. Set to the average of the 
        final 10 data points of population_2 if None.

    Returns
    -------
    k : float
        The transfer rate.
    """
    if q2_hardcode is None:
        q2 = np.mean(population_2[-10:])  # Average of last 10 points
    else:
        q2 = q2_hardcode

    R = (q2 ** -2) * np.trapz(np.abs(q2 - population_2), x=t_axis)
    return 1 / R


def _rate_fit(t_axis, k, q2):
    """
    Finds the population dynamics of an acceptor state from a rate equation. Assumes
    all population begins in the donor state.

    Parameters
    ----------
    t_axis : np.array(float)
        The time points of the analyzed population dynamics.
    k : float
        The transport rate (inverse time units).
    q2 : float
        Equilibrium population of the acceptor state.

    Returns
    -------
    population_2 : np.array(float)
        The population of the acceptor state at all times.
    """
    return q2 * (1 - np.exp(-(1 / q2) * k * t_axis))


def set_publication_style():
    """Set publication-quality plot styling."""
    # Use LaTeX for all text rendering
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "text.latex.preamble": r"\usepackage{amsmath}"
    })

    # Font sizes
    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    LARGE_SIZE = 24
    BIGGER_SIZE = 28

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=LARGE_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=LARGE_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    # Set figure DPI for display and saving
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['savefig.dpi'] = 300

    # Set line widths and other parameters
    plt.rcParams['axes.linewidth'] = 3     # Width of the axes lines
    plt.rcParams['xtick.major.width'] = 3  # Width of the major x tick marks
    plt.rcParams['ytick.major.width'] = 3  # Width of the major y tick marks
    plt.rcParams['xtick.major.size'] = 6     # Size of the major x tick marks
    plt.rcParams['ytick.major.size'] = 6     # Size of the major y tick marks

    # Style settings
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['grid.linestyle'] = '--'
    plt.rcParams['grid.linewidth'] = 0.5
    plt.rcParams['grid.alpha'] = 0.7
    
    # Create a colorblind-friendly color cycle
    colors = [
        '#000000',  # Black
        '#1b9e77',  # Dark teal
        '#d95f02',  # Dark orange
        '#7570b3',  # Purple
        '#e7298a',  # Pink
        '#66a61e',  # Green
        '#e6ab02',  # Yellow
        '#a6761d',  # Brown
        '#666666',  # Gray
    ]
    
    plt.rcParams['axes.prop_cycle'] = cycler(color=colors)
    
    return colors

def plot_spectra(params):
    """
    General function to plot spectral data with multiple curves.
    
    Parameters
    ----------
    params : dict
        Dictionary containing plot parameters:
        - 'curves' : list of tuples, each containing:
                    (x_data, y_data, label, [optional kwargs])
        - 'spectra_type' : string, 'Fluorescence' or 'Absorption'
        - 'temp' : temperature in K
        - 'N_ens' : ensemble size
        - 'dir_path' : directory path for saving figures
        - 'panel_label' : panel label (e.g., 'A', 'B')
        - 'x_lim' : tuple with x-axis limits
        - 'y_lim' : tuple with y-axis limits (default (0, 1.05))
        - 'save_pdf' : whether to save PDF (default False)
        - 'show_plot' : whether to show plot (default True)
        - 'normalize' : whether to normalize each curve (default True)
        - 'filename' : custom filename (default based on spectra_type)
        - 'figsize' : figure size as tuple (default (10, 7.5))
        - 'dpi' : DPI for saving (default 300)
        - 'legend_loc' : position for legend (default 'upper right')
        - 'legend_bbox_to_anchor' : bbox_to_anchor for legend positioning
        - 'show_panel_label' : whether to show panel label (default True)
    """
    # Set default values for optional parameters
    spectra_type = params.get('spectra_type', 'Fluorescence')
    panel_label = params.get('panel_label', 'A')
    x_lim = params.get('x_lim', (630, 750) if spectra_type == 'Fluorescence' else (600, 700))
    y_lim = params.get('y_lim', (0, 1.05))
    save_pdf = params.get('save_pdf', False)
    show_plot = params.get('show_plot', True)
    dir_path = params.get('dir_path', './')
    temp = params.get('temp', 300)
    N_ens = params.get('N_ens', 1000)
    normalize = params.get('normalize', True)
    figsize = params.get('figsize', (10, 7.5))
    dpi = params.get('dpi', 300)
    curves = params.get('curves', [])
    legend_loc = params.get('legend_loc', 'upper right')
    legend_bbox_to_anchor = params.get('legend_bbox_to_anchor', None)
    show_panel_label = params.get('show_panel_label', True)
    
    # Print debugging info
    print(f"Plotting {spectra_type} with {len(curves)} curves")
    for i, curve in enumerate(curves):
        print(f"  Curve {i+1}: {curve[2]}")
    
    # Generate default filename if not provided
    default_filename = f"{spectra_type.lower()}_{temp}K_N_ens_{N_ens}"
    filename = params.get('filename', default_filename)
    
    # Ensure output directory exists
    os.makedirs(f'{dir_path}Spectra_figures/', exist_ok=True)
    
    # Set publication style and get color cycle
    colors = set_publication_style()
    
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot each curve from the provided data
    for i, curve_data in enumerate(curves):
        # Extract curve data
        if len(curve_data) >= 3:
            x_data, y_data, label = curve_data[0:3]
            
            # Get any additional kwargs if provided
            kwargs = {}
            if len(curve_data) >= 4 and curve_data[3] is not None:
                kwargs = curve_data[3]
        else:
            raise ValueError("Each curve must be a tuple of (x_data, y_data, label, [optional kwargs])")
            
        # Set default styling parameters if not in kwargs
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 4 if i > 0 else 6
        if 'zorder' not in kwargs:
            kwargs['zorder'] = i + 1
        if 'alpha' not in kwargs:
            kwargs['alpha'] = 0.85 if i > 0 else 1.0
            
        # Create a copy of y_data to avoid modifying the original
        y_data_plot = np.copy(y_data)
            
        # Normalize data if requested
        if normalize:
            y_max = np.max(y_data_plot)
            if y_max > 0:  # Avoid division by zero
                y_data_plot = y_data_plot / y_max
        
        # Plot the curve
        ax.plot(x_data, y_data_plot, label=label, **kwargs)
    
    # Set limits and ticks
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(['0', '0.5', '1'])

    # Add minor ticks for more precision
    ax.minorticks_on()
    ax.tick_params(which='minor', bottom=True, left=True)

    # Labels with LaTeX formatting
    ax.set_xlabel('Wavelength (nm)', fontweight='bold')
    ax.set_ylabel(f'{spectra_type} (a.u.)', fontweight='bold')

    # Legend with improved styling
    if len(curves) > 0:
        if legend_bbox_to_anchor:
            legend = ax.legend(frameon=True, loc=legend_loc, 
                              bbox_to_anchor=legend_bbox_to_anchor, 
                              framealpha=0.9)
        else:
            legend = ax.legend(frameon=True, loc=legend_loc, framealpha=0.9)
            
        legend.get_frame().set_linewidth(1)
        legend.get_frame().set_edgecolor('black')

    # Remove grid for cleaner look
    ax.grid(False)

    # Add subtle background grid
    ax.grid(True, which='major', linestyle='--', linewidth=0.5, color='lightgray', alpha=0.5)

    # Add subtle box around the plot
    for spine in ax.spines.values():
        spine.set_linewidth(1.5)
        
    # Adjust layout
    plt.tight_layout()

    # Add panel label for multi-panel figures
    if show_panel_label:
        ax.text(0.03, 0.95, panel_label, transform=ax.transAxes, 
                fontsize=28, fontweight='bold', va='top')

    # Save with transparent background for publication
    plt.savefig(f"{dir_path}Spectra_figures/{filename}.png", 
               dpi=dpi, bbox_inches='tight', transparent=False)
               
    # Save PDF if requested
    if save_pdf:
        plt.savefig(f"{dir_path}Spectra_figures/{filename}.pdf", 
                   format='pdf', bbox_inches='tight', transparent=True)
    
    # Show plot if requested
    if show_plot:
        plt.show()
    else:
        plt.close()
        
    return fig, ax

# Alias functions for backwards compatibility and convenience
def plot_fluorescence(params):
    """Wrapper for plot_spectra with Fluorescence as spectra_type."""
    # Make a copy of the params dict to avoid modifying the original
    params_copy = params.copy()
    params_copy['spectra_type'] = 'Fluorescence'
    return plot_spectra(params_copy)

def plot_absorption(params):
    """Wrapper for plot_spectra with Absorption as spectra_type."""
    # Make a copy of the params dict to avoid modifying the original
    params_copy = params.copy()
    params_copy['spectra_type'] = 'Absorption'
    return plot_spectra(params_copy)
