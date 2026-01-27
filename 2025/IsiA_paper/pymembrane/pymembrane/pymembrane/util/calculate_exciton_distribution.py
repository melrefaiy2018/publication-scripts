import numpy as np
from typing import List, Dict, Tuple, Union, Optional
import random
from pymembrane.util.physical_constants import hbar, c
from matplotlib import pyplot as plt
from matplotlib.pyplot import rcParams
import seaborn as sns
import pandas as pd
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde
import os 


# Reset Matplotlib settings to default values
plt.rcdefaults()
# Set display options to print the complete DataFrame without truncation
pd.set_option('display.max_columns', None)  # Display all columns
pd.set_option('display.expand_frame_repr', False)  # Disable wrapping

# Plotting parameters
params = {
    'axes.labelsize': 30,
    'axes.spines.right': False,
    'axes.spines.top': False,
    'axes.spines.left': True,
    'font.size': 30,
    'legend.fontsize': 16,
    'xtick.labelsize': 30,
    'xtick.major.top': False,
    'xtick.minor.top': False,
    'xtick.minor.bottom': False,
    'ytick.labelsize': 30,
    'ytick.major.right': False,
    'ytick.minor.right': False,
    'xtick.major.pad': 8,
    'ytick.major.pad': 8,
    'legend.handletextpad': 0.4,
    'legend.labelspacing': 0.4,
    'figure.figsize': [7, 6]
}

rcParams.update(params)

def calculate_exciton_distribution(
    H_data: pd.DataFrame, 
    list_pigment_domains: List[List[str]], 
    N_ens: int, 
    sigma_e: float, 
    protein_atomic: object, 
    calc_type: str
) -> Union[
    Tuple[Dict[str, Tuple[List[float], List[float]]], List[str]],
    Tuple[Dict[str, Tuple[List[float], List[float]]], List[str], Dict[str, dict]]
]:
    """
    Calculate the exciton state distribution for each pigment site.

    Parameters:
    ----------
    H_data : pd.DataFrame
        DataFrame containing the Hamiltonian matrix for the system.
    list_pigment_domains : list of lists
        List of pigment domains where each domain is a list of pigment labels.
    N_ens : int
        Number of disorder realizations to average over.
    sigma_e : float
        Standard deviation of the Gaussian disorder applied to the site energies.
    protein_atomic : object
        Object containing dipole moment information for each pigment.
    calc_type : str
        Type of calculation ('abs' for absorption or 'ld' for linear dichroism).

    Returns:
    -------
    dict_distribution_by_site : dict
        Dictionary with pigment labels as keys and a tuple of lists containing
        eigenvalues and corresponding probabilities.
    list_site_label : list
        List of pigment labels ordered as in the Hamiltonian matrix.
    dict_Data : dict, optional
        Dictionary with dipole moments and related factors for each pigment.
        This is only returned if calc_type is 'ld'.
    """
    # Get list_site_label from Hamiltonian matrix columns
    list_site_label = H_data.columns.tolist()

    # Create a mapping from pigment labels to their indices in the Hamiltonian matrix
    label_to_index = {label: idx for idx, label in enumerate(list_site_label)}

    # Convert DataFrame to numpy array
    H = H_data.values

    # Initialize the distribution dictionary
    dict_distribution_by_site = {label: ([], []) for label in label_to_index.keys()}

    dict_Data = {}

    for index in range(N_ens):
        rng = np.random.RandomState(seed=index)

        # Create a random energy perturbation vector for each pigment
        energy_perturbation = rng.normal(0, scale=sigma_e, size=len(label_to_index))
        H_perturbed = H + np.diag(energy_perturbation)

        # Calculate the distribution for the perturbed Hamiltonian
        if calc_type == 'abs':
            calculate_abs_distribution(H_perturbed, label_to_index, list_pigment_domains, dict_distribution_by_site)
        else:
            calculate_ld_distribution(H_perturbed, label_to_index, list_pigment_domains, dict_distribution_by_site, protein_atomic, dict_Data)

    if calc_type == 'abs':
        return dict_distribution_by_site, list_site_label
    else:
        return dict_distribution_by_site, list_site_label, dict_Data
    
def calculate_abs_distribution(
    H2_site: np.ndarray, 
    label_to_index: Dict[str, int], 
    list_pigment_domains: List[List[str]], 
    dict_distribution_by_site: Dict[str, Tuple[List[float], List[float]]]
) -> None:
    """
    Calculate the exciton state distribution for a single perturbed Hamiltonian.

    Parameters:
    ----------
    H2_site : np.ndarray
        Perturbed Hamiltonian matrix for the system.
    label_to_index : dict
        Mapping from pigment labels to their indices in the Hamiltonian matrix.
    list_pigment_domains : list of lists
        List of pigment domains where each domain is a list of pigment labels.
    dict_distribution_by_site : dict
        Dictionary to store the distribution results, with pigment labels as keys
        and a tuple of two lists: eigenvalues and corresponding squared eigenvector
        coefficients (probabilities).
    
    Returns:
    -------
    None
    """
    # Construct the domain coefficients and energies
    for domain in list_pigment_domains:
        # Get the indices of pigments in the domain
        domain_indices = [label_to_index[pigment] for pigment in domain]

        # Solve the eigenproblem for the domain sub-Hamiltonian
        output = np.linalg.eigh(H2_site[np.ix_(domain_indices, domain_indices)])
        domain_energies = output[0]  # Eigenvalues (energies)
        domain_coefficients = output[1]  # Eigenvectors (coefficients)

        # Store the eigenvectors and eigenvalues for each site in the domain
        for site in domain:
            site_index = label_to_index[site]
            index_rel = domain_indices.index(site_index)

            # Append energies and squared eigenvector coefficients to the distribution dictionary
            dict_distribution_by_site[site][0].extend(domain_energies)
            dict_distribution_by_site[site][1].extend(np.abs(domain_coefficients[index_rel, :]) ** 2)

def calculate_angle(
    dipole_exciton: np.ndarray, 
    reference_direction: np.ndarray
) -> float:
    """
    Calculate the angle between the dipole_exciton and the reference_direction.

    Parameters:
    ----------
    dipole_exciton : np.ndarray
        Dipole moment vector of the exciton.
    reference_direction : np.ndarray
        Reference direction vector (e.g., membrane normal).

    Returns:
    -------
    theta : float
        Angle between the dipole_exciton and the reference_direction in radians.
    """
    # Calculate the dot product
    dot_product = np.dot(dipole_exciton, reference_direction)
    # Calculate the magnitudes
    vector_mag = np.linalg.norm(dipole_exciton) * np.linalg.norm(reference_direction)
    # Calculate cos(theta)
    cos_theta = dot_product / vector_mag
    # Ensure cos_theta is in the valid range for arccos due to floating point precision
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    # Return the angle in radians
    theta = np.arccos(cos_theta)
    
    return theta

def calculate_ld_distribution(
    H2_site: np.ndarray, 
    label_to_index: Dict[str, int], 
    list_pigment_domains: List[List[str]], 
    dict_distribution_by_site: Dict[str, Tuple[List[float], List[float]]], 
    protein_atomic: object, 
    dict_Data: Dict[str, Tuple[np.ndarray, float, float]]
) -> None:
    """
    Calculate the linear dichroism (LD) exciton state distribution for a single perturbed Hamiltonian.

    Parameters:
    ----------
    H2_site : np.ndarray
        Perturbed Hamiltonian matrix for the system.
    label_to_index : dict
        Mapping from pigment labels to their indices in the Hamiltonian matrix.
    list_pigment_domains : list of lists
        List of pigment domains where each domain is a list of pigment labels.
    dict_distribution_by_site : dict
        Dictionary to store the distribution results. The key is the pigment label, and the value is a tuple 
        of two lists: one for eigenvalues (energies) and one for the modified probabilities.
    protein_atomic : object
        Object containing dipole moment information for each pigment.
    dict_Data : dict
        Dictionary to store dipole moments and related factors (cosine angle and LD factor) for each pigment.
    
    Returns:
    -------
    None
    """
    membrane_normal = np.array([0, 0, 1])  # Define the membrane normal vector.

    # Construct the domain coefficients and energies.
    for domain in list_pigment_domains:
        print(domain)  # Debug print to show the current domain.
        
        # Get the indices of pigments in the domain.
        domain_indices = [label_to_index[pigment] for pigment in domain]
        
        # Solve the eigenproblem for the domain sub-Hamiltonian.
        output = np.linalg.eigh(H2_site[np.ix_(domain_indices, domain_indices)])
        domain_energies = output[0]  # Eigenvalues (energies).
        domain_coefficients = output[1]  # Eigenvectors (coefficients).

        # Process each site in the domain.
        for site in domain:
            # Get the dipole vector of the current site.
            dipole_vector = protein_atomic.dict_pigments[f'{site}'].dipole
            
            # Calculate the angle between the dipole vector and the membrane normal.
            cos_theta = calculate_angle(dipole_vector, membrane_normal)
            factor = 1 - 3 * cos_theta ** 2  # LD factor based on the angle.

            # Get the index of the site in the Hamiltonian matrix.
            site_index = label_to_index[site]
            index_rel = domain_indices.index(site_index)
            
            # Store dipole moments, cos_theta, and LD factor for each pigment in dict_Data.
            dict_Data[site] = (dipole_vector, cos_theta, factor)
            
            # Modify the probabilities using the LD factor.
            modified_probabilities = np.abs(domain_coefficients[index_rel, :]) ** 2 * factor
            
            # Store the eigenvalues (energies) and modified probabilities in the distribution dictionary.
            dict_distribution_by_site[site][0].extend(domain_energies)
            dict_distribution_by_site[site][1].extend(modified_probabilities)

def plot_ld_distributions(
    data: Dict[str, Tuple[np.ndarray, np.ndarray]], 
    LD_exp: Optional[np.ndarray] = None, 
    LD_sim: Optional[np.ndarray] = None, 
    x_min: int = 600, 
    x_max: int = 720, 
    x_interval: int = 20, 
    chain_name: Optional[str] = None, 
    path_saving: str = f'{os.getcwd()}/'
) -> None:
    """
    Plot the energy distribution for all keys in the data dictionary, separating
    positive and negative linear dichroism (LD) probabilities.

    Parameters:
    ----------
    data : dict
        Dictionary containing energy values and corresponding probabilities.
        The keys are labels (e.g., pigments), and values are tuples of energies 
        and corresponding probabilities.
    LD_exp : np.ndarray, optional
        Experimental linear dichroism data for comparison. The first column should 
        be the wavelengths, and the second column the corresponding experimental values.
    LD_sim : np.ndarray, optional
        Simulated linear dichroism data for comparison. The first column should 
        be the wavelengths, and the second column the corresponding simulated values.
    x_min : int, optional
        Minimum wavelength (nm) for the x-axis. Default is 600.
    x_max : int, optional
        Maximum wavelength (nm) for the x-axis. Default is 720.
    x_interval : int, optional
        Interval between x-ticks on the x-axis. Default is 20.
    chain_name : str, optional
        If provided, the plot will be saved with this name included in the filename.
    path_saving : str, optional
        Path where the plot image will be saved. Default is the current working directory.

    Returns:
    -------
    None
    """
    # Get unique colors and markers for each key
    keys = list(data.keys())
    unique_colors, unique_markers = get_unique_colors_and_markers(keys)
    
    # Initialize the figure and axes
    fig, ax = plt.subplots(figsize=(13, 10))

    # Iterate over each key in the dictionary
    for i, (key, (energies, probabilities)) in enumerate(data.items()):
        # Convert energies to wavelengths
        wavelengths = [1e7 / energy for energy in energies]

        # Convert to numpy arrays for easier manipulation
        wavelengths = np.array(wavelengths)
        probabilities = np.array(probabilities)

        # Filter out NaNs and Infs
        valid_indices = np.isfinite(wavelengths) & np.isfinite(probabilities)
        wavelengths = wavelengths[valid_indices]
        probabilities = probabilities[valid_indices]

        # Separate positive and negative probabilities for LD plotting
        pos_probabilities = np.array([w if w >= 0 else 0 for w in probabilities], dtype=float)
        neg_probabilities = np.array([w if w < 0 else 0 for w in probabilities], dtype=float)

        # Normalize positive and negative probabilities for visualization
        if np.sum(np.abs(pos_probabilities)) > 0:
            pos_probabilities /= np.sum(np.abs(pos_probabilities))
        if np.sum(np.abs(neg_probabilities)) > 0:
            neg_probabilities /= np.sum(np.abs(neg_probabilities))

        # Plot histograms for positive and negative probabilities
        if np.any(pos_probabilities > 0):
            ax.hist(wavelengths, bins=30, weights=pos_probabilities, alpha=0.5, label=f'{key} Positive', color=unique_colors[i])

        if np.any(neg_probabilities < 0):
            ax.hist(wavelengths, bins=30, weights=neg_probabilities, alpha=0.5, label=f'{key} Negative', color=unique_colors[i], linestyle='--')

        # Compute and plot KDE for positive probabilities
        if np.any(pos_probabilities > 0):
            kde_pos = gaussian_kde(wavelengths, weights=pos_probabilities)
            x_pos = np.linspace(min(wavelengths), max(wavelengths), 1000)
            y_pos = kde_pos(x_pos)
            ax.plot(x_pos, y_pos, color=unique_colors[i], lw=2)

        # Compute and plot KDE for negative probabilities
        if np.any(neg_probabilities < 0):
            kde_neg = gaussian_kde(wavelengths, weights=np.abs(neg_probabilities))  # Use absolute values for KDE
            x_neg = np.linspace(min(wavelengths), max(wavelengths), 1000)
            y_neg = kde_neg(x_neg)
            ax.plot(x_neg, -y_neg, color=unique_colors[i], lw=2, linestyle='--')  # Flip the KDE curve to the negative side

    # Plot experimental and simulated LD data if provided
    if LD_exp is not None:
        ax.plot(LD_exp[:, 0], LD_exp[:, 1], linewidth=8, color='black', label='Experiment')
    if LD_sim is not None:
        ax.plot(LD_sim[:, 0], LD_sim[:, 1], color='red', label='Sim.', linewidth=8, alpha=0.7)

    # Customize the plot
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Density')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(width=3, length=7)
    ax.legend(frameon=False)
    for spine in ['bottom', 'left']:
        ax.spines[spine].set_linewidth(3)
        ax.tick_params(width=3, length=7)

    # Set x-axis limits and intervals
    ax.set_xlim(x_min, x_max)
    ax.set_xticks(np.arange(x_min, x_max, x_interval))

    # Save and show the plot
    plt.tight_layout()
    if chain_name:
        plt.savefig(f'{path_saving}LD_distribution_{chain_name}.png', dpi=600)
    else:
        plt.savefig(f'{path_saving}LD_distribution.png', dpi=600)
    plt.show()

def get_unique_colors_and_markers(
    list_site_label: List[str], 
    color_mapping: Optional[Dict[str, str]] = None
) -> Tuple[List[str], List[str]]:
    """
    Generate unique colors and markers for plotting based on site labels.

    Parameters:
    ----------
    list_site_label : List[str]
        List of pigment labels or other identifiers for which colors and markers are needed.
    color_mapping : Optional[Dict[str, str]], optional
        Optional mapping to override specific label colors.

    Returns:
    -------
    unique_colors : List[str]
        List of unique color codes for each site label.
    unique_markers : List[str]
        List of unique markers for each site label.
    """
    custom_color_palette = ['#E63946', '#457B9D', '#1D3557', '#2A9D8F', '#7209B7', '#EA4C89',
                            '#6A0572', '#9A031E', '#4A90E2', '#4ECDC4', '#1E90FF', '#FF6B6B',
                            '#4287f5', '#a51c1c', '#5438dc', '#218c74', '#161D6F', '#2C363F',
                            '#6D326D', '#101820']
    custom_markers = ['o', 's', '^', 'v', 'D', 'p', 'H', 'X', '*', '+', '.', '|', '_', '1', '2', '3', '4', '<', '>',
                      '8', 's', 'p', 'P', 'X', 'D', 'd', 'H', 'h', '^', 'v', '<', '>']

    num_unique_colors = min(len(list_site_label), len(custom_color_palette))
    num_unique_markers = min(len(list_site_label), len(custom_markers))

    unique_colors = random.sample(custom_color_palette, num_unique_colors)
    unique_markers = random.sample(custom_markers, num_unique_markers)

    if color_mapping:
        for label in list_site_label:
            if label in color_mapping:
                unique_colors[list_site_label.index(label)] = color_mapping[label]

    return unique_colors, unique_markers

def setup_plot(ax: plt.Axes, xlabel: str, ylabel: str) -> None:
    """
    Set up the plot aesthetics.

    Parameters:
    ----------
    ax : plt.Axes
        Matplotlib axes object to customize.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    """
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(width=3, length=7)
    for spine in ['bottom', 'left', 'right', 'top']:
        ax.spines[spine].set_linewidth(3)
    plt.grid(axis='x')
    plt.tight_layout()


def annotate_highest_peak(
    ax: plt.Axes, 
    x: np.ndarray, 
    y: np.ndarray, 
    label: str, 
    color: str, 
    existing_labels: List[Tuple[float, float]],
    fontsize: int = None,         # New parameters for individual font properties
    fontweight: str = None,
    family: str = None
) -> None:
    """
    Annotate the highest peak in a given plot.

    Parameters:
    ----------
    ax : plt.Axes
        Matplotlib axes object where the peak will be annotated.
    x : np.ndarray
        Array of x-values (wavelengths).
    y : np.ndarray
        Array of y-values (probabilities or intensities).
    label : str
        The label to annotate (typically a site label).
    color : str
        Color for the annotation text and arrow.
    existing_labels : List[Tuple[float, float]]
        List of previously placed labels' coordinates to avoid overlap.
    fontsize : int, optional
        Font size for the annotation.
    fontweight : str, optional
        Font weight for the annotation (e.g., 'bold').
    family : str, optional
        Font family for the annotation (e.g., 'serif', 'sans-serif').
    """
    label = label.split('_')[-1]  # Example: 4_CLA_601 -> 601
    peaks, _ = find_peaks(y)
    if peaks.size > 0:
        highest_peak = peaks[np.argmax(y[peaks])]
        x_peak = x[highest_peak]
        y_peak = y[highest_peak]

        # Adjust the label position to avoid overlap
        offset = 0.1 * max(y)  # Initial offset
        for (existing_x, existing_y) in existing_labels:
            if abs(existing_x - x_peak) < 5 and abs(existing_y - y_peak) < offset:
                y_peak += offset
                break

        ax.annotate(f'{label}',
                    xy=(x_peak, y_peak),
                    xytext=(x_peak, y_peak + offset),
                    arrowprops=dict(facecolor=color, shrink=0.05),
                    ha='center', color=color,
                    fontsize=fontsize,  # Apply fontsize
                    fontweight=fontweight,  # Apply fontweight
                    family=family)  # Apply font family

        existing_labels.append((x_peak, y_peak + offset))
def plot_exciton_distribution(
    exciton_distribution: Dict[str, Tuple[List[float], List[float]]], 
    list_site_label: List[str], 
    path_saving: str, 
    show_labels: bool = True, 
    color_mapping: Optional[Dict[str, str]] = None, 
    x_min: int = 600, 
    x_max: int = 720, 
    x_interval: int = 20, 
    legend: bool = False,
    fontsize=16, 
    fontweight='bold', 
    font_style='serif'
) -> None:
    """
    Plot the exciton distribution for each site label.

    Parameters:
    ----------
    exciton_distribution : dict
        Dictionary with site labels as keys and tuples containing energy values 
        and corresponding probabilities as values.
    list_site_label : list of str
        List of pigment site labels to be plotted.
    path_saving : str
        Directory path where the plot image will be saved.
    show_labels : bool, optional
        If True, annotate the highest peaks for each site label. Default is True.
    color_mapping : dict, optional
        Optional mapping to override specific label colors.
    x_min : int, optional
        Minimum value for the x-axis (wavelength). Default is 600.
    x_max : int, optional
        Maximum value for the x-axis (wavelength). Default is 720.
    x_interval : int, optional
        Interval for the x-axis ticks. Default is 20.
    legend : bool, optional
        If True, display the legend. Default is False.

    Returns:
    -------
    None
    """
    fig, ax = plt.subplots(figsize=(13, 10))
    unique_colors, unique_markers = get_unique_colors_and_markers(list_site_label, color_mapping)
    existing_labels = []

    for i, site in enumerate(list_site_label):
        list_energy = exciton_distribution[site][0]
        list_nm = [c * 2 * np.pi * hbar / energy for energy in list_energy]
        list_weight = exciton_distribution[site][1]

        color = unique_colors[i]
        marker = unique_markers[i % len(unique_markers)]

        kde = sns.kdeplot(x=list_nm, weights=list_weight, cut=2, linewidth=4, alpha=0.7, ax=ax, color=color, label=f'{site}')

        if show_labels:
            # Get data for the current KDE plot
            line = kde.get_lines()[-1]
            x_data, y_data = line.get_data()

            # Annotate the highest peak
            annotate_highest_peak(ax, x_data, y_data, site, color, existing_labels, fontsize=fontsize, fontweight=fontweight, family=font_style)

    setup_plot(ax, 'Wavelength (nm)', 'Density')
    plt.xticks(np.arange(x_min, x_max + 1, x_interval))
    if legend:
        ax.legend(frameon=False)
    plt.savefig(f'{path_saving}exciton_distribution.png', dpi=600)
    plt.show()

def plot_combined_absorption_and_exciton(
    sliced_lambda_axis_a: np.ndarray, 
    absorption_data: np.ndarray,
    exp_absorption: Optional[np.ndarray], 
    exciton_distribution: Dict[str, Tuple[List[float], List[float]]], 
    list_site_label: List[str], 
    temp: int, 
    path_saving: str,
    show_labels: bool = True, 
    color_mapping: Optional[Dict[str, str]] = None, 
    x_min: Optional[int] = None, 
    x_max: Optional[int] = None,
    fontsize=16, 
    fontweight='bold', 
    font_style='serif'
) -> None:
    """
    Plot the combined absorption data and exciton distribution.

    Parameters:
    ----------
    sliced_lambda_axis_a : np.ndarray
        Wavelength axis data for absorption.
    absorption_data : np.ndarray
        Absorption data to be plotted.
    exp_absorption : np.ndarray, optional
        Experimental absorption data for comparison.
    exciton_distribution : dict
        Dictionary with site labels as keys and tuples containing energy values 
        and corresponding probabilities as values.
    list_site_label : list of str
        List of pigment site labels to be plotted.
    temp : int
        Temperature at which the simulation or experiment is conducted.
    path_saving : str
        Directory path where the plot image will be saved.
    show_labels : bool, optional
        If True, annotate the highest peaks for each site label. Default is True.
    color_mapping : dict, optional
        Optional mapping to override specific label colors.
    x_min : int, optional
        Minimum value for the x-axis (wavelength).
    x_max : int, optional
        Maximum value for the x-axis (wavelength).
    fontsize : int, optional
        Font size for the annotation.
    fontweight : str, optional
        Font weight for the annotation (e.g., 'bold').
    font_style : str, optional
        Font family for the annotation (e.g., 'serif', 'sans')

    Returns:
    -------
    None
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    ax1.plot(sliced_lambda_axis_a, absorption_data / np.max(absorption_data), linewidth=8, color='green', label=f'A_mohamed_{temp}', alpha=0.5)
    if exp_absorption is not None:
        ax1.plot(exp_absorption[:, 0], exp_absorption[:, 1] / np.max(exp_absorption[:, 1]), color='black', alpha=0.95, label='Experiment', linewidth=8)
    
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylabel(f'Absorption (a.u.)')
    ax1.set_yticks([0, 0.5, 1])
    ax1.set_yticklabels(labels=[0, 0.5, 1])
    # ax1.legend()

    unique_colors, unique_markers = get_unique_colors_and_markers(list_site_label, color_mapping)
    existing_labels = []

    for i, site in enumerate(list_site_label):
        list_energy = exciton_distribution[site][0]
        list_nm = [c * 2 * np.pi * hbar / energy for energy in list_energy]
        list_weight = exciton_distribution[site][1]

        color = unique_colors[i % len(unique_colors)]
        # marker = unique_markers[i % len(unique_markers)]

        kde = sns.kdeplot(x=list_nm, weights=list_weight, cut=2, linewidth=4, alpha=0.7, ax=ax2, color=color, label=f'{site}')

        if show_labels:
            # Get data for the current KDE plot
            line = kde.get_lines()[-1]
            x_data, y_data = line.get_data()

            # Annotate the highest peak
            annotate_highest_peak(ax2, x_data, y_data, site, color, existing_labels, fontsize=fontsize, fontweight=fontweight, family=font_style)

    setup_plot(ax2, 'Wavelength (nm)', 'Density')
    # ax2.legend(frameon=False)
    plt.savefig(f'{path_saving}exciton_distribution_with_absorption.png', dpi=600)
    plt.show()