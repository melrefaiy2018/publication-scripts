import os
import matplotlib.pyplot as plt
import numpy as np
import re
from tabulate import tabulate
import pandas as pd

# Plotting
SMALL_SIZE = 16
MEDIUM_SIZE = 25
BIGGER_SIZE = 40

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
import matplotlib.pyplot as plt

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

# Call this function in your scripts before generating plots
apply_default_plot_settings()

from typing import List , Dict , Optional


class ProtonationProfile :
    def __init__(self , filepath: str) -> None :
        """
        Initialize with the path to the fort.38 file.

        Parameters:
        - filepath: str, the path to the fort.38 file.
        """
        self.filepath = filepath

    def load_data(self , protein: "Protein") -> None :
        """
        Load data from the fort.38 file and populate the Protein structure.

        Parameters:
        - protein: Protein object to fill with data.
        """
        with open(self.filepath , 'r') as file :
            try :
                # Try to read the first line, which should contain pH values
                pH_values: List[float] = list(map(float , next(file).split()[1 :]))
            except StopIteration :
                # If the file is empty, simply return without processing
                return

            # Process each line representing a residue-conformer
            for line in file :
                line = line.strip()
                if not line :
                    continue

                parts = line.split()
                residue_conformer = parts[0]  # e.g., "ASN01A0002_001" or "GLU-1A0004_005"
                probabilities: List[float] = list(
                    map(float , parts[1 :]))  # List of probabilities at different pH levels

                # Ensure the residue-conformer string is 14 characters long
                if len(residue_conformer) != 14 :
                    raise ValueError(
                        f"Invalid residue-conformer length: {residue_conformer} (length: {len(residue_conformer)})")

                # Extract the protonation state
                protonation_state: str = residue_conformer[3 :5]

                # Validate protonation state: either '±1' or '0' followed by a digit
                if not (protonation_state in ['-1' , '+1', 'DM'] or (
                        protonation_state.startswith('0') and protonation_state[1].isdigit())) :
                    raise ValueError(f"Invalid protonation state: {protonation_state} in {residue_conformer}")

                # Parse residue details
                residue_name: str = residue_conformer[:3]
                chain_name: str = residue_conformer[5]
                residue_index: str = residue_conformer[6 :10]
                conformer_id: str = residue_conformer[11 :]

                # Add data to the hierarchical structure
                chain = protein.get_chain(chain_name)
                if not chain :
                    chain = Chain(chain_name)
                    protein.add_chain(chain)

                residue = chain.get_residue(residue_name , protonation_state , residue_index)
                if not residue :
                    residue = Residue(residue_name , protonation_state , residue_index)
                    chain.add_residue(residue)

                conformer = Conformer(conformer_id , probabilities , pH_values)
                residue.add_conformer(conformer)


class Conformer :
    def __init__(self , conformer_id: str , probabilities: List[float] , pH_values: List[float]) -> None :
        """
        Conformer represents an individual conformer of a residue with its specific probabilities across pH levels.

        Parameters:
        - conformer_id: str, the unique conformer ID (e.g., "001", "002").
        - probabilities: List of float, probabilities at different pH levels.
        - pH_values: List of float, corresponding pH values.
        """
        self.conformer_id = conformer_id
        self.probabilities = probabilities
        self.pH_values = pH_values

    def get_probability_at_pH(self , pH: float) -> Optional[float] :
        """
        Get the probability at a specific pH value.

        Parameters:
        - pH: float, the pH value to get the probability for.

        Returns:
        - probability: float, the probability at the given pH, or None if pH not found.
        """
        if pH in self.pH_values :
            index = self.pH_values.index(pH)
            return self.probabilities[index]
        return None

    def __repr__(self) -> str :
        return f"Conformer {self.conformer_id}, Probabilities: {self.probabilities}"


class Residue :
    def __init__(self , residue_name: str , protonation_state: str , residue_index: str) -> None :
        """
        Residue represents a specific residue within a chain that can have multiple conformers.

        Parameters:
        - residue_name: str, the name of the residue (e.g., "ASP").
        - protonation_state: str, the protonation state (e.g., "01", "02").
        - residue_index: str, the residue index (e.g., "0165").
        """
        self.residue_name = residue_name
        self.protonation_state = protonation_state
        self.residue_index = residue_index
        self.conformers: List[Conformer] = []  # List to hold all conformers for this residue
        self.total_probability = 0.0  # Sum of probabilities for all conformers at pH 7 (or another specific pH)

    def add_conformer(self , conformer: Conformer) -> None :
        """
        Add a conformer to the residue and update the total probability for pH 7.

        Parameters:
        - conformer: Conformer object.
        """
        self.conformers.append(conformer)

        # Get the probability at pH 7 (assuming index 7 corresponds to pH 7)
        probability_at_pH7 = conformer.get_probability_at_pH(7.0)
        if probability_at_pH7 is not None :
            self.total_probability += probability_at_pH7

    def __iter__(self) -> iter :
        """
        Allow iteration over the conformers of the residue.
        """
        return iter(self.conformers)

    def __repr__(self) -> str :
        return f"Residue {self.residue_name}{self.protonation_state}{self.residue_index}"


class Chain :
    def __init__(self , chain_name: str) -> None :
        """
        Chain represents a protein chain that contains multiple residues.

        Parameters:
        - chain_name: str, the name of the chain (e.g., "A", "B", "D").
        """
        self.chain_name = chain_name
        self.residues: Dict[tuple , Residue] = { }  # Dictionary to hold residues by residue index

    def add_residue(self , residue: Residue) -> None :
        """
        Add a residue to the chain.

        Parameters:
        - residue: Residue object.
        """
        residue_key = (residue.residue_name , residue.protonation_state , residue.residue_index)
        if residue_key not in self.residues :
            self.residues[residue_key] = residue

    def get_residue(self , residue_name: str , protonation_state: str , residue_index: str) -> Optional[Residue] :
        """
        Retrieve a residue from the chain based on residue name, protonation state, and index.

        Parameters:
        - residue_name: str
        - protonation_state: str
        - residue_index: str

        Returns:
        - Residue object if found, None otherwise.
        """
        residue_key = (residue_name , protonation_state , residue_index)
        return self.residues.get(residue_key , None)

    def __iter__(self) -> iter :
        """
        Allow iteration over the residues in the chain.
        """
        return iter(self.residues.values())

    def __repr__(self) -> str :
        return f"Chain {self.chain_name}"


class Protein :
    def __init__(self , protein_name: str) -> None :
        """
        Protein represents a complete protein, consisting of multiple chains.

        Parameters:
        - protein_name: str, the name of the protein.
        """
        self.protein_name = protein_name
        self.chains: Dict[str , Chain] = { }  # Dictionary to hold chains by chain name

    def add_chain(self , chain: Chain) -> None :
        """
        Add a chain to the protein.

        Parameters:
        - chain: Chain object.
        """
        if chain.chain_name not in self.chains :
            self.chains[chain.chain_name] = chain

    def get_chain(self , chain_name: str) -> Optional[Chain] :
        """
        Retrieve a chain from the protein based on chain name.

        Parameters:
        - chain_name: str

        Returns:
        - Chain object if found, None otherwise.
        """
        return self.chains.get(chain_name , None)

    def __iter__(self) -> iter :
        """
        Allow iteration over the chains in the protein.
        """
        return iter(self.chains.values())

    def __repr__(self) -> str :
        return f"Protein {self.protein_name}"


def residue_titration_curve_across_pH(protein, target_residue, residue_index, chains_list):
    """
    Visualize the distribution of protonation states for a given residue across all pH levels.
    
    Parameters:
    - protein: Protein object containing chains and residues.
    - target_residue: str, the name of the residue to visualize (e.g., "LYS").
    - residue_index: str, the residue index to analyze (e.g., "0003").
    - chains_list: list of str, the names of the chains to visualize (e.g., ["A", "B", "C"])
    """
    # Store pH values and probabilities for each protonation state
    pH_values = None
    protonation_probabilities = {}

    # Iterate over the specified chains
    for chain_name in chains_list:
        chain = protein.get_chain(chain_name)
        if chain:
            # Find all matching residues
            for residue_key, residue in chain.residues.items():
                if residue.residue_name == target_residue and residue.residue_index == residue_index:
                    # Process all conformers for this residue
                    for conformer in residue.conformers:
                        # Set pH_values if not already set (assuming all conformers have the same pH values)
                        if pH_values is None:
                            pH_values = conformer.pH_values
                        # Collect probabilities for each protonation state
                        key = f"{residue.protonation_state}-{conformer.conformer_id}"
                        protonation_probabilities[key] = conformer.probabilities

    # Plot the protonation distribution across pH levels
    if pH_values:
        plt.figure(figsize=(10, 6))
        
        # Plot each protonation state across pH levels
        for protonation_state, probabilities in protonation_probabilities.items():
            plt.plot(pH_values, probabilities, label=protonation_state)
        
        # Labeling the plot
        plt.xlabel('pH Level')
        plt.ylabel('Probability')
        plt.title(f'Protonation States for {target_residue} {residue_index} across Chains {", ".join(chains_list)}')
        plt.legend(title="Protonation States")
        plt.grid(True)
        plt.tight_layout()

        # Show the plot
        plt.show()
    else:
        print(f"No data found for residue {target_residue}{residue_index} in the specified chains.")



def plot_protonation_distribution(protein, residue_name, residue_index, chains_list, pH_value):
    """
    Plot the distribution of protonation states for a specific residue in a list of chains at a specified pH value.
    
    Parameters:
    - protein: Protein object containing chains and residues.
    - residue_name: str, the name of the residue to analyze (e.g., "ASP").
    - residue_index: str, the residue index to analyze (e.g., "0003").
    - chains_list: list of str, a list of chain names to include in the analysis (e.g., ["A", "B", "D"]).
    - pH_value: float, the specific pH at which to visualize the protonation distribution.
    """
    protonation_states = []
    probabilities = []
    
    # Iterate over the specified chains
    for chain_name in chains_list:
        chain = protein.get_chain(chain_name)
        if chain:
            # Check all residues in the chain for the specified residue name and index
            for residue_key, residue in chain.residues.items():
                if residue.residue_name == residue_name and residue.residue_index == residue_index:
                    # Get the total probability for each conformer at the specified pH value
                    for conformer in residue.conformers:
                        probability_at_pH = conformer.get_probability_at_pH(pH_value)
                        if probability_at_pH is not None:
                            protonation_states.append(f"{residue.protonation_state}-{conformer.conformer_id}")
                            probabilities.append(probability_at_pH)
    
    # Plot the distribution if there are any protonation states found
    if protonation_states:
        plt.figure(figsize=(8, 6))
        plt.bar(protonation_states, probabilities, color='skyblue')
        plt.xlabel("Protonation States")
        plt.ylabel(f"Total Probability at pH {pH_value}")
        plt.title(f"Distribution of Protonation States for {residue_name}{residue_index} in Chains {chains_list} at pH {pH_value}")
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.show()
    else:
        print(f"No data found for residue {residue_name}{residue_index} in the specified chains at pH {pH_value}.")

def visualize_protonation(self, target_residue, chains_list, pH_value):
    """
    Visualize the protonation states for all residues of the specified type across a list of chains at a specified pH value.
    Residues with protonation state 01 or 02 will be colored blue, and others will be colored red.
    """
    protonation_sums = {}
    x_labels = []
    residue_indices = []
    
    # Iterate over the specified chains
    for chain_name in chains_list:
        chain = self.get_chain(chain_name)  # Retrieve the chain from the protein
        if chain:
            for residue_key, residue in chain.residues.items():
                if residue.residue_name == target_residue:
                    for conformer in residue.conformers:
                        probability_at_pH = conformer.get_probability_at_pH(pH_value)
                        if probability_at_pH is not None:
                            key = f"{residue.protonation_state}-{residue.residue_index}-{chain_name}"
                            protonation_sums[key] = protonation_sums.get(key, 0) + probability_at_pH
                            if key not in x_labels:
                                x_labels.append(f"{chain_name}-{residue.protonation_state}-{residue.residue_index}")
                                residue_indices.append(residue.residue_index)  # Track residue index

    residues = list(protonation_sums.keys())
    probabilities = list(protonation_sums.values())
    
    # Create a dictionary to group bars by residue index
    grouped_indices = {}
    for i, residue in enumerate(residues):
        residue_index = residue.split('-')[1]
        if residue_index not in grouped_indices:
            grouped_indices[residue_index] = []
        grouped_indices[residue_index].append(i)
    
    # Determine x positions for each group of residue index
    x_positions = []
    current_position = 0
    for residue_index in sorted(grouped_indices.keys()):  # Sort residue indices
        for i in grouped_indices[residue_index]:
            x_positions.append(current_position)
            current_position += 0.5  # Keep bars close for the same residue index
        current_position += 1.5  # Add gap between different residue indices
    
    # Create a list to hold colors: blue for protonation states "01" or "02", red for others
    colors = ['blue' if res.split('-')[0] in ['01', '02'] else 'red' for res in residues]
    
    # Create figure and axis
    plt.figure(figsize=(16, 6))
    
    # Bar plot
    plt.bar(x_positions, probabilities, color=colors)
    
    # Set x-ticks at the same positions as the bars and label them with chain name, protonation state, and residue index
    plt.xticks(x_positions, residues, rotation=90, ha='center')  # Center-align the labels

    # Add labels and title
    plt.ylabel(f"Probability at pH {pH_value}")
    plt.title(f'{target_residue} Protonation State Probabilities in Chains {", ".join(chains_list)} at pH {pH_value}')
    
    # Adjust layout to prevent overlap of x-axis labels
    plt.tight_layout()
    
    # Save the figure to a directory named after the residue
    directory = target_residue
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Save the figure as a PNG file inside the directory
    filename = f"{directory}/{target_residue}_protonation_chains_{'_'.join(chains_list)}_pH{pH_value}.png"
    plt.savefig(filename)
    print(f"Figure saved to {filename}")
    
    # Optionally display the plot
    plt.show()

def visualize_stacked_protonation(self, target_residue, chains_list, pH_value, legend_box='T'):
    protonation_sums = {}
    residue_indices = []

    # Regex to capture protonation state correctly, including negative states
    protonation_state_pattern = re.compile(r'^([+-]?\d+)-(\d+)-(.+)$')

    for chain_name in chains_list:
        chain = self.get_chain(chain_name)
        if chain:
            for residue_key, residue in chain.residues.items():
                if residue.residue_name == target_residue:
                    for conformer in residue.conformers:
                        probability_at_ph = conformer.get_probability_at_pH(pH_value)
                        if probability_at_ph is not None:
                            key = f"{residue.protonation_state}-{residue.residue_index}-{chain_name}"
                            protonation_sums[key] = protonation_sums.get(key, 0) + probability_at_ph
                            if (chain_name, residue.residue_index) not in residue_indices:
                                residue_indices.append((chain_name, residue.residue_index))

    # Sort residue indices by chain name and residue index (dropping residue name)
    residue_indices = sorted(residue_indices, key=lambda x: (x[0], x[1]))

    if not protonation_sums:
        print("No data in protonation_sums, nothing to plot.")
        return  # Exit early if there's no data

    # Extract conformers using the regex to handle negative protonation states
    conformers = sorted(list(set([protonation_state_pattern.match(key).group(1) for key in protonation_sums.keys()])))

    # Create a dictionary to accumulate values
    bar_data = {f"{chain}-{res_index}": {conf: 0 for conf in conformers}
                for chain, res_index in residue_indices}

    # Accumulate probabilities
    for key, prob in protonation_sums.items():
        match = protonation_state_pattern.match(key)
        if match:
            conf = match.group(1)  # Protonation state (e.g., '-1', '01', '02')
            res_index = match.group(2)  # Residue index as string
            chain_name = match.group(3)  # Chain name
            combined_key = f"{chain_name}-{res_index}"  # Combine chain name and residue index
            if combined_key in bar_data:
                bar_data[combined_key][conf] += prob

    # Create lists for the stacked bar chart (using chain and residue index, dropping residue name)
    labels = [f"{chain}-{res_index}" for chain, res_index in residue_indices]
    data = np.array([[bar_data[f"{chain}-{res_index}"][conf] for conf in conformers]
                     for chain, res_index in residue_indices])

    if data.size == 0:
        print("No data to plot.")
        return  # Exit early if there's no data

    # Create the figure with a more professional aspect ratio and increased size
    fig, ax = plt.subplots(figsize=(12, 8))

    # Set bottom to zero for stacking bars
    bottom = np.zeros(len(labels))

    # Define a colorblind-friendly color palette
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

    # Plot the stacked bars
    for i, conf in enumerate(conformers):
        ax.bar(labels, data[:, i], bottom=bottom, color=colors[i % len(colors)], label=conf)
        bottom += data[:, i]

    # Improve axis labels and spacing
    ax.set_xlabel("Chain-Residue", labelpad=15, fontsize=22)
    ax.set_ylabel("Conformer Probability", labelpad=15, fontsize=22)

    # Rotate x-axis labels and adjust alignment
    plt.xticks(rotation=90, ha='right')

    # Add a more descriptive title with better padding
    plt.title(f'Conformer Probability Distribution for Each {target_residue} at pH {pH_value}', pad=20, fontsize=22)

    if legend_box == 'T':
        # Move the legend outside the plot
        plt.legend(title='Protonation States', bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    elif legend_box == 'F':
        # turn off legend
        plt.legend().set_visible(False)
    else:
        pass


    # Ensure layout doesn't overlap and is properly spaced
    plt.tight_layout()

    # Save the figure to a directory named after the residue
    directory = target_residue
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Save the figure as a high-resolution PNG file inside the directory
    filename = f"{directory}/{target_residue}_protonation_chains_{'_'.join(chains_list)}_pH{pH_value}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')  # Save at 300 dpi for publication quality
    print(f"Figure saved to {filename}")

    plt.show()
    plt.clf()

    import pandas as pd


def show_stacked_protonation_table(self, target_residue, chains_list, pH_value):
    """
    Show and save the conformer probabilities for each residue of a specified type across chains at a given pH value in a formatted table.
    
    Parameters:
    - target_residue: str, the residue to display (e.g., "ASP").
    - chains_list: list of str, chains to include in the table.
    - pH_value: float, the pH level at which to show the protonation distribution.
    """
    protonation_sums = {}
    residue_indices = []

    # Regex to capture protonation state correctly, including negative states
    protonation_state_pattern = re.compile(r'^([+-]?\d+)-(\d+)-(.+)$')

    for chain_name in chains_list:
        chain = self.get_chain(chain_name)
        if chain:
            for residue_key, residue in chain.residues.items():
                if residue.residue_name == target_residue:
                    for conformer in residue.conformers:
                        probability_at_ph = conformer.get_probability_at_pH(pH_value)
                        if probability_at_ph is not None:
                            key = f"{residue.protonation_state}-{residue.residue_index}-{chain_name}"
                            protonation_sums[key] = protonation_sums.get(key, 0) + probability_at_ph
                            if (chain_name, residue.residue_index) not in residue_indices:
                                residue_indices.append((chain_name, residue.residue_index))

    # Sort residue indices by chain name and residue index
    residue_indices = sorted(residue_indices, key=lambda x: (x[0], x[1]))

    if not protonation_sums:
        print("No data in protonation_sums, nothing to display.")
        return  # Exit early if there's no data

    # Extract conformers using the regex to handle negative protonation states
    conformers = sorted(list(set([protonation_state_pattern.match(key).group(1) for key in protonation_sums.keys()])))

    # Create a dictionary to accumulate values
    table_data = []
    for chain_name, res_index in residue_indices:
        row = {'Chain': chain_name, 'Residue Index': res_index}
        for conf in conformers:
            row[conf] = protonation_sums.get(f"{conf}-{res_index}-{chain_name}", 0)
        table_data.append(row)

    # Convert to DataFrame
    df = pd.DataFrame(table_data)

    # Format the table for output
    table_str = tabulate(df, headers="keys", tablefmt="pretty", showindex=False)
    
    # Display the table
    print(f"\nProtonation State Probabilities for {target_residue} at pH {pH_value}")
    print(table_str)

    # Prepare filename based on parameters
    chains_str = "_".join(chains_list)
    filename = f"{target_residue}_pH{pH_value}_chains_{chains_str}.txt"
    
    # Save the table to a text file
    with open(filename, "w") as file:
        file.write(f"Protonation State Probabilities for {target_residue} at pH {pH_value}\n")
        file.write(table_str)

    print(f"Table saved to {filename}")




#==========================================================================================================================================
# Now we can integrate the loader with the previously defined hierarchy
# Example usage:
# from pymembrane.util.mcce_neededFiles.ProtonationProfileVisualizer import *

# # Create a protein object
# protein = Protein('MyProtein')

# # Load data from fort.38
# loader = ProtonationProfile('fort.38')
# loader.load_data(protein)

# # # You can now query the protein structure, e.g., retrieve a residue and check its total probability
# # retrieved_chain = protein.get_chain('A')
# # if retrieved_chain:
# #     retrieved_residue = retrieved_chain.get_residue('LYS', '+1', '0096')
# #     if retrieved_residue:
# #         print(f"Total probability for residue ASP01A0003: {retrieved_residue.total_probability}")
# #     else:
# #         print("Residue not found.")
# # else:
# #     print("Chain not found.")


# # 2. Print Structure:
# # Example: Iterate through protein, chains, residues, and conformers
# for chain in protein:
#     print(f"Chain: {chain.chain_name}")
#     for residue in chain:
#         print(f"Residue: {residue.residue_name}{residue.protonation_state}{residue.residue_index}")
#         print(f"  Protonation state: {residue.protonation_state}")
#         print(f"  Total probability: {residue.total_probability}")
#         for conformer in residue:
#             print(f"    {conformer}")

# #3. Plot Protonation Distribution for a Specific Residue:
# # plot_protonation_distribution(protein, "ASP", "0003", ["A", "B", "D"])

# # 4.Visualize Protonation States for a Specific Residue Type:
# # visualize_protonation(protein, "LYS", ["A", "B", "D"])
# # visualize_protonation(protein, "ASP", ["A", "B", "D"])  
# # visualize_protonation(protein, "GLU", ["A", "B", "D"])
# # visualize_protonation(protein, "HIS", ["A", "B", "D"])
# # visualize_protonation(protein, "ARG", ["A", "B", "D"])




# # Plot protonation distribution for a specific residue at pH 6.5
# plot_protonation_distribution(protein, "ASP", "0003", ["A"], pH_value=6.0)

# # Visualize protonation states for LYS across chains at pH 7.0
# visualize_protonation(protein, "LYS", ["A",'B',], pH_value=7.0)
# visualize_protonation(protein, "ASP", ['C','D'], pH_value=3.0)


# # residue_titration_curve_across_pH(protein, "LYS", "0096", ["A", "B", "D"])
# # residue_titration_curve_across_pH(protein, "HIS", "0043", ["A"])

# for ph in range(14):
#     visualize_protonation(protein, "ASP", ['C','D'], pH_value=ph)

# residue_titration_curve_across_pH(protein, "ASP", "0003", ["A"])


# visualize_stacked_protonation(protein, "ASP", chain, pH_value=7.0)
