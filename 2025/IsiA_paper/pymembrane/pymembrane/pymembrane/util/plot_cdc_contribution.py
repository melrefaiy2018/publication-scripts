import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

def create_directories(paths):
    for path in paths:
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            print(f'Directory {path} already exists and may not be empty!')

def save_contributions_to_file(file_path, index, keys, values):
    with open(file_path, 'w') as dict_file:
        for key, value in zip(keys, values): 
            dict_file.write(f'{index}:     {key} : {value} \n')

def plot_residue_contribution(dict_total_contribution, cutoff, main_path):
    residue_contribution_path = os.path.join(main_path, f'residue_contribution_{cutoff}')
    selected_contribution_path = os.path.join(main_path, f'contribution_{cutoff}')
    
    create_directories([residue_contribution_path, selected_contribution_path])

    for index, contributions in dict_total_contribution.items():
        sorted_contributions = dict(sorted(contributions.items(), key=lambda item: item[1]))
        
        selected_keys = [k for k, v in sorted_contributions.items() if abs(v) > cutoff]
        selected_values = [v for k, v in sorted_contributions.items() if abs(v) > cutoff]

        file_path = os.path.join(selected_contribution_path, f'selected_contribution_to__{index}.txt')
        save_contributions_to_file(file_path, index, selected_keys, selected_values)

        if selected_values:
            fig, ax = plt.subplots( tight_layout=True)
            sns.barplot(x=selected_keys, y=selected_values)
            plt.ylabel("Energy Contribution (cm-1)", fontweight='bold')
            plt.xticks(rotation=90)
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=12)
            plt.title(f"Site Energy Shift Contribution of {index}, cutoff={cutoff}", fontweight='bold')
            plt.tight_layout()
            plt.savefig(os.path.join(residue_contribution_path, f'pigment_{index}.png'))
            plt.close()

def parse_interaction_data(interaction_path):
    if not os.path.exists(interaction_path):
        raise FileNotFoundError(f"The file {interaction_path} does not exist.")
    
    interaction_data = []
    with open(interaction_path, "r") as f:
        for line in f:
            parts = line.strip().split(":")
            residue_pair = parts[0].strip()
            residue_a = parts[1].strip()
            interaction_value = float(parts[2])
            interaction_data.append((residue_pair, residue_a, interaction_value))
    
    return interaction_data

def calculate_residue_centers(pdb_path):
    dict_pigment_location = {}
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    residue_centers = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                center_of_mass = np.mean([atom.get_coord() for atom in residue], axis=0)
                residue_centers[f"{chain.id}_{residue.resname}_{residue.id[1]}"] = center_of_mass
    return residue_centers

def plot_circle_value(center_of_mass, residue, interaction_value, color, circle_size, fontsize, show_values, cutoff, shape='o'):
    plt.scatter(center_of_mass[0], center_of_mass[1], c=[color], s=circle_size, alpha=0.4, edgecolors="black", marker=shape)
    # Only display the label if the absolute interaction value is above the cutoff
    if show_values and abs(interaction_value) > cutoff:
        plt.text(center_of_mass[0], center_of_mass[1], f"{residue}\n{round(interaction_value):d}", color="black", 
                 ha="center", va="center", fontsize=fontsize)

def plot_spatial_interaction(residue_centers, interaction_data, interacting_pigment, font_size, show_values, main_path, directory_path, cutoff=0):
    residue_interactions = {entry[1]: entry[2] for entry in interaction_data}

    plt.figure(figsize=(10, 10))
    
    # Plot the interacting pigment location in green with a large size and distinct shape
    if '_'.join(interacting_pigment.split('_')[1:]) in residue_centers:
        pigment = '_'.join(interacting_pigment.split('_')[1:])  # Extracting the relevant part for pymembrane compatibility
        plot_circle_value(residue_centers[pigment], interacting_pigment, 0, 'green', 800, font_size, show_values=False, cutoff=cutoff, shape='h')
    
    # Plot the residues that contribute to the site energy of the interacting pigment
    for residue, center in residue_centers.items():
        # Skip the interacting pigment itself as we have already plotted it above
        if residue == interacting_pigment:
            continue
        
        interaction_value = residue_interactions.get(residue, 0)
        if interaction_value > cutoff:
            # Identify the color and shape based on the type of residue
            if 'NA' in residue:
                # Color Na and Cl ions differently and use a different shape (e.g., square)
                plot_circle_value(center, residue, interaction_value, 'purple', 300, font_size, show_values, cutoff, shape='s')
            elif 'CL' in residue:
                # Color Na and Cl ions differently and use a different shape (e.g., square)
                plot_circle_value(center, residue, interaction_value, 'cyan', 300, font_size, show_values, cutoff, shape='s')

            else:
                # Default color coding based on interaction value
                color = 'blue' if interaction_value > 0 else 'red'
                circle_size = 100 * np.abs(interaction_value)
                plot_circle_value(center, residue, interaction_value, color, circle_size, font_size, show_values, cutoff, shape='o')

    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    plt.title(f"{interacting_pigment}")
    plt.tight_layout()

    save_path_fig = os.path.join(main_path, '2d_interaction_fig')
    save_path_data = os.path.join(main_path, '2d_interaction_data')

    create_directories([save_path_fig, save_path_data])
    plt.savefig(os.path.join(save_path_fig, f'spatial_interaction_{interacting_pigment}.png'), dpi=600)
    plt.close()

    # Save interaction data to a file
    with open(os.path.join(save_path_data, f'residue_interactions_{interacting_pigment}.txt'), 'w') as txt_file:
        txt_file.write(f'Residues contributing to the energy shift of {interacting_pigment}:\n')
        txt_file.write('========================================================\n\n')
        for residue, interaction_value in residue_interactions.items():
            txt_file.write(f"Residue: {residue}, Interaction Value: {interaction_value}\n")



def plot_residue_contribution_and_spatial_interactions(dict_total_contribution, cutoff, pdb_path, save_path, font_size=10, show_values=False):
    main_save_path = os.path.join(save_path, 'CDC_analysis')
    create_directories([main_save_path])
    
    plot_residue_contribution(dict_total_contribution, cutoff, main_save_path)

    for pigment in dict_total_contribution.keys():
        interaction_path = os.path.join(main_save_path, f'contribution_{cutoff}', f'selected_contribution_to__{pigment}.txt')
        interaction_data = parse_interaction_data(interaction_path)
        residue_centers = calculate_residue_centers(pdb_path)
        plot_spatial_interaction(residue_centers, interaction_data, pigment, font_size, show_values, main_save_path, os.path.dirname(pdb_path))

# Example usage:
# plot_residue_contribution_and_spatial_interactions(
#     dict_total_contribution=dict_total_contribution_site_energy, 
#     cutoff=20, 
#     pdb_path=pdb_file, 
#     save_path=dir_path,
#     font_size=10, 
#     show_values=True
# )
