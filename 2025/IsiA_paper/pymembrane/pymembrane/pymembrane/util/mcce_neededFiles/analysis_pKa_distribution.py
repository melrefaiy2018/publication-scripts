#!/usr/bin/env python3
"""
Script to analyze and plot pKa values by residue type from MCCE pK.out files.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from collections import defaultdict
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 14,
    'axes.labelsize': 16,
    'axes.titlesize': 18,
    'legend.fontsize': 12,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14
})
REF_PKA = {
    'ASP': 3.9,
    'GLU': 4.3,
    'HIS': 6.5,
    'LYS': 10.5,
    'ARG': 12.0
}

class pKaAnalyzer:
    def __init__(self, file_path):
        self.file_path = Path(file_path)
        self.data = self.parse_pk_file()

    def parse_pk_file(self):
        data_rows = []

        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.strip() and not line.startswith('pH'):
                parts = line.split()
                if len(parts) >= 2:
                    residue_full = parts[0]
                    pka_raw = parts[1]

                    if len(residue_full) < 9:
                        print(f"Skipping residue with unexpected format: {residue_full}")
                        continue

                    residue_type = residue_full[0:3]
                    chain_name = residue_full[4]
                    residue_number = int(residue_full[5:9])

                    try:
                        if pka_raw == '>14.0':
                            pka_value = 14.0
                            is_extreme = True
                        elif pka_raw == '<0.0':
                            pka_value = 0.0
                            is_extreme = True
                        else:
                            pka_value = float(pka_raw)
                            is_extreme = False

                        data_rows.append({
                            'residue_full': residue_full,
                            'residue_type': residue_type,
                            'chain_name': chain_name,
                            'residue_number': residue_number,
                            'pka_value': pka_value,
                            'is_extreme': is_extreme
                        })
                    except ValueError:
                        print(f"Warning: Could not parse pKa value '{pka_raw}' for residue {residue_full}")

        if not data_rows:
            print("Warning: No valid data parsed from the file.")

        return pd.DataFrame(data_rows)


def parse_residue_index(residue_number):
    residue_str = str(residue_number)
    return int(residue_str[1:]) if len(residue_str) > 1 else residue_number


def generate_protonation_tables_by_residue_type(analyzer, save_path=None):
    if analyzer.data is None:
        raise ValueError("Data not loaded. Call parse_pk_file() first.")

    data_filtered = analyzer.data[~analyzer.data['is_extreme']].copy()
    data_filtered['protonation_state_at_pH7'] = np.where(
        data_filtered['pka_value'] > 7.0, 'Protonated', 'Deprotonated'
    )

    residue_types = data_filtered['residue_type'].unique()
    tables = {}

    for r_type in residue_types:
        subset = data_filtered[data_filtered['residue_type'] == r_type][[
            'residue_full', 'residue_number', 'pka_value', 'protonation_state_at_pH7'
        ]]
        tables[r_type] = subset

        if save_path:
            file_path = Path(save_path) / f'{r_type}_protonation_table.csv'
            subset.to_csv(file_path, index=False)

    return tables


def compare_to_reference(analyzer, save_path=None):
    if analyzer.data is None:
        raise ValueError("Data not loaded. Call parse_pk_file() first.")

    df = analyzer.data[~analyzer.data['is_extreme']].copy()
    df['protonation_state_at_pH7'] = np.where(df['pka_value'] > 7, 'Protonated', 'Deprotonated')
    df['expected_pKa'] = df['residue_type'].map(REF_PKA)
    df['expected_state'] = np.where(df['expected_pKa'] > 7, 'Protonated', 'Deprotonated')
    df['matches_ref'] = df['protonation_state_at_pH7'] == df['expected_state']

    summary = df.groupby('residue_type').agg(
        total=('residue_full', 'count'),
        matches=('matches_ref', 'sum')
    )
    summary['pct_agree'] = summary['matches'] / summary['total'] * 100

    if save_path:
        comparison_file = Path(save_path) / 'protonation_reference_comparison.csv'
        summary.to_csv(comparison_file)

    return df, summary


def plot_pka_vs_residue_number(analyzer,pH_effective, title="pKa Values per Residue Index", save_path=None):
    data = analyzer.data.copy()
    if data.empty:
        print("Warning: No data available in analyzer for plotting.")
        return

    # Filter to only residues defined in REF_PKA
    data = data[data['residue_type'].isin(REF_PKA.keys())]
    if data.empty:
        print("Warning: No residues with defined REF_PKA.")
        return

    data['residue_index'] = data['residue_number']  # Plot actual residue numbers on x-axis

    plt.figure(figsize=(8.27, 11.69))  # A4 portrait size in inches

    markers = {'ASP': 'o', 'GLU': 's', 'HIS': 'D', 'LYS': '^', 'ARG': 'v'}
    palette = {'ASP': '#E64B35', 'GLU': '#4DBBD5', 'HIS': '#00A087', 'LYS': '#3C5488', 'ARG': '#F39B7F'}

    for residue_type in data['residue_type'].unique():
        subset = data[data['residue_type'] == residue_type]
        plt.scatter(
            subset['residue_index'],
            subset['pka_value'],
            label=residue_type,
            marker=markers.get(residue_type, 'o'),
            color=palette.get(residue_type, 'gray'),
            s=100,
            edgecolor='black',
            linewidth=0.5
        )
    min_pH_effective = pH_effective - 1
    max_pH_effective = pH_effective + 1
    plt.axhspan(min_pH_effective, max_pH_effective, color='gray', alpha=0.2)
    labile = data[(data['pka_value'] >= min_pH_effective) & (data['pka_value'] <= max_pH_effective)]
    for _, row in labile.iterrows():
        plt.text(row['residue_index'], row['pka_value'] + 0.3, row['residue_full'], fontsize=10, ha='center', va='bottom')

    plt.xlabel('Residue Index')
    plt.ylabel('pKa Value')
    plt.title(title, weight='bold')
    plt.legend(title='Residue Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    sns.despine()
    plt.tight_layout()

    if save_path:
        plt.savefig(Path(save_path) / 'enhanced_pka_vs_resindex_A4.png', dpi=600, bbox_inches='tight')
        plt.savefig(Path(save_path) / 'enhanced_pka_vs_resindex_A4.svg', format='svg', dpi=600, bbox_inches='tight')

    plt.show()


def generate_pka_summary_table(analyzer, save_path=None):
    data = analyzer.data.copy()
    if data.empty:
        print("Warning: No data available for summary table.")
        return pd.DataFrame()

    # Filter to only residues defined in REF_PKA
    data = data[data['residue_type'].isin(REF_PKA.keys())]
    if data.empty:
        print("Warning: No residues with defined REF_PKA.")
        return pd.DataFrame()

    summary_table = data[['residue_full', 'residue_type', 'chain_name', 'residue_number', 'pka_value']].copy()
    summary_table = summary_table.sort_values(by='pka_value')

    if save_path:
        summary_file = Path(save_path) / 'pka_summary_table.csv'
        summary_table.to_csv(summary_file, index=False)
        print(f"Summary table saved to {summary_file}")

    return summary_table


def get_residues_around_pH_effective(analyzer, pH_effective=7.0, tolerance=1.0, save_path=None):
    """
    Find and save residues with pKa values around the effective pH.
    
    Parameters:
    - analyzer: pKaAnalyzer object
    - pH_effective: The effective pH value (default 7.0)
    - tolerance: The range around pH_effective to consider (default 1.0)
    - save_path: Directory to save the results
    
    Returns:
    - DataFrame of residues around pH_effective
    """
    data = analyzer.data.copy()
    if data.empty:
        print("Warning: No data available in analyzer.")
        return pd.DataFrame()

    # Filter to only residues defined in REF_PKA and non-extreme values
    data = data[data['residue_type'].isin(REF_PKA.keys()) & ~data['is_extreme']]
    if data.empty:
        print("Warning: No residues with defined REF_PKA.")
        return pd.DataFrame()

    # Find residues around pH_effective
    min_pH = pH_effective - tolerance
    max_pH = pH_effective + tolerance
    
    labile_residues = data[(data['pka_value'] >= min_pH) & (data['pka_value'] <= max_pH)].copy()
    
    # Sort by pKa value
    labile_residues = labile_residues.sort_values(by='pka_value')
    
    # Add additional information
    labile_residues['pH_effective'] = pH_effective
    labile_residues['tolerance'] = tolerance
    labile_residues['distance_from_pH_eff'] = abs(labile_residues['pka_value'] - pH_effective)
    
    print(f"\nResidues around pH_effective = {pH_effective} (±{tolerance}):")
    print(f"Found {len(labile_residues)} residues in the range [{min_pH:.1f}, {max_pH:.1f}]")
    print("=" * 80)
    
    if not labile_residues.empty:
        # Print residues categorized by type and sorted by residue index
        residue_types_order = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG']
        
        for residue_type in residue_types_order:
            type_residues = labile_residues[labile_residues['residue_type'] == residue_type]
            if not type_residues.empty:
                # Sort by residue number within each type
                type_residues = type_residues.sort_values(by='residue_number')
                
                print(f"\n{residue_type} residues ({len(type_residues)} found):")
                print("-" * 50)
                for _, row in type_residues.iterrows():
                    print(f"  {row['residue_full']:12} | pKa: {row['pka_value']:6.2f} | Distance from pH_eff: {row['distance_from_pH_eff']:5.2f}")
    
    if save_path and not labile_residues.empty:
        os.makedirs(save_path, exist_ok=True)
        output_file = Path(save_path) / f'labile_residues_pH{pH_effective:.1f}_tol{tolerance:.1f}.csv'
        # Sort the dataframe by residue type and then by residue number before saving
        labile_residues_sorted = labile_residues.copy()
        labile_residues_sorted['residue_type_order'] = labile_residues_sorted['residue_type'].map({
            'ASP': 1, 'GLU': 2, 'HIS': 3, 'LYS': 4, 'ARG': 5
        })
        labile_residues_sorted = labile_residues_sorted.sort_values(['residue_type_order', 'residue_number'])
        labile_residues_sorted = labile_residues_sorted.drop('residue_type_order', axis=1)
        labile_residues_sorted.to_csv(output_file, index=False)
        print(f"\nLabile residues saved to: {output_file}")
    
    return labile_residues


def save_categorized_residues_to_text(analyzer, pH_effective=7.0, tolerance=1.0, save_path=None):
    """
    Save the categorized residue output to a text file with the same format as console output.
    
    Parameters:
    - analyzer: pKaAnalyzer object
    - pH_effective: The effective pH value (default 7.0)
    - tolerance: The range around pH_effective to consider (default 1.0)
    - save_path: Directory to save the results
    
    Returns:
    - Path to the saved text file
    """
    data = analyzer.data.copy()
    if data.empty:
        print("Warning: No data available in analyzer.")
        return None
    
    # Filter to only residues defined in REF_PKA and non-extreme values
    data = data[data['residue_type'].isin(REF_PKA.keys()) & ~data['is_extreme']]
    if data.empty:
        print("Warning: No residues with defined REF_PKA.")
        return None
    
    # Find residues around pH_effective
    min_pH = pH_effective - tolerance
    max_pH = pH_effective + tolerance
    
    labile_residues = data[(data['pka_value'] >= min_pH) & (data['pka_value'] <= max_pH)].copy()
    
    if labile_residues.empty:
        print("No residues found around the effective pH.")
        return None
    
    if not save_path:
        print("No save path provided.")
        return None
    
    # Sort by pKa value
    labile_residues = labile_residues.sort_values(by='pka_value')
    
    # Add additional information
    labile_residues['distance_from_pH_eff'] = abs(labile_residues['pka_value'] - pH_effective)
    
    # Create the output directory if it doesn't exist
    os.makedirs(save_path, exist_ok=True)
    
    # Create the text file
    output_file = Path(save_path) / f'labile_residues_categorized_pH{pH_effective:.1f}_tol{tolerance:.1f}.txt'
    
    with open(output_file, 'w') as f:
        f.write(f"Residues around pH_effective = {pH_effective} (±{tolerance}):\n")
        f.write(f"Found {len(labile_residues)} residues in the range [{min_pH:.1f}, {max_pH:.1f}]\n")
        f.write("=" * 80 + "\n")
        
        if not labile_residues.empty:
            # Print residues categorized by type and sorted by residue index
            residue_types_order = ['ASP', 'GLU', 'HIS', 'LYS', 'ARG']
            
            for residue_type in residue_types_order:
                type_residues = labile_residues[labile_residues['residue_type'] == residue_type]
                if not type_residues.empty:
                    # Sort by residue number within each type
                    type_residues = type_residues.sort_values(by='residue_number')
                    
                    f.write(f"\n{residue_type} residues ({len(type_residues)} found):\n")
                    f.write("-" * 50 + "\n")
                    for _, row in type_residues.iterrows():
                        f.write(f"  {row['residue_full']:12} | pKa: {row['pka_value']:6.2f} | Distance from pH_eff: {row['distance_from_pH_eff']:5.2f}\n")
        
        f.write(f"\n\nFile generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Analysis parameters: pH_effective = {pH_effective}, tolerance = ±{tolerance}\n")
        f.write(f"Input file: {analyzer.file_path}\n")
    
    print(f"Categorized residue output saved to: {output_file}")
    return output_file

# # Example usage:
# # CP24:
# # ===============================================================================
# input_file = "/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP24/CP24_LYS_gold/fit_HIS_183/pK.out"
# save_dir = "/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP24/CP24_LYS_gold/fit_HIS_183/output_dir"
# title = "pKa Values per Ionizable Residue (CP24 Protein)"
# analyzer = pKaAnalyzer(input_file)
# analyzer.parse_pk_file()
# tables = generate_protonation_tables_by_residue_type(analyzer, save_path=save_dir)
# df, summary = compare_to_reference(analyzer, save_path=save_dir)
# plot_pka_vs_residue_number(analyzer,title, save_path=save_dir)
# generate_pka_summary_table(analyzer, save_path=save_dir)
# # Get and save residues around pH effective (8.45 for CP24)
# labile_residues = get_residues_around_pH_effective(analyzer, pH_effective=8.45, tolerance=1.0, save_path=save_dir)

# # # CP26:
# # # ===============================================================================
# input_file = "/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP26/gold_GLU/pK.out"
# save_dir = "/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP26/gold_GLU/output_dir"

# title = "pKa Values per Ionizable Residue (CP26 Protein)"
# analyzer = pKaAnalyzer(input_file)
# analyzer.parse_pk_file()
# tables = generate_protonation_tables_by_residue_type(analyzer, save_path=save_dir)
# df, summary = compare_to_reference(analyzer, save_path=save_dir)
# plot_pka_vs_residue_number(analyzer,title, pH_effective=5.3, save_path=save_dir)
# generate_pka_summary_table(analyzer, save_path=save_dir)

# # Get and save residues around pH effective (5.3 for CP26)
# labile_residues = get_residues_around_pH_effective(analyzer, pH_effective=5.3, tolerance=1.0, save_path=save_dir)
# # Save the categorized output to a text file
# text_output_file = save_categorized_residues_to_text(analyzer, pH_effective=5.3, tolerance=1.0, save_path=save_dir)

# Display the results in a formatted table if residues are found
# # ===============================================================================
# # CP29:
# # ===============================================================================
# # input_file = "/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP29/mcce/pK.out"
# # save_dir = "/Users/mohamed/Documents/Research/Projects/SuperComplex_Kajwal/correct_calculation/CP29/mcce/output_dir"
# # os.makedirs(save_dir, exist_ok=True)
# # title = "pKa Values per Ionizable Residue (CP29 Protein)"
# # analyzer = pKaAnalyzer(input_file)
# # analyzer.parse_pk_file()
# # tables = generate_protonation_tables_by_residue_type(analyzer, save_path=save_dir)
# # df, summary = compare_to_reference(analyzer, save_path=save_dir)
# # plot_pka_vs_residue_number(analyzer,title, save_path=save_dir)
# # generate_pka_summary_table(analyzer, save_path=save_dir)
# Get and save residues around pH effective (5.3 for CP26)
# labile_residues = get_residues_around_pH_effective(analyzer, pH_effective=8.4, tolerance=1.0, save_path=save_dir)

# # ===============================================================================