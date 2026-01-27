"""
Author: Mohamed A. A. Elrefaiy
Date: November 2024
===========================================================================
IsiA Monomer Hamiltonian and Spectra Calculation:
Calculates excitonic Hamiltonian, absorption, and fluorescence spectra
for the IsiA protein complex at pH 7.
"""

import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)
from pymembrane.structure.atomic_pigment import ChlorophyllAtomic, PheophytinAtomic
from pymembrane.parameters.thylakoid.plants.hamiltonian.pigment_renger_hfcis import coupling_data_by_type_mcce as coupling_renger_by_type_mcce
from pymembrane.parameters.thylakoid.plants.hamiltonian.pigment_renger_hfcis import q0011_charges_by_type_mcce as q0011_renger_by_type_mcce
from pymembrane.parameters.thylakoid.plants.hamiltonian.pigment_renger_lineshape import kubo_renger_cla
from pymembrane.structure.atomic_protein import ElectrostaticProteinAtomic
from pymembrane.parameters.thylakoid.protein_default import define_domain
from pymembrane.exciton.linear_spectra import LinearSpectra
from pymembrane.util.plot_cdc_contribution import plot_residue_contribution_and_spatial_interactions
from pymembrane.util.helper_functions import *
import numpy as np
from Bio import PDB
from scipy.stats import norm
from pymembrane.util.physical_constants import hbar, c
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
import json
from pathlib import Path
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import matplotlib as mpl
import seaborn as sns
import re
import logging
from datetime import datetime

# ============================================================================
# PATH SETUP
# ============================================================================

# Get the script's location and navigate to project root
script_file = Path(__file__).resolve()
hamiltonian_dir = script_file.parent           # scripts/hamiltonian/
scripts_dir = hamiltonian_dir.parent           # scripts/
project_root = scripts_dir.parent              # IsiA_paper/

# Define input paths
pdb_file = project_root / 'NeededData' / 'structure' / 'extended_most_occ_pH7.pdb'

# Define output directory - all results in one place
results_base_dir = hamiltonian_dir / 'results' / 'hamiltonian_results'

# Create subdirectories for this script's outputs
site_energy_dir = results_base_dir / 'SiteEnergy_Data'
hamiltonian_dir_out = results_base_dir / 'Hamiltonian'
spectra_data_dir = results_base_dir / 'Spectra_Data'
spectra_figures_dir = results_base_dir / 'Spectra_figures'
domain_viz_dir = results_base_dir / 'Domain_visualization'
pka_analysis_dir = results_base_dir / 'pKa_analysis'

# Legacy aliases for compatibility
data_dir = results_base_dir
results_dir = project_root / 'results'
# Remove figures_dir alias to avoid confusion - use specific subdirectories instead

# Verify PDB file exists
if not pdb_file.exists():
    raise FileNotFoundError(
        f"PDB file not found at: {pdb_file}\n"
        f"Please ensure 'extended_most_occ_pH7.pdb' is in the structure/ directory"
    )

# Convert to string for compatibility
pdb_file_str = str(pdb_file)

# ============================================================================
# SETUP LOGGING
# ============================================================================

# Create the results directory first
create_directory(results_base_dir)

# Generate timestamp for log file
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = results_base_dir / f"calculation_log_{timestamp}.log"

# Configure logging to write to both file and console
# Create custom logger to avoid capturing matplotlib/other library logs
logger = logging.getLogger('IsiA_Calculation')
logger.setLevel(logging.INFO)

# File handler
file_handler = logging.FileHandler(log_file)
file_handler.setLevel(logging.INFO)
file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(file_formatter)

# Console handler
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)

# Add handlers to logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# Suppress verbose logging from matplotlib and other libraries
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)

# Print status
logger.info("="*80)
logger.info("IsiA Monomer Hamiltonian and Spectra Calculation")
logger.info("="*80)
logger.info(f"Project root: {project_root}")
logger.info(f"PDB file: {pdb_file.name}")
logger.info(f"Results output: {results_base_dir}")
logger.info(f"Log file: {log_file.name}")
logger.info("="*80)
logger.info("")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def create_directory(path):
    """Create directory if it doesn't exist"""
    Path(path).mkdir(parents=True, exist_ok=True)

def save_dict_to_json(dict_to_save, file_path):
    """
    Save dictionary to a JSON file, handling NumPy arrays and other scientific data types.
    """
    class NumpyEncoder(json.JSONEncoder):
        def default(self, obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            if isinstance(obj, np.integer):
                return int(obj)
            if isinstance(obj, np.floating):
                return float(obj)
            if isinstance(obj, np.bool_):
                return bool(obj)
            if isinstance(obj, complex):
                return {"real": obj.real, "imag": obj.imag}
            if hasattr(obj, '__dict__'):
                return obj.__dict__
            return super().default(obj)
    
    try:
        with open(file_path, 'w') as f:
            json.dump(dict_to_save, f, cls=NumpyEncoder, indent=2)
        logger.info(f"Successfully saved data to: {file_path}")
    except Exception as e:
        logger.error(f"Error saving to {file_path}: {str(e)}")

def load_dict_from_json(file_path):
    """Load dictionary from JSON file."""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        logger.info(f"Successfully loaded data from: {file_path}")
        return data
    except Exception as e:
        logger.error(f"Error loading from {file_path}: {str(e)}")
        return None

# ============================================================================
# BUILD ATOMIC PROTEIN
# ============================================================================

logger.info("Loading PDB structure and building atomic protein...")
pdb_atomic = PDB.PDBParser().get_structure('Isia', pdb_file_str)
IsiA_atomic = ElectrostaticProteinAtomic(pdb_atomic, path_extended_pdb=pdb_file_str, name='Isia')

# Define protein parameters
IsiA_atomic.fval_tresp = 0.72
IsiA_atomic.dielectric_cdc = 2
IsiA_atomic.dielectric_tresp = 1

# Validate required parameters are set
required_params = {
    'fval_tresp': IsiA_atomic.fval_tresp,
    'dielectric_cdc': IsiA_atomic.dielectric_cdc,
    'dielectric_tresp': IsiA_atomic.dielectric_tresp
}

for param_name, param_value in required_params.items():
    if param_value is None or not hasattr(IsiA_atomic, param_name):
        raise ValueError(
            f"Required parameter '{param_name}' is not defined. "
            f"Please set IsiA_atomic.{param_name} before line 160."
        )

logger.info(f"✓ Parameters validated: fval_tresp={IsiA_atomic.fval_tresp}, "
      f"dielectric_cdc={IsiA_atomic.dielectric_cdc}, "
      f"dielectric_tresp={IsiA_atomic.dielectric_tresp}")

IsiA_atomic.prepare_pigments('CLA', ChlorophyllAtomic,
                        list_chain=None,
                        lineshape=kubo_renger_cla,
                        dict_coupling_data=coupling_renger_by_type_mcce['CLA_IPPC'],
                        dict_q0011_charges=q0011_renger_by_type_mcce['CLA'],
                        disorder=70)
# ============================================================================
# CALCULATE SITE ENERGY SHIFT
# ============================================================================

logger.info("\nCalculating site energy shifts...")
e0_a = 14900
e0_b = 0 
N_ens = 100
temp = 300

dict_site_energy_shift, dict_total_contribution_site_energy = IsiA_atomic.calculate_cdc_site_shift(verbose=False)
dict_site_energy_raw = IsiA_atomic.calculate_total_site_energy(dict_site_energy_shift, E_0a=e0_a, E_0b=e0_b)
logger.info(dict_site_energy_raw)
dict_site_raw_sorted = sort_dict_site_energy(dict_site_energy_raw)
logger.info(dict_site_raw_sorted)

# ============================================================================
# SAVE SITE ENERGY DATA
# ============================================================================

logger.info("\nSaving site energy data...")
create_directory(site_energy_dir)

# Define base filename
base_filename = f'Isia_chain_A_dielc_{IsiA_atomic.dielectric_cdc}_a{e0_a}_b{e0_b}_N_ens_{N_ens}'

# Save site energy shift data
save_path_shift = site_energy_dir / f'{base_filename}_shift.json'
save_dict_to_json(dict_site_energy_shift, save_path_shift)

# Save raw site energy data
save_path_raw = site_energy_dir / f'{base_filename}_raw.json'
save_dict_to_json(dict_site_energy_raw, save_path_raw)

# Save total contribution data
save_path_contrib = site_energy_dir / f'{base_filename}_contrib.json'
save_dict_to_json(dict_total_contribution_site_energy, save_path_contrib)

# Save metadata
metadata = {
    "parameters": {
        "dielectric_cdc": IsiA_atomic.dielectric_cdc,
        "E_0a": e0_a,
        "E_0b": e0_b,
        "N_ens": N_ens,
        "chain": "A"
    },
    "data_description": {
        "shift": "Site energy shifts from CDC calculations",
        "raw": "Total site energies including vacuum contributions",
        "contrib": "Total contribution to site energy calculations"
    },
    "units": {
        "energy": "cm^-1",
        "note": "All energies in wavenumber units"
    }
}

save_path_meta = site_energy_dir / f'{base_filename}_metadata.json'
save_dict_to_json(metadata, save_path_meta)

# Save summary statistics
energies = list(dict_site_energy_raw.values())
summary_stats = {
    "n_sites": len(energies),
    "min_energy": float(min(energies)),
    "max_energy": float(max(energies)),
    "mean_energy": float(np.mean(energies)),
    "std_energy": float(np.std(energies)),
    "energy_range": float(max(energies) - min(energies))
}

save_path_summary = site_energy_dir / f'{base_filename}_summary.json'
save_dict_to_json(summary_stats, save_path_summary)

# Print summary
logger.info("\n" + "="*60)
logger.info("SITE ENERGY DATA SAVED")
logger.info("="*60)
logger.info(f"Directory: {site_energy_dir}")
logger.info(f"Base filename: {base_filename}")
logger.info("\nFiles saved:")
logger.info("  - *_shift.json     : Site energy shifts")
logger.info("  - *_raw.json       : Raw total site energies") 
logger.info("  - *_contrib.json   : Total contribution data")
logger.info("  - *_metadata.json  : Calculation parameters and descriptions")
logger.info("  - *_summary.json   : Summary statistics")

logger.info(f"\nSummary statistics:")
logger.info(f"  Sites: {summary_stats['n_sites']}")
logger.info(f"  Energy range: {summary_stats['energy_range']:.1f} cm⁻¹")
logger.info(f"  Mean ± SD: {summary_stats['mean_energy']:.1f} ± {summary_stats['std_energy']:.1f} cm⁻¹")

# Verification
logger.info("\n" + "="*60)
logger.info("VERIFICATION")
logger.info("="*60)
loaded_data = load_dict_from_json(save_path_raw)
if loaded_data:
    logger.info(f"Verification successful: {len(loaded_data)} sites loaded")
    logger.info("Sample data:")
    for i, (key, value) in enumerate(loaded_data.items()):
        if i < 3:
            logger.info(f"  {key}: {value}")
        elif i == 3:
            logger.info("  ...")
            break

# ============================================================================
# PLOT RESIDUE CONTRIBUTION
# ============================================================================

logger.info("\nPlotting residue contribution and spatial interactions...")
# Note: This section is commented out - uncomment if needed
# create_directory(spectra_figures_dir)
# plot_residue_contribution_and_spatial_interactions(
#     dict_total_contribution=dict_total_contribution_site_energy,
#     cutoff=5,
#     pdb_path=pdb_file_str,
#     save_path=str(spectra_figures_dir),
#     font_size=10,
#     show_values=True)

# ============================================================================
# CONSTRUCT AND SAVE HAMILTONIAN
# ============================================================================

logger.info("\nConstructing Hamiltonian...")
create_directory(hamiltonian_dir_out)

df_hamiltonian = IsiA_atomic.construct_hamiltonian(e0_a=e0_a, e0_b=e0_b)
path_hamiltonian = hamiltonian_dir_out / f'raw_IsiA_chain_A_dielc_{IsiA_atomic.dielectric_cdc}_a{e0_a}_b{e0_b}.csv'
IsiA_atomic._save_hamiltonian_to_csv(df_hamiltonian, str(path_hamiltonian))
df_hamiltonian = pd.read_csv(path_hamiltonian, index_col=0)
# ============================================================================
# BUILD HAMILTONIAN DOMAINS
# ============================================================================

logger.info("\nBuilding Hamiltonian domains...")
IsiA_atomic.load_hamiltonian(str(path_hamiltonian))
list_pigments_by_domain = IsiA_atomic.build_hamiltonian_domains('coupling', coupling_cutoff=20.0)

# Save list_pigments_by_domain as text file
domain_list_file = hamiltonian_dir_out / 'list_pigments_by_domain_coupling.txt'
with open(domain_list_file, 'w') as f:
    for domain_idx, pigments in enumerate(list_pigments_by_domain):
        f.write(f"Domain: {domain_idx}\n")
        f.write("Pigments:\n")
        for pigment in pigments:
            f.write(f"  - {pigment}\n")
        f.write("\n")

IsiA_atomic.define_domains(list_pigments_by_domain)

# ============================================================================
# CALCULATE SPECTRA
# ============================================================================

logger.info("\nCalculating absorption and fluorescence spectra...")
logger.info(f"Running {N_ens} ensemble calculations...")

for index in tqdm(np.arange(N_ens)):
    IsiA_atomic.H2_hamiltonian.add_disorder(index)
    linear_isia = LinearSpectra(IsiA_atomic)
    if index == 0:
        w_axis, final_abs_isia = linear_isia.calc_absorption(temp)
        w_axis, final_flo_isia = linear_isia.calc_fluorescence(temp)
    else:
        _, abs_isia = linear_isia.calc_absorption(temp)
        _, flo_isia = linear_isia.calc_fluorescence(temp)
        final_abs_isia += abs_isia
        final_flo_isia += flo_isia

lambda_axis_a = c / ((w_axis) / (2 * np.pi * hbar))
lambda_axis_f = c / ((w_axis) / (2 * np.pi * hbar))

A_axis, A_intensity = get_sliced_data(final_abs_isia, lambda_axis_a, 100)
F_axis, F_intensity = get_sliced_data(final_flo_isia, lambda_axis_f, 100)

A_normalized = A_intensity / np.max(A_intensity)
F_normalized = F_intensity / np.max(F_intensity)

# ============================================================================
# SAVE SPECTRA DATA
# ============================================================================

logger.info("\nSaving spectra data...")
create_directory(spectra_data_dir)

absorption_file = spectra_data_dir / f'A_normalized_{temp}K_N_ens_{N_ens}.npy'
fluorescence_file = spectra_data_dir / f'F_normalized_{temp}K_N_ens_{N_ens}.npy'

np.save(absorption_file, np.column_stack((A_axis, A_normalized)))
np.save(fluorescence_file, np.column_stack((F_axis, F_normalized)))

# Also save as CSV for easier access
absorption_csv = spectra_data_dir / 'absorption_spectrum.csv'
fluorescence_csv = spectra_data_dir / 'fluorescence_spectrum.csv'

np.savetxt(absorption_csv, np.column_stack((A_axis, A_normalized)), 
           delimiter=',', header='Wavelength (nm), Normalized Intensity')
np.savetxt(fluorescence_csv, np.column_stack((F_axis, F_normalized)), 
           delimiter=',', header='Wavelength (nm), Normalized Intensity')
df_hamiltonian.to_csv(path_hamiltonian)

logger.info(f"  → Absorption: {absorption_csv.name}")
logger.info(f"  → Fluorescence: {fluorescence_csv.name}")
logger.info(f"  → Hamiltonian: {path_hamiltonian.name}")
# ============================================================================
# LOAD EXPERIMENTAL DATA (OPTIONAL - Comment out if not available)
# ============================================================================

logger.info("\nLoading experimental data...")
# Experimental data path within project structure
exp_path = project_root / 'NeededData' / 'Experimental_Data'

try:
    Fl_exp = np.load(exp_path / 'Fl_IsiA_monomer_300k_Gabriela.npy')
    Abs_exp = np.load(exp_path / 'Abs_IsiA_monomer_300K_Gabriela.npy')
    has_exp_data = True
    logger.info("  ✓ Experimental data loaded")
except FileNotFoundError:
    logger.info("  ⚠ Experimental data not found - will plot simulation only")
    has_exp_data = False

# ============================================================================
# PLOT SPECTRA (Publication Quality)
# ============================================================================

logger.info("\nGenerating publication-quality figures...")
create_directory(spectra_figures_dir)

# Publication-style settings
mpl.rcParams.update({
    'font.family': 'Arial',
    'font.size': 12,
    'axes.linewidth': 1.5,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.minor.width': 1.0,
    'ytick.minor.width': 1.0,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.minor.size': 3,
    'ytick.minor.size': 3,
    'legend.fontsize': 12,
    'legend.frameon': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'svg.fonttype': 'none',
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

exp_color = '#2E3440'
sim_color = "#268D19"

# ===============================
# ====== Plot Fluorescence ======
# ===============================
fig, ax = plt.subplots(figsize=(4.5, 3.5))

if has_exp_data:
    ax.plot(Fl_exp[:, 0], Fl_exp[:, 1] / np.max(Fl_exp[:, 1]), color=exp_color, label='Experimental', linewidth=4, zorder=3)

ax.plot(F_axis, F_normalized, color=sim_color, label='Simulation', linewidth=4, alpha=0.9, zorder=2)

ax.set_xlim(640, 750)
ax.set_ylim(0, 1.05)
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.set_xlabel('Wavelength (nm)', fontweight='normal')
ax.set_ylabel('Normalized fluorescence', fontweight='normal')
ax.legend(loc='upper right', frameon=False, handlelength=1.5)
ax.tick_params(direction='in', which='both')

plt.tight_layout(pad=0.3)
plt.savefig(spectra_figures_dir / f"fluorescence_{temp}K_N_ens_{N_ens}.pdf", 
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()
plt.close()

# =============================
# ====== Plot Absorption ======
# =============================
fig, ax = plt.subplots(figsize=(4.5, 3.5))

if has_exp_data:
    ax.plot(Abs_exp[:, 0], Abs_exp[:, 1] / np.max(Abs_exp[:, 1]), 
            color=exp_color, label='Experimental', linewidth=4, zorder=3)

ax.plot(A_axis, A_normalized, 
        color=sim_color, label='Simulation', linewidth=4, 
        alpha=0.9, zorder=2)

ax.set_xlim(600, 730)
ax.set_ylim(0, 1.05)
ax.xaxis.set_major_locator(MultipleLocator(25))
ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.yaxis.set_major_locator(MultipleLocator(0.5))
ax.yaxis.set_minor_locator(MultipleLocator(0.1))
ax.set_xlabel('Wavelength (nm)', fontweight='normal')
ax.set_ylabel('Normalized absorption', fontweight='normal')
ax.legend(loc='upper right', frameon=False, handlelength=1.5)
ax.tick_params(direction='in', which='both')

plt.tight_layout(pad=0.3)
plt.savefig(spectra_figures_dir / f"absorption_{temp}K_N_ens_{N_ens}.pdf", 
            dpi=300, bbox_inches='tight', transparent=True)
plt.show()
plt.close()

logger.info(f"  ✓ Spectra figures saved to: {spectra_figures_dir}")

# ============================================================================
# CALCULATE AND PLOT EXCITON DISTRIBUTION
# ============================================================================

logger.info("\nCalculating exciton distribution...")
from pymembrane.util.calculate_exciton_distribution import *

H_data = pd.read_csv(path_hamiltonian, index_col=0)

exciton_distribution_by_site, site_labels = calculate_exciton_distribution(
    H_data=H_data,
    list_pigment_domains=IsiA_atomic.list_pigments_by_domain,
    protein_atomic=IsiA_atomic,
    calc_type='abs',
    N_ens=N_ens,
    sigma_e=70 # Broadening parameter in cm^-1
)

logger.info("  ✓ Exciton distribution calculated")
logger.info(f"  → Plotting exciton distribution to: {spectra_figures_dir}")

plot_exciton_distribution(
    exciton_distribution_by_site,
    site_labels,
    str(spectra_figures_dir),
    x_min=650,
    x_max=695,
    x_interval=10
)

logger.info("  ✓ Exciton distribution plot saved")

if has_exp_data:
    logger.info("  → Plotting combined absorption and exciton distribution")
    plot_combined_absorption_and_exciton(
        sliced_lambda_axis_a=A_axis,  
        absorption_data=A_intensity,  
        exp_absorption=Abs_exp,
        exciton_distribution=exciton_distribution_by_site,
        list_site_label=site_labels,
        color_mapping=False,
        show_labels=True,
        temp=temp,
        path_saving=str(spectra_figures_dir),
        x_min=650,
        x_max=695,
        fontsize=20
    )
    logger.info("  ✓ Combined absorption and exciton plot saved")

# ============================================================================
# PLOT ENERGETIC DOMAINS
# ============================================================================

logger.info("\nPlotting energetic domains...")
# Save domain energy plot to spectra_figures_dir
domain_energy_path = spectra_figures_dir / "energetic_domains.png"
plot_domain_energy(linear_protein=linear_isia, fig_size=(30,15), save_path=str(domain_energy_path))
logger.info(f"  ✓ Energetic domains plot saved to: {domain_energy_path}")

# ============================================================================
# VISUALIZE DOMAINS
# ============================================================================

logger.info("\nVisualizing pigment networks...")
dict_domain_info, dict_pigment_info = build_domains_data(IsiA_atomic, list_pigments_by_domain)

create_directory(domain_viz_dir)

network_viz_path = domain_viz_dir / 'Pigments_Network_Visualization.png'
visualize_pigments_network(
    list_pigments_by_domain, 
    dict_pigment_info, 
    str(network_viz_path)
)
logger.info(f"  ✓ Network visualization saved to: {network_viz_path}")

pyvis_network_path = domain_viz_dir / 'pyvis_pigments_network.html'
visualize_pyvis_pigments_network(
    list_pigments_by_domain, 
    dict_pigment_info, 
    str(pyvis_network_path)
)
logger.info(f"  ✓ Interactive network visualization saved to: {pyvis_network_path}")

# ============================================================================
# PLOT SITE ENERGY SHIFT (Gradient Bar Chart)
# ============================================================================

logger.info("\nPlotting site energy shift visualization...")

# Sort by shift value
sorted_items = sorted(dict_site_energy_shift.items(), key=lambda x: x[1])
labels_full, shifts = zip(*sorted_items)
labels = [re.search(r'\d+$', name).group() for name in labels_full]

# Set color normalization
vmin, vmax = min(shifts), max(shifts)
norm = mpl.colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
cmap = plt.cm.seismic_r
colors = [cmap(norm(s)) for s in shifts]

# Plot
fig, ax = plt.subplots(figsize=(11, 5))
bars = ax.bar(labels, shifts, color=colors, edgecolor='black', width=0.65)

# Annotate values
for bar, shift in zip(bars, shifts):
    offset = 3 if shift >= 0 else -10
    va = 'bottom' if shift >= 0 else 'top'
    ax.text(bar.get_x() + bar.get_width() / 2, shift + offset, f"{shift:.1f}",
            ha='center', va=va, fontsize=11)

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, pad=0.02)
cbar.set_label('Site Energy Shift (cm$^{-1}$)', fontsize=12)

# Style
ax.axhline(0, color='black', linewidth=1)
ax.set_ylabel("Site Energy Shift (cm$^{-1}$)")
ax.set_title("IsiA Chlorophyll Site Energy Shifts", pad=10)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, rotation=90)
ax.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()

site_energy_shift_path = spectra_figures_dir / "IsiA_SiteEnergy_Shift_Gradient.pdf"
plt.savefig(site_energy_shift_path, dpi=300, bbox_inches='tight', transparent=True)
logger.info(f"  ✓ Site energy shift plot saved to: {site_energy_shift_path}")
plt.show()
plt.close()

# ============================================================================
# OPTIONAL: pKa ANALYSIS (Comment out if not needed)
# ============================================================================

# Uncomment if you have pKa analysis data
logger.info("\nAnalyzing pKa distribution...")
from pymembrane.util.mcce_neededFiles.analysis_pKa_distribution import *

input_file = project_root / "NeededData" / "mcce" / "pK.out"
create_directory(pka_analysis_dir)

if input_file.exists():
    analyzer = pKaAnalyzer(str(input_file))
    analyzer.parse_pk_file()
    tables = generate_protonation_tables_by_residue_type(analyzer, save_path=str(pka_analysis_dir))
    df, summary = compare_to_reference(analyzer, save_path=str(pka_analysis_dir))
    plot_pka_vs_residue_number(analyzer, pH_effective=6.7, title="pKa Values per Ionizable Residue (IsiA Protein)", 
                               save_path=str(pka_analysis_dir))
    generate_pka_summary_table(analyzer, save_path=str(pka_analysis_dir))
    labile_residues = get_residues_around_pH_effective(analyzer, pH_effective=6.7, 
                                                       tolerance=1.0, save_path=str(pka_analysis_dir))
else:
    logger.info("  ⚠ pKa input file not found, skipping pKa analysis")

# ============================================================================
# COMPLETION MESSAGE
# ============================================================================

logger.info("\n" + "="*80)
logger.info("✓ CALCULATION COMPLETE!")
logger.info("="*80)
logger.info(f"All results saved to: {results_base_dir}")
logger.info(f"Log file saved to: {log_file}")
logger.info("\nKey outputs:")
logger.info(f"  - Hamiltonian data: {hamiltonian_dir_out}")
logger.info(f"  - Site energy data: {site_energy_dir}")
logger.info(f"  - Spectra data: {spectra_data_dir}")
logger.info(f"  - Figures: {spectra_figures_dir}")
logger.info(f"  - Domain visualization: {domain_viz_dir}")
logger.info(f"  - pKa analysis: {pka_analysis_dir}")
logger.info("="*80)
logger.info(f"Calculation completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Clean up any stray figure files that may have been saved to root directory
logger.info("\nOrganizing output files...")
import glob
import shutil

# Move any PNG/PDF files from root results dir to Spectra_figures
for pattern in ['*.png', '*.pdf']:
    for file in glob.glob(str(results_base_dir / pattern)):
        file_path = Path(file)
        if file_path.parent == results_base_dir:
            dest_path = spectra_figures_dir / file_path.name
            shutil.move(str(file_path), str(dest_path))
            logger.info(f"  → Moved {file_path.name} to Spectra_figures/")

logger.info("✓ All files organized in proper subdirectories")