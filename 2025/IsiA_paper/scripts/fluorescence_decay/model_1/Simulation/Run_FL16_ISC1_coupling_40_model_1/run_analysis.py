#!/usr/bin/env python3
import os 

# These files need to be in the current directory "Simulation" to run the analysis:
os.system(f'cp -r ../class_save_ensamble.py .')
os.system(f'cp -r ../Fluoresence_analysis.py .')
os.system(f'cp -r ../load_param.py .')

from Fluoresence_analysis import FluorescenceData
from load_param import load_parameter_scan_definition
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Base data directory
data_dir = 'Data'

# Path to experimental data
"""
Experimental data directory resolution
-------------------------------------
Uses Path to discover the repository root by walking upward from this file
until it finds the standard folder 'NeededData/Experimental_Data'. This keeps
the script portable across machines. You may override the autodetected path by
setting the environment variable ISIA_EXP_DATA_DIR to a custom directory.
"""

def _resolve_experimental_data_dir() -> Path:
    # 1) Environment variable override (portable and explicit)
    exp_dir_env = os.environ.get('ISIA_EXP_DATA_DIR')
    if exp_dir_env:
        return Path(exp_dir_env).expanduser().resolve()

    # 2) Autodiscover by walking up parents looking for NeededData/Experimental_Data
    here = Path(__file__).resolve()
    for parent in [here.parent, *here.parents]:
        candidate = parent / 'NeededData' / 'Experimental_Data'
        if candidate.exists() and candidate.is_dir():
            return candidate.resolve()

    # 3) If not found, raise a helpful error
    raise FileNotFoundError(
        "Could not locate 'NeededData/Experimental_Data' by walking up from this script.\n"
        "Please ensure you are running within the repository checkout, or set\n"
        "the environment variable ISIA_EXP_DATA_DIR to the experimental data directory."
    )

exp_dir = _resolve_experimental_data_dir()
print(f"Looking for experimental data in: {exp_dir}")

# Run ensemble calculation if requested
run_ensemble = False  # turn this on if you are running for the first time to generate seed files
if run_ensemble:
    print(f"\nRunning ensemble calculation for {data_dir}...")
    from class_save_ensamble import run_ensemble_calculation
    success = run_ensemble_calculation(data_dir, cleanup=False)
    if not success:
        print("Warning: Ensemble calculation failed or found no data.")
        print("Make sure the data directory contains valid seed files.")

correct_run_directory = '.'
dict_loaded_params = load_parameter_scan_definition(correct_run_directory)
print(dict_loaded_params)

# Initialize with all available components
fl_data = FluorescenceData(
    main_dir_name=data_dir,
    components=None,  # Use all available
    exp_dir=exp_dir,
    start_time=0.0,  # Include all time points
    params=dict_loaded_params,
)

# Create a dict of monomer/CT ratios to analyze (0% to 100% monomer in steps of 10%)
ratios = {}

# Calculate complementary percentages
mon_percent = 100
ct_percent = 0 # model2
rp1_ct_percent = 0 # model_4
rp1_percent = 0 # model_3
# This includes all percentages needed for correct file naming
label = f"Monomer {mon_percent}%: {mon_percent}% M, {ct_percent}% CT, {rp1_ct_percent}% RP1+CT, {rp1_percent}% RP1"
ratios[label] = (mon_percent, ct_percent, rp1_ct_percent, rp1_percent)

# Plot all ratio combinations and save the figures and CSV data
fl_data.plot_multiple_ratios(ratios, save=True, show=False)

# Example: Create specific mixture (70% monomer, 30% CT)
# fl_data.create_specific_mixture(
#     monomer_pct=70,
#     ct_pct=30,
#     rp1_pct=0,
#     rp1_only_pct=0,
#     save=True,
#     show=False,
# )

print("\nAnalysis complete. Biexponential fit results have been saved to CSV files.")