"""
Module: Model Comparison and Best-Fit Analysis
Purpose: Compare fluorescence decay simulations with experimental data
Author: Mohamed A. A. Elrefaiy
Date: 2025
Dependencies: numpy, scipy, matplotlib, pathlib, Best_model_pick

Usage:
    python run_best_model_pick.py

Description:
    This script analyzes fluorescence decay simulations from all models in the parent directories
    and compares them against experimental data. It:
    1. Loads simulation results from each model variant
    2. Calculates goodness-of-fit metrics (RMSE, χ², etc.)
    3. Identifies best-matching model parameters
    4. Generates comparison plots and statistical analysis
    5. Exports results for machine learning applications

Output:
    - Comparison tables (CSV format)
    - Fit metrics and residuals
    - Visual comparisons (PNG/PDF)
    - Machine learning dataset (NumPy format)

Example:
    >>> python run_best_model_pick.py
    >>> # Outputs saved to outputs/ subdirectory

Notes:
    - Requires Simulation/ directory with model outputs
    - Experimental data must be present in NeededData/Experimental_Data/
    - Results are reproducible (no randomness in comparison)
"""

import os
import time
from pathlib import Path
from Best_model_pick import SimulationAnalyzer
import matplotlib.pyplot as plt
plt.style.use('default')


# Analysis parameters
base_dir = 'Simulation'
exp_wavelength_file = Path('../../../NeededData/Experimental_Data') / 'Fl_IsiA_monomer_300K_Gabriela_2025.npy'
exp_time_params = {'a1_pct': 23.7, 'tau1': 0.21, 'a2_pct': 76.3, 'tau2': 3.89, 'tau_avg': 3.02}
wavelength_range = (650, 760)
N_best = 10

# Create the analyzer
analyzer = SimulationAnalyzer(
    base_dir=base_dir,
    exp_wavelength_file=exp_wavelength_file,
    exp_time_params=exp_time_params,
    wavelength_range=wavelength_range,
    N_best=N_best
)

# Run the analysis
print(f"Starting simulation analysis at {time.strftime('%Y-%m-%d %H:%M:%S')}")
print("-" * 80)
analyzer.run_analysis()

# Export ALL simulation data for machine learning
print("\nExporting all simulation data for machine learning...")
all_data_path = analyzer.export_all_simulation_data()
print(f"All simulation data exported to: {all_data_path}")

# Generate additional visualizations
print("\nGenerating error space visualization...")
analyzer.plot_error_space()

print("\nGenerating visual model comparison...")
analyzer.visual_model_comparison(top_n=5)

print("-" * 80)
print(f"Analysis completed at {time.strftime('%Y-%m-%d %H:%M:%S')}")