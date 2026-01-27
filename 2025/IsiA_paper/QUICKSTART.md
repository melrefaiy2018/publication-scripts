# Quick Start Guide

Get up and running with the IsiA spectroscopy analysis in 5 minutes.

## Prerequisites

- Python 3.8 or higher
- Git
- 2-4 GB free disk space
- 8+ GB RAM (recommended for running simulations)

## Installation (2 minutes)

```bash
# 1. Clone the repository
git clone https://github.com/melrefaiy2018/IsiA_paper.git
cd IsiA_paper

# 2. Create and activate virtual environment
python3 -m venv env
source env/bin/activate  # On Windows: env\Scripts\activate

# 3. Install dependencies
pip install -r requirements.txt

# 4. Install pymembrane library
pip install -e pymembrane
```

## Quick Test Run (3 minutes)

### Option 1: Calculate Hamiltonian and Spectra

Calculate the excitonic Hamiltonian for the IsiA monomer:

```bash
cd scripts/hamiltonian
python cdc_IsiA_monomer_average_pH7.py
```

**Expected output:** Hamiltonian matrices and spectra saved to `NeededData/hamiltonian/`

**Runtime:** ~30-60 seconds

**What it does:**
- Reads protein structure from PDB file
- Calculates site energies and couplings
- Generates absorption and fluorescence spectra

### Option 2: Run Fluorescence Decay Simulation

Simulate fluorescence decay using Model 1:

```bash
cd scripts/fluorescence_decay/model_1
python run_model_1.py
```

**Expected output:** Time-resolved fluorescence data in `Simulation/` subdirectory

**Runtime:** ~2-5 minutes (depending on system)

**What it does:**
- Loads Hamiltonian and experimental data
- Runs ensemble-averaged decay simulations
- Generates wavelength-resolved fluorescence traces

## Verify Installation

Check that outputs were created successfully:

```bash
# Check Hamiltonian outputs
ls -lh NeededData/hamiltonian/*.csv

# Check fluorescence simulation outputs
ls -lh scripts/fluorescence_decay/model_1/Simulation/
```

You should see:
- CSV files with Hamiltonian matrices and spectra
- NPY/NPZ files with time-resolved fluorescence data
- Analysis figures (if matplotlib is installed)

## Example Outputs

### Hamiltonian Calculation
- `hamiltonian_matrix.csv`: 17×17 excitonic Hamiltonian (cm⁻¹)
- `absorption_spectrum.csv`: Calculated absorption spectrum
- `fluorescence_spectrum.csv`: Calculated fluorescence spectrum
- `SiteEnergy_Data/`: Site energy distributions for each chromophore

### Fluorescence Decay Simulation
- `ensemble_*.npy`: Ensemble-averaged fluorescence traces
- `time_resolved_*.npz`: Full time-resolved data
- `wavelength_resolved_*.npz`: Wavelength-resolved decay traces
- `analysis_*.csv`: Statistical analysis of decay components

## Next Steps

Once you've verified the basic installation works:

1. **Explore different models:** Try models 2, 3, 4, etc. in `scripts/fluorescence_decay/`
2. **Parameter optimization:** Use `find_best_parameter_model/` scripts to compare models
3. **Read the full documentation:** See [README.md](README.md) for detailed usage
4. **Customize parameters:** Edit `unified_parameters.py` in each model directory

## Common Issues

### Import Error: No module named 'pymembrane'
**Solution:** Make sure you ran `pip install -e pymembrane` from the repository root

### FileNotFoundError: Cannot find experimental data
**Solution:** Verify you're running scripts from the correct directory (see paths above)

### MemoryError during simulation
**Solution:** Reduce ensemble size (`N_ens`) in `unified_parameters.py`

### Slow performance
**Solution:**
- Use fewer time points or wavelengths
- Run on a system with more CPU cores
- Consider using the SLURM parallel scripts for cluster execution

## Getting Help

- **Full documentation:** [README.md](README.md)
- **Data dictionary:** [NeededData/DATA_DICTIONARY.md](NeededData/DATA_DICTIONARY.md)
- **File formats:** [scripts/fluorescence_decay/NPZ_FORMAT.md](scripts/fluorescence_decay/NPZ_FORMAT.md)
- **Issues:** Open an issue on GitHub
- **Citation:** See [CITATION.cff](CITATION.cff)

## System Requirements

**Minimum:**
- Python 3.8+
- 4 GB RAM
- 2 CPU cores
- 2 GB disk space

**Recommended:**
- Python 3.9+
- 16 GB RAM
- 8+ CPU cores
- 4 GB disk space
- Linux or macOS (Windows also supported)

## Quick Validation

To verify your installation is working correctly, run:

```bash
# Test Python imports
python -c "import numpy; import scipy; import pymembrane; print('All imports successful!')"

# Check data files exist
python -c "import os; assert os.path.exists('NeededData/structure/6kig_structure_prepared_mcce_input.pdb'); print('Data files found!')"
```

If both commands succeed, you're ready to go!
