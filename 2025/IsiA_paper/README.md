# IsiA Protein Spectroscopy Analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/downloads/)
[![Zenodo](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.XXXXXX-blue)](https://zenodo.org/)

## Description

This repository contains computational data and analysis code accompanying the manuscript:

> **Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States**
>
> Mohamed A. A. Elrefaiy, Dvir Harris, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, Gabriela S. Schlau-Cohen
>
> *Manuscript in preparation, 2025*

This work presents a comprehensive computational analysis of the IsiA protein complex, including:
- Excitonic Hamiltonian calculations from cryo-EM structures
- Absorption and fluorescence spectra predictions
- Time-resolved fluorescence decay simulations with multiple kinetic models
- Model comparison and parameter optimization

---

## Quick Start (5 minutes)

**Want to get running immediately?** See [QUICKSTART.md](QUICKSTART.md) for a fast setup guide with example commands.

### Minimal Installation

```bash
# Clone and setup
git clone https://github.com/melrefaiy2018/IsiA_paper.git
cd IsiA_paper

# Create environment
python3 -m venv env
source env/bin/activate

# Install dependencies
pip install -r requirements.txt
pip install -e pymembrane

# Run first example (takes ~1 minute)
cd scripts/hamiltonian
python cdc_IsiA_monomer_average_pH7.py
```

---

## System Requirements

### Minimum Requirements
- **OS:** Linux, macOS, or Windows (WSL2 recommended)
- **Python:** 3.8 or higher
- **RAM:** 8 GB (16 GB recommended)
- **Disk space:** 2-4 GB free
- **Processor:** 2+ CPU cores (4+ cores recommended for parallel simulations)

### Recommended Configuration
- **Python:** 3.9 or 3.10
- **RAM:** 16+ GB (for ensemble simulations with 10,000+ trajectories)
- **Disk space:** 4+ GB (for storing simulation outputs)
- **Processor:** 8+ CPU cores (for SLURM cluster execution)
- **Compiler:** GCC/Clang (for building C extensions)

### Tested Environments
- Ubuntu 20.04 LTS, Python 3.8, numpy 1.20
- macOS 12.x, Python 3.9, numpy 1.21
- Windows 11 WSL2 (Ubuntu 20.04), Python 3.10, numpy 1.23

---

## Repository Structure

```
IsiA_paper/
├── README.md                  # This file
├── QUICKSTART.md              # 5-minute setup guide
├── CITATION.cff               # Citation metadata
├── CHANGELOG.md               # Version history
├── AUTHORS.md                 # Contributor information
├── LICENSE                    # MIT license
├── requirements.txt           # Python dependencies
├── .gitignore                 # Git ignore patterns
├── .zenodo.json               # Zenodo metadata
│
├── NeededData/                # Input data and generated outputs
│   ├── Experimental_Data/     # Experimental fluorescence/absorption
│   │   ├── archive/           # Older versions (for reference)
│   │   ├── csv/               # CSV format data files
│   │   └── README.md          # Data specifications
│   ├── hamiltonian/           # Hamiltonian and spectral outputs
│   │   ├── Hamiltonian/       # Detailed matrices
│   │   ├── SiteEnergy_Data/   # Per-pigment energies
│   │   ├── Spectra_Data/      # Raw spectral calculations
│   │   └── README.md          # Format documentation
│   ├── mcce/                  # MCCE electrostatics output
│   │   └── README.md          # pKa and protonation data
│   ├── structure/             # Protein PDB files
│   │   └── README.md          # Structure descriptions
│   └── DATA_DICTIONARY.md     # Parameter definitions and units
│
├── scripts/                   # Executable analysis scripts
│   ├── hamiltonian/           # Hamiltonian calculation
│   │   ├── cdc_IsiA_monomer_average_pH7.py
│   │   └── README.md
│   └── fluorescence_decay/    # Fluorescence decay simulations
│       ├── model_1/ to model_123/  # 8 different model variants
│       │   ├── run_model_X.py
│       │   ├── unified_parameters.py
│       │   ├── Fluoresence_analysis.py
│       │   ├── find_best_parameter_model/  # Model comparison scripts
│       │   ├── Simulation/   # Output directory
│       │   └── README.md
│       ├── NPZ_FORMAT.md      # File format documentation
│       └── README.md
│
└── pymembrane/                # Local pymembrane library
    ├── setup.py
    └── README.md
```

### Key Directories and Files

- **`NeededData/`**: All input data and calculated outputs
  - `Experimental_Data/`: Reference spectroscopy measurements
  - `hamiltonian/`: Pre-calculated Hamiltonian matrices and spectra
  - `mcce/`: Electrostatic and protonation state data
  - `structure/`: Protein atomic coordinates (PDB format)
  - `DATA_DICTIONARY.md`: Complete parameter definitions

- **`scripts/`**: Reproducible analysis workflows
  - `hamiltonian/`: Calculate excitonic Hamiltonian from structure
  - `fluorescence_decay/`: Simulate time-resolved fluorescence
  - 8 different models with parameter optimization
  - SLURM scripts for cluster execution

- **Documentation files:**
  - `QUICKSTART.md`: 5-minute setup
  - `CITATION.cff`: Proper citation format
  - `CHANGELOG.md`: Version history
  - `AUTHORS.md`: Contributor roles

---

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/melrefaiy2018/IsiA_paper.git
cd IsiA_paper
```

### 2. Create Python Virtual Environment

```bash
# Using venv (recommended)
python3 -m venv env
source env/bin/activate  # On Windows: env\Scripts\activate

# Or using conda
conda create -n isia python=3.9
conda activate isia
```

### 3. Install Dependencies

```bash
# Core scientific packages
pip install -r requirements.txt

# Install local pymembrane library in editable mode
pip install -e pymembrane
```

### 4. Verify Installation

```bash
# Test imports
python -c "import numpy; import scipy; import pymembrane; print('✓ All imports successful')"

# Check data files
python -c "import os; assert os.path.exists('NeededData/structure/6kig_structure_prepared_mcce_input.pdb'); print('✓ Data files found')"
```

### Troubleshooting Installation

**Issue:** `ModuleNotFoundError: No module named 'pymembrane'`
- **Solution:** Ensure you ran `pip install -e pymembrane` from the repository root

**Issue:** ImportError for numpy/scipy on macOS
- **Solution:** Install compiler tools: `xcode-select --install`

**Issue:** Permission denied on structure files
- **Solution:** `chmod 644 NeededData/structure/*.pdb`

See [QUICKSTART.md](QUICKSTART.md#common-issues) for more troubleshooting.

---

## Computational Requirements

### Hamiltonian Calculation
- **Execution time:** 30 seconds - 2 minutes
- **Memory:** ~500 MB
- **Output size:** ~25 KB (matrices) + ~20 KB (spectra)
- **Parallelization:** Single-threaded (no speedup from multi-core)

### Fluorescence Decay Simulations (per model)
- **Execution time:** 2-30 minutes (varies by ensemble size and parameters)
- **Memory:** 2-8 GB (for default 1000-10000 trajectory ensemble)
- **Output size:** 10-500 MB per model (depends on time/wavelength resolution)
- **Parallelization:** Fully parallelizable (SLURM scripts provided)
- **SLURM job request:** Typical 4-16 cores, 30 minutes - 2 hours

### Model Comparison
- **Execution time:** 5-10 minutes per model
- **Memory:** ~1 GB
- **Output:** Goodness-of-fit metrics and comparison plots

### Full Workflow (All 8 Models)
- **Total execution time:** 3-4 hours (single CPU) or 30-60 minutes (8-core parallel)
- **Total disk space:** ~1-2 GB for all outputs
- **Recommended:** Run on HPC cluster or powerful workstation

---

## How to Reproduce the Results

### 1. Calculate the Excitonic Hamiltonian

Calculate site energies and electronic couplings to generate the excitonic Hamiltonian:

```bash
cd scripts/hamiltonian
python cdc_IsiA_monomer_average_pH7.py
```

**What it does:**
- Reads protein structure (PDB format)
- Extracts MCCE-derived electrostatic energies
- Calculates transition dipole couplings (point-dipole approximation)
- Constructs 17×17 Hamiltonian matrix
- Calculates absorption and fluorescence spectra

**Outputs generated:**
```
NeededData/hamiltonian/
├── hamiltonian_matrix.csv          # 17×17 Hamiltonian (cm⁻¹)
├── absorption_spectrum.csv         # Calculated absorption
├── fluorescence_spectrum.csv       # Calculated fluorescence
├── Hamiltonian/                    # Detailed matrices
├── SiteEnergy_Data/                # Per-pigment site energies
└── Spectra_Data/                   # Spectral components
```

**Expected results:**
- Absorption peak at ~680 nm
- Hamiltonian values: diagonal 14500-15500 cm⁻¹, couplings ±100 cm⁻¹

### 2. Simulate Fluorescence Decay

Simulate time-resolved fluorescence using multiple kinetic models:

```bash
# Run a single model
cd scripts/fluorescence_decay/model_1
python run_model_1.py

# Or run all models in sequence
for model in model_{1,2,3,4,12,13,14,123}; do
    cd scripts/fluorescence_decay/$model
    python run_model_${model/model_/}.py
done
```

**Parameters to customize** (edit `unified_parameters.py`):
- `N_ens`: Number of trajectories (100-10000)
- `t_max`: Simulation time (5-20 ns)
- `dt`: Time step (0.001-0.01 ns)
- Model-specific parameters (see comments in each `unified_parameters.py`)

**For SLURM clusters:**
```bash
sbatch run_model1_slurm_parallel.py
```

**Outputs generated:**
```
scripts/fluorescence_decay/model_1/Simulation/
├── ensemble_average.npy            # Mean fluorescence trace
├── time_resolved_fluorescence.npz  # Time-resolved data
├── wavelength_resolved_*.npz       # Wavelength-specific traces
├── analysis_*.csv                  # Statistical metrics
└── figures/                         # Analysis plots
```

### 3. Compare Models and Find Best Fit

Compare all models against experimental data:

```bash
cd scripts/fluorescence_decay/model_1/find_best_parameter_model
python run_best_model_pick.py
```

**What it does:**
- Loads simulation results from all models
- Compares with experimental fluorescence data
- Calculates RMSE and χ² goodness-of-fit
- Generates comparison plots
- Ranks models by fit quality

**Outputs:**
```
find_best_parameter_model/outputs/
├── model_comparison.csv            # Fit metrics for all models
├── best_fit_parameters.csv         # Best-fit model parameters
├── comparison_plots.png            # Visual comparison
└── residuals_analysis.csv          # Fit residuals
```

---

## Data Availability and Formats

### Experimental Data
- **Primary file:** `NeededData/Experimental_Data/Fl_IsiA_monomer_300K_Gabriela_2025.npy`
- **Format:** NumPy binary (1D or 2D array)
- **Temperature:** 300 K (room temperature)
- **Time range:** 0-10 ns
- **Wavelength range:** 600-750 nm

See [NeededData/Experimental_Data/README.md](NeededData/Experimental_Data/README.md) for full specifications.

### Hamiltonian Data
- **Hamiltonian matrix:** 17×17 symmetric, units cm⁻¹
- **Spectra:** Wavelength vs. intensity, normalized to 1.0
- **Format:** CSV files (human-readable)

See [NeededData/hamiltonian/README.md](NeededData/hamiltonian/README.md) for details.

### Simulation Outputs
- **Time-resolved:** NPZ format with time, wavelength, and intensity arrays
- **Ensemble:** Multiple trajectories with statistics
- **Format:** Documented in [scripts/fluorescence_decay/NPZ_FORMAT.md](scripts/fluorescence_decay/NPZ_FORMAT.md)

### Parameter Dictionary
All parameter abbreviations and units documented in [NeededData/DATA_DICTIONARY.md](NeededData/DATA_DICTIONARY.md):
- `E_0a`, `E_0b`: Reference energies (cm⁻¹)
- `dielc`: Dielectric constant
- `N_ens`: Ensemble size
- And many more...

---

## Expected Results and Outputs

### Example: Hamiltonian Calculation Output
```
Hamiltonian diagonal (site energies):
  [14823.4, 14921.3, 14756.8, ..., 15102.1] cm⁻¹

Absorption spectrum peak: 680.3 nm
Fluorescence spectrum peak: 698.7 nm
```

### Example: Fluorescence Decay Simulation
```
Time-resolved fluorescence (model_1):
  Time [ns]    Intensity [a.u.]
  0.0          1.00
  0.1          0.87
  0.5          0.45
  1.0          0.18
  5.0          0.02
  10.0         < 0.001
```

### Example: Model Comparison Results
```
Model Rankings (by RMSE):
  1. model_1:   RMSE = 0.023
  2. model_2:   RMSE = 0.031
  3. model_123: RMSE = 0.045
  ...
```

---

## Troubleshooting

### Common Issues and Solutions

**Problem:** `FileNotFoundError: Cannot find 'NeededData/structure/...'`
- **Check:** Are you running from the repository root?
- **Solution:** `cd` to repository root and use absolute paths if needed

**Problem:** Memory error during simulation
- **Cause:** Too many ensemble members or high time/wavelength resolution
- **Solution:** Reduce `N_ens` in `unified_parameters.py`, or increase system RAM

**Problem:** SLURM job fails with module not found
- **Solution:** Load Python module on cluster: `module load python/3.9`
- Ensure pymembrane is installed in correct environment

**Problem:** Slow simulation (>10 min for small ensemble)
- **Cause:** Inefficient parameter choices or system load
- **Solution:** Reduce output frequency or wavelength resolution

**Problem:** Matplotlib cannot display figures
- **Solution:** Use non-interactive backend: `matplotlib.use('Agg')`

For more help, see [QUICKSTART.md#common-issues](QUICKSTART.md#common-issues).

---

## Citation

### Cite the Manuscript

If you use this code or data, please cite:

```bibtex
@article{elrefaiy2025isia,
  author = {Elrefaiy, Mohamed A. A. and Harris, Dvir and Toporik, Hila and
            Gisriel, Christopher J. and Mazor, Yuval and Raccah, Doran I. G. B. and
            Schlau-Cohen, Gabriela S.},
  title = {Quenching of the Photosynthetic Antenna {IsiA} is Facilitated by its Red-Emitting States},
  journal = {[Journal Name]},
  year = {2025},
  note = {In preparation}
}
```

### Cite the Code/Data Repository

```bibtex
@software{elrefaiy2025isia_zenodo,
  author = {Elrefaiy, Mohamed A. A. and others},
  title = {IsiA Protein Spectroscopy Analysis: Code and Data},
  year = {2025},
  url = {https://github.com/melrefaiy2018/IsiA_paper},
  doi = {10.5281/zenodo.XXXXXX}
}
```

Or use the [CITATION.cff](CITATION.cff) file for automatic citation formatting.

---

## Funding and Acknowledgments

This research was supported by:
- [Funding agency 1] - Grant [number]
- [Funding agency 2] - Grant [number]
- MIT Department of Chemistry

We thank:
- Gabriela S. Schlau-Cohen lab members for experimental guidance
- MIT Center for Excitonics for computational resources
- The open-source scientific Python community

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

All code is provided "as is" for research and educational purposes.

---

## Documentation Guide

- **[QUICKSTART.md](QUICKSTART.md)** - 5-minute setup guide
- **[NeededData/DATA_DICTIONARY.md](NeededData/DATA_DICTIONARY.md)** - Parameter definitions
- **[scripts/fluorescence_decay/NPZ_FORMAT.md](scripts/fluorescence_decay/NPZ_FORMAT.md)** - File format details
- **[NeededData/Experimental_Data/README.md](NeededData/Experimental_Data/README.md)** - Experimental data specs
- **[NeededData/hamiltonian/README.md](NeededData/hamiltonian/README.md)** - Hamiltonian data
- **[NeededData/structure/README.md](NeededData/structure/README.md)** - PDB structure files
- **[NeededData/mcce/README.md](NeededData/mcce/README.md)** - Electrostatics data
- **[CHANGELOG.md](CHANGELOG.md)** - Version history
- **[AUTHORS.md](AUTHORS.md)** - Contributor information

---

## Contributing

This repository documents published research. For questions or issues:

1. **Report bugs:** Open an issue on GitHub
2. **Discuss methods:** Use GitHub discussions
3. **Suggest improvements:** Submit pull requests (for documentation)

Note: Major changes to analysis code may not be accepted as they could affect reproducibility.

---

## Contact

**Principal Investigator:** Gabriela S. Schlau-Cohen (MIT Chemistry)

**Lead Developer:** Mohamed A. A. Elrefaiy (MIT Chemistry)

**Questions?** Open an issue on GitHub or contact the authors.

---

*Last updated: January 26, 2025*

*Repository version: 1.0.0*

*Status: Preparing for Zenodo release*
