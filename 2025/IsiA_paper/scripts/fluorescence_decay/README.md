# Fluorescence Decay Simulation Scripts

This directory contains scripts for simulating time-resolved fluorescence decay dynamics of the IsiA protein using multiple kinetic models.

## Overview

These scripts simulate how excited states of chlorophyll pigments relax through radiative (fluorescence) and non-radiative (quenching) channels. Eight different models explore how red-emitting states facilitate the quenching process in the IsiA photosynthetic antenna.

---

## Models

### Model Variants

| Model | Description | Key Mechanism | Typical Runtime |
|-------|-------------|---|---|
| **model_1** | Single red-state quenching | One low-energy state enhances quenching | 5-15 min |
| **model_2** | Two independent red states | Two separate quenching pathways | 8-20 min |
| **model_3** | Three red-state populations | Energy-dependent quenching rates | 10-25 min |
| **model_4** | Four red-state ensemble | Full state space coverage | 12-30 min |
| **model_12** | Combined 1+2 analysis | Comparative study | 8-20 min |
| **model_13** | Combined 1+3 analysis | Comparative study | 10-25 min |
| **model_14** | Combined 1+4 analysis | Comparative study | 12-30 min |
| **model_123** | All three models | Comprehensive comparison | 15-40 min |

### Model Architecture

Each model directory contains identical structure:

```
model_X/
├── run_model_X.py                      # Main simulation script
├── unified_parameters.py               # Model parameters (customize here)
├── Fluoresence_analysis.py             # Post-simulation analysis
├── run_model_X_slurm_parallel.py       # SLURM cluster submission
├── find_best_parameter_model/          # Model comparison tools
│   ├── run_best_model_pick.py          # Compare models vs. experiment
│   └── plot_best_models.py             # Visualization
└── Simulation/                         # Output directory (auto-created)
    ├── ensemble_average.npy            # Averaged fluorescence trace
    ├── time_resolved_fluorescence.npz  # Full time-resolved data
    ├── wavelength_resolved_*.npz       # Wavelength-specific traces
    ├── analysis_*.csv                  # Statistical metrics
    └── figures/                        # Generated plots
```

---

## Directory Structure

```
scripts/fluorescence_decay/
├── README.md                    # This file
├── NPZ_FORMAT.md                # Output file format specifications
├── model_1/                     # Model 1: Single red-state
├── model_2/                     # Model 2: Two red-states
├── model_3/                     # Model 3: Three red-states
├── model_4/                     # Model 4: Four red-states
├── model_12/                    # Combined 1+2
├── model_13/                    # Combined 1+3
├── model_14/                    # Combined 1+4
└── model_123/                   # Combined 1+2+3
```

---

## Installation and Setup

### Prerequisites

```bash
# Install dependencies (if not already done)
pip install -r ../../requirements.txt
pip install -e ../../pymembrane

# Ensure data files are present
ls NeededData/hamiltonian/hamiltonian_matrix.csv
ls NeededData/Experimental_Data/Fl_IsiA_monomer_300K_Gabriela_2025.npy
```

### Verify Setup

```bash
# Test imports
python -c "import numpy; import scipy; print('✓ Ready to run simulations')"

# Check Hamiltonian data exists
python -c "import numpy as np; H = np.loadtxt('../../NeededData/hamiltonian/hamiltonian_matrix.csv', delimiter=','); print(f'✓ Hamiltonian loaded: {H.shape}')"
```

---

## Usage

### Single Model Execution

```bash
# Run a single model
cd model_1
python run_model_1.py

# Expected output: Generates Simulation/ directory with results
```

### Run All Models Sequentially

```bash
# Run all models one after another
for model in model_{1,2,3,4,12,13,14,123}; do
    echo "Running $model..."
    cd scripts/fluorescence_decay/$model
    python run_model_${model/model_/}.py
    cd ../../../
done
```

### SLURM Cluster Execution

```bash
# Submit single job to cluster
cd model_1
sbatch run_model_1_slurm_parallel.py

# Check job status
squeue -u $(whoami)

# Monitor output
tail -f slurm-*.out
```

### Batch SLURM Submission

```bash
#!/bin/bash
# Submit all models as array job

for model in model_{1,2,3,4,12,13,14,123}; do
    cd scripts/fluorescence_decay/$model
    sbatch run_model_${model/model_/}.py
    cd ../../../
done
```

---

## Parameters

### Key Simulation Parameters

Edit `unified_parameters.py` in each model directory:

| Parameter | Type | Range | Default | Description | Units |
|-----------|------|-------|---------|-------------|-------|
| `N_ens` | int | 100-10000 | 1000 | Number of trajectories | - |
| `t_max` | float | 5-20 | 10 | Max simulation time | ns |
| `dt` | float | 0.001-0.01 | 0.005 | Time step | ns |
| `T` | float | 273-350 | 300 | Temperature | K |
| `wl_min` | float | 600-650 | 600 | Min wavelength | nm |
| `wl_max` | float | 700-750 | 750 | Max wavelength | nm |
| `wl_step` | float | 0.5-5.0 | 1.0 | Wavelength resolution | nm |

### Model-Specific Parameters

**model_1 (Single red state):**
```python
f_red = 0.3           # Red-state population fraction
k_red = 0.8           # Red-state quenching rate (ns⁻¹)
k_blue = 0.1          # Blue-state quenching rate (ns⁻¹)
```

**model_2 (Two red states):**
```python
f_red1 = 0.2          # First red-state fraction
f_red2 = 0.15         # Second red-state fraction
k_red1 = 0.9          # First red-state rate
k_red2 = 0.5          # Second red-state rate
```

**model_3 (Three red states):**
```python
f_red1, f_red2, f_red3 = 0.15, 0.12, 0.08
k_red1, k_red2, k_red3 = 1.0, 0.7, 0.3
```

**model_4 (Four red states):**
```python
f_red = [0.12, 0.10, 0.08, 0.05]   # State fractions
k_red = [1.1, 0.8, 0.5, 0.2]       # State quenching rates
```

---

## Computational Requirements

### Per-Model Estimates

**Default Parameters (N_ens=1000, t_max=10 ns):**

| Metric | Value |
|--------|-------|
| **Runtime** | 5-15 minutes (single CPU) |
| **Memory** | 2-4 GB |
| **Disk I/O** | ~50-100 MB output |
| **CPU Cores** | 1 (can parallelize trajectories) |

**Large Simulations (N_ens=10000, t_max=20 ns, high resolution):**

| Metric | Value |
|--------|-------|
| **Runtime** | 30-60 minutes (single CPU) |
| **Memory** | 8-16 GB |
| **Disk I/O** | ~500-1000 MB output |
| **CPU Cores** | 8+ recommended |

### Cluster Scheduling

**Typical SLURM request:**
```bash
#SBATCH --job-name=isia_model
#SBATCH --time=00:60:00          # 60 minutes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4        # Use 4 cores if parallelized
#SBATCH --mem=16G                # 16 GB RAM
#SBATCH --partition=short        # Adjust for your cluster
```

---

## Output Files

### File Types

**NPY Files** (NumPy binary arrays):
- `ensemble_average.npy` - Mean fluorescence trace
- Individual trajectory data (if saved)

**NPZ Files** (Compressed archives):
- `time_resolved_fluorescence.npz` - Complete time-resolved data
- `wavelength_resolved_*.npz` - Wavelength-specific traces
- Multiple arrays with named keys

**CSV Files** (Text data):
- `analysis_statistics.csv` - Mean, std, quantiles
- `decay_metrics.csv` - Lifetime, amplitude
- Comparison metrics

**Figures** (PNG/PDF):
- Time decay curves
- 2D heatmaps (time vs wavelength)
- Model comparisons

### File Format Details

See [NPZ_FORMAT.md](NPZ_FORMAT.md) for complete specifications including:
- Array dimensions and shapes
- Data types and units
- How to load and process files
- Example Python code

### Output Location

```
model_X/Simulation/
├── ensemble_average.npy                    # Primary output
├── time_resolved_fluorescence.npz          # Complete data
├── wavelength_resolved_680nm.npz           # Example wavelength
├── wavelength_resolved_700nm.npz           # Example wavelength
├── analysis_ensemble_statistics.csv        # Statistics
├── analysis_decay_components.csv           # Fitted decays
└── figures/
    ├── time_resolved_heatmap.png
    ├── wavelength_traces.png
    ├── ensemble_comparison.png
    └── decay_curves.png
```

---

## Analysis and Comparison

### Post-Simulation Analysis

```bash
# Within model directory
python Fluoresence_analysis.py
# Generates analysis plots and CSV metrics
```

### Compare Models

```bash
# Compare specific models against experimental data
cd model_1/find_best_parameter_model
python run_best_model_pick.py

# This compares all 8 models and ranks by fit quality
```

### Expected Outputs

```
Model Comparison Results:
  Model 1:   RMSE = 0.0234  ✓ Best fit
  Model 2:   RMSE = 0.0312
  Model 3:   RMSE = 0.0445
  ...
  Model 123: RMSE = 0.0856
```

---

## Expected Results

### Validation Checks

For each model, verify:

1. **Simulation completes without errors**
2. **Output files generated** (ensemble_average.npy, etc.)
3. **Fluorescence decays from 1.0 to ~0.001 within 10 ns**
4. **Decay is smooth (no discontinuities)**
5. **Multiple wavelengths show similar decay patterns**

### Example Output Values

```
Time [ns]  Intensity [a.u.]  Wavelength [nm]
0.0        1.000             680
0.5        0.450             680
1.0        0.180             680
5.0        0.015             680
10.0       0.001             680

Decay Parameters:
  τ₁ = 0.5 ns (amplitude = 0.4)
  τ₂ = 2.5 ns (amplitude = 0.6)
  Average lifetime = 1.7 ns
```

---

## Troubleshooting

### Common Issues

**Issue:** `ModuleNotFoundError: No module named 'pymembrane'`
- **Solution:** `pip install -e ../../pymembrane` from repo root

**Issue:** `FileNotFoundError: Cannot find hamiltonian_matrix.csv`
- **Solution:** Run from model directory or use absolute paths
- Check: Hamiltonian calculation completed? Run `scripts/hamiltonian/cdc_IsiA_monomer_average_pH7.py`

**Issue:** `MemoryError: Unable to allocate X.XX GiB`
- **Solution:** Reduce N_ens or use smaller time/wavelength ranges
- Recommended: N_ens ≤ 5000 for 16GB RAM

**Issue:** Simulation runs very slowly (>1 hour for small ensemble)**
- **Solution:** Check system load, close other applications
- Consider reducing wl_step or n_time_points

**Issue:** NaN or Inf values in output**
- **Cause:** Invalid parameters (e.g., negative rates, invalid temperatures)
- **Solution:** Check unified_parameters.py for physically reasonable values

**Issue:** Different results each run (stochastic variation)**
- **Expected:** Ensemble simulations include randomness
- **Solution:** Increase N_ens for more stable averages
- Use same random seed for reproducibility

### Debugging

Enable verbose output by adding to script:

```python
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
logger.debug(f"Starting simulation with N_ens={N_ens}, t_max={t_max}")
```

---

## Performance Optimization

### Speed Up Simulations

1. **Reduce ensemble size:** N_ens = 1000 → 500 (faster, noisier)
2. **Coarsen time resolution:** dt = 0.005 → 0.01 ns
3. **Coarsen wavelength resolution:** wl_step = 1.0 → 5.0 nm
4. **Use SSD storage:** Faster I/O
5. **Parallelize:** Use multi-threading (if implemented)

### Memory Optimization

1. **Store every Nth trajectory:** Save disk space
2. **Use float32 instead of float64:** Reduces size by 50%
3. **Compress NPZ files:** np.savez_compressed()
4. **Clear temporary arrays:** del large_array after use

---

## Citation

If you use these fluorescence decay simulations, please cite:

Mohamed A. A. Elrefaiy, Dvir Harris, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, and Gabriela S. Schlau-Cohen. "Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States." *Manuscript in preparation*, 2025.

---

## Related Documentation

- **NPZ output format:** [NPZ_FORMAT.md](NPZ_FORMAT.md)
- **Hamiltonian input:** [../hamiltonian/README.md](../hamiltonian/README.md)
- **Experimental data:** [../../NeededData/Experimental_Data/README.md](../../NeededData/Experimental_Data/README.md)
- **Parameter definitions:** [../../NeededData/DATA_DICTIONARY.md](../../NeededData/DATA_DICTIONARY.md)
- **Main README:** [../../README.md](../../README.md)

---

*Last updated: January 26, 2025*
