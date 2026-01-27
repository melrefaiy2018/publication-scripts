# NPZ and NPY File Format Documentation

This document describes the structure and contents of all NumPy binary files (.npy and .npz) used in the fluorescence decay simulations.

## Table of Contents
- [Overview](#overview)
- [File Types](#file-types)
- [Common File Formats](#common-file-formats)
- [Loading and Accessing Data](#loading-and-accessing-data)
- [Examples](#examples)

---

## Overview

NumPy binary files are used throughout this project for efficient storage of numerical data:
- **.npy files:** Single array storage
- **.npz files:** Multiple arrays stored together (like a dictionary)

**Advantages:**
- Fast loading/saving
- Preserves numerical precision
- Efficient compression
- Native Python/NumPy support

---

## File Types

### 1. Experimental Data Files

#### `Fl_IsiA_monomer_300K_Gabriela_2025.npy`
**Location:** `NeededData/Experimental_Data/`

**Description:** Experimental fluorescence decay data at 300K

**Structure:**
- **Format:** Single 1D or 2D NumPy array
- **Dimensions:**
  - If 1D: `(n_time_points,)` - time-resolved trace
  - If 2D: `(n_wavelengths, n_time_points)` - wavelength-resolved data
- **Data type:** `float64`
- **Units:** Normalized fluorescence intensity (0-1, arbitrary units)

**Loading:**
```python
import numpy as np
fl_data = np.load('NeededData/Experimental_Data/Fl_IsiA_monomer_300K_Gabriela_2025.npy')
```

---

### 2. Hamiltonian Data Files

#### `hamiltonian_matrix.npy` or `.csv`
**Location:** `NeededData/hamiltonian/`

**Description:** Excitonic Hamiltonian matrix

**Structure:**
- **Format:** 2D square matrix
- **Dimensions:** `(17, 17)` for IsiA monomer
- **Data type:** `float64`
- **Units:** cm⁻¹ (wavenumbers)
- **Properties:**
  - Symmetric matrix
  - Diagonal: site energies
  - Off-diagonal: excitonic couplings

**Loading:**
```python
# If NPY format
H = np.load('NeededData/hamiltonian/hamiltonian_matrix.npy')

# If CSV format
H = np.loadtxt('NeededData/hamiltonian/hamiltonian_matrix.csv', delimiter=',')
```

---

### 3. Simulation Output Files

#### Time-Resolved Fluorescence: `time_resolved_*.npz`
**Location:** `scripts/fluorescence_decay/model_*/Simulation/`

**Description:** Complete time-resolved fluorescence data with metadata

**Keys and Structure:**
```python
{
    'time': array,          # Time axis, shape: (n_time,), units: ns
    'wavelength': array,    # Wavelength axis, shape: (n_wl,), units: nm
    'fluorescence': array,  # Fluorescence data, shape: (n_time, n_wl)
    'parameters': dict,     # Simulation parameters (stored as structured array)
    'temperature': float,   # Temperature in Kelvin
}
```

**Dimensions:**
- `time`: Typically 0-10 ns with variable resolution (100-1000 points)
- `wavelength`: Typically 600-750 nm with 1 nm resolution
- `fluorescence`: Time × wavelength array of intensity values

**Loading:**
```python
data = np.load('time_resolved_output.npz', allow_pickle=True)
time = data['time']
wavelength = data['wavelength']
fluorescence = data['fluorescence']
params = data['parameters'].item()  # Convert to dict if stored as object
```

---

#### Wavelength-Resolved Data: `wavelength_resolved_*.npz`
**Location:** `scripts/fluorescence_decay/model_*/Simulation/`

**Description:** Fluorescence decay traces at specific wavelengths

**Keys and Structure:**
```python
{
    'time': array,              # Time axis, shape: (n_time,)
    'wavelengths': array,       # Selected wavelengths, shape: (n_selected_wl,)
    'traces': array,            # Decay traces, shape: (n_selected_wl, n_time)
    'wavelength_indices': array # Indices mapping to full wavelength array
}
```

**Example wavelengths:** [680, 690, 700, 710, 720] nm (red region of spectrum)

**Loading:**
```python
data = np.load('wavelength_resolved_output.npz')
time = data['time']
wls = data['wavelengths']
traces = data['traces']  # Shape: (n_wl, n_time)

# Plot decay at first wavelength
import matplotlib.pyplot as plt
plt.plot(time, traces[0])
plt.xlabel('Time (ns)')
plt.ylabel('Fluorescence (a.u.)')
plt.title(f'Decay at {wls[0]} nm')
```

---

#### Ensemble Data: `ensemble_*.npy` or `ensemble_*.npz`
**Location:** `scripts/fluorescence_decay/model_*/Simulation/`

**Description:** Ensemble-averaged data from multiple trajectories

**Format Options:**

**Option 1: Single Array (.npy)**
```python
# Shape: (n_time,) - averaged over all ensemble members and wavelengths
ensemble_avg = np.load('ensemble_average.npy')
```

**Option 2: Multiple Arrays (.npz)**
```python
{
    'ensemble_traces': array,  # Shape: (n_ensemble, n_time, n_wl)
    'ensemble_avg': array,     # Shape: (n_time, n_wl) - mean over ensemble
    'ensemble_std': array,     # Shape: (n_time, n_wl) - std over ensemble
    'n_ensemble': int          # Number of ensemble members
}
```

**Loading:**
```python
data = np.load('ensemble_output.npz')
traces = data['ensemble_traces']      # Individual trajectories
avg = data['ensemble_avg']            # Mean
std = data['ensemble_std']            # Standard deviation
n_ens = data['n_ensemble'].item()     # Convert to Python int
```

---

### 4. Analysis Output Files

#### Best-Fit Results: `best_fit_results_*.npz`
**Location:** `scripts/fluorescence_decay/model_*/find_best_parameter_model/outputs/`

**Description:** Results from parameter optimization

**Keys and Structure:**
```python
{
    'parameters': array,        # Best-fit parameter values, structured array
    'chi_squared': float,       # Chi-squared goodness of fit
    'rmse': float,             # Root mean square error
    'fitted_trace': array,     # Best-fit fluorescence trace
    'experimental_trace': array, # Experimental data for comparison
    'residuals': array,        # Fitted - experimental
    'param_grid': array        # Parameter grid searched (if grid search)
}
```

**Loading:**
```python
results = np.load('best_fit_results.npz', allow_pickle=True)
best_params = results['parameters'].item()
chi2 = results['chi_squared']
fitted = results['fitted_trace']
experimental = results['experimental_trace']
```

---

## Common File Formats

### Structured Arrays (Parameter Storage)

Parameters are often stored as structured NumPy arrays to preserve names and types:

```python
# Creating a structured array
params = np.array([(14800.0, 2.5, 100, 0.3)],
                  dtype=[('E_0a', 'f8'), ('dielc', 'f8'),
                         ('N_ens', 'i4'), ('f_red', 'f8')])

# Saving
np.savez('parameters.npz', parameters=params)

# Loading and accessing
data = np.load('parameters.npz')
params = data['parameters']
E_0a = params['E_0a'][0]
```

---

## Loading and Accessing Data

### Basic Loading

```python
import numpy as np

# Load single array (.npy)
data = np.load('file.npy')

# Load multiple arrays (.npz)
data = np.load('file.npz')
print(data.files)  # List all keys
array1 = data['key1']
array2 = data['key2']
```

### Loading with Pickle (for complex objects)

Some files contain Python objects (dicts, lists) that require `allow_pickle=True`:

```python
data = np.load('file.npz', allow_pickle=True)
params_dict = data['parameters'].item()  # Convert 0-d array to dict
```

### Memory-Mapped Loading (for large files)

For very large files, use memory mapping to avoid loading entire file:

```python
data = np.load('large_file.npy', mmap_mode='r')  # Read-only memory map
subset = data[0:1000]  # Only loads requested portion
```

---

## Examples

### Example 1: Load and Plot Time-Resolved Data

```python
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.load('time_resolved_fluorescence.npz')
time = data['time']
wavelength = data['wavelength']
fl = data['fluorescence']

# Plot heatmap
plt.figure(figsize=(10, 6))
plt.pcolormesh(time, wavelength, fl.T, shading='auto', cmap='viridis')
plt.xlabel('Time (ns)')
plt.ylabel('Wavelength (nm)')
plt.colorbar(label='Fluorescence (a.u.)')
plt.title('Time-Resolved Fluorescence')
plt.show()

# Plot decay at specific wavelength
wl_idx = np.argmin(np.abs(wavelength - 700))  # Find 700 nm
plt.figure()
plt.plot(time, fl[:, wl_idx])
plt.xlabel('Time (ns)')
plt.ylabel('Fluorescence (a.u.)')
plt.title('Decay at 700 nm')
plt.show()
```

### Example 2: Compare Simulation with Experiment

```python
import numpy as np
import matplotlib.pyplot as plt

# Load experimental data
exp_data = np.load('NeededData/Experimental_Data/Fl_IsiA_monomer_300K_Gabriela_2025.npy')

# Load simulation results
sim_data = np.load('scripts/fluorescence_decay/model_1/Simulation/ensemble_average.npy')

# Load time axis (from separate file or create)
time = np.linspace(0, 10, len(sim_data))  # Assume 0-10 ns

# Plot comparison
plt.figure(figsize=(10, 6))
plt.plot(time, exp_data, 'o', label='Experimental', alpha=0.6)
plt.plot(time, sim_data, '-', label='Simulation', linewidth=2)
plt.xlabel('Time (ns)')
plt.ylabel('Fluorescence (a.u.)')
plt.legend()
plt.title('Simulation vs Experiment')
plt.show()

# Calculate RMSE
rmse = np.sqrt(np.mean((sim_data - exp_data)**2))
print(f'RMSE: {rmse:.4f}')
```

### Example 3: Analyze Ensemble Statistics

```python
import numpy as np
import matplotlib.pyplot as plt

# Load ensemble data
data = np.load('ensemble_output.npz')
traces = data['ensemble_traces']  # Shape: (n_ens, n_time, n_wl)
avg = data['ensemble_avg']
std = data['ensemble_std']

# Calculate statistics
mean_trace = np.mean(traces, axis=(0, 2))  # Average over ensemble and wavelength
std_trace = np.std(traces, axis=(0, 2))

# Plot with error bars
time = np.arange(len(mean_trace)) * 0.01  # Assume 0.01 ns steps
plt.figure(figsize=(10, 6))
plt.plot(time, mean_trace, 'k-', label='Mean')
plt.fill_between(time, mean_trace - std_trace, mean_trace + std_trace,
                 alpha=0.3, label='±1 σ')
plt.xlabel('Time (ns)')
plt.ylabel('Fluorescence (a.u.)')
plt.legend()
plt.title(f'Ensemble Average (N={traces.shape[0]})')
plt.show()
```

### Example 4: Inspect File Contents

```python
import numpy as np

def inspect_npz(filename):
    """Print structure of an NPZ file."""
    data = np.load(filename, allow_pickle=True)
    print(f"File: {filename}")
    print(f"Keys: {list(data.files)}")
    print()

    for key in data.files:
        arr = data[key]
        print(f"Key: {key}")
        print(f"  Type: {type(arr)}")
        print(f"  Dtype: {arr.dtype}")
        print(f"  Shape: {arr.shape}")
        if arr.size < 10:
            print(f"  Values: {arr}")
        print()

# Usage
inspect_npz('time_resolved_fluorescence.npz')
```

---

## Data Validation

### Checking Data Integrity

```python
import numpy as np

def validate_time_resolved(filename):
    """Validate time-resolved fluorescence data."""
    data = np.load(filename)

    # Check required keys
    required_keys = ['time', 'wavelength', 'fluorescence']
    for key in required_keys:
        assert key in data.files, f"Missing key: {key}"

    time = data['time']
    wl = data['wavelength']
    fl = data['fluorescence']

    # Check dimensions
    assert fl.shape[0] == len(time), "Time dimension mismatch"
    assert fl.shape[1] == len(wl), "Wavelength dimension mismatch"

    # Check for NaN or Inf
    assert not np.any(np.isnan(fl)), "Data contains NaN"
    assert not np.any(np.isinf(fl)), "Data contains Inf"

    # Check physical ranges
    assert np.all(time >= 0), "Negative time values"
    assert np.all(wl > 0), "Non-positive wavelengths"
    assert np.all(fl >= 0), "Negative fluorescence values"

    print(f"✓ File validated: {filename}")
    print(f"  Time range: {time[0]:.3f} - {time[-1]:.3f} ns")
    print(f"  Wavelength range: {wl[0]:.1f} - {wl[-1]:.1f} nm")
    print(f"  Fluorescence range: {fl.min():.3e} - {fl.max():.3e}")

# Usage
validate_time_resolved('time_resolved_fluorescence.npz')
```

---

## File Size Considerations

### Typical File Sizes

| File Type | Dimensions | Typical Size | Compressed Size |
|-----------|------------|--------------|-----------------|
| Experimental 1D | (1000,) | ~8 KB | ~2 KB |
| Hamiltonian 17×17 | (17, 17) | ~2 KB | ~1 KB |
| Time-resolved | (1000, 150) | ~1.2 MB | ~200 KB |
| Ensemble (1000) | (1000, 1000, 150) | ~1.2 GB | ~200 MB |

### Compression

Use `np.savez_compressed()` for better compression:

```python
# Save with compression
np.savez_compressed('output_compressed.npz',
                   time=time, wavelength=wl, fluorescence=fl)

# File will be significantly smaller with minimal speed penalty
```

---

## Best Practices

1. **Always use descriptive keys** in NPZ files
2. **Include metadata** (units, parameters, timestamps)
3. **Validate data** after loading
4. **Use compression** for large files
5. **Document array shapes** in comments
6. **Check data types** match expectations
7. **Use `allow_pickle=True`** cautiously (security risk with untrusted files)

---

## Troubleshooting

### Common Issues

**Issue:** `ValueError: Object arrays cannot be loaded when allow_pickle=False`
**Solution:** Add `allow_pickle=True` when loading:
```python
data = np.load('file.npz', allow_pickle=True)
```

**Issue:** `KeyError: 'array_name'`
**Solution:** Check available keys:
```python
data = np.load('file.npz')
print(data.files)  # List all available keys
```

**Issue:** File is very large and causes memory error
**Solution:** Use memory mapping:
```python
data = np.load('large_file.npy', mmap_mode='r')
```

---

*For additional questions or issues with file formats, please open a GitHub issue or consult the main README.*

*Last updated: January 26, 2025*
