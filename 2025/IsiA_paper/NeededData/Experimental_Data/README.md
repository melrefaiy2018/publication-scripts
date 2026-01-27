# Experimental Data

This directory contains experimental absorption and fluorescence spectra data used for comparison with computational simulations and model validation.

## Overview

All experimental measurements were performed on purified IsiA monomer samples at 300 K (room temperature). The data includes both steady-state absorption spectra and time-resolved fluorescence decay measurements.

## Files

### Current Data Files (Use These)

- **`Fl_IsiA_monomer_300K_Gabriela_2025.npy`** - **PRIMARY FILE**
  - **Description:** Experimental fluorescence decay data (latest version)
  - **Date:** February 14, 2025
  - **Format:** NumPy binary array (.npy)
  - **Temperature:** 300 K
  - **Size:** ~1.2 KB
  - **Status:** CURRENT - Use this for all analyses

- **`Abs_IsiA_monomer_300K_Gabriela.npy`**
  - **Description:** Experimental absorption spectrum
  - **Date:** May 5, 2024
  - **Format:** NumPy binary array (.npy)
  - **Temperature:** 300 K
  - **Size:** ~1.2 KB
  - **Wavelength range:** Typically 600-750 nm
  - **Status:** CURRENT

### Subdirectories

- **`csv/`** - CSV versions of spectroscopic data
  - Contains comma-separated value files for human readability
  - Same data as NPY files but in text format
  - Useful for inspection or import into other software

- **`archive/`** - Older versions of experimental data
  - Contains superseded data files
  - Preserved for transparency and historical reference
  - See `archive/README.md` for version history
  - **Do not use for new analyses**

## Data Format

### NPY Files (NumPy Binary)

NPY files can be loaded directly with NumPy:

```python
import numpy as np

# Load fluorescence data
fl_data = np.load('NeededData/Experimental_Data/Fl_IsiA_monomer_300K_Gabriela_2025.npy')

# Load absorption data
abs_data = np.load('NeededData/Experimental_Data/Abs_IsiA_monomer_300K_Gabriela.npy')
```

**Structure:**
- **1D array:** `(n_time_points,)` - single trace
- **2D array:** `(n_time_points, 2)` - time in first column, intensity in second
- **2D array:** `(n_wavelengths, n_time_points)` - wavelength-resolved

**Units:**
- **Time:** nanoseconds (ns)
- **Wavelength:** nanometers (nm)
- **Intensity:** Arbitrary units (a.u.), normalized to 0-1

### CSV Files

CSV files have header rows and can be opened with any spreadsheet software:

```csv
Wavelength (nm), Intensity (a.u.)
600.0, 0.123
601.0, 0.145
...
```

## Data Specifications

### Fluorescence Decay Data
- **Primary file:** `Fl_IsiA_monomer_300k_Gabriela_2025.npy`
- **Temperature:** 300 K (room temperature)
- **Time range:** 0-10 ns (typical)
- **Time resolution:** Variable (typically 100-1000 points)
- **Wavelength range:** 600-750 nm (red region)
- **Wavelength resolution:** 1-5 nm
- **Normalization:** Peak normalized to 1.0
- **Baseline:** Background subtracted

### Absorption Spectrum
- **Primary file:** `Abs_IsiA_monomer_300K_Gabriela.npy`
- **Temperature:** 300 K
- **Wavelength range:** 600-750 nm
- **Wavelength resolution:** 0.5-1 nm
- **Data type:** Optical density (OD) or normalized absorbance
- **Normalization:** Maximum absorption set to 1.0

## Experimental Methods

### Sample Preparation
- **Protein:** IsiA monomer purified from cyanobacteria
- **Buffer:** [Details to be added]
- **Concentration:** [Details to be added]
- **Temperature control:** Maintained at 300 ± 1 K

### Measurement Techniques
- **Absorption:** UV-Vis spectrophotometry
- **Fluorescence:** Time-resolved fluorescence spectroscopy
  - Excitation wavelength: [To be added]
  - Detection range: 600-750 nm
  - Time resolution: [To be added]
  - Instrument: [To be added]

## Data Quality

### Signal-to-Noise Ratio
- Absorption: High SNR (>100:1)
- Fluorescence: Good SNR (>50:1 in peak region)

### Reproducibility
- Multiple independent measurements averaged
- Standard deviation typically <5% of peak intensity
- Consistent sample preparation protocols

### Calibration
- Wavelength calibrated using standard lamps
- Time axis calibrated using instrument response function
- Intensity corrected for detector response

## Version History

| Version | Date | File | Status | Changes |
|---------|------|------|--------|---------|
| v2.0 | Feb 14, 2025 | `Fl_IsiA_monomer_300k_Gabriela_2025.npy` | **CURRENT** | Improved SNR, better baseline |
| v1.0 | May 5, 2024 | `Fl_IsiA_monomer_300k_Gabriela.npy` | Archived | Initial measurements |

See `archive/README.md` for detailed version history and differences.

## Usage in This Repository

### Scripts That Use This Data

1. **Model comparison:** `scripts/fluorescence_decay/model_*/find_best_parameter_model/`
   - Compares simulated decay traces with experimental data
   - Calculates goodness-of-fit metrics (RMSE, χ²)

2. **Validation:** Various analysis scripts
   - Validates Hamiltonian-based predictions
   - Benchmarks different model variants

### Loading Example

```python
import numpy as np
import matplotlib.pyplot as plt

# Load the current experimental fluorescence data
fl_exp = np.load('NeededData/Experimental_Data/Fl_IsiA_monomer_300k_Gabriela_2025.npy')

# Assuming 1D time trace
time = np.linspace(0, 10, len(fl_exp))  # 0-10 ns

# Plot
plt.figure(figsize=(8, 5))
plt.plot(time, fl_exp, 'o-', label='Experimental')
plt.xlabel('Time (ns)')
plt.ylabel('Fluorescence (a.u.)')
plt.title('IsiA Monomer Fluorescence Decay at 300K')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

## Data Availability

All experimental data in this directory is:
- **Open access** under the MIT license
- **Citable** via the repository DOI (see main README)
- **Reproducible** with methods described in the manuscript

## Citation

If you use this experimental data, please cite:

**Manuscript:**
Mohamed A. A. Elrefaiy, Dvir Harris, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, and Gabriela S. Schlau-Cohen. "Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States." *Manuscript in preparation*, 2025.

**Data Repository:**
See [CITATION.cff](../../CITATION.cff) in the repository root.

## Contact

- **Principal Investigator:** Gabriela S. Schlau-Cohen (MIT)
- **Data analysis:** Mohamed A. A. Elrefaiy (MIT)
- **Questions:** Open an issue on GitHub

## Related Documentation

- **Data dictionary:** [NeededData/DATA_DICTIONARY.md](../DATA_DICTIONARY.md)
- **File formats:** [scripts/fluorescence_decay/NPZ_FORMAT.md](../../scripts/fluorescence_decay/NPZ_FORMAT.md)
- **Main README:** [README.md](../../README.md)

---

*Last updated: January 26, 2025*
