# Archived Experimental Data

This directory contains older versions of the experimental fluorescence data files that have been superseded by more recent versions.

## Files in Archive

### Fl_IsiA_monomer_300k_Gabriela.npy
- **Date:** May 5, 2024
- **Status:** DEPRECATED - Use the 2025 version instead
- **Size:** ~1.0 KB
- **Description:** Original fluorescence decay data for IsiA monomer at 300K
- **Reason for archiving:** Replaced by updated measurements with improved calibration

### Fl_IsiA_monomer_300K_Gabriela_new.csv
- **Date:** February 14, 2025
- **Status:** DEPRECATED - CSV format superseded by NPY format
- **Size:** ~1.7 KB
- **Format:** CSV (plain text)
- **Description:** Intermediate version of fluorescence data in CSV format
- **Reason for archiving:** NPY format provides better numerical precision and is more efficient

## Current Version

The current, canonical version of the experimental data is:
- **File:** `Fl_IsiA_monomer_300k_Gabriela_2025.npy` (in parent directory)
- **Date:** February 14, 2025
- **Format:** NumPy binary (.npy)
- **Description:** Latest fluorescence decay measurements with improved experimental setup

## Version History

1. **v1.0 (May 2024):** Initial experimental measurements
2. **v2.0 (February 2025):** Updated measurements with:
   - Improved signal-to-noise ratio
   - Better baseline correction
   - Enhanced calibration procedures
   - Switch to binary NPY format for precision

## Usage Notes

- **For reproducing published results:** Use `Fl_IsiA_monomer_300k_Gabriela_2025.npy` from the parent directory
- **For historical reference:** These archived files are preserved for transparency and reproducibility of earlier analyses
- **Data format:** All .npy files can be loaded with `numpy.load()`

## Contact

For questions about data versions or differences between files, please contact:
- Gabriela S. Schlau-Cohen (MIT)
- Mohamed A. A. Elrefaiy (MIT)
