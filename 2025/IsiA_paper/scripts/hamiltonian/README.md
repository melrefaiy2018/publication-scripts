# Hamiltonian Calculation Scripts

This directory contains scripts for calculating the excitonic Hamiltonian, site energies, and spectral predictions for the IsiA protein complex.

## Overview

The main script calculates the electronic coupling between chlorophyll pigments to construct the excitonic Hamiltonian matrix. This matrix describes all excitonic interactions and forms the foundation for fluorescence decay simulations.

## Files

### Main Script

- **`cdc_IsiA_monomer_average_pH7.py`**
  - **Purpose:** Calculate excitonic Hamiltonian and absorption/fluorescence spectra
  - **Type:** Sequential calculation (single-threaded)
  - **Input:** Protein structure (PDB) and MCCE electrostatics data
  - **Output:** Hamiltonian matrix and spectra (CSV format)
  - **Runtime:** 30 seconds - 2 minutes
  - **Memory:** ~500 MB
  - **Parallelization:** None (I/O bound)

## Input Requirements

### Data Files Required

| File | Location | Description | Required |
|------|----------|-------------|----------|
| PDB structure | `NeededData/structure/extended_most_occ_pH7.pdb` | Protein atomic coordinates | Yes |
| MCCE output | `NeededData/mcce/pK.out` | Site energies and pKa values | Yes |

### Parameter Settings

Inside the script, key parameters can be adjusted:

| Parameter | Type | Default | Description | Units |
|-----------|------|---------|-------------|-------|
| `E_0a` | float | 14800 | Reference energy A | cm⁻¹ |
| `E_0b` | float | 14800 | Reference energy B | cm⁻¹ |
| `dielc` | float | 3.0 | Effective dielectric constant | dimensionless |
| `T` | float | 300 | Temperature | K |
| `linewidth` | float | 100 | Spectral broadening (FWHM) | cm⁻¹ |

## Usage

### Basic Execution

```bash
cd scripts/hamiltonian
python cdc_IsiA_monomer_average_pH7.py
```

### With Custom Parameters

Edit `cdc_IsiA_monomer_average_pH7.py` to modify parameters, then run:

```bash
python cdc_IsiA_monomer_average_pH7.py
```

### Parallel Execution (if modified for multi-threading)

```bash
# Note: Current implementation is single-threaded
# Future versions may support OpenMP parallelization
python cdc_IsiA_monomer_average_pH7.py
```

## Output Files

### Generated Files

All outputs are saved to `NeededData/hamiltonian/`:

```
NeededData/hamiltonian/
├── hamiltonian_matrix.csv           # 17×17 Hamiltonian matrix (cm⁻¹)
├── absorption_spectrum.csv          # Calculated absorption spectrum
├── fluorescence_spectrum.csv        # Calculated fluorescence emission
├── Hamiltonian/                     # Detailed Hamiltonian data
│   ├── hamiltonian_eigenvalues.csv  # Exciton state energies
│   ├── hamiltonian_eigenvectors.csv # Exciton state compositions
│   └── [other matrices]
├── SiteEnergy_Data/                 # Per-pigment analysis
│   ├── site_energies.csv            # Site energy values
│   ├── site_variations.csv          # Energy distributions
│   └── [other energy data]
└── Spectra_Data/                    # Raw spectral calculations
    ├── oscillator_strengths.csv     # Transition strengths
    ├── transitions.csv              # Electronic transitions
    └── [spectral components]
```

### File Format Details

**Hamiltonian Matrix (hamiltonian_matrix.csv):**
- Square 17×17 symmetric matrix
- Units: cm⁻¹ (wavenumbers)
- Diagonal: Site energies (~14500-15500 cm⁻¹)
- Off-diagonal: Excitonic couplings (~±100 cm⁻¹)

**Absorption Spectrum (absorption_spectrum.csv):**
- Two columns: Wavelength (nm) and Intensity (a.u.)
- Wavelength range: ~600-750 nm
- Peak position: ~680 nm (typical)
- Normalized to 1.0

**Fluorescence Spectrum (fluorescence_spectrum.csv):**
- Same format as absorption
- Red-shifted relative to absorption (Stokes shift)
- Peak position: ~698 nm (typical)

## Computational Requirements

### Runtime Estimates

- **Typical execution:** 1-2 minutes
- **Data loading:** ~10 seconds
- **Matrix construction:** ~30 seconds
- **Eigenvalue decomposition:** ~10 seconds
- **Spectral calculation:** ~20 seconds
- **I/O and output:** ~30 seconds

### Memory Usage

- **Peak memory:** ~500-800 MB
- **Hamiltonian storage:** ~2 KB
- **Full matrices:** ~20 KB
- **Spectra data:** ~40 KB

### Processor Requirements

- **CPU cores:** 1 (single-threaded)
- **CPU speed:** Standard (no special requirements)
- **GPU acceleration:** Not implemented

## Expected Results

### Validation Checklist

Run the script and verify:

1. **Hamiltonian Matrix**
   - [ ] Shape is 17×17
   - [ ] Symmetric (H = H^T)
   - [ ] Diagonal values between 14500-15500 cm⁻¹
   - [ ] Off-diagonal couplings between ±150 cm⁻¹
   - [ ] No NaN or Inf values

2. **Absorption Spectrum**
   - [ ] Peak wavelength near 680 nm
   - [ ] Bandwidth ~50-100 nm
   - [ ] Smooth, continuous function
   - [ ] Normalized between 0-1

3. **Fluorescence Spectrum**
   - [ ] Peak wavelength near 700 nm (red-shifted)
   - [ ] Stokes shift ~20 nm
   - [ ] Similar bandwidth to absorption
   - [ ] Normalized between 0-1

### Example Output Values

```
Hamiltonian diagonal (site energies) [cm⁻¹]:
  [14823.4, 14921.3, 14756.8, 15102.1, 14834.2, ...]

Hamiltonian off-diagonal samples (couplings) [cm⁻¹]:
  [-45.3, 12.1, -23.4, 67.8, -12.5, ...]

Absorption spectrum:
  Peak: 680.3 nm
  Bandwidth: 62 nm (FWHM)
  Maximum intensity: 1.000

Fluorescence spectrum:
  Peak: 698.7 nm
  Stokes shift: 18.4 nm
  Bandwidth: 65 nm (FWHM)
```

## Troubleshooting

### Common Issues

**Issue:** `FileNotFoundError: Cannot find 'NeededData/structure/...'`
- **Cause:** Running script from wrong directory
- **Solution:** Run from repository root or use absolute paths

**Issue:** `ValueError: Array shapes don't match`
- **Cause:** Inconsistent structure or MCCE data
- **Solution:** Verify PDB and pK.out files are not corrupted

**Issue:** `NaN values in output`
- **Cause:** Invalid parameter values (e.g., dielc=0)
- **Solution:** Check parameter values in script, use physically reasonable ranges

**Issue:** Slow performance (>5 minutes)`
- **Cause:** Slow disk I/O or system load
- **Solution:** Check disk space, close other applications, ensure SSD if possible

**Issue:** Absorption peak at wrong wavelength
- **Cause:** Incorrect E_0a or dielc parameters
- **Solution:** Adjust parameters based on experimental data

### Debug Mode

To enable verbose output, modify the script to add debug prints:

```python
print(f"DEBUG: Loading structure from {pdb_file}")
print(f"DEBUG: Hamiltonian shape: {H.shape}")
print(f"DEBUG: Site energy range: {H.diag().min():.1f}-{H.diag().max():.1f} cm⁻¹")
print(f"DEBUG: Coupling range: {H[H!=np.diag(np.diag(H))].min():.1f} to {H[H!=np.diag(np.diag(H))].max():.1f} cm⁻¹")
```

## Performance Optimization

### Speeding Up Calculations

1. **Use SSD storage:** I/O is significant bottleneck
2. **Reduce spectral resolution:** Fewer wavelength points = faster
3. **Compile with optimizations:** Use MKL-linked numpy/scipy

### Memory Optimization

- Current implementation uses minimal memory
- Not suitable for GPU acceleration (matrix too small)
- Parallelization not beneficial due to small problem size

## Relationship to Fluorescence Simulations

### Hamiltonian → Fluorescence Decay

This script's outputs feed into fluorescence decay simulations:

1. **Hamiltonian matrix** → Energy landscape for dynamics
2. **Site energies** → Individual pigment contributions
3. **Spectra** → Validation of calculated parameters

### Integration

```
hamiltonian_matrix.csv
         ↓
scripts/fluorescence_decay/model_*/run_model_X.py
         ↓
Time-resolved fluorescence predictions
```

## Citation

If you use this Hamiltonian calculation, please cite:

Mohamed A. A. Elrefaiy, Dvir Harris, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, and Gabriela S. Schlau-Cohen. "Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States." *Manuscript in preparation*, 2025.

## Related Documentation

- **Hamiltonian data:** [../../NeededData/hamiltonian/README.md](../../NeededData/hamiltonian/README.md)
- **Structure files:** [../../NeededData/structure/README.md](../../NeededData/structure/README.md)
- **MCCE data:** [../../NeededData/mcce/README.md](../../NeededData/mcce/README.md)
- **Parameter dictionary:** [../../NeededData/DATA_DICTIONARY.md](../../NeededData/DATA_DICTIONARY.md)
- **Main README:** [../../README.md](../../README.md)

---

*Last updated: January 26, 2025*
