# Data Dictionary

This document provides definitions and explanations for all parameters, abbreviations, and naming conventions used throughout the IsiA spectroscopy analysis project.

## Table of Contents
- [Parameter Abbreviations](#parameter-abbreviations)
- [File Naming Conventions](#file-naming-conventions)
- [Units and Standards](#units-and-standards)
- [Physical Constants](#physical-constants)
- [Data Structure Conventions](#data-structure-conventions)

---

## Parameter Abbreviations

### Hamiltonian and Energetics

| Parameter | Full Name | Description | Typical Range | Units |
|-----------|-----------|-------------|---------------|-------|
| `E_0a` | Reference Energy A | Reference energy for site energy calculations | 14500-15500 | cm⁻¹ |
| `E_0b` | Reference Energy B | Alternative reference energy for specific sites | 14500-15500 | cm⁻¹ |
| `dielc` | Dielectric Constant | Effective dielectric constant for coupling calculations | 2-10 | dimensionless |
| `V_ij` | Excitonic Coupling | Electronic coupling between sites i and j | -100 to +100 | cm⁻¹ |
| `H_ij` | Hamiltonian Element | Element of the excitonic Hamiltonian matrix | varies | cm⁻¹ |
| `eps_i` | Site Energy | Energy of excitation localized on site i | 14000-16000 | cm⁻¹ |

### Spectroscopy

| Parameter | Full Name | Description | Typical Range | Units |
|-----------|-----------|-------------|---------------|-------|
| `lambda` / `wl` | Wavelength | Wavelength of light | 600-750 | nm |
| `omega` / `freq` | Frequency | Frequency of light | 13000-17000 | cm⁻¹ |
| `sigma` | Linewidth | Full width at half maximum of spectral line | 50-200 | cm⁻¹ |
| `mu` | Transition Dipole | Transition dipole moment magnitude | 5-10 | Debye |
| `OD` | Optical Density | Absorbance (dimensionless) | 0-1 | a.u. |
| `Fl` | Fluorescence | Fluorescence intensity | 0-1 | a.u. |

### Dynamics and Kinetics

| Parameter | Full Name | Description | Typical Range | Units |
|-----------|-----------|-------------|---------------|-------|
| `t` | Time | Time variable for dynamics | 0-10 | ns |
| `dt` | Time Step | Time step for numerical integration | 0.001-0.1 | ns |
| `k_ET` | Energy Transfer Rate | Rate constant for energy transfer | 0.1-100 | ns⁻¹ |
| `k_decay` | Decay Rate | Fluorescence decay rate constant | 0.01-10 | ns⁻¹ |
| `tau` | Lifetime | Fluorescence lifetime (1/k_decay) | 0.1-100 | ns |
| `N_ens` | Ensemble Size | Number of trajectories in ensemble average | 100-10000 | dimensionless |

### Structural Parameters

| Parameter | Full Name | Description | Units |
|-----------|-----------|-------------|-------|
| `r_ij` | Distance | Distance between chromophores i and j | Å |
| `R` | Center-to-Center Distance | Distance between chromophore centers | Å |
| `theta` | Angle | Orientation angle between transition dipoles | degrees or radians |
| `kappa` | Orientation Factor | Geometric factor in Förster theory | dimensionless |

### Model-Specific Parameters

| Parameter | Full Name | Description | Model | Units |
|-----------|-----------|-------------|-------|-------|
| `f_red` | Red State Fraction | Fraction of population in red-emitting states | Model 1-4 | 0-1 |
| `k_q` | Quenching Rate | Rate constant for non-radiative quenching | All models | ns⁻¹ |
| `alpha` | Mixing Coefficient | Linear combination coefficient | Model 123 | dimensionless |
| `beta` | Heterogeneity Parameter | Disorder parameter for Gaussian distribution | All models | cm⁻¹ |

---

## File Naming Conventions

### Temperature Notation
- **Standard:** Always use capital K (Kelvin)
- **Examples:** `300K`, `77K`, `4K`
- **Incorrect:** `300k`, `300_k`, `300kelvin`

### Version Indicators
- **Date-based:** Use ISO format YYYY-MM-DD
  - Example: `data_2025-02-14.npy`
- **Sequential:** Use v1, v2, v3, etc.
  - Example: `analysis_v2.csv`
- **Avoid:** "new", "final", "latest" (ambiguous)

### Model Naming
- **Pattern:** `model_X` where X indicates the model type
  - `model_1`: Single red-state model
  - `model_2`: Two red-state model
  - `model_3`: Three red-state model
  - `model_4`: Four red-state model
  - `model_12`: Combined models 1+2
  - `model_13`: Combined models 1+3
  - `model_14`: Combined models 1+4
  - `model_123`: Combined models 1+2+3

### Data Type Suffixes
- `.npy` - NumPy binary array (single array)
- `.npz` - NumPy compressed archive (multiple arrays)
- `.csv` - Comma-separated values (human-readable)
- `.pdb` - Protein Data Bank structure file
- `.out` - MCCE output file

### Descriptor Components
Files typically follow the pattern: `[system]_[condition]_[temperature]_[source]_[version].[ext]`

**Example:** `Fl_IsiA_monomer_300K_Gabriela_2025.npy`
- `Fl` = Fluorescence
- `IsiA_monomer` = System being studied
- `300K` = Temperature
- `Gabriela` = Data source/collector
- `2025` = Version year
- `.npy` = NumPy binary format

---

## Units and Standards

### Energy Units
- **Primary:** cm⁻¹ (wavenumber)
- **Conversion:** 1 eV = 8065.54 cm⁻¹
- **Conversion:** E(cm⁻¹) = 10⁷ / λ(nm)

### Time Units
- **Primary:** nanoseconds (ns)
- **Conversion:** 1 ns = 10⁻⁹ s
- **Conversion:** 1 ps = 0.001 ns

### Length Units
- **Structural:** Ångström (Å)
- **Spectroscopy:** nanometers (nm)
- **Conversion:** 1 nm = 10 Å

### Temperature Units
- **Always:** Kelvin (K)
- **Conversion:** K = °C + 273.15

### Angle Units
- **Computational:** Radians
- **Output/Display:** Degrees
- **Conversion:** π radians = 180 degrees

---

## Physical Constants

Constants used throughout the calculations:

| Symbol | Name | Value | Units |
|--------|------|-------|-------|
| `h` | Planck's constant | 6.62607015 × 10⁻³⁴ | J·s |
| `c` | Speed of light | 2.99792458 × 10⁸ | m/s |
| `k_B` | Boltzmann constant | 1.380649 × 10⁻²³ | J/K |
| `N_A` | Avogadro's number | 6.02214076 × 10²³ | mol⁻¹ |

---

## Data Structure Conventions

### Hamiltonian Matrices
- **Format:** Square symmetric matrix
- **Dimensions:** N×N where N = number of chromophores
- **IsiA monomer:** 17×17 matrix (17 chlorophylls)
- **Diagonal elements:** Site energies (cm⁻¹)
- **Off-diagonal elements:** Excitonic couplings (cm⁻¹)

### Spectral Data Arrays
- **Wavelength axis:** Typically 600-750 nm with 1 nm resolution
- **Frequency axis:** Typically 13000-17000 cm⁻¹
- **Intensity:** Normalized to 0-1 (arbitrary units)

### Time-Resolved Data
- **Structure:** 2D or 3D arrays
- **Dimensions:** (time, wavelength) or (ensemble, time, wavelength)
- **Time axis:** Typically 0-10 ns with variable resolution
- **Storage:** Usually NPZ format with labeled arrays

### Ensemble Data
- **N_ens:** Number of independent trajectories
- **Storage:** Each trajectory stored separately or averaged
- **Typical size:** 1000-10000 trajectories for good statistics

---

## Special Notations

### pH Nomenclature
- `pH7` or `pH_7` = pH 7.0 (neutral)
- Used in structure files to indicate protonation state

### Chromophore Numbering
- **Convention:** 1-indexed (first chromophore is #1, not #0)
- **IsiA monomer:** Chromophores 1-17
- **Matching:** Corresponds to PDB residue numbering

### State Labels
- `S0` = Ground state
- `S1` = First excited state
- `S2` = Second excited state
- `red1`, `red2`, etc. = Red-shifted states (low energy)
- `blue` = Blue-shifted states (high energy)

---

## Common Variable Names in Code

### Arrays and Matrices
- `H` = Hamiltonian matrix
- `E` = Eigenvalues (energies)
- `U` = Eigenvectors (transformation matrix)
- `rho` = Density matrix
- `pop` = Population array

### Indices and Counters
- `i`, `j`, `k` = Site indices
- `t_idx` = Time index
- `wl_idx` = Wavelength index
- `n_ens` = Ensemble member index

### File Paths
- `base_dir` = Base directory path
- `data_dir` = Data directory path
- `output_dir` = Output directory path
- `pdb_file` = PDB structure file path

---

## Version Control

- **Last updated:** 2025-01-26
- **Applies to:** All data files and scripts in the repository
- **Maintained by:** Mohamed A. A. Elrefaiy

For questions or clarifications, please open an issue on GitHub or contact the authors.
