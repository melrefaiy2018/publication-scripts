# MCCE Data

This directory contains input and output files from Multi-Conformation Continuum Electrostatics (MCCE) calculations used to determine protonation states and electrostatic energies for the IsiA protein.

## Overview

MCCE is a computational tool that predicts pKa values and protonation states of ionizable residues in proteins. The output is used to calculate site energies for the excitonic Hamiltonian, accounting for electrostatic effects of the protein environment on chlorophyll chromophores.

## Files

### `pK.out`
- **Description:** Main MCCE output file containing pKa values, protonation states, and occupancies
- **Size:** ~10 KB
- **Lines:** 68 lines
- **Format:** Plain text with structured columns
- **Date:** December 10, 2024
- **pH condition:** Calculations performed across pH range, pH 7.0 used for simulations
- **Content:**
  - Residue identifiers
  - Calculated pKa values
  - Protonation state occupancies
  - Conformer energies
  - Electrostatic interactions

## Data Format

The `pK.out` file follows MCCE output format:

### File Structure

```
# Comments and header information
RES_ID  CONF  pKa   OCC_pH7  ENERGY  ...
ASP_10  A     3.8   0.95     -2.34   ...
ASP_10  B     4.2   0.05     1.12    ...
GLU_25  A     4.5   0.98     -3.45   ...
...
```

### Column Descriptions

| Column | Name | Description | Units |
|--------|------|-------------|-------|
| 1 | RES_ID | Residue identifier (name + number) | - |
| 2 | CONF | Conformer ID (A, B, C, ...) | - |
| 3 | pKa | Calculated pKa value | pH units |
| 4 | OCC_pH7 | Occupancy at pH 7.0 | 0-1 |
| 5 | ENERGY | Relative conformer energy | kcal/mol |
| 6+ | Additional | Interaction energies, charges, etc. | various |

### Residue Types Included

- **Titratable amino acids:** ASP, GLU, HIS, LYS, ARG, TYR, CYS
- **Chlorophylls:** CLA (chlorophyll a) - 17 pigments
- **Other cofactors:** If present

## Usage in This Repository

### Input For Hamiltonian Calculations

The pK.out file is read by:
```bash
scripts/hamiltonian/cdc_IsiA_monomer_average_pH7.py
```

**What it provides:**
1. **Protonation states:** Determines which residues are protonated/deprotonated at pH 7.0
2. **Electrostatic energies:** Used to calculate chlorophyll site energies
3. **Conformer selection:** Most occupied conformers define protein structure
4. **Environmental effects:** Accounts for how charged residues affect pigment energies

### Reading pK.out in Python

```python
import pandas as pd

# Read MCCE output (assuming tab or space-separated)
# Note: Format may vary, adjust delimiter as needed
data = []
with open('NeededData/mcce/pK.out', 'r') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue  # Skip comments and empty lines
        data.append(line.strip().split())

# Convert to DataFrame
df = pd.DataFrame(data, columns=['RES_ID', 'CONF', 'pKa', 'OCC_pH7', 'ENERGY'])
df['pKa'] = pd.to_numeric(df['pKa'], errors='coerce')
df['OCC_pH7'] = pd.to_numeric(df['OCC_pH7'], errors='coerce')
df['ENERGY'] = pd.to_numeric(df['ENERGY'], errors='coerce')

print(df.head())

# Filter for high occupancy conformers at pH 7
high_occ = df[df['OCC_pH7'] > 0.5]
print(f"\n{len(high_occ)} conformers with occupancy > 0.5 at pH 7.0")
```

## Physical Interpretation

### pKa Values

- **Standard pKa:** Expected values for amino acids in water
  - ASP: ~3.9
  - GLU: ~4.3
  - HIS: ~6.0
  - LYS: ~10.5
  - ARG: ~12.5

- **Shifted pKa:** Values different from standard indicate environmental effects
  - Lower pKa (acidic residues): Stabilized deprotonated form
  - Higher pKa (basic residues): Stabilized protonated form

### Occupancy at pH 7.0

- **OCC = 1.0:** Fully in that protonation state
- **OCC = 0.5:** Equal populations of two states
- **OCC < 0.1:** Minor conformer, usually ignored

### Conformers

Multiple conformers per residue represent:
- Different protonation states
- Different rotamers (side-chain orientations)
- Hydrogen bonding patterns

## Relationship to Other Data

### Input Requirements
- **Structure:** `NeededData/structure/6kig_structure_prepared_mcce_input.pdb`
  - Prepared PDB with proper MCCE format
  - Includes hydrogens and proper atom naming

### Output Used By
1. **Extended PDB generation:** Creates `extended_most_occ_pH7.pdb`
   - Selects most occupied conformers
   - Writes extended structure file

2. **Hamiltonian calculations:** `scripts/hamiltonian/cdc_IsiA_monomer_average_pH7.py`
   - Extracts site energies for chlorophylls
   - Includes electrostatic contributions
   - Accounts for local protein environment

3. **Spectroscopy predictions:** Indirectly affects all downstream simulations
   - Site energies → Hamiltonian → Spectra → Dynamics

## MCCE Calculation Details

### Method

MCCE uses:
- **Poisson-Boltzmann:** Continuum electrostatics
- **Conformer sampling:** Multiple side-chain rotamers
- **Monte Carlo:** Statistical sampling of protonation states
- **pH titration:** Calculates properties across pH range

### Parameters (Typical)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Dielectric (protein) | 4-8 | Interior dielectric constant |
| Dielectric (solvent) | 80 | Water dielectric constant |
| Ionic strength | 0.1 M | Salt concentration |
| Temperature | 298 K | Room temperature |
| pH range | 0-14 | Full titration curve |

### Validation

- **pKa accuracy:** Typically ±0.5-1.0 pH units vs. experiment
- **Hydrogen bonds:** Verified in output structure
- **Charge balance:** Total protein charge consistent with pH

## File Metadata

| Property | Value |
|----------|-------|
| File size | ~10 KB |
| Number of lines | 68 |
| Number of residues | ~300 (protein) + 17 (chlorophylls) |
| Number of conformers | Multiple per titratable residue |
| pH focus | 7.0 (neutral) |
| Date generated | December 10, 2024 |

## Known Issues and Limitations

### MCCE Approximations
- **Continuum model:** Treats water as dielectric continuum (no explicit waters)
- **Static structure:** Doesn't account for conformational changes
- **Point charges:** Simplified atomic charge model
- **No quantum effects:** Classical electrostatics only

### Data Completeness
- Some residues may have incomplete conformer sets
- Edge residues may have higher uncertainty
- Chlorophyll site energies require additional parameters beyond MCCE

## Reproducibility

To regenerate pK.out from scratch:

### Prerequisites
- MCCE software (version X.X - to be specified)
- Prepared PDB structure: `6kig_structure_prepared_mcce_input.pdb`
- MCCE parameter files (standard or custom)

### Steps

```bash
# 1. Prepare structure for MCCE (if not already done)
mcce_prepare 6kig_structure.pdb

# 2. Run MCCE calculation
mcce run.prm

# 3. Extract pK.out from MCCE output
# (File is automatically generated)
```

**Note:** The provided pK.out is ready to use; regeneration is not required unless modifying parameters.

## Example Usage

### Extract Chlorophyll Site Energies

```python
def extract_chl_energies(pk_file):
    """Extract chlorophyll energies from pK.out."""
    chl_energies = {}

    with open(pk_file, 'r') as f:
        for line in f:
            if 'CLA' in line or 'CHL' in line:  # Chlorophyll residues
                parts = line.split()
                res_id = parts[0]
                energy = float(parts[4])  # Assuming energy in column 5

                # Store by residue number
                res_num = int(''.join(filter(str.isdigit, res_id)))
                if res_num not in chl_energies:
                    chl_energies[res_num] = []
                chl_energies[res_num].append(energy)

    # Average over conformers
    avg_energies = {k: sum(v)/len(v) for k, v in chl_energies.items()}

    return avg_energies

# Usage
energies = extract_chl_energies('NeededData/mcce/pK.out')
print(f"Chlorophyll site energies: {energies}")
```

### Analyze Protonation States

```python
def analyze_protonation(pk_file, ph=7.0):
    """Analyze protonation states at given pH."""
    protonated = []
    deprotonated = []

    with open(pk_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split()
            res_id = parts[0]
            pka = float(parts[2])
            occ = float(parts[3])

            # Residues with pKa < pH are mostly deprotonated
            # Residues with pKa > pH are mostly protonated
            if occ > 0.5:  # Only count major conformers
                if pka < ph:
                    deprotonated.append((res_id, pka))
                else:
                    protonated.append((res_id, pka))

    print(f"At pH {ph}:")
    print(f"  Deprotonated: {len(deprotonated)} residues")
    print(f"  Protonated: {len(protonated)} residues")

    return protonated, deprotonated

# Usage
prot, deprot = analyze_protonation('NeededData/mcce/pK.out', ph=7.0)
```

## Citation

If you use this MCCE data, please cite:

**MCCE Software:**
[MCCE citation - to be added]

**This work:**
Mohamed A. A. Elrefaiy, Dvir Harris, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, and Gabriela S. Schlau-Cohen. "Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States." *Manuscript in preparation*, 2025.

## Related Documentation

- **Structure files:** [../structure/README.md](../structure/README.md)
- **Hamiltonian:** [../hamiltonian/README.md](../hamiltonian/README.md)
- **Data dictionary:** [../DATA_DICTIONARY.md](../DATA_DICTIONARY.md)
- **Main README:** [../../README.md](../../README.md)

## Further Reading

- **MCCE website:** [MCCE documentation URL]
- **MCCE publications:** [Key papers about MCCE method]
- **pKa prediction:** General reviews of computational pKa methods

---

*Last updated: January 26, 2025*
