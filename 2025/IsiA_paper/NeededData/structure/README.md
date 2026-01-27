# Protein Structure Files

This directory contains PDB (Protein Data Bank) structure files for the IsiA protein used in all computational simulations.

## Overview

The IsiA protein structures serve as the foundation for Hamiltonian calculations and fluorescence decay simulations. These files contain atomic coordinates, residue information, and chromophore positions.

## Files

### `6kig_structure_prepared_mcce_input.pdb`
- **Description:** Original PDB structure prepared as input for MCCE calculations
- **PDB ID:** 6KIG (cryoEM structure)
- **Resolution:** Cryo-EM map
- **Source:** Protein Data Bank (www.rcsb.org)
- **Preparation:** Converted to MCCE-compatible format
- **Size:** ~332 KB
- **Date:** December 10, 2024
- **Content:**
  - Full protein structure with hydrogen atoms
  - All 17 chlorophyll chromophores
  - Protein backbone and sidechains
  - Formatted for MCCE electrostatics calculations

### `extended_most_occ_pH7.pdb`
- **Description:** Extended structure with most probable protonation states at pH 7.0
- **Derivation:** Generated from MCCE calculations
- **Size:** ~708 KB
- **Date:** December 10, 2024
- **Status:** **PRIMARY FILE** - Use this for all simulations
- **Content:**
  - Complete atomic coordinates
  - Hydrogens added based on pH 7.0 protonation
  - Most occupied conformers selected
  - Includes all chromophores with proper geometry
- **Usage:** Input for Hamiltonian and fluorescence decay calculations

## Data Format

Both files follow standard PDB format with extensions for MCCE compatibility:

### Standard PDB Records

```
ATOM      1  N   MET A   1      10.123  20.456  30.789  1.00 50.00           N
ATOM      2  CA  MET A   1      11.234  21.567  31.890  1.00 50.00           C
...
HETATM 1234  MG  CLA A 101     45.678  12.345  67.890  1.00 40.00          MG
```

**Standard PDB columns:**
- Columns 1-6: Record type (ATOM, HETATM, etc.)
- Columns 7-11: Atom serial number
- Columns 13-16: Atom name
- Column 17: Alternate location indicator
- Columns 18-20: Residue name
- Column 22: Chain identifier
- Columns 23-26: Residue sequence number
- Columns 31-38: X coordinate (Å)
- Columns 39-46: Y coordinate (Å)
- Columns 47-54: Z coordinate (Å)
- Columns 55-60: Occupancy
- Columns 61-66: Temperature factor (B-factor)
- Columns 77-78: Element symbol

### IsiA Monomer Details

- **Protein chains:** Single chain (typically chain A)
- **Residues:** ~300 amino acids
- **Chromophores:** 17 chlorophyll a molecules
  - Residue names: CLA, CHL, or similar MCCE nomenclature
  - Include Mg²⁺ central ion
  - Include phytol tail
- **Total atoms:** ~10,000 (including hydrogens)

## Chlorophyll Numbering

The 17 chlorophyll molecules in IsiA monomer are numbered:
- **Chl 1-17:** Sequential numbering from N- to C-terminus
- **Coordinates:** Extracted from cryo-EM structure 6KIG
- **Geometry:** Optimized for MCCE calculations

## Relationship to Other Data

### Input For
1. **MCCE calculations:** `6kig_structure_prepared_mcce_input.pdb` → MCCE software
2. **Hamiltonian:** `extended_most_occ_pH7.pdb` → `scripts/hamiltonian/cdc_IsiA_monomer_average_pH7.py`
3. **Chromophore extraction:** Both files used to identify pigment positions

### Derived From
- **Original source:** PDB entry 6KIG
- **Reference:** [Citation to original structure paper]
- **Modifications:**
  - Hydrogens added using MCCE
  - Protonation states optimized for pH 7.0
  - Missing atoms/residues modeled
  - Conformers selected based on occupancy

### Related Files
- **MCCE output:** `NeededData/mcce/pK.out` (derived from this structure)
- **Hamiltonian:** `NeededData/hamiltonian/` (uses chromophore coordinates)

## Usage

### Loading PDB in Python

```python
# Using BioPython
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('IsiA', 'NeededData/structure/extended_most_occ_pH7.pdb')

# Access atoms
for model in structure:
    for chain in model:
        for residue in chain:
            print(f"Residue: {residue.get_resname()} {residue.id[1]}")
            for atom in residue:
                print(f"  Atom: {atom.name} at {atom.coord}")
```

### Extracting Chlorophyll Coordinates

```python
# Extract Mg positions of chlorophylls
import numpy as np
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('IsiA', 'NeededData/structure/extended_most_occ_pH7.pdb')

mg_coords = []
for residue in structure.get_residues():
    if residue.get_resname() in ['CLA', 'CHL']:  # Chlorophyll residues
        for atom in residue:
            if atom.name == 'MG':  # Magnesium central ion
                mg_coords.append(atom.coord)

mg_coords = np.array(mg_coords)
print(f"Found {len(mg_coords)} chlorophylls")
print(f"Coordinates shape: {mg_coords.shape}")
```

### Visualizing Structure

**PyMOL:**
```python
# In PyMOL
load NeededData/structure/extended_most_occ_pH7.pdb
show cartoon
select chlorophylls, resn CLA+CHL
show sticks, chlorophylls
color green, chlorophylls
```

**VMD:**
```tcl
# In VMD
mol new NeededData/structure/extended_most_occ_pH7.pdb
mol modstyle 0 0 NewCartoon
mol representation Licorice
mol selection "resname CLA CHL"
mol addrep 0
```

## Structure Quality

### Validation Metrics
- **Cryo-EM resolution:** ~3-4 Å (original 6KIG)
- **Ramachandran:** >95% in favored regions (expected)
- **Clashes:** Minimal after MCCE optimization
- **Geometry:** Standard bond lengths and angles

### Known Issues
- Hydrogen positions are modeled (not in original cryo-EM data)
- Some loop regions may have higher uncertainty
- Phytol tails may have multiple conformations

## File Metadata

| File | Size | Atoms | Residues | Chromophores | Date |
|------|------|-------|----------|--------------|------|
| `6kig_structure_prepared_mcce_input.pdb` | 332 KB | ~10,000 | ~300 | 17 | Dec 10, 2024 |
| `extended_most_occ_pH7.pdb` | 708 KB | ~15,000 | ~300 | 17 | Dec 10, 2024 |

## Protonation States

The `extended_most_occ_pH7.pdb` file includes:
- **pH condition:** 7.0 (neutral)
- **Titratable residues:** ASP, GLU (deprotonated); LYS, ARG (protonated); HIS (varies)
- **Selection criterion:** Most occupied conformer from MCCE
- **Hydrogen placement:** Optimized by MCCE electrostatics

## Reproducibility

To regenerate the extended PDB from scratch:

1. Start with original PDB 6KIG from RCSB
2. Prepare for MCCE (add hydrogens, format)
3. Run MCCE calculations at pH 7.0
4. Extract most occupied conformers
5. Generate extended PDB output

**Note:** The provided files are ready to use; regeneration is not required.

## Citation

If you use these structures, please cite:

**Original structure:**
[6KIG PDB entry citation - to be added]

**This work:**
Mohamed A. A. Elrefaiy, Dvir Harris, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, and Gabriela S. Schlau-Cohen. "Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States." *Manuscript in preparation*, 2025.

## Related Documentation

- **MCCE data:** [../mcce/README.md](../mcce/README.md)
- **Hamiltonian:** [../hamiltonian/README.md](../hamiltonian/README.md)
- **Data dictionary:** [../DATA_DICTIONARY.md](../DATA_DICTIONARY.md)
- **Main README:** [../../README.md](../../README.md)

---

*Last updated: January 26, 2025*
