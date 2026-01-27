# pymembrane

A Python toolkit for **mesoscale modeling** of photosynthetic membrane structures and **exciton dynamics** simulations. Build equilibrated thylakoid membrane architectures using Monte Carlo methods and simulate energy transfer in pigment-protein complexes.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## What does pymembrane do?

`pymembrane` addresses the challenge of modeling photosynthetic systems across multiple scales:

1. **Membrane Structure Building**: Construct 2D thylakoid membrane architectures from atomic-resolution protein structures (PDB files)
2. **Structural Equilibration**: Use Monte Carlo simulations (translation, rotation, swapping) with simulated annealing to find energetically favorable protein arrangements
3. **Exciton Calculations**: Build electronic Hamiltonians for pigment aggregates and compute exciton states
4. **Spectroscopy**: Calculate linear optical spectra (absorption, fluorescence) using Kubo lineshape theory
5. **Energy Transfer Dynamics**: Simulate excitation energy transfer using kinetic Monte Carlo with generalized Förster rates

## Package Architecture

```
pymembrane/
│
├── structure/              # Structural modeling components
│   ├── atomic_protein.py   # Protein classes (Static/Dynamic/Pigment/Electrostatic)
│   ├── atomic_pigment.py   # Pigment molecules (Chlorophyll, Pheophytin)
│   ├── cg_particle.py      # Coarse-grained protein wrappers
│   ├── membrane.py         # Membrane container & Monte Carlo simulation
│   └── electrostatics.py   # Electronic coupling calculations
│
├── exciton/                # Spectroscopy & dynamics
│   ├── hamiltonian.py      # Electronic Hamiltonian management
│   ├── exciton_aggregate.py # Exciton domain handling
│   ├── linear_spectra.py   # Absorption/fluorescence calculations
│   ├── forster.py          # FRET rate calculations
│   ├── kmc_dynamics.py     # Kinetic Monte Carlo simulations
│   ├── rate_matrix.py      # Inter-domain rate matrix
│   └── lineshape.py        # Spectral lineshape functions
│
├── parameters/             # System-specific data
│   └── thylakoid/
│       ├── plants/         # LHCII, PSII parameters
│       ├── cyano/          # IsiA parameters
│       └── chlamy/         # Chlamydomonas systems
│
└── util/                   # Utilities
    ├── membrane_factory.py # Factory functions
    ├── load_hamiltonian.py # Hamiltonian I/O
    ├── coupling_data.py    # TrEsp coupling management
    └── physical_constants.py # Physical constants
```

## Class Hierarchy

### Structure Module

```
Protein Representations:
┌─────────────────────────────────────────────────┐
│ StaticProteinAtomic (Base)                      │  ← BioPython PDB wrapper
│  ├─ Read-only atomic coordinates                │
│  ├─ Geometric properties (hull, radius, COM)    │
│  └─ Atom/residue/chain access methods           │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ DynamicProteinAtomic                            │  ← Adds modifications
│  ├─ translate(), rotate_by_axis()               │
│  ├─ rotate_by_matrix()                          │
│  └─ _flatten(), _flip_up_down()                 │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ PigmentProteinAtomic                            │  ← Pigment management
│  ├─ prepare_pigments()                          │
│  ├─ load_hamiltonian()                          │
│  ├─ build_hamiltonian_domains()                 │
│  └─ H2_hamiltonian property                     │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ ElectrostaticProteinAtomic                      │  ← Charge calculations
│  ├─ set_atomic_charges()                        │
│  └─ calculate_cdc_site_shift()                  │
└─────────────────────────────────────────────────┘

Coarse-Grained Wrapper:
┌─────────────────────────────────────────────────┐
│ CGProtein                                       │  ← Efficient representation
│  ├─ References atomic protein                   │
│  ├─ CG position & orientation                   │
│  └─ Propagates transformations to atomic level  │
└─────────────────────────────────────────────────┘

Membrane System:
┌─────────────────────────────────────────────────┐
│ StaticAggregate (Base)                          │  ← Protein collection
│  └─ Visualization methods                       │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ PigmentAggregate                                │  ← Pigment aggregation
│  ├─ prepare_hamiltonian()                       │
│  └─ define_domain()                             │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ ProteinAggregate                                │  ← Force field
│  └─ calculate_total_potential()                 │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ Membrane                                        │  ← Boundary conditions
│  ├─ rho_of_r(), rdf(), nndf()                   │
│  └─ check_contained()                           │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ MonteCarloMembrane                              │  ← MC sampling
│  ├─ sample_translation/rotation/swapping        │
│  ├─ run_monte_carlo()                           │
│  └─ simulated_annealing()                       │
└─────────────────────────────────────────────────┘
```

### Exciton Module

```
Hamiltonian Management:
┌─────────────────────────────────────────────────┐
│ Hamiltonian                                     │  ← Electronic structure
│  ├─ update_coupling()                           │
│  ├─ subset_by_name/index()                      │
│  └─ build_energetic/coupling_domains()          │
└─────────────────────────────────────────────────┘

Exciton States:
┌─────────────────────────────────────────────────┐
│ ExcitonDomain                                   │  ← Single domain
│  ├─ Diagonalizes Hamiltonian                    │
│  ├─ thermal_pop(temp)                           │
│  ├─ abs_and_fl_exc_by_temp()                    │
│  └─ calculate_lifetime_by_temp()                │
└─────────────────────────────────────────────────┘

Spectroscopy:
┌─────────────────────────────────────────────────┐
│ LinearSpectra                                   │  ← Optical properties
│  ├─ calc_absorption(temp)                       │
│  ├─ calc_fluorescence(temp)                     │
│  └─ calc_linear_dichroism(temp, ref_dir)        │
└─────────────────────────────────────────────────┘
           ↓ extends
┌─────────────────────────────────────────────────┐
│ RateMatrix                                      │  ← Inter-domain transfer
│  └─ Sparse rate matrix construction             │
└─────────────────────────────────────────────────┘

Dynamics:
┌─────────────────────────────────────────────────┐
│ KMCMap                                          │  ← Rate matrix preprocessing
│  ├─ Connection indices                          │
│  └─ Transition probabilities                    │
└─────────────────────────────────────────────────┘
           ↓
┌─────────────────────────────────────────────────┐
│ KMCTrajectory                                   │  ← Single trajectory
│  └─ propagate(t_final, initial_state)           │
└─────────────────────────────────────────────────┘
           ↓
┌─────────────────────────────────────────────────┐
│ KMCEnsemble                                     │  ← Multiple trajectories
│  ├─ propagate_ensemble()                        │
│  └─ population_vs_time()                        │
└─────────────────────────────────────────────────┘
```

## Typical Workflow

```
1. STRUCTURE BUILDING
   │
   ├─► Load atomic structures
   │    prepare_protein('lhcii') → StaticProteinAtomic
   │
   ├─► Add pigments
   │    protein.prepare_pigments(resname='CLA', ...)
   │
   ├─► Create coarse-grained representations
   │    CGProtein(protein, location=[x, y, z])
   │
   └─► Build membrane
        Membrane(list_proteins, membrane_size, boundary='square')

2. STRUCTURAL EQUILIBRATION
   │
   ├─► Setup Monte Carlo
   │    MonteCarloMembrane(proteins, size=1200, boundary='square')
   │
   └─► Run simulated annealing
        membrane.simulated_annealing(
            start_temp=5000, end_temp=100,
            n_annealing_points=50,
            n_steps_per_block=1000
        )

3. HAMILTONIAN CONSTRUCTION
   │
   ├─► Load/build Hamiltonian for each protein
   │    protein.load_hamiltonian('ham.npz')
   │
   ├─► Partition into exciton domains
   │    protein.build_hamiltonian_domains(calc_method='tresp', r_max=20)
   │
   └─► Aggregate inter-protein couplings
        membrane.prepare_hamiltonian(r_max=100)

4. SPECTROSCOPY CALCULATIONS
   │
   ├─► Create spectra object
   │    LinearSpectra(pigmented_protein)
   │
   ├─► Calculate absorption
   │    freq, abs = spectra.calc_absorption(temp=77)
   │
   └─► Calculate fluorescence
        freq, fl = spectra.calc_fluorescence(temp=77)

5. DYNAMICS SIMULATION
   │
   ├─► Build rate matrix
   │    RateMatrix(pigmented_protein)
   │
   ├─► Calculate FRET rates
   │    generalized_forster_exciton(domain_a, domain_b, temp)
   │
   ├─► Setup KMC ensemble
   │    KMCEnsemble(rate_matrix)
   │
   └─► Run trajectories
        ensemble.propagate_ensemble(
            n_trajectories=100,
            t_final=1e-12,  # 1 picosecond
            initial_state=0
        )
```

## Installation

### Requirements

- Python >= 3.10
- Core dependencies:
  - `biopython` - PDB file handling
  - `numba` <= 0.57 - JIT compilation for performance
  - `numpy` < 1.24.0 - Numerical arrays
  - `scipy` - Scientific computing
  - `matplotlib` - Visualization
  - `interpolation` - Spectral interpolation
  - `shapely` - 2D geometric operations

### Install from source

```bash
git clone https://github.com/MesoscienceLab/pymembrane.git
cd pymembrane
pip install -e .
```

## Quick Start

### Example 1: Build and equilibrate a membrane

```python
import pymembrane as pm

# Load pre-configured protein structures
atomic_lhcii = pm.prepare_protein('lhcii')
atomic_psii = pm.prepare_protein('psii')

# Create coarse-grained representations for membrane simulation
cg_lhcii = pm.CGProtein(atomic_lhcii, name='lhcii', location=[100, 50, 0])
cg_psii = pm.CGProtein(atomic_psii, name='psii', location=[-100, 0, 0])

# Build membrane with multiple copies
proteins = [cg_lhcii] * 20 + [cg_psii] * 10

membrane = pm.MonteCarloMembrane(
    proteins,
    membrane_size=1200,  # 1200 Å x 1200 Å
    boundary='square'
)

# Equilibrate structure using simulated annealing
membrane.simulated_annealing(
    start_temp=5000,     # Initial temperature (kB*T in cm⁻¹)
    end_temp=100,        # Final temperature
    start_r=1000,        # Initial max translation (Å)
    end_r=10,            # Final max translation (Å)
    phi_max=40,          # Max rotation angle (degrees)
    n_annealing_points=50,
    n_steps_per_block=1000
)

# Analyze structure
print(f"Final energy: {membrane.calculate_total_potential()}")
print(f"Acceptance rates: {membrane.acceptance_rate('translate')}")
```

### Example 2: Calculate absorption spectrum

```python
from pymembrane.exciton import LinearSpectra

# Load protein with pigments
protein = pm.prepare_protein('lhcii')

# The protein comes pre-configured with pigments and Hamiltonian
# Create spectroscopy object
spectra = LinearSpectra(protein)

# Calculate absorption at 77K
freq, absorption = spectra.calc_absorption(temp=77)

# Plot
import matplotlib.pyplot as plt
plt.plot(freq, absorption)
plt.xlabel('Frequency (cm⁻¹)')
plt.ylabel('Absorption')
plt.title('LHCII Absorption at 77K')
plt.show()
```

### Example 3: Simulate energy transfer dynamics



## Included Scripts & Utilities

Helper scripts shipped inside `pymembrane/` that you can call directly or import in your own workflows:

- **Membrane setup**
  - `util/membrane_factory.py` — `calculate_membrane_stoichiometry()` estimates per-protein counts for a target membrane size/density and can instantly build a `Membrane` object; run with your own parameters inside a Python session or `python -m pymembrane.util.membrane_factory` to see the baked-in demo.
  - `util/parameterize.py` — `parameterize_default()` sweeps distances/angles between proteins to generate clash grids for force-field parameterization (supports multiprocessing).
  - `util/membrane_snapshot_loader.py` — `load_membrane_snapshot(path, index=None)` safely locates `.npz` snapshots with `list_pos`, `list_R2_orient`, `list_protein_type` keys.

- **Structure I/O**
  - `util/save_atomic_proteins_to_pdb.py` — `save_proteins_to_pdb([cg_proteins], "out.pdb")` exports CG proteins with preserved chain IDs/charges for Chimera/PyMOL.
  - `util/extract_pdb_chains.py` — `extract_chain_from_pdb("3ARC", "3ARC.pdb", ["A","D"])` splits specific chains to standalone PDBs.
  - `util/pdb_structure_alignment_and_mapping.py` — alignment/mapping helpers for matching atoms between PDB structures.
  - `util/writeExtendedPDB.py` — writer used by the exporters to keep charge annotations intact.

- **Hamiltonian & exciton analysis**
  - `util/load_hamiltonian.py` — load CSV Hamiltonians and remap pigment names with `load_hamiltonian()` / `match_pigment_names()`.
  - `util/calculate_exciton_distribution.py` — generate exciton energy/population distributions over disorder (`calculate_exciton_distribution(...)` for absorption/LD).
  - `util/helper_functions.py` — plotting presets plus wavelength/energy converters (`convert_to_nm`, `convert_nm_to_wavenumber`, etc.).
  - `util/plot_cdc_contribution.py` — build CDC residue-contribution bar charts and 2D spatial maps from interaction text files.

- **Specialized pipelines & datasets**
  - `util/mcce_neededFiles/run_mcce2Extended_conversion_terminal.py` — CLI to convert MCCE PDBs to extended PDBs: `python pymembrane/util/mcce_neededFiles/run_mcce2Extended_conversion_terminal.py -i input.pdb -o extended.pdb -d MEM [-c]`.
  - Additional MCCE helpers live in `util/mcce_neededFiles/` (pKa fitting, protonation profile visualization, renaming/conversion utilities).
  - Publication-focused scripts and bundled inputs live in `util/NeededFiles_PSII_paper/` and `util/NeededFiles_WSCP/` for reproducing spectra/CDC analyses.

## Key Classes Reference

### Core Entry Points (from `pymembrane`)

| Import | Purpose | Key Methods |
|--------|---------|-------------|
| `prepare_protein()` | Load pre-configured protein structures | Returns `StaticProteinAtomic` |
| `Membrane` | Membrane container with boundary conditions | `rdf()`, `nndf()`, `rho_of_r()` |
| `MonteCarloMembrane` | MC membrane equilibration | `run_monte_carlo()`, `simulated_annealing()` |
| `CGProtein` | Coarse-grained protein wrapper | `translate()`, `rotate_by_axis()` |
| `StaticProteinAtomic` | Read-only atomic protein | Property access: `mass`, `center_of_mass`, `hull3d` |
| `DynamicProteinAtomic` | Modifiable atomic protein | `translate()`, `rotate_by_matrix()` |

### Spectroscopy Classes (from `pymembrane.exciton`)

| Class | Purpose | Key Methods |
|-------|---------|-------------|
| `Hamiltonian` | Electronic Hamiltonian management | `update_coupling()`, `build_coupling_domains()` |
| `ExcitonDomain` | Single exciton domain | `thermal_pop()`, `abs_and_fl_exc_by_temp()` |
| `LinearSpectra` | Linear optical spectra | `calc_absorption()`, `calc_fluorescence()` |
| `RateMatrix` | Inter-domain transfer rates | Inherits `LinearSpectra` methods |
| `KMCEnsemble` | Kinetic Monte Carlo ensemble | `propagate_ensemble()`, `population_vs_time()` |

## Available Pre-Configured Systems

The `parameters/thylakoid/` directory contains parameters for several photosynthetic complexes:

- **Plant systems** (`plants/`):
  - LHCII (Light-Harvesting Complex II) - major antenna
  - PSII (Photosystem II) - reaction center

- **Cyanobacterial systems** (`cyano/`):
  - IsiA (iron-stress-induced A) - antenna ring

- **Chlamydomonas systems** (`chlamy/`):
  - Various algal complexes

Each system includes:
- PDB structure files
- Hamiltonian coupling matrices
- Lineshape parameters (Renger model)
- Disorder models (diagonal and off-diagonal)
- Force field parameters

## Feature Completeness

### ✓ Fully Implemented & Tested

- **Membrane structure modeling**: All protein classes, membrane boundaries, visualization
- **Monte Carlo simulations**: Translation, rotation, swapping, simulated annealing
- **Distribution functions**: Radial distribution (RDF), nearest neighbor (NNDF), density profiles
- **Pigment management**: Automatic detection, Hamiltonian construction (dipole/TrEsp coupling)
- **Exciton states**: Hamiltonian diagonalization, thermal populations
- **Linear spectroscopy**: Absorption, fluorescence, linear dichroism with Kubo lineshape
- **FRET calculations**: Generalized Förster rates with excitonic corrections
- **Kinetic Monte Carlo**: Single and ensemble trajectories with Gillespie algorithm
- **Visualization**: 2D/3D convex/concave hulls, atomic structures, dipole moments

### ⚠ Partially Implemented

- **Electrostatic calculations**: Charge loading works; charge-density coupling (CDC) present but incomplete
- **2D electronic spectroscopy**: Class structure exists; core functionality under development

### ✗ Not Included

- **Molecular dynamics**: No all-atom MD force fields (only clash detection)
- **Quantum dynamics**: No coherent exciton propagation (KMC only)
- **Advanced force fields**: Limited to geometric clash detection

## Use Cases

`pymembrane` is designed for computational researchers studying:

- **Photosynthetic light-harvesting**: Energy transfer pathways in photosynthetic proteins
- **Membrane organization**: Protein arrangement effects on energy transfer efficiency
- **Structure-function relationships**: How membrane architecture impacts optical properties
- **Spectroscopy prediction**: Calculating absorption/fluorescence for comparison with experiments
- **Mesoscale modeling**: Bridging atomic and cellular scales (nm to μm)

**Not suitable for:**
- All-atom molecular dynamics simulations
- Quantum coherent dynamics
- Non-photosynthetic membrane systems (no parameterization)

## Documentation

Detailed documentation is available in the source code docstrings. To explore:

```python
import pymembrane as pm
help(pm.Membrane)               # Membrane container
help(pm.CGProtein)              # Coarse-grained proteins
help(pm.MonteCarloMembrane)     # Monte Carlo simulation

from pymembrane.exciton import LinearSpectra, KMCEnsemble
help(LinearSpectra)             # Spectroscopy calculations
help(KMCEnsemble)               # Kinetic Monte Carlo
```

Key files to read:
- `structure/membrane.py` - Membrane systems and Monte Carlo
- `structure/atomic_protein.py` - Protein class hierarchy
- `exciton/linear_spectra.py` - Spectroscopy methods
- `exciton/kmc_dynamics.py` - Kinetic Monte Carlo implementation

## Contributing

This package is maintained by the MesoScience Lab at the University of Texas at Austin. Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

## Authors & Maintainers

- **Andrea Tagliabue** (andrea.tagliabue@edu.unige.it)
- **Mohamed El Refaiy** (moelrefaiy@gmail.com)
- **Bailey Raber** (baysmoot12@gmail.com)
- **Kajwal Patra** (kajwal.patra@austin.utexas.edu)
- **Doran I. G. B. Raccah** (doran.raccah@utexas.edu) - *Primary Maintainer*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2025 MesoScienceLab

## Citation

If you use `pymembrane` in your research, please cite:

```bibtex
@software{pymembrane2025,
  author = {Tagliabue, Andrea and Elrefaiy, Mohamed and Raber, Bailey and Patra, Kajwal and Raccah, Doran I. G. B.},
  title = {pymembrane: Mesoscale modeling of photosynthetic membranes},
  year = {2025},
  url = {https://github.com/MesoscienceLab/pymembrane}
}
```

## Repository

**GitHub**: https://github.com/MesoscienceLab/pymembrane

## Keywords

photosynthesis • thylakoid membrane • exciton dynamics • light-harvesting • LHCII • PSII • FRET • Monte Carlo • spectroscopy • mesoscale modeling • kinetic Monte Carlo
