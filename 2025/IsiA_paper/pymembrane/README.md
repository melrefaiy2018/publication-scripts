# PyMembrane Library

This directory contains a local copy of the `pymembrane` library, used for membrane protein simulations and excitonic Hamiltonian calculations in the IsiA project.

## Overview

The pymembrane library provides Python utilities for:
- Membrane protein structure handling
- Excitonic coupling calculations
- Chromophore positioning and analysis
- Spectral predictions from quantum chemistry

## Current Status

⚠️ **Note:** This directory currently contains placeholder files only. To complete the installation, the full pymembrane source code needs to be added.

## Installation

### Option 1: From Local Directory (Recommended for this project)

```bash
# From repository root
pip install -e pymembrane
```

### Option 2: From Remote Repository

```bash
# If pymembrane is hosted on GitHub (URL to be filled in):
pip install -e git+https://github.com/[username]/pymembrane.git#egg=pymembrane

# Or with specific commit hash:
pip install -e git+https://github.com/[username]/pymembrane.git@[commit-hash]#egg=pymembrane
```

## Expected Structure

The complete pymembrane package should have the following structure:

```
pymembrane/
├── setup.py                          # Package setup script
├── setup.cfg                         # Setup configuration
├── pyproject.toml                    # Build configuration (PEP 517)
├── README.md                         # This file
├── LICENSE                           # License file (should match MIT)
├── requirements.txt                  # Dependencies
├── pymembrane/                       # Main package directory
│   ├── __init__.py                   # Package initialization
│   ├── __version__.py                # Version information
│   ├── parameters/                   # Parameter definitions
│   │   ├── __init__.py
│   │   ├── hamiltonian.py            # Hamiltonian parameters
│   │   ├── chromophores.py           # Chromophore definitions
│   │   └── [other parameter modules]
│   ├── exciton/                      # Excitonic calculations
│   │   ├── __init__.py
│   │   ├── coupling.py               # Excitonic coupling
│   │   ├── hamiltonian.py            # Hamiltonian construction
│   │   └── [other exciton modules]
│   ├── structure/                    # Structure handling
│   │   ├── __init__.py
│   │   ├── pdb.py                    # PDB file processing
│   │   ├── chromophores.py           # Chromophore extraction
│   │   └── [other structure modules]
│   ├── spectra/                      # Spectral calculations
│   │   ├── __init__.py
│   │   ├── absorption.py             # Absorption spectrum
│   │   ├── fluorescence.py           # Fluorescence spectrum
│   │   └── [other spectral modules]
│   └── utils/                        # Utility functions
│       ├── __init__.py
│       ├── io.py                     # Input/output
│       ├── constants.py              # Physical constants
│       └── [other utility modules]
├── tests/                            # Unit tests
│   ├── __init__.py
│   ├── test_hamiltonian.py
│   ├── test_spectra.py
│   └── [other test files]
├── docs/                             # Documentation (optional)
│   └── [documentation files]
└── examples/                         # Usage examples (optional)
    └── [example scripts]
```

## Obtaining pymembrane

### From Developers

Contact the pymembrane developers for the latest source code:
- [Developer name/email - to be filled]
- [GitHub profile - to be filled]
- [Lab website - to be filled]

### From Repository

If pymembrane is public:
```bash
git clone https://github.com/[username]/pymembrane.git pymembrane_source
cp -r pymembrane_source/pymembrane pymembrane/
```

## Verification After Installation

Once pymembrane source is added, verify installation:

```bash
# Test imports
python -c "import pymembrane; print(pymembrane.__version__)"

# Check modules
python -c "from pymembrane import parameters, exciton, structure, spectra; print('✓ All modules loaded')"

# Run tests (if available)
cd pymembrane
python -m pytest tests/
```

## Integration with IsiA Analysis

The pymembrane library is used by:
- `scripts/hamiltonian/cdc_IsiA_monomer_average_pH7.py` - Hamiltonian calculations
- `scripts/fluorescence_decay/model_*/run_model_*.py` - Fluorescence simulations

These scripts import pymembrane modules for:
- Loading and processing protein structures
- Calculating excitonic couplings
- Managing chromophore data
- Generating spectral predictions

## Requirements

pymembrane requires:
- numpy >= 1.20.0
- scipy >= 1.7.0
- pandas >= 1.3.0
- (Other requirements to be specified in pymembrane/requirements.txt)

## License

The pymembrane library should use a license compatible with MIT:
- MIT License (preferred for compatibility)
- BSD License (acceptable)
- Apache 2.0 (acceptable)
- GPL (not compatible - would restrict this project)

## Attribution

When using pymembrane, ensure proper citation:

```bibtex
@software{pymembrane,
  author = {[Authors]},
  title = {PyMembrane: Python Library for Membrane Protein Simulations},
  year = {[Year]},
  url = {[GitHub URL]}
}
```

## Troubleshooting

### ModuleNotFoundError: No module named 'pymembrane'

This error means pymembrane source needs to be added to this directory.

**Solution:**
1. Obtain pymembrane source code
2. Place source files in this directory maintaining structure
3. Run: `pip install -e .` from repository root
4. Verify with: `python -c "import pymembrane"`

### ImportError: cannot import name 'X' from 'pymembrane'

The required module is missing or has different structure.

**Solution:**
1. Check pymembrane version matches requirements
2. Verify all required modules are present
3. Review pymembrane package structure

## TODO - Zenodo Release

Before final release, ensure:
- [ ] pymembrane source code is added to this directory
- [ ] All dependencies are specified in requirements.txt
- [ ] License is compatible (MIT or BSD preferred)
- [ ] Attribution and citation are correct
- [ ] Installation can be verified successfully

## Related Documentation

- **Main README:** [../README.md](../README.md)
- **Installation guide:** [../README.md#installation](../README.md#installation)
- **Requirements:** [../requirements.txt](../requirements.txt)
- **Quick start:** [../QUICKSTART.md](../QUICKSTART.md)

---

*Last updated: January 26, 2025*

*Status: Awaiting pymembrane source code*
