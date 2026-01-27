# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Add ORCID identifiers to author metadata
- Include additional model variants for extended analysis
- Add interactive visualization notebooks

---

## [1.0.0] - 2025-01-26

### Added - Zenodo Release Preparation
- Created comprehensive CITATION.cff file for proper citation metadata
- Created .zenodo.json for Zenodo repository upload configuration
- Created QUICKSTART.md with 5-minute setup instructions
- Created DATA_DICTIONARY.md with comprehensive parameter definitions
- Created AUTHORS.md listing all contributors and their roles
- Created CHANGELOG.md (this file) for version tracking
- Created .gitignore for better repository hygiene
- Added archive/ subdirectory in Experimental_Data with version history
- Added detailed README files for all data subdirectories

### Changed
- Renamed all `find_best_paramter_model/` directories to `find_best_parameter_model/` (fixed typo)
- Updated README references from "paramter" to "parameter"
- Fixed file timestamps for future-dated files (pK.out, PDB structures) to December 2024
- Archived older experimental data files to Experimental_Data/archive/
- Standardized experimental data: `Fl_IsiA_monomer_300k_Gabriela_2025.npy` is now the canonical version
- Enhanced main README.md with expanded installation instructions
- Improved documentation for all NeededData subdirectories

### Removed
- Deleted all 42 .DS_Store files (macOS system metadata)
- Deleted all 12 .pyc files (Python bytecode)
- Removed calculation_log_20251110_112447.log (outdated log file)

### Fixed
- Corrected directory naming convention (parameter vs paramter)
- Fixed inconsistent file timestamps
- Resolved duplicate data file confusion with clear versioning

---

## [0.9.0] - 2024-12-10

### Added - Initial Repository Structure
- Initial commit of all analysis scripts and data
- Hamiltonian calculation scripts (cdc_IsiA_monomer_average_pH7.py)
- Eight fluorescence decay models (model_1 through model_123)
- MCCE structure and output files
- Experimental fluorescence data
- Basic README documentation
- MIT License

### Changed
- Migrated from local development to GitHub repository
- Organized scripts into hamiltonian/ and fluorescence_decay/ directories
- Structured data into NeededData/ with subdirectories

---

## Data Version History

### Experimental Data

#### v2.0 (February 14, 2025)
- **File:** `Fl_IsiA_monomer_300k_Gabriela_2025.npy`
- **Changes:**
  - Improved signal-to-noise ratio
  - Enhanced baseline correction
  - Better calibration procedures
  - Switched to binary NPY format for numerical precision
- **Status:** Current canonical version

#### v1.0 (May 5, 2024)
- **File:** `Fl_IsiA_monomer_300k_Gabriela.npy`
- **Changes:** Initial experimental measurements
- **Status:** Archived

### Hamiltonian Data

#### December 2024
- Initial Hamiltonian calculations using MCCE pH 7 protonation states
- 17×17 excitonic Hamiltonian for IsiA monomer
- Site energies and coupling matrices

---

## Model Development History

### Model Evolution

**Model 1** (Original): Single red-state quenching model
**Model 2**: Two red-state model with independent quenching
**Model 3**: Three red-state model with energy-dependent rates
**Model 4**: Four red-state model with complete state space
**Model 12**: Combined analysis of models 1 and 2
**Model 13**: Combined analysis of models 1 and 3
**Model 14**: Combined analysis of models 1 and 4
**Model 123**: Comprehensive three-model comparison

### Parameter Optimization
- Implemented parameter sweep functionality for all models
- Added best-fit model selection based on experimental comparison
- Created analysis pipelines for model comparison

---

## Code Changes by Category

### Phase 1: Initial Development (May-November 2024)
- Core Hamiltonian calculation engine
- Fluorescence decay simulation framework
- Basic SLURM parallelization scripts

### Phase 2: Model Expansion (November-December 2024)
- Extended from 3 models to 8 model variants
- Improved ensemble averaging statistics
- Enhanced output data structures

### Phase 3: Documentation and Release Prep (January 2025)
- Comprehensive documentation overhaul
- Zenodo release preparation
- Code cleanup and standardization

---

## Manuscript Milestones

- **2024-05:** Initial data collection
- **2024-11:** Computational analysis completion
- **2024-12:** Initial manuscript draft
- **2025-01:** Zenodo release preparation
- **2025-XX:** Manuscript submission (planned)
- **2025-XX:** Publication (planned)

---

## Breaking Changes

### Version 1.0.0
- **Directory renaming:** Scripts referencing `find_best_paramter_model/` will need to update paths to `find_best_parameter_model/`
- **Data file changes:** Older experimental data moved to archive/; scripts should use `Fl_IsiA_monomer_300k_Gabriela_2025.npy`

---

## Technical Debt and Future Work

### Current Limitations
- Some models have long execution times (hours to days)
- Limited error handling in parallel execution scripts
- Memory usage could be optimized for large ensembles

### Planned Improvements
- Add progress bars for long-running simulations
- Implement checkpoint/resume functionality
- Add automated unit tests for core functions
- Create Docker container for reproducibility
- Add GPU acceleration for large ensemble calculations

---

## Contributors

See [AUTHORS.md](AUTHORS.md) for a complete list of contributors.

---

## How to Report Issues

Found a bug or have a suggestion? Please:
1. Check existing GitHub Issues
2. Create a new issue with:
   - Clear description of the problem
   - Steps to reproduce
   - Expected vs actual behavior
   - System information (OS, Python version)

---

*This changelog will be updated with each significant release and before major milestones.*
