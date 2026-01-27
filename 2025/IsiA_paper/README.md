# IsiA Protein Spectroscopy Analysis

## Description

This repository contains data and code accompanying the manuscript:

**Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States**
*Mohamed A. A. Elrefaiy, Dvir Harris,, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, Gabriela S. Schlau-Cohen*


This work presents a computational analysis of the IsiA protein complex, including:
* Excitonic Hamiltonian calculations
* Absorption and fluorescence spectra predictions
* Time-resolved fluorescence decay simulations

## Repository Structure

The repository is organized as follows:

```
IsiA_data_and_code/
├── README.md                  # This file
├── LICENSE                    # License information
├── requirements.txt           # Python dependencies
├── NeededData/                # Input data and generated outputs
│   ├── Experimental_Data/     # Raw & processed experimental spectra
│   ├── fluorescence_decay/    # Ensemble-averaged fluorescence decay data
│   ├── hamiltonian/           # Hamiltonian & spectra outputs
│   ├── mcce/                  # MCCE input and output files
│   └── structure/             # Protein structure files
├── pymembrane/                # Local copy of the pymembrane library
└── scripts/                   # Executable analysis scripts
    ├── hamiltonian/
    └── fluorescence_decay/
```

### Key Directories and Files

*   **`NeededData/`**: This directory contains all the necessary input data for the calculations and stores the results.
    *   `Experimental_Data/`: Contains the experimental absorption and fluorescence spectra for comparison.
    *   `fluorescence_decay/`: Stores the output of the fluorescence decay simulations.
    *   `hamiltonian/`: Contains the calculated Hamiltonian, site energies, and spectra.
    *   `mcce/`: Contains input and output files for MCCE (Multi-Conformation Continuum Electrostatics) calculations.
    *   `structure/`: Contains the PDB structures used in the simulations.
*   **`scripts/`**: This directory contains the Python scripts used to perform the calculations and analysis.
    *   `hamiltonian/`: Scripts for calculating the excitonic Hamiltonian and spectra.
    *   `fluorescence_decay/`: Scripts for simulating the fluorescence decay dynamics, organized into different models.
*   **`pymembrane/`**: A local copy of the `pymembrane` library, used for membrane protein simulations.
*   **`requirements.txt`**: A list of the Python packages required to run the scripts in this repository.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/IsiA_data_and_code.git
    cd IsiA_data_and_code
    ```

2.  **Create a Python virtual environment (recommended):**
    ```bash
    python3 -m venv env
    source env/bin/activate
    ```

3.  **Install the required Python packages:**
    ```bash
    pip install -r requirements.txt
    ```

4.  **Install `pymembrane`:**
    The `pymembrane` library is included as a local copy. Install it in editable mode:
    ```bash
    pip install -e pymembrane
    ```

## How to Reproduce the Results

### 1. Calculate the Excitonic Hamiltonian and Spectra

This step calculates the excitonic Hamiltonian, site energies, and spectra for the IsiA monomer.

1.  Navigate to the `scripts/hamiltonian` script directory:
    ```bash
    cd scripts/hamiltonian
    ```

2.  Run the main script:
    ```bash
    python cdc_IsiA_monomer_average_pH7.py
    ```

    This script will:
    *   Read the protein structure and MCCE data.
    *   Calculate the site energies and excitonic couplings.
    *   Generate the Hamiltonian matrix.
    *   Calculate the absorption and fluorescence spectra.

    **Outputs:**
    The results will be saved in the `NeededData/hamiltonian/` directory, including:
    *   `hamiltonian_matrix.csv`: The calculated excitonic Hamiltonian.
    *   `absorption_spectrum.csv`: The calculated absorption spectrum.
    *   `fluorescence_spectrum.csv`: The calculated fluorescence spectrum.
    *   `SiteEnergy_Data/`: Detailed site energy information.
    *   `Spectra_Data/`: Raw spectral data.

### 2. Simulate Fluorescence Decay

This step simulates the fluorescence decay dynamics using different models. The scripts are organized into subdirectories within `scripts/fluorescence_decay/`. Each `model_*` directory represents a different set of parameters or assumptions.

1.  Navigate to the desired model directory, for example, `model_1`:
    ```bash
    cd scripts/fluorescence_decay/model_1
    ```

2.  Run the simulation script:
    ```bash
    python run_model_1.py
    ```
    This will run the simulation for the parameters defined in `unified_parameters.py`.

    For parallel execution on a SLURM cluster, you can use the provided SLURM script:
    ```bash
    sbatch run_model1_slurm_parallel.py
    ```

    **Outputs:**
    The simulation results, including time-resolved and wavelength-resolved data, will be saved in the `simulation_data` directory at the root of the project.

### 3. Find the Best-Fit Model

After running the simulations for different models, you can use the `find_best_paramter_model` scripts to compare the simulation results with experimental data and identify the best-fitting model.

1.  Navigate to the `find_best_paramter_model` directory within a model's folder, for example:
    ```bash
    cd scripts/fluorescence_decay/model_1/find_best_paramter_model
    ```

2.  Run the analysis script:
    ```bash
    python run_best_model_pick.py
    ```

    This script will:
    *   Load the simulation and experimental data.
    *   Calculate the error between the simulated and experimental decay traces.
    *   Generate plots comparing the different models.

    **Outputs:**
    The analysis results, including plots and CSV files with error metrics, will be saved in the `outputs` subdirectory.

## Data Availability

All data used in this study, including the final results, are available in the `NeededData` directory. The raw simulation data is available in the `simulation_data` directory.

## Citation

### Manuscript (In Preparation)

If you use this code or data in your research, please cite the following manuscript:

**Quenching of the Photosynthetic Antenna IsiA is Facilitated by its Red-Emitting States**

Dvir Harris, Mohamed A. A. Elrefaiy, Hila Toporik, Christopher J. Gisriel, Yuval Mazor, Doran I. G. B. Raccah, Gabriela S. Schlau-Cohen

*Manuscript in preparation. 2025.*

### GitHub Repository

This repository is also citable. You can cite it using the DOI provided by Zenodo:

```
Harris, D., Elrefaiy, M. A. A., Toporik, H., Gisriel, C. J., Mazor, Y., Raccah, D. I. G. B., & Schlau-Cohen, G. (2025). 
IsiA protein spectroscopy analysis: Code and data [Computer software]. GitHub. 
https://github.com/melrefaiy2018/Papers
```

*Note: Once the manuscript is published, the citation information will be updated with the journal details and DOI.*

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.