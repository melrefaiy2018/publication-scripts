# Standard library imports
import argparse
import logging
import os
import sys
import warnings
from itertools import product
import yaml  # Added for YAML saving

from pathlib import Path
# Get the current working directory
current_dir = Path(__file__).resolve().parent.parent / 'NeededFiles'
print(f"Current directory: {current_dir}")

# Third-party imports
import matplotlib.pyplot as plt
import numpy as np
from Bio import BiopythonWarning
from scipy import sparse
from tqdm import tqdm

# Local pymembrane imports
from pymembrane.exciton.nonlinear_spectra import NonlinearSpectra
from pymembrane.util.physical_constants import c, hbar
from construct_monomer_with_CT import construct_monomer

# Import unified parameters
from unified_parameters import get_parameter_set, print_configuration, list_configurations

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)
# use default matplotlib
plt.style.use('default')  # Use the default style sheet (white background)

def parse_args():
    """Parse command line arguments for the script."""
    parser = argparse.ArgumentParser(description="Run Fl decay model (Model 12) with specified parameters.")
    parser.add_argument('--seed', type=int, default=1, help="Random seed for the ensemble.")
    parser.add_argument('--param_index', type=int, default=None,
                        help="Index of parameter combination to run (for parallel execution).")
    parser.add_argument('--total_params', type=int, default=None,
                        help="Total number of parameter combinations (for validation).")
    parser.add_argument('--config', type=str, default='model_1',
                        help="Parameter configuration to use (see unified_parameters.py)")
    parser.add_argument('--list_configs', action='store_true',
                        help="List available parameter configurations and exit")
    parser.add_argument('--show_config', type=str, default=None,
                        help="Show details of a specific configuration and exit")
    parser.add_argument('--list_params', action='store_true',
                        help="Just list all parameter combinations and exit.")

    return parser.parse_args()

def setup_logging(run_dir, seed):
    """Set up logging for the script.

    Parameters:
    -----------
    run_dir : str
        Directory where logs will be stored
    seed : int
        Random seed for this run

    Returns:
    --------
    logger : logging.Logger
        Configured logger object
    """
    # Create log directory within the run directory
    log_dir = os.path.join(run_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    # Get the script name dynamically and use it to create the log file name
    # Use a fixed name if __file__ is not defined (e.g., in interactive environments)
    try:
        script_name = os.path.splitext(os.path.basename(__file__))[0]
    except NameError:
        script_name = "run_model_1"  # Fallback name
    log_file = os.path.join(log_dir, f'{script_name}_seed_{seed}.log')

    # Set up logging
    # Check if handlers already exist to avoid duplicates if called multiple times
    logger = logging.getLogger()
    # Clear existing handlers to avoid duplicate logs if function is called multiple times
    # in the same process (e.g., during sequential runs in main)
    if logger.hasHandlers():
        logger.handlers.clear()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    logger.setLevel(logging.INFO)  # Ensure level is set

    return logger

def construct_ratematrix(protein_atomic, seed, run_dir, time_fl, time_isc):
    """
    Construct the rate matrix for the given protein_atomic object.

    Parameters:
    -----------
    protein_atomic : object
        An instance of the ElectrostaticProteinAtomic class.
    seed : int
        Random seed for disorder
    run_dir : str
        Directory to save the rate matrix
    time_fl : float
        Fluorescence time constant
    time_isc : float
        Intersystem crossing time constant

    Returns:
    --------
    nonlinear_isia : NonlinearSpectra
        The constructed NonlinearSpectra object
    """
    logger = logging.getLogger()

    protein_atomic.H2_hamiltonian.add_disorder(seed)
    nonlinear_isia = NonlinearSpectra(protein_atomic, temp=300, k_ij=None)

    nonlinear_isia.add_fluorescence(1 / (time_fl * 10 ** 6))
    nonlinear_isia.add_isc(1 / (time_isc * 10 ** 6))
    nonlinear_isia._correct_diagonal()
    logger.info(f'Saving the rate matrix with seed: {seed}')

    # Create the rate matrix directory inside the run directory
    rate_matrix_dir = os.path.join(run_dir, 'Rate_matrix')
    os.makedirs(rate_matrix_dir, exist_ok=True)

    sparse.save_npz(os.path.join(rate_matrix_dir, f'K_ij_isia_monomer_{seed}.npz'), nonlinear_isia.k_ij)
    return nonlinear_isia

def save_compressed_data(file_path, **kwargs):
    """
    Save data with compression to reduce file size.

    Parameters:
    -----------
    file_path : str
        Path where the compressed data will be saved
    **kwargs : dict
        Key-value pairs to be saved in the compressed file
    """
    logger = logging.getLogger()
    logger.info(f'Saving compressed data to {file_path}')
    np.savez_compressed(file_path, **kwargs)

def calculate_fluorescence(nonlinear_isia, prob_i0, run_dir, seed, suffix=""):
    """
    Calculate time-resolved fluorescence and save the results.

    Parameters:
    -----------
    nonlinear_isia : NonlinearSpectra
        NonlinearSpectra object
    prob_i0 : ndarray
        Initial probability distribution
    run_dir : str
        Directory to save the results
    seed : int
        Random seed
    suffix : str
        Suffix to add to the output file name

    Returns:
    --------
    tuple
        list_t, lambda_axis, list_fl_t
    """
    logger = logging.getLogger()
    logger.info(f'Calculating the time resolved fluorescence{suffix}')  # Added suffix to log

    N_t = 800000
    dt = 10
    N_save = 1000
    list_t, w_axis, list_fl_t = nonlinear_isia.calculate_time_resolved_fluorescence(
        P1_0=prob_i0, N_t=N_t, dt=dt, N_save=N_save
    )
    
    # Ensure w_axis is array and handle potential division by zero
    w_axis = np.asarray(w_axis)
    lambda_axis = np.full_like(w_axis, np.nan, dtype=float)  # Initialize with NaN
    valid_mask = w_axis != 0
    if np.any(valid_mask):
        lambda_axis[valid_mask] = c / ((w_axis[valid_mask]) / (2 * np.pi * hbar))
    else:
        logger.warning(f"w_axis contains only zeros or is empty for suffix '{suffix}'. lambda_axis will contain NaNs.")

    # Save the data
    logger.info(f'Saving fluorescence decay data{suffix} with compression')  # Added suffix to log
    data_dir = os.path.join(run_dir, 'Data')
    os.makedirs(data_dir, exist_ok=True)

    # Save data with compression
    data_file = os.path.join(data_dir, f'FL_decay_monomer{suffix}_{seed}.npz')
    # Better compression by downsampling and converting to float32
    # Reduce time points (take every 10th point)
    list_t_reduced = list_t[::10]
    list_fl_t_reduced = list_fl_t[::10]

    # Convert to float32 for further compression
    list_t_reduced = np.array(list_t_reduced, dtype=np.float32)
    lambda_axis_reduced = np.array(lambda_axis, dtype=np.float32)  # Convert potentially NaN-filled axis
    list_fl_t_reduced = np.array(list_fl_t_reduced, dtype=np.float32)

    save_compressed_data(data_file,
                         list_t=list_t_reduced,
                         lambda_axis=lambda_axis_reduced,
                         list_fl_t=list_fl_t_reduced)

    # Check file existence before getting size
    if os.path.exists(data_file):
        try:
            file_size_mb = os.path.getsize(data_file) / (1024*1024)
            logger.info(f'Data{suffix} saved with reduced size: {file_size_mb:.2f} MB')  # Added suffix to log
        except OSError as e:
            logger.error(f"Could not get size of saved file {data_file}: {e}")
    else:
        logger.error(f"Failed to save data file: {data_file}")

    return list_t, lambda_axis, list_fl_t

def generate_parameter_combinations(config_name='model_1'):
    """Generate all combinations of parameters to test."""
    PARAMETER_SETS = get_parameter_set(config_name)
    param_keys = list(PARAMETER_SETS.keys())
    param_values = [PARAMETER_SETS[key] for key in param_keys]

    combinations = list(product(*param_values))
    return combinations, param_keys

def save_parameters_to_yaml(run_dir, params):
    """
    Saves the given parameters dictionary to a YAML file in the specified directory.

    Parameters:
    -----------
    run_dir : str
        The main directory for the simulation run.
    params : dict
        Dictionary containing the parameters to save.
    """
    logger = logging.getLogger()  # Get logger instance
    param_dir = os.path.join(run_dir, 'param')
    os.makedirs(param_dir, exist_ok=True)
    param_file = os.path.join(param_dir, 'run_parameters.yaml')

    try:
        with open(param_file, 'w') as f:
            # Use sort_keys=False to maintain order if desired (requires PyYAML 5.1+)
            # Use default_flow_style=False for block style (more readable)
            yaml.dump(params, f, default_flow_style=False, sort_keys=False)
        logger.info(f"Run parameters saved to {param_file}")
    except Exception as e:
        logger.error(f"Failed to save parameters to {param_file}: {e}")

def run_simulation(time_fl, time_isc, coupling_value, seed=1, config_name='model_1'):
    """
    Run a complete simulation with the given parameters.

    Parameters:
    -----------
    time_fl : float
        Fluorescence time constant (in us)
    time_isc : float
        Intersystem crossing time constant (in us)
    coupling_value : int
        Coupling value for CT state (in cm-1)
    seed : int
        Random seed for the ensemble
    config_name : str
        Configuration name used for this run
    """
    # Create a main run directory that includes all parameters in the name
    run_dir = f'Simulation/Run_FL{time_fl}_ISC{time_isc}_coupling_{coupling_value}_{config_name}'
    os.makedirs(run_dir, exist_ok=True)

    # Setup logging for this specific run
    logger = setup_logging(run_dir, seed)
    logger.info(f"Starting analysis with seed = {seed}")
    logger.info(f"Run directory: {run_dir}")
    logger.info(f"Configuration: {config_name}")
    # Updated log message to include units
    logger.info(f"Parameters: FL={time_fl} us, ISC={time_isc} us, Coupling={coupling_value} cm-1")

    # --- Save Parameters to YAML ---
    pdb_file = 'extended_most_occ_pH7_CT.pdb'  # Define pdb_file used (no CT version needed)
    
    # Build only the *specific* parameter set used
    used_parameters = {
        'TIME_FL': time_fl,
        'TIME_ISC': time_isc,
        'COUPLING_value': coupling_value
    }

    # Define the parameters dictionary to be saved, including all relevant info
    run_params_to_save = {
        'configuration_used': config_name,
        'run_specific_parameters': {
            'seed': seed,
            'time_fl_us': time_fl,
            'time_isc_us': time_isc,
            'coupling_value_cm-1': coupling_value,
        },
        'fixed_parameters': {
             'pdb_file': pdb_file,
             'e0_a_cm-1': 14900,  # Parameter for construct_monomer*
             'e0_b_cm-1': 0,      # Parameter for construct_monomer*
             'N_ens': 100,        # Parameter for construct_monomer*
             'temp_K': 300,       # Parameter for construct_monomer*
             'dielectric_tresp': 1,  # Parameter for construct_monomer*
             'dielectric_dipole': 1.5,  # Parameter for construct_monomer*
             'dielectric_cdc': 2,       # Parameter for construct_monomer*
        },
        'calculation_settings': {  # Parameters from calculate_fluorescence
            'N_t': 800000,
            'dt_fs': 10,  # Assuming dt is in femtoseconds
            'N_save': 1000
        },
        'used_parameter_set': used_parameters
    }
    save_parameters_to_yaml(run_dir, run_params_to_save)  # Call the saving function
    # --- End Save Parameters ---

    # Copy the needed files to the run directory
    # Check if files exist before copying
    analysis_file_path = current_dir / 'run_analysis.py'
    if analysis_file_path.exists():
        os.system(f'cp "{analysis_file_path}" {run_dir}')
    else:
        logger.warning(f"Could not find '{analysis_file_path}' to copy.")

    # Copy the script itself and unified parameters
    try:
        script_name = os.path.basename(__file__)
        os.system(f'cp {script_name} {run_dir}')
        os.system(f'cp unified_parameters.py {run_dir}')
    except NameError:
        logger.warning("Could not copy script file (__file__ not defined).")

    # Prepare protein structures
    logger.info("Preparing protein structures")
    # Using fixed parameters as in the original call
    isia_monomer = construct_monomer(
        pdb_file, dir_path=run_dir, e0_a=14900, e0_b=0, N_ens=100,
        temp=300, dielectric_tresp=1, dielectric_dipole=1, dielectric_cdc=2, seed=seed
    )

    # Build the rate matrices
    logger.info(f'Calculating the rate matrices with seed: {seed}')
    nonlinear_isia = construct_ratematrix(isia_monomer, seed, run_dir, time_fl, time_isc)

    # Calculate initial probability distributions
    logger.info('Calculating initial probability distributions')
    # Add checks for list_domain existence and content before calculating probabilities
    prob_i0 = None
    if hasattr(nonlinear_isia, 'list_domain') and nonlinear_isia.list_domain:
        try:
            sums = [np.sum(['A_CLA' in name for name in domain.list_names])
                    for domain in nonlinear_isia.list_domain]
            prob_i0 = np.array(sums, dtype=float)
            total_prob = np.sum(prob_i0)
            if total_prob > 0:
                prob_i0 = prob_i0 / total_prob
            else:
                logger.warning("Sum of initial probabilities is zero. Check 'A_CLA' names in list_domain.")
                prob_i0 = None  # Indicate invalid probability
        except (AttributeError, TypeError) as e:
            logger.error(f"Error accessing list_names for prob_i0 calculation: {e}")
            prob_i0 = None

    # Calculate fluorescence only if probabilities are valid
    if prob_i0 is not None:
        calculate_fluorescence(nonlinear_isia, prob_i0, run_dir, seed)
    else:
        logger.error("Skipping fluorescence calculation due to invalid initial probability.")

    logger.info(f"Analysis complete. Results saved in {run_dir}")

def main():
    """Main function to run the script, supporting both single and parallel execution."""
    args = parse_args()
    
    # Handle configuration queries
    if args.list_configs:
        print("Available parameter configurations:")
        for config in list_configurations():
            print(f"  - {config}")
        sys.exit(0)
    
    if args.show_config:
        print_configuration(args.show_config)
        sys.exit(0)
    
    seed = args.seed
    config_name = args.config

    # Generate all parameter combinations for the specified configuration
    all_combinations, param_keys = generate_parameter_combinations(config_name)
    total_combinations = len(all_combinations)
    
    # If just listing parameters, print them and exit
    if args.list_params:
        print(f"Total parameter combinations for config '{config_name}': {total_combinations}")
        for i, combo in enumerate(all_combinations):
            param_dict = dict(zip(param_keys, combo))
            print(f"Combination {i}: {param_dict}")
        return

    # Validate total_params if provided
    if args.total_params is not None and args.total_params != total_combinations:
        print(f"ERROR: Provided total_params ({args.total_params}) doesn't match actual combinations ({total_combinations})")
        sys.exit(1)

    # Determine which parameter combination to run
    if args.param_index is not None:
        # --- Parallel Execution ---
        if args.param_index < 0 or args.param_index >= total_combinations:
            print(f"ERROR: param_index ({args.param_index}) out of range (0-{total_combinations-1})")
            sys.exit(1)

        # Run only the specified combination
        combo = all_combinations[args.param_index]
        param_dict = dict(zip(param_keys, combo))
        print(f"Running combination {args.param_index+1}/{total_combinations}: {param_dict} with seed {seed}")  # Added seed info

        # Extract only the parameters needed by run_simulation signature
        run_simulation(
            time_fl=param_dict["TIME_FL"],
            time_isc=param_dict["TIME_ISC"],
            coupling_value=param_dict["COUPLING_value"],
            seed=seed,  # Pass the seed from args
            config_name=config_name
        )
    else:
        # --- Sequential Execution ---
        print(f"Running all {total_combinations} parameter combinations for config '{config_name}' sequentially with seed {seed}")  # Added seed info
        # Use tqdm for progress bar
        for i, combo in enumerate(tqdm(all_combinations, desc="Running Simulations")):
            param_dict = dict(zip(param_keys, combo))
            # Add separator for clarity in logs/output when running sequentially
            print(f"\n--- Running combination {i+1}/{total_combinations}: {param_dict} with seed {seed} ---")

            run_simulation(
                time_fl=param_dict["TIME_FL"],
                time_isc=param_dict["TIME_ISC"],
                coupling_value=param_dict["COUPLING_value"],
                seed=seed,  # Pass the seed from args
                config_name=config_name
            )

if __name__ == "__main__":
    main()