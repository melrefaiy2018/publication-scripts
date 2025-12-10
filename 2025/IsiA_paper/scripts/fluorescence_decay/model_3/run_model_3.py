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
os.system(f'cp -r "{current_dir}/construct_monomer_with_CT.py" .')
os.system(f'cp -r "{current_dir}/extended_most_occ_pH7_CT.pdb" .')


import matplotlib.pyplot as plt
import numpy as np
from Bio import BiopythonWarning
from scipy import sparse
from tqdm import tqdm

from pymembrane.exciton.nonlinear_spectra import NonlinearSpectra
from pymembrane.util.physical_constants import c, hbar
from construct_monomer_with_CT import construct_monomer_with_CT, construct_monomer
from unified_parameters import get_parameter_set, PARAMETER_SETS, CT_ENERGY, print_configuration, list_configurations

# Suppress Biopython warnings
warnings.simplefilter('ignore', BiopythonWarning)
plt.style.use('default')


def parse_args():
    """Parse command line arguments for the script."""
    parser = argparse.ArgumentParser(description="Run Fl decay model with specified parameters.")
    parser.add_argument('--config', type=str, default='model_3',
                        help="Configuration name to use (default: model_3)")
    parser.add_argument('--seed', type=int, default=1, help="Random seed for the ensemble.")
    parser.add_argument('--param_index', type=int, default=None, 
                        help="Index of parameter combination to run (for parallel execution).")
    parser.add_argument('--total_params', type=int, default=None,
                        help="Total number of parameter combinations (for validation).")
    parser.add_argument('--list_params', action='store_true',
                        help="Just list all parameter combinations and exit.")
    
    return parser.parse_args()

def setup_logging(run_dir, seed):
    """Set up logging for the script."""
    log_dir = os.path.join(run_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)

    script_name = os.path.splitext(os.path.basename(__file__))[0]
    log_file = os.path.join(log_dir, f'{script_name}_seed_{seed}.log')

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )

    return logging.getLogger()

# --- Added function for saving parameters ---
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
    logger = logging.getLogger() # Get logger instance
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

def run_simulation(time_fl, time_isc, time_to_rp1, coupling_value, seed=1):
    """
    Run a complete simulation with the given parameters.
    """
    # Create directory for this parameter combination (without seed in name)
    run_dir = f'Simulation/Run_FL{time_fl}_ISC{time_isc}_RP{time_to_rp1}_coupling_{coupling_value}_CT_energy_{CT_ENERGY}'
    os.makedirs(run_dir, exist_ok=True)
    
    # Copy needed files to the run directory (only once)
    if not os.path.exists(os.path.join(run_dir, 'run_analysis.py')):
        os.system(f'cp -r ../NeededFiles/run_analysis.py {run_dir}')
        os.system(f'cp -r run_model3.py {run_dir}')

    # Setup logging with seed-specific log file
    logger = setup_logging(run_dir, seed)
    logger.info(f"Starting analysis with seed = {seed}")
    logger.info(f"Run directory: {run_dir}")
    logger.info(f"Parameters: FL={time_fl}, ISC={time_isc}, RP1={time_to_rp1}, Coupling={coupling_value}")


    # --- Save Parameters to YAML ---
    pdb_file = 'extended_most_occ_pH7_CT.pdb' # Define pdb_file used
    # Build only the *specific* parameter set used
    used_parameters = {
        'TIME_FL': time_fl,
        'TIME_ISC': time_isc,
        'TIME_to_CT': 0,
        'TIME_to_RP1': time_to_rp1,
        'COUPLING_value': coupling_value
    }

    # --- Save Parameters to YAML ---

    # Define the parameters dictionary to be saved, including all relevant info
    run_params_to_save = {
        'run_specific_parameters': {
            'seed': seed,
            'time_fl_us': time_fl,
            'time_isc_us': time_isc,
            'coupling_value_cm-1': coupling_value,
        },
        'fixed_parameters': {
             'ct_energy_cm-1': CT_ENERGY,
             'pdb_file': pdb_file,
             'e0_a_cm-1': 14900,
             'e0_b_cm-1': 0,
             'N_ens': 100,
             'temp_K': 300,
             'dielectric_tresp': 1,
             'dielectric_dipole': 1.5,
             'dielectric_cdc': 2,
             'ct_pigment': '502'
        },
        'calculation_settings': {
            'N_t': 800000,
            'dt_fs': 10,
            'N_save': 1000
        },
        'used_parameter_set': used_parameters
    }
    save_parameters_to_yaml(run_dir, run_params_to_save) # Call the saving function
    # --- End Save Parameters ---


    # Build structures with and without CT
    isia_monomer_rp1 = construct_monomer_with_CT(
        pdb_file, dir_path=run_dir, e0_a=14900, ct_energy=CT_ENERGY, e0_b=0, N_ens=100, 
        temp=300, dielectric_tresp=1, dielectric_dipole=1.5, dielectric_cdc=2, 
        coupling_ct=coupling_value, ct_pigment='502'
    )

    def construct_ratematrix_with_rp1(protein_atomic, seed, run_dir, time_fl, time_isc, time_to_rp1):
        """Construct rate matrix with RP1 state"""
        protein_atomic.H2_hamiltonian.add_disorder(seed)
        nonlinear_isia = NonlinearSpectra(protein_atomic, temp=300, k_ij=None)

        # Add RP1 coupled to CT
        nonlinear_isia.add_rp1_to_ct(
            ct_pigment_name='Isia_Z_CLA_600', 
            k_ct=0,  # forward rate (CT → RP1)
            k_rp1=1/(time_to_rp1*10**6) if time_to_rp1 > 0 else 0  # backward rate (RP1 → CT)
        )
        
        nonlinear_isia.add_fluorescence(1 / (time_fl * 10 ** 6))
        nonlinear_isia.add_isc(1 / (time_isc * 10 ** 6))
        nonlinear_isia._correct_diagonal()
        
        rate_matrix_dir = os.path.join(run_dir, 'Rate_matrix')
        os.makedirs(rate_matrix_dir, exist_ok=True)
        
        # Include seed in the rate matrix filename
        sparse.save_npz(os.path.join(rate_matrix_dir, f'K_ij_isia_monomer_RP1_{seed}.npz'), nonlinear_isia.k_ij)
        return nonlinear_isia

    # Define function for calculating fluorescence
    def calculate_fluorescence(nonlinear_isia, prob_i0, run_dir, seed, suffix=""):
        """Calculate time-resolved fluorescence"""
        logger.info('Calculating the time resolved fluorescence')
        
        N_t = 800000
        dt = 10
        N_save = 1000
        list_t, w_axis, list_fl_t = nonlinear_isia.calculate_time_resolved_fluorescence(
            P1_0=prob_i0, N_t=N_t, dt=dt, N_save=N_save
        )
        lambda_axis = c / ((w_axis) / (2 * np.pi * hbar))
        
        # Save the data
        data_dir = os.path.join(run_dir, 'Data')
        os.makedirs(data_dir, exist_ok=True)

        # Include seed in the data filename
        data_file = os.path.join(data_dir, f'FL_decay_monomer{suffix}_{seed}.npz')
        
        # Reduce time points and convert to float32 for better compression
        list_t_reduced = np.array(list_t[::10], dtype=np.float32)
        lambda_axis_reduced = np.array(lambda_axis, dtype=np.float32)
        list_fl_t_reduced = np.array(list_fl_t[::10], dtype=np.float32)
        
        np.savez_compressed(data_file, 
                            list_t=list_t_reduced, 
                            lambda_axis=lambda_axis_reduced, 
                            list_fl_t=list_fl_t_reduced)
        
        logger.info(f'Data saved to {data_file}')
        
        return list_t, lambda_axis, list_fl_t

    # Build the rate matrices
    logger.info(f'Calculating the rate matrices with seed: {seed}')
    nonlinear_isia_with_rp1 = construct_ratematrix_with_rp1(
        isia_monomer_rp1, seed, run_dir, time_fl, time_isc, time_to_rp1
    )

    # Calculate initial probability distributions
    logger.info('Calculating initial probability distributions')
    prob_i0_rp1 = np.array([np.sum(['A_CLA' in name for name in domain.list_names]) 
                           for domain in nonlinear_isia_with_rp1.list_domain])
    prob_i0_rp1 = prob_i0_rp1 / np.sum(prob_i0_rp1)

    # Calculate fluorescence
    calculate_fluorescence(nonlinear_isia_with_rp1, prob_i0_rp1, run_dir, seed, suffix="_RP1")

    logger.info(f"Analysis complete for seed {seed}. Results saved in {run_dir}")

def main():
    """Main function to run the script, supporting both single and parallel execution."""
    args = parse_args()
    seed = args.seed
    
    # Load the specified configuration
    try:
        config_params = get_parameter_set(args.config)
    except ValueError as e:
        print(f"ERROR: {e}")
        print(f"Available configurations: {list_configurations()}")
        sys.exit(1)
    
    print(f"Using configuration: {args.config}")
    
    # Generate all parameter combinations from the config
    param_keys = list(config_params.keys())
    param_values = [config_params[key] for key in param_keys]
    all_combinations = list(product(*param_values))
    total_combinations = len(all_combinations)
    
    # If just listing parameters, print them and exit
    if args.list_params:
        print(f"Total parameter combinations: {total_combinations}")
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
        if args.param_index < 0 or args.param_index >= total_combinations:
            print(f"ERROR: param_index ({args.param_index}) out of range (0-{total_combinations-1})")
            sys.exit(1)
        
        # Run only the specified combination (for parallel execution)
        combo = all_combinations[args.param_index]
        param_dict = dict(zip(param_keys, combo))
        print(f"Running combination {args.param_index+1}/{total_combinations}: {param_dict}")
        
        run_simulation(
            param_dict["TIME_FL"],
            param_dict["TIME_ISC"],
            param_dict["TIME_to_RP1"],
            param_dict["COUPLING_value"],
            seed
        )
    else:
        # Run all combinations sequentially (for single execution)
        print(f"Running all {total_combinations} parameter combinations")
        for i, combo in enumerate(all_combinations):
            param_dict = dict(zip(param_keys, combo))
            print(f"Running combination {i+1}/{total_combinations}: {param_dict}")
            
            run_simulation(
                param_dict["TIME_FL"],
                param_dict["TIME_ISC"],
                param_dict["TIME_to_RP1"],
                param_dict["COUPLING_value"],
                seed
            )

if __name__ == "__main__":
    main()
