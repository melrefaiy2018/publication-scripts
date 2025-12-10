import yaml
import os
import logging

# Set up a basic logger for the function
# In a larger application, you might pass a logger instance or use a configured root logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_parameter_scan_definition(run_dir):
    """
    Loads the parameter_scan_definition dictionary from the run_parameters.yaml
    file located within the 'param' subdirectory of the specified run directory.

    Parameters:
    -----------
    run_dir : str
        The path to the specific simulation run directory (e.g.,
        'Simulation/Run_FL16_ISC5_coupling_250_CTenergy_14450_seed_1').

    Returns:
    --------
    dict or None:
        The dictionary stored under the 'parameter_scan_definition' key in the
        YAML file, or None if the file is not found, cannot be read, or the
        key is missing.
    """
    param_file_path = os.path.join(run_dir, 'param', 'run_parameters.yaml')
    logger.info(f"Attempting to load parameters from: {param_file_path}")

    if not os.path.exists(param_file_path):
        logger.error(f"Parameter file not found at: {param_file_path}")
        return None

    try:
        with open(param_file_path, 'r') as f:
            all_params = yaml.safe_load(f) # Use safe_load for security

        if all_params is None:
            logger.error(f"YAML file is empty or invalid: {param_file_path}")
            return None

        # Check if the expected key exists
        if 'used_parameter_set' in all_params:
            param_scan_def = all_params['used_parameter_set']
            logger.info(f"Successfully loaded parameter_scan_definition from {param_file_path}")
            return param_scan_def
        else:
            logger.error(f"'used_parameter_set' key not found in {param_file_path}")
            return None

    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML file {param_file_path}: {e}")
        return None
    except IOError as e:
        logger.error(f"Error reading file {param_file_path}: {e}")
        return None
    except Exception as e:
        logger.error(f"An unexpected error occurred while loading parameters: {e}")
        return None

# # Correct path: Ends with the run directory name
# correct_run_directory = '/Users/mohamed/Documents/Research/Projects/IsiA/IsiA_2025/IsiA_monomer/CT_state/sorted/Apr23/scaning_models/model_2/Simulation/Run_FL16_ISC5_coupling_250_CTenergy_14450/'

# # Call the function with the correct directory path
# dict_loaded_params = load_parameter_scan_definition(correct_run_directory)

# if dict_loaded_params:
#     print("Successfully loaded parameters:")
#     print(dict_loaded_params)
# else:
#     print("Failed to load parameters. Check logs for errors.")
