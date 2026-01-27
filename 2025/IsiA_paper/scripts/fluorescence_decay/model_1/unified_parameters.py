# unified_parameters.py
"""
Unified parameter configuration for all simulation runs.
This file serves as the single source of truth for parameter sets.
"""

# Global constant for CT energy
CT_ENERGY = 14800

# Define different parameter sets for different types of simulations
PARAMETER_SETS = {
    # Configuration 1: Basic CT model (used in run_model_1.py)
    "model_1" : {
    "TIME_FL": [16],
    "TIME_ISC": [1, 2, 3, 5, 10],
    "COUPLING_value": [20, 40, 60, 80, 100, 150, 200, 250, 300]
    },

    # Configuration 2: Basic CT model (used in run_model_2.py)
    "model_2": {
        "TIME_FL": [16],
        "TIME_ISC": [1, 2, 3, 5, 10],
        "COUPLING_value": [20, 40, 60, 80, 100, 150, 200, 250, 300]
    },
    # Configuration 3: RP1 model (used in run_model_3.py)
    "model_3" : {
    "TIME_FL": [16],
    "TIME_ISC": [3, 5],
    "TIME_to_RP1": [0.2, 0.4, 0.6, 0.8, 1, 2, 3],
    "COUPLING_value": [20, 50, 100, 150, 200, 250, 300]
    },

    # Configuration 4: CT/RP1 model (used in run_model_4.py)
    "model_4" : {
        "TIME_FL": [16],
        "TIME_ISC": [3, 5],
        "TIME_to_CT": [0.2, 0.4, 0.6, 0.8, 1, 2, 3, 8],
        "TIME_to_RP1": [0.2, 0.4, 0.6, 0.8, 1, 2, 3],
        "COUPLING_value": [20, 50, 100, 150, 200, 250, 300]
    },

    "model_12" : {
    "TIME_FL": [16],
    "TIME_ISC": [1, 2, 2.5, 2.7, 2.8, 3, 4, 5],
    "COUPLING_value": [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 150, 180, 200, 250, 300]
    },
    # Configuration 13: CT/RP1 model (used in run_model_13.py)

    "model_13" : {
        "TIME_FL": [16],
        "TIME_ISC": [3, 5],
        "TIME_to_RP1": [0.2, 0.4, 0.6, 0.8, 1, 2, 3],
        "COUPLING_value": [20, 50, 100, 150, 200, 250, 300]
    },   
    # Configuration 4: Test/Debug (small parameter set for quick testing)
    "test": {
        "TIME_FL": [16],
        "TIME_ISC": [5],
        "COUPLING_value": [100, 200]
    }

}
# Default configuration to use if none specified
DEFAULT_CONFIG = "model_1"

def get_parameter_set(config_name=None):
    """
    Get parameter set by configuration name.
    
    Parameters:
    -----------
    config_name : str, optional
        Name of the configuration to get. If None, returns default.
    
    Returns:
    --------
    dict
        Parameter set dictionary
    """
    if config_name is None:
        config_name = DEFAULT_CONFIG
    
    if config_name not in PARAMETER_SETS:
        raise ValueError(f"Configuration '{config_name}' not found. Available: {list(PARAMETER_SETS.keys())}")
    
    return PARAMETER_SETS[config_name]

def list_configurations():
    """List all available configurations."""
    return list(PARAMETER_SETS.keys())

def print_configuration(config_name):
    """Print the details of a specific configuration."""
    params = get_parameter_set(config_name)
    print(f"\nConfiguration: {config_name}")
    print("-" * 40)
    for key, values in params.items():
        print(f"{key}: {values}")
    print("-" * 40)
