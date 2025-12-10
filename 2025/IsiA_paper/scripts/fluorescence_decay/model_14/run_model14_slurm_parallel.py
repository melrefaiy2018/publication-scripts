import os
import sys
import multiprocessing.pool as pool
import logging

"""
IsiA Model 2 - SLURM-based Parallel Ensemble Runner

This script is designed to work with SLURM job submission.
Each node runs this script with a unique node_id.

For single node: python3 run_model3_slurm_parallel.py 0
For multi-node: srun python3 run_model3_slurm_parallel.py $SLURM_PROCID

Following guidelines from Tarun Gera's pyhops hpcthulhu setup guide:
- Each core restricted to 1 thread (OMP_NUM_THREADS=1)
- Multiprocessing pool for parallelization
- All 96 cores utilized per node
"""

# Configuration
script_name = "run_model_14.py"
config_name = os.environ.get('ENSEMBLE_CONFIG', 'model_14')
ensembles_per_node = int(os.environ.get('ENSEMBLE_SIZE', '20'))

# Determine node index
if 'SLURM_PROCID' in os.environ:
    node_index = int(os.environ['SLURM_PROCID'])
else:
    node_index = int(sys.argv[1]) if len(sys.argv) > 1 else 0

# Setup logging
log_dir = 'logs'
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f'node_{node_index}.log')

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - Node {0} - %(levelname)s - %(message)s'.format(node_index),
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger()

def task(seed):
    """
    Run a single ensemble member (one seed) on one CPU core.
    
    Following hpcthulhu best practices:
    - Each task restricted to 1 core
    - OMP_NUM_THREADS=1 to prevent oversubscription
    
    Parameters:
    -----------
    seed : int
        Unique seed identifier for disorder generation
    
    Returns:
    --------
    int : Exit code (0 = success, non-zero = failure)
    """
    # Restrict this task to 1 core (critical for efficiency)
    os.environ['OMP_NUM_THREADS'] = '1'
    
    # Run the simulation
    cmd = f"python3 {script_name} --seed {seed} --config {config_name}"
    exit_code = os.system(cmd)
    
    return exit_code

if __name__ == '__main__':
    print("\n" + "="*70)
    print("IsiA Model 2 - SLURM Parallel Ensemble Runner")
    print("="*70)
    print(f"Node Index: {node_index}")
    print(f"Configuration: {config_name}")
    print(f"Ensemble members per node: {ensembles_per_node}")
    print(f"Script: {script_name}")
    print("="*70 + "\n")
    
    logger.info(f"Starting Node {node_index} execution")
    logger.info(f"Configuration: {config_name}")
    logger.info(f"Ensemble members: {ensembles_per_node}")
    
    # Get total parameter combinations
    try:
        from itertools import product
        from unified_parameters import get_parameter_set
        
        params = get_parameter_set(config_name)
        keys = list(params.keys())
        values = [params[key] for key in keys]
        combinations = list(product(*values))
        num_combos = len(combinations)
        
        logger.info(f"Parameter combinations to run: {num_combos}")
    except Exception as e:
        logger.error(f"Failed to load parameter combinations: {e}")
        sys.exit(1)
    
    # Calculate seed range for this node
    start_seed = node_index * ensembles_per_node + 1
    end_seed = start_seed + ensembles_per_node
    
    total_jobs = num_combos * ensembles_per_node
    
    print(f"Node {node_index} will run seeds {start_seed} to {end_seed-1}")
    print(f"Per seed: {num_combos} parameter combinations")
    print(f"Total jobs this node: {total_jobs}")
    print(f"Total cores available: 96")
    print("="*70 + "\n")
    
    logger.info(f"Running seeds {start_seed} to {end_seed-1}")
    logger.info(f"Total simulations this node: {total_jobs}")
    
    # Create process pool with 96 workers (one per core)
    num_cores = os.cpu_count() or 96
    logger.info(f"Creating process pool with {num_cores} workers")
    
    successful = 0
    failed = 0
    
    with pool.Pool(num_cores) as process_pool:
        # Use imap to run tasks and maintain order
        for i, result in enumerate(process_pool.imap(task, range(start_seed, end_seed))):
            current_seed = start_seed + i
            
            if result == 0:
                successful += 1
                logger.info(f"Seed {current_seed}: COMPLETED")
                print(f"[Node {node_index}] Seed {current_seed}: ✓")
            else:
                failed += 1
                logger.error(f"Seed {current_seed}: FAILED (exit code: {result})")
                print(f"[Node {node_index}] Seed {current_seed}: ✗ (exit code {result})")
    
    print("\n" + "="*70)
    print(f"Node {node_index} - Execution Summary")
    print("="*70)
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total: {successful + failed}")
    print("="*70 + "\n")
    
    logger.info(f"Node {node_index} complete. Successful: {successful}, Failed: {failed}")
    
    if failed > 0:
        sys.exit(1)
    else:
        sys.exit(0)
