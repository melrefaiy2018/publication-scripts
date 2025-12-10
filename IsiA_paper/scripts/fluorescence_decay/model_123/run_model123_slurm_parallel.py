import os
import sys
import multiprocessing.pool as pool
import logging
from itertools import product

"""
IsiA Model 123 - SLURM-based Parallel Ensemble Runner - CORRECTED VERSION

This script is designed to work with SLURM job submission.
Each node runs this script with a unique node_id.

FIXED: Now properly generates ALL parameter combinations and distributes them across cores

For single node: python3 run_model123_slurm_parallel.py 0
For multi-node: srun python3 run_model123_slurm_parallel.py $SLURM_PROCID

Following guidelines from Tarun Gera's pyhops hpcthulhu setup guide:
- Each core restricted to 1 thread (OMP_NUM_THREADS=1)
- Multiprocessing pool for parallelization
- All 96 cores utilized per node
"""

# Configuration
script_name = "run_model_123.py"
config_name = os.environ.get('ENSEMBLE_CONFIG', 'model_123_diff_Vct')
ensembles_per_param_set = int(os.environ.get('ENSEMBLE_SIZE', '20'))  # Changed: one seed per param set

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

def task(param_index_and_seed):
    """
    Run a single parameter combination with a specific seed on one CPU core.
    
    Following hpcthulhu best practices:
    - Each task restricted to 1 core
    - OMP_NUM_THREADS=1 to prevent oversubscription
    
    Parameters:
    -----------
    param_index_and_seed : tuple
        (param_index, seed) - parameter combination index and seed value
    
    Returns:
    --------
    tuple : (success: bool, param_index: int, seed: int)
    """
    param_index, seed = param_index_and_seed
    
    # Restrict this task to 1 core (critical for efficiency)
    os.environ['OMP_NUM_THREADS'] = '1'
    
    # Run the simulation with both param_index and seed
    # THIS IS THE KEY FIX: Now includes --param_index!
    cmd = f"python3 {script_name} --seed {seed} --config {config_name} --param_index {param_index}"
    exit_code = os.system(cmd)
    
    return (exit_code == 0, param_index, seed)

if __name__ == '__main__':
    print("\n" + "="*80)
    print("IsiA Model 123 - SLURM Parallel Ensemble Runner - CORRECTED VERSION")
    print("="*80)
    print(f"Node Index: {node_index}")
    print(f"Configuration: {config_name}")
    print(f"Ensemble members per parameter set: {ensembles_per_param_set}")
    print(f"Script: {script_name}")
    print("="*80 + "\n")
    
    logger.info(f"Starting Node {node_index} execution")
    logger.info(f"Configuration: {config_name}")
    logger.info(f"Ensemble size: {ensembles_per_param_set}")
    
    # Get total parameter combinations - THIS IS KEY!
    try:
        from unified_parameters import get_parameter_set
        
        params = get_parameter_set(config_name)
        keys = list(params.keys())
        values = [params[key] for key in keys]
        combinations = list(product(*values))
        num_combos = len(combinations)
        
        logger.info(f"Parameter combinations to run: {num_combos}")
        
        print(f"\n✓ Parameter Information:")
        print(f"  Configuration: {config_name}")
        print(f"  Total parameter combinations: {num_combos}")
        for key in keys:
            print(f"    - {key}: {len(values[keys.index(key)])} values")
        
    except Exception as e:
        logger.error(f"Failed to load parameter combinations: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
    
    # CORRECTED: Distribute parameter combinations across nodes and cores
    # Each node gets a portion of ALL parameter combinations
    num_nodes = int(os.environ.get('SLURM_JOB_NUM_NODES', 1))
    total_cores_available = os.cpu_count() or 96
    
    # For this node, determine which parameter combinations to run
    # Simple distribution: split combinations across nodes
    combos_per_node = max(1, (num_combos + num_nodes - 1) // num_nodes)
    start_combo_idx = node_index * combos_per_node
    end_combo_idx = min(start_combo_idx + combos_per_node, num_combos)
    
    # Create list of (param_index, seed) tuples to run
    tasks_to_run = []
    for param_idx in range(start_combo_idx, end_combo_idx):
        for seed_idx in range(1, ensembles_per_param_set + 1):
            tasks_to_run.append((param_idx, seed_idx))
    
    total_jobs = len(tasks_to_run)
    
    print(f"\n✓ Node Distribution:")
    print(f"  Total nodes: {num_nodes}")
    print(f"  This node (Node {node_index}): parameter combinations {start_combo_idx} to {end_combo_idx-1}")
    print(f"  Combinations this node: {end_combo_idx - start_combo_idx}")
    print(f"  Ensemble members per combo: {ensembles_per_param_set}")
    print(f"  Total jobs this node: {total_jobs}")
    print(f"  CPU cores available: {total_cores_available}")
    print("="*80 + "\n")
    
    logger.info(f"Node {node_index} running parameter combinations {start_combo_idx}-{end_combo_idx-1}")
    logger.info(f"Total jobs this node: {total_jobs}")
    logger.info(f"Available cores: {total_cores_available}")
    
    successful = 0
    failed = 0
    failed_jobs = []
    
    # Create process pool with available cores
    num_workers = min(total_cores_available, total_jobs)
    logger.info(f"Creating process pool with {num_workers} workers")
    
    with pool.Pool(num_workers) as process_pool:
        # Use imap_unordered to run tasks in parallel (order doesn't matter)
        for i, (success, param_idx, seed) in enumerate(
            process_pool.imap_unordered(task, tasks_to_run), 1
        ):
            if success:
                successful += 1
                logger.info(f"Job {i}/{total_jobs}: param_index={param_idx}, seed={seed} - COMPLETED")
                if i % 10 == 0 or i == total_jobs:
                    print(f"[Node {node_index}] Progress: {i}/{total_jobs} jobs completed")
            else:
                failed += 1
                failed_jobs.append((param_idx, seed))
                logger.error(f"Job {i}/{total_jobs}: param_index={param_idx}, seed={seed} - FAILED")
                print(f"[Node {node_index}] Job {i}/{total_jobs}: param_index={param_idx}, seed={seed} - ✗")
    
    print("\n" + "="*80)
    print(f"Node {node_index} - Execution Summary")
    print("="*80)
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Total: {successful + failed}")
    
    if failed > 0:
        print(f"\nFailed jobs (first 10):")
        for param_idx, seed in failed_jobs[:10]:
            print(f"  - param_index={param_idx}, seed={seed}")
        if len(failed_jobs) > 10:
            print(f"  ... and {len(failed_jobs) - 10} more")
    
    print("="*80 + "\n")
    
    logger.info(f"Node {node_index} complete. Successful: {successful}, Failed: {failed}")
    
    if failed > 0:
        sys.exit(1)
    else:
        sys.exit(0)