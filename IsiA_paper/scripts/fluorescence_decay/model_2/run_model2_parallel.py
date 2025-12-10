import os
import sys
import multiprocessing.pool as pool

"""
Parallel Ensemble Runner for IsiA Model 2
Distributes multiple ensemble seeds across CPU cores on a single node.
Each node runs this script with a unique node_id.

Usage:
  python3 run_model2_parallel.py <node_id>
  
Where node_id is provided by SLURM_PROCID in multi-node runs.
"""

# Configuration
node_index = int(sys.argv[1]) if len(sys.argv) > 1 else 0
ensembles_per_node = 20        # Number of ensemble seeds per node
script_name = "run_model_2.py"
config_name = "model_2"

# Model 2 has: 1 FL × 5 ISC × 9 Coupling = 45 parameter combinations
# So: 45 combos × 20 ensembles = 900 simulations per node
# 4 nodes × 900 = 3600 total jobs

def task(seed):
    """
    Run a single ensemble member (one seed) using one CPU core.
    
    Parameters:
    -----------
    seed : int
        Unique seed identifier for disorder generation
    """
    # Restrict each task to 1 core (critical for cluster efficiency)
    os.environ['OMP_NUM_THREADS'] = '1'
    
    # Run the model_2 script with this seed
    cmd = f"python3 {script_name} --seed {seed} --config {config_name}"
    exit_code = os.system(cmd)
    
    if exit_code != 0:
        print(f"WARNING: Seed {seed} exited with code {exit_code}")
    
    return exit_code

if __name__ == '__main__':
    print(f"\n{'='*60}")
    print(f"IsiA Model 2 - Parallel Ensemble Runner")
    print(f"{'='*60}")
    print(f"Node Index: {node_index}")
    print(f"Ensembles per node: {ensembles_per_node}")
    print(f"Configuration: {config_name}")
    print(f"Script: {script_name}")
    print(f"{'='*60}\n")
    
    # Calculate seed range for this node
    start_seed = node_index * ensembles_per_node + 1
    end_seed = start_seed + ensembles_per_node
    
    print(f"Running seeds {start_seed} to {end_seed-1}")
    print(f"Total simulations this node: {ensembles_per_node} × 45 parameter combos = 900 jobs\n")
    
    # Create a process pool of 96 workers (one per core)
    with pool.Pool(96) as process_pool:
        # Use imap to run tasks in parallel and maintain order
        for i, result in enumerate(process_pool.imap(task, range(start_seed, end_seed))):
            current_seed = start_seed + i
            if result == 0:
                print(f"[Node {node_index}] Seed {current_seed}: COMPLETED")
            else:
                print(f"[Node {node_index}] Seed {current_seed}: FAILED (exit code: {result})")
    
    print(f"\n{'='*60}")
    print(f"Node {node_index} processing complete!")
    print(f"Results saved to: ./Simulation/Run_FL*_ISC*_coupling_*_CT_energy_*_{config_name}/")
    print(f"{'='*60}\n")
