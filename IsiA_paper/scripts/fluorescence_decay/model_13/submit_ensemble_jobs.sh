#!/bin/bash
# bash submit_ensemble_jobs.sh --config model_13 --ensemble-size 20

# SLURM-based Ensemble Run for IsiA Model 13
# Follows hpcthulhu best practices from Tarun Gera's guide
# Usage: bash submit_ensemble_jobs.sh [OPTIONS]

set -e

echo "╔════════════════════════════════════════════════════════════╗"
echo "║  IsiA Model13 - SLURM Ensemble Submission                  ║"
echo "║  All Parameters × 20 Ensemble Members                      ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Configuration
CONFIG="model_13"
ENSEMBLE_SIZE=20
NUM_NODES=1
CPUS_PER_TASK=96
TIME_LIMIT="72:00:00"
EMAIL="melrefaiy@utexas.edu"

# Parse arguments
while [ $# -gt 0 ]; do
    case $1 in
        --config)
            CONFIG="$2"
            shift 2
            ;;
        --ensemble-size)
            ENSEMBLE_SIZE="$2"
            shift 2
            ;;
        --nodes)
            NUM_NODES="$2"
            shift 2
            ;;
        --cpus)
            CPUS_PER_TASK="$2"
            shift 2
            ;;
        --time)
            TIME_LIMIT="$2"
            shift 2
            ;;
        --email)
            EMAIL="$2"
            shift 2
            ;;
        --help)
            echo "Usage: bash submit_ensemble_jobs.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --config CONFIG          Configuration (default: model_13)"
            echo "  --ensemble-size N        Ensemble members (default: 20)"
            echo "  --nodes N                Number of nodes (default: 1)"
            echo "  --cpus N                 CPUs per task (default: 96)"
            echo "  --time HH:MM:SS          Time limit (default: 72:00:00)"
            echo "  --email EMAIL            Email for SLURM notifications"
            echo "  --help                   Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Calculate expected statistics
echo "Configuration: $CONFIG"
echo "Ensemble size: $ENSEMBLE_SIZE"
echo "Nodes: $NUM_NODES"
echo "CPUs per node: $CPUS_PER_TASK"
echo ""

# Get parameter count
PARAM_COUNT=$(python3 -c "
from unified_parameters import get_parameter_set
from itertools import product

params = get_parameter_set('$CONFIG')
keys = list(params.keys())
values = [params[key] for key in keys]
combinations = list(product(*values))
print(len(combinations))
")

TOTAL_JOBS=$((PARAM_COUNT * ENSEMBLE_SIZE))

echo "📊 Job Statistics:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Parameter combinations: $PARAM_COUNT"
echo "  Ensemble members each: $ENSEMBLE_SIZE"
echo "  Total simulations: $TOTAL_JOBS"
echo "  Expected output: ~$((PARAM_COUNT * ENSEMBLE_SIZE * 70 / 1000)) GB"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Verify required files
echo "🔍 Verifying files..."
for file in run_model_13.py extended_most_occ_pH7_CT.pdb unified_parameters.py construct_monomer_with_CT.py; do
    if [ ! -f "$file" ]; then
        echo "  ✗ $file NOT FOUND"
        exit 1
    else
        echo "  ✓ $file"
    fi
done

echo ""

# Create directories
mkdir -p logs
mkdir -p slurm_jobs

# Generate SLURM job script
SLURM_SCRIPT="slurm_jobs/ensemble_run_${CONFIG}_${ENSEMBLE_SIZE}.sub"

cat > "$SLURM_SCRIPT" <<'SLURM_EOF'
#!/bin/bash

# SLURM Ensemble Run - Generated Script
# Based on hpcthulhu guidelines from Tarun Gera's pyhops guide

#SBATCH -J model13_@CONFIG@_@ENSEMBLE_SIZE@  # Job name
#SBATCH -o logs/ensemble_%j.out    # Output file
#SBATCH -e logs/ensemble_%j.err    # Error file
#SBATCH -N @NUM_NODES@             # Number of nodes
#SBATCH -n @NUM_NODES@             # Number of MPI tasks (one per node)
#SBATCH --cpus-per-task=@CPUS_PER_TASK@  # CPUs per task
#SBATCH -t @TIME_LIMIT@            # Time limit
#SBATCH --mail-type=all            # Email notifications
#SBATCH --mail-user=@EMAIL@        # Email address

# System information
echo "╔════════════════════════════════════════════════════════════╗"
echo "║          IsiA Model13 - SLURM Ensemble Run                 ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

echo "SLURM Job Information:"
echo "  Job ID: $SLURM_JOB_ID"
echo "  Node(s): $SLURM_JOB_NODELIST"
echo "  Num Nodes: $SLURM_JOB_NUM_NODES"
echo "  CPUs per task: $SLURM_CPUS_PER_TASK"
echo "  Total CPUs: $((SLURM_JOB_NUM_NODES * SLURM_CPUS_PER_TASK))"
echo ""

module list
pwd
date
echo ""

# Activate environment
source ~/miniconda3/bin/activate pymembrane

# Change to project directory
cd @PROJECT_DIR@

echo "Project directory: $(pwd)"
echo ""

# Configuration
CONFIG="@CONFIG@"
ENSEMBLE_SIZE=@ENSEMBLE_SIZE@
NUM_NODES=@NUM_NODES@

echo "═══════════════════════════════════════════════════════════"
echo "Simulation Parameters:"
echo "  Configuration: $CONFIG"
echo "  Ensemble size: $ENSEMBLE_SIZE"
echo "  Nodes: $NUM_NODES"
echo "═══════════════════════════════════════════════════════════"
echo ""

# Get parameter combinations
python3 << 'PYTHON_BLOCK'
import sys
import itertools
from unified_parameters import get_parameter_set

params = get_parameter_set("@CONFIG@")
keys = list(params.keys())
values = [params[key] for key in keys]
combinations = list(itertools.product(*values))

print(f"Generated {len(combinations)} parameter combinations")
print(f"Ensemble members per combination: @ENSEMBLE_SIZE@")
print(f"Total simulations: {len(combinations) * @ENSEMBLE_SIZE@}")
PYTHON_BLOCK

echo ""
echo "Starting parallel simulations..."
echo ""

# Run using srun (SLURM's parallel job launcher)
# This distributes the work across nodes

if [ @NUM_NODES@ -eq 1 ]; then
    # Single node - use run_13_slurm_parallel.py
    python3 run_model13_slurm_parallel.py 0
else
    # Multiple nodes - use srun to distribute
    srun bash -c '
        NODE_ID=$((SLURM_PROCID))
        echo "Node $NODE_ID starting on $HOSTNAME"
        python3 run_model13_slurm_parallel.py $NODE_ID
        echo "Node $NODE_ID complete"
    '
fi

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "Simulation Results Summary"
echo "═══════════════════════════════════════════════════════════"

# Count output files
TOTAL_FILES=$(find Simulation -name "*.npz" 2>/dev/null | wc -l)
echo "Total output files generated: $TOTAL_FILES"
echo "Expected: $(($(python3 -c "
from unified_parameters import get_parameter_set
from itertools import product
params = get_parameter_set('@CONFIG@')
keys = list(params.keys())
values = [params[key] for key in keys]
combinations = list(product(*values))
print(len(combinations))
") * @ENSEMBLE_SIZE@))"

# Summary
echo ""
echo "═══════════════════════════════════════════════════════════"
date
echo "Job complete!"
echo "Results in: Simulation/"
echo "Logs in: logs/"
echo "═══════════════════════════════════════════════════════════"

SLURM_EOF

# Replace placeholders
sed -i "s|@NUM_NODES@|$NUM_NODES|g" "$SLURM_SCRIPT"
sed -i "s|@CPUS_PER_TASK@|$CPUS_PER_TASK|g" "$SLURM_SCRIPT"
sed -i "s|@TIME_LIMIT@|$TIME_LIMIT|g" "$SLURM_SCRIPT"
sed -i "s|@EMAIL@|$EMAIL|g" "$SLURM_SCRIPT"
sed -i "s|@CONFIG@|$CONFIG|g" "$SLURM_SCRIPT"
sed -i "s|@ENSEMBLE_SIZE@|$ENSEMBLE_SIZE|g" "$SLURM_SCRIPT"
sed -i "s|@PROJECT_DIR@|$(pwd)|g" "$SLURM_SCRIPT"

chmod +x "$SLURM_SCRIPT"

echo "✅ SLURM script created: $SLURM_SCRIPT"
echo ""

# Submit job
echo "🚀 Submitting to SLURM queue..."
echo ""

JOB_ID=$(sbatch "$SLURM_SCRIPT" | awk '{print $4}')

echo "✓ Job submitted successfully!"
echo "  Job ID: $JOB_ID"
echo ""

echo "📊 Monitoring Commands:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Check job status:"
echo "    squeue -j $JOB_ID"
echo ""
echo "  Watch output:"
echo "    tail -f logs/ensemble_${JOB_ID}.out"
echo ""
echo "  Count progress:"
echo "    watch -n 10 'find Simulation -name \"*.npz\" | wc -l'"
echo ""
echo "  Check for errors:"
echo "    tail logs/ensemble_${JOB_ID}.err"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "📝 Job Details:"
echo "  Configuration: $CONFIG"
echo "  Parameter combinations: $PARAM_COUNT"
echo "  Ensemble size: $ENSEMBLE_SIZE"
echo "  Total jobs: $TOTAL_JOBS"
echo "  Expected time: $TIME_LIMIT"
echo "  Expected output: ~$((PARAM_COUNT * ENSEMBLE_SIZE * 70 / 1000)) GB"
echo ""

echo "✅ Ready to run!"
