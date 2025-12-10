#!/bin/bash
# bash submit_analysis_job.sh --time 24:00:00
# SLURM-based Job Submission for Analysis
# Usage: bash submit_analysis_job.sh [OPTIONS]

set -e

echo "╔════════════════════════════════════════════════════════════╗"
echo "║  Analysis Job - SLURM Submission                           ║"
echo "║  run_best_model_pick.py                              ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

# Configuration
NUM_NODES=1
CPUS_PER_TASK=96
TIME_LIMIT="24:00:00"
EMAIL="melrefaiy@utexas.edu"
JOB_NAME="model_14_analysis"

# Parse arguments
while [ $# -gt 0 ]; do
    case $1 in
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
        --name)
            JOB_NAME="$2"
            shift 2
            ;;
        --help)
            echo "Usage: bash submit_analysis_job.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --nodes N                Number of nodes (default: 1)"
            echo "  --cpus N                 CPUs per task (default: 96)"
            echo "  --time HH:MM:SS          Time limit (default: 24:00:00)"
            echo "  --email EMAIL            Email for SLURM notifications"
            echo "  --name NAME              Job name (default: analysis_job)"
            echo "  --help                   Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Display configuration
echo "Configuration:"
echo "  Nodes: $NUM_NODES"
echo "  CPUs per task: $CPUS_PER_TASK"
echo "  Time limit: $TIME_LIMIT"
echo "  Job name: $JOB_NAME"
echo ""

# Verify required files
echo "🔍 Verifying files..."
if [ ! -f "run_best_model_pick.py" ]; then
    echo "  ✗ run_best_model_pick.py NOT FOUND"
    exit 1
else
    echo "  ✓ run_best_model_pick.py"
fi

echo ""

# Create directories
mkdir -p logs
mkdir -p slurm_jobs

# Generate SLURM job script
SLURM_SCRIPT="slurm_jobs/analysis_job_${JOB_NAME}.sub"

cat > "$SLURM_SCRIPT" <<'SLURM_EOF'
#!/bin/bash

# SLURM Analysis Job - Generated Script
# Runs run_best_model_pick.py

#SBATCH -J @JOB_NAME@           # Job name
#SBATCH -o logs/analysis_%j.out # Output file
#SBATCH -e logs/analysis_%j.err # Error file
#SBATCH -N @NUM_NODES@          # Number of nodes
#SBATCH -n @NUM_NODES@          # Number of MPI tasks (one per node)
#SBATCH --cpus-per-task=@CPUS_PER_TASK@  # CPUs per task
#SBATCH -t @TIME_LIMIT@         # Time limit
#SBATCH --mail-type=all         # Email notifications
#SBATCH --mail-user=@EMAIL@     # Email address

# System information
echo "╔════════════════════════════════════════════════════════════╗"
echo "║          Analysis Job - SLURM Run                          ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""

echo "SLURM Job Information:"
echo "  Job ID: $SLURM_JOB_ID"
echo "  Job Name: $SLURM_JOB_NAME"
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

echo "═══════════════════════════════════════════════════════════"
echo "Starting Analysis..."
echo "═══════════════════════════════════════════════════════════"
echo ""

# Run the analysis script
python3 run_best_model_pick.py

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "Analysis Complete"
echo "═══════════════════════════════════════════════════════════"
date
echo ""

SLURM_EOF

# Replace placeholders
sed -i "s|@NUM_NODES@|$NUM_NODES|g" "$SLURM_SCRIPT"
sed -i "s|@CPUS_PER_TASK@|$CPUS_PER_TASK|g" "$SLURM_SCRIPT"
sed -i "s|@TIME_LIMIT@|$TIME_LIMIT|g" "$SLURM_SCRIPT"
sed -i "s|@EMAIL@|$EMAIL|g" "$SLURM_SCRIPT"
sed -i "s|@JOB_NAME@|$JOB_NAME|g" "$SLURM_SCRIPT"
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
echo "    tail -f logs/analysis_${JOB_ID}.out"
echo ""
echo "  Check for errors:"
echo "    tail -f logs/analysis_${JOB_ID}.err"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

echo "📝 Job Details:"
echo "  Job Name: $JOB_NAME"
echo "  Nodes: $NUM_NODES"
echo "  CPUs: $CPUS_PER_TASK"
echo "  Time limit: $TIME_LIMIT"
echo "  Log file: logs/analysis_${JOB_ID}.out"
echo ""

echo "✅ Ready to run!"
