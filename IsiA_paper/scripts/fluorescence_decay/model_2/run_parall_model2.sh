#!/bin/bash

# Parse command line arguments
CONFIG_NAME="model_2"
ENSEMBLE_SIZE=20
RESOURCE_PERCENT=0.7

while [ $# -gt 0 ]; do
    case $1 in
        --config)
            CONFIG_NAME="$2"
            shift 2
            ;;
        --ensemble-size)
            ENSEMBLE_SIZE="$2"
            shift 2
            ;;
        --resource-percent)
            RESOURCE_PERCENT="$2"
            shift 2
            ;;
        --help)
            echo "Usage: $0 [--config CONFIG_NAME] [--ensemble-size N] [--resource-percent PCT]"
            echo "  --config: Parameter configuration to use (default: rp1_model)"
            echo "  --ensemble-size: Number of ensemble members per parameter set (default: 20)"
            echo "  --resource-percent: Fraction of CPU cores to use (default: 0.7)"
            echo ""
            echo "Available configurations can be seen with:"
            echo "  python3 run_model_2.py --list_configs"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Create directories for logs and output
mkdir -p log
mkdir -p debug

# Calculate system resources and target usage
TOTAL_CORES=$(nproc)
MAX_PROCESSES=$(printf "%.0f" $(echo "$TOTAL_CORES * $RESOURCE_PERCENT" | bc))

echo "System has $TOTAL_CORES CPU cores, using $MAX_PROCESSES processes (${RESOURCE_PERCENT}%)"
echo "Using parameter configuration: $CONFIG_NAME"
echo "Ensemble size per parameter set: $ENSEMBLE_SIZE"

# Generate parameter combinations using the Python script
echo "Generating parameter combinations for configuration: $CONFIG_NAME..."
python3 -c "
import sys
import itertools

# Import unified parameters
try:
    from unified_parameters import get_parameter_set
    params = get_parameter_set('$CONFIG_NAME')
except ImportError:
    print('ERROR: Could not import unified_parameters.py')
    print('Make sure unified_parameters.py is in the current directory')
    sys.exit(1)
except ValueError as e:
    print(f'ERROR: {e}')
    print('Available configurations:')
    from unified_parameters import list_configurations
    for config in list_configurations():
        print(f'  - {config}')
    sys.exit(1)

keys = list(params.keys())
values = [params[key] for key in keys]
combinations = list(itertools.product(*values))

# Write parameter combinations to file
with open('parameter_combinations.txt', 'w') as f:
    # Write header with parameter names
    f.write('|'.join(['idx'] + keys) + '\n')
    for i, combo in enumerate(combinations):
        f.write(f'{i}|' + '|'.join(map(str, combo)) + '\n')
        
print(f'Generated {len(combinations)} parameter combinations')
print(f'Parameter names: {keys}')
"

# Check if parameter generation was successful
if [ ! -f "parameter_combinations.txt" ]; then
    echo "ERROR: Failed to generate parameter combinations"
    exit 1
fi

# Count total parameter combinations (subtract 1 for header)
TOTAL_COMBOS=$(($(wc -l < parameter_combinations.txt) - 1))
echo "Total parameter combinations to run: $TOTAL_COMBOS"
echo "Each combination will run $ENSEMBLE_SIZE ensemble members"
echo "Total jobs to run: $((TOTAL_COMBOS * ENSEMBLE_SIZE))"

# Function to check number of running processes
count_running_jobs() {
    ps aux | grep "python3.*run_model_2.py" | grep -v grep | wc -l
}

# Initialize job tracking files
echo "" > log/completed_jobs.txt
echo "" > log/active_jobs.txt
echo "" > log/failed_jobs.txt

# Read parameter names from header
PARAM_NAMES=$(head -1 parameter_combinations.txt | cut -d'|' -f2-)

# Process each parameter combination (skip header)
# Process each parameter combination (skip header)
tail -n +2 parameter_combinations.txt | while IFS='|' read -r idx param_values; do
    # Create a parameter string for directory naming using the raw param_values
    # Replace | with _ and remove any leading/trailing spaces
    PARAM_STR=$(echo "$param_values" | sed 's/|/_/g' | sed 's/^ *//' | sed 's/ *$//')
    
    echo "Processing parameter set $((idx+1))/$TOTAL_COMBOS: $PARAM_STR"
    
    # Run ensemble members for this parameter set
    for seed in $(seq 1 $ENSEMBLE_SIZE); do
        # Wait until we have resources available
        while [ $(count_running_jobs) -ge $MAX_PROCESSES ]; do
            echo "Waiting for resources... ($(count_running_jobs)/$MAX_PROCESSES processes running)"
            sleep 5
        done
        
        # Create debug directory for this job
        debug_dir="debug/job_${idx}_seed_${seed}"
        mkdir -p "$debug_dir"
        
        echo "Starting job $((idx+1))/$TOTAL_COMBOS, ensemble member $seed/$ENSEMBLE_SIZE: $PARAM_STR"
        
        # Run with the specific parameters and seed
        (
            # Set the working directory explicitly
            cd "$(pwd)"
            
            # Run the job with parameter index, seed, and configuration
            python3 run_model_2.py --param_index $idx --total_params $TOTAL_COMBOS --seed $seed --config $CONFIG_NAME > "$debug_dir/output.log" 2> "$debug_dir/error.log"
            
            # Check status
            STATUS=$?
            if [ $STATUS -eq 0 ]; then
                echo "$idx|$seed|$PARAM_STR|SUCCESS" >> log/completed_jobs.txt
            else
                echo "$idx|$seed|$PARAM_STR|FAILED" >> log/failed_jobs.txt
                echo "Exit code: $STATUS" >> "$debug_dir/status.txt"
            fi
        ) &
        
        # Record the running job
        PID=$!
        echo "$idx|$seed|$PARAM_STR|$PID" >> log/active_jobs.txt
        
        # Sleep briefly to avoid overwhelming file system with job starts
        sleep 0.5
    done
    
    # Short pause between parameter sets
    sleep 1
done

# Wait for all jobs to complete
echo "All jobs submitted, waiting for completion..."
wait
echo "All jobs completed!"

# Summarize results
SUCCESS_COUNT=$(wc -l < log/completed_jobs.txt)
FAILED_COUNT=$(wc -l < log/failed_jobs.txt)
TOTAL_EXPECTED=$((TOTAL_COMBOS * ENSEMBLE_SIZE))

echo "Summary:"
echo "- Configuration: $CONFIG_NAME"
echo "- Total expected: $TOTAL_EXPECTED"
echo "- Successful: $SUCCESS_COUNT"
echo "- Failed: $FAILED_COUNT"

if [ $FAILED_COUNT -gt 0 ]; then
    echo "First 10 failed jobs:"
    head -10 log/failed_jobs.txt
    echo "Check debug directories for detailed error logs"
    
    # Provide command to examine a sample error log
    FIRST_FAILED=$(head -1 log/failed_jobs.txt | cut -d'|' -f1,2)
    if [ ! -z "$FIRST_FAILED" ]; then
        IDX=$(echo $FIRST_FAILED | cut -d'|' -f1)
        SEED=$(echo $FIRST_FAILED | cut -d'|' -f2)
        echo "To view error for first failed job, run: cat debug/job_${IDX}_seed_${SEED}/error.log"
    fi
fi

# Clean up temporary files
rm -f parameter_combinations.txt

echo "Run complete!"