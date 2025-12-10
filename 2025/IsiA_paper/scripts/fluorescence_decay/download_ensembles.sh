#!/bin/bash

#=============================================================================
# Configuration - Edit these variables for your setup
#=============================================================================
REMOTE_USER="mae3742"
REMOTE_HOST="hpcthulhu-head.cm.utexas.edu"
REMOTE_BASE_PATH="/home/mae3742/projects/IsiA/fl_decay"
LOCAL_DESTINATION="/Users/mohamed/Library/CloudStorage/Box-Box/MesoscienceLab/ActiveProjects/IsiA/IsiA_manuscript/Data/IsiA_data_and_code/simulation_data"

#=============================================================================
# Main download function
#=============================================================================
echo "=========================================="
echo "Starting download of Ensemble directories"
echo "=========================================="
echo "Remote: ${REMOTE_USER}@${REMOTE_HOST}:${REMOTE_BASE_PATH}"
echo "Local: ${LOCAL_DESTINATION}"
echo ""

# Create local destination if it doesn't exist
mkdir -p "${LOCAL_DESTINATION}"

# Discover all model directories on the remote server
echo "Discovering all model directories on server..."

# Store models in an array (bash 3.2 compatible way)
MODELS=()
while IFS= read -r model_name; do
    [ -n "$model_name" ] && MODELS+=("$model_name")
done < <(ssh "${REMOTE_USER}@${REMOTE_HOST}" "cd ${REMOTE_BASE_PATH} && ls -d model_* 2>/dev/null" | sort)

if [ ${#MODELS[@]} -eq 0 ]; then
    echo "ERROR: No model directories found at ${REMOTE_BASE_PATH}"
    exit 1
fi

echo "Found ${#MODELS[@]} model(s) to process"
echo "Models: ${MODELS[*]}"
echo ""

# Counter for total downloads
total_downloads=0
successful_downloads=0
failed_downloads=0

# Loop through each model
for model in "${MODELS[@]}"; do
    echo ""
    echo "============================================"
    echo "Processing ${model}..."
    echo "============================================"
    
    # Check if Simulation directory exists
    if ! ssh "${REMOTE_USER}@${REMOTE_HOST}" "[ -d '${REMOTE_BASE_PATH}/${model}/Simulation' ]" 2>/dev/null; then
        echo "  ⚠ No Simulation directory found for ${model}"
        continue
    fi
    
    # Get list of simulation directories for this model
    echo "  Searching for Run_* directories..."
    
    # Store simulations in an array
    SIMULATIONS=()
    while IFS= read -r sim_name; do
        [ -n "$sim_name" ] && SIMULATIONS+=("$sim_name")
    done < <(ssh "${REMOTE_USER}@${REMOTE_HOST}" "cd ${REMOTE_BASE_PATH}/${model}/Simulation && ls -d Run_* 2>/dev/null" | sort)
    
    if [ ${#SIMULATIONS[@]} -eq 0 ]; then
        echo "  ⚠ No Run_* directories found for ${model}"
        continue
    fi
    
    echo "  Found ${#SIMULATIONS[@]} simulation(s)"
    echo ""
    
    # Loop through each simulation directory
    for sim_name in "${SIMULATIONS[@]}"; do
        # Check if Ensemble directory exists
        ensemble_path="${REMOTE_BASE_PATH}/${model}/Simulation/${sim_name}/Data/Ensemble"
        
        echo "  [${model}] Checking: ${sim_name}"
        
        # Verify the Ensemble directory exists before attempting download
        if ssh "${REMOTE_USER}@${REMOTE_HOST}" "[ -d '${ensemble_path}' ]" 2>/dev/null; then
            echo "    → Ensemble found! Downloading..."
            
            # Create local directory structure
            local_model_dir="${LOCAL_DESTINATION}/${model}/${sim_name}/Data"
            mkdir -p "$local_model_dir"
            
            # Download the Ensemble directory using SCP with recursive option
            # Remove -q flag to see progress
            scp -r "${REMOTE_USER}@${REMOTE_HOST}:${ensemble_path}" "${local_model_dir}/" > /dev/null 2>&1
            
            if [ $? -eq 0 ]; then
                echo "    ✓ Success"
                successful_downloads=$((successful_downloads + 1))
            else
                echo "    ✗ Failed"
                failed_downloads=$((failed_downloads + 1))
            fi
            total_downloads=$((total_downloads + 1))
        else
            echo "    ✗ Ensemble directory not found, skipping"
        fi
    done
    
    echo ""
    echo "  ${model} complete: ${#SIMULATIONS[@]} simulations processed"
done

echo ""
echo "=========================================="
echo "Download Summary"
echo "=========================================="
echo "Total downloads attempted: ${total_downloads}"
echo "Successful: ${successful_downloads}"
echo "Failed: ${failed_downloads}"
echo "Data saved to: ${LOCAL_DESTINATION}"
echo "=========================================="