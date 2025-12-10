#!/bin/bash

# Set your Simulation directory
base_dir="."  # <-- Change this to your actual path

# Loop through each 'Run_' directory
for run_dir in "$base_dir"/Run_*; do
  if [ -d "$run_dir" ]; then
    echo "Checking $run_dir"

    # If Analysis_Data exists, delete it
    if [ -d "$run_dir/Analysis_Data" ]; then
      echo "Deleting $run_dir/Analysis_Data"
      rm -rf "$run_dir/Analysis_Data"
    fi

    # If Analysis_Figures exists, delete it
    if [ -d "$run_dir/Analysis_Figures" ]; then
      echo "Deleting $run_dir/Analysis_Figures"
      rm -rf "$run_dir/Analysis_Figures"
    fi

  fi
done

echo "Done!"
