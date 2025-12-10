import os
import subprocess
import sys
import shutil
from pathlib import Path
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, as_completed

def process_directory(directory, source_script_path):
    """
    Process a single Run_ directory by copying the analysis script and running it.
    
    Args:
        directory (Path): The Run_ directory to process
        source_script_path (Path): Path to the source run_analysis.py script
        
    Returns:
        tuple: (directory name, success status, output, error)
    """
    dir_path = directory
    dir_name = dir_path.name
    analysis_figures_path = dir_path / 'Analysis_Figures'
    
    # Skip if Analysis_Figures already exists
    if analysis_figures_path.is_dir():
        return (dir_name, False, "", f"Skipping: 'Analysis_Figures' directory already exists.")
    
    # Copy run_analysis.py into the directory
    try:
        shutil.copy2(source_script_path, dir_path)
    except Exception as e:
        return (dir_name, False, "", f"Error copying 'run_analysis.py': {e}")
    
    # Check if the script exists after copying
    script_path = dir_path / 'run_analysis.py'
    if not script_path.exists():
        return (dir_name, False, "", f"Script not found after copying.")
    
    # Run the analysis script
    current_dir = Path.cwd()
    try:
        os.chdir(dir_path)
        result = subprocess.run([sys.executable, 'run_analysis.py'],
                               capture_output=True, text=True, check=False)
        
        if result.returncode == 0:
            return (dir_name, True, result.stdout, result.stderr)
        else:
            return (dir_name, False, result.stdout, f"Analysis finished with errors (return code: {result.returncode})\n{result.stderr}")
            
    except Exception as e:
        return (dir_name, False, "", f"An unexpected error occurred: {e}")
    finally:
        os.chdir(current_dir)

def run_analysis_in_parallel(base_path='.', max_workers=None):
    """
    Finds directories starting with 'Run_', and processes them in parallel.
    If an 'Analysis_Figures' subdirectory exists, the directory is skipped.
    Otherwise, it copies 'run_analysis.py' into the directory and executes it.
    
    Args:
        base_path (str): The base directory to search for 'Run_' directories.
        max_workers (int, optional): Maximum number of worker processes to use.
            If None, uses the number of CPU cores.
    """
    # If max_workers is None, use the number of CPU cores
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    # Save the directory where the main script was launched
    main_script_dir = Path.cwd()
    source_script_path = main_script_dir / 'run_analysis.py'

    # Check if the source script exists
    if not source_script_path.exists():
        print(f"Error: Source script '{source_script_path}' not found. Exiting.")
        return

    # Get all items in the base directory
    try:
        base_path = Path(base_path)
        items = list(base_path.iterdir())
    except FileNotFoundError:
        print(f"Error: Base directory '{base_path}' not found.")
        return
    except Exception as e:
        print(f"Error listing directory '{base_path}': {e}")
        return

    # Filter for directories that start with 'Run_'
    run_directories = [item for item in items if item.is_dir() and item.name.startswith('Run_')]

    if not run_directories:
        print(f"No directories starting with 'Run_' found in '{base_path}'.")
        return

    print(f"Found {len(run_directories)} 'Run_' directories in '{base_path}'.")
    print(f"Processing with {max_workers} parallel workers...")

    processed_count = 0
    skipped_count = 0
    
    # Use ProcessPoolExecutor to parallelize the analysis
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all directories for processing
        future_to_dir = {
            executor.submit(process_directory, directory, source_script_path): directory
            for directory in run_directories
        }
        
        # Process results as they complete
        for future in as_completed(future_to_dir):
            dir_name, success, output, error = future.result()
            
            print(f"\n{'=' * 30}")
            print(f"Results for '{dir_name}':")
            
            if success:
                print(f"Successfully finished analysis.")
                processed_count += 1
                
                # Print stdout if available
                if output:
                    print(f"--- Output from {dir_name} ---")
                    print(output.strip())
                    print(f"--- End Output from {dir_name} ---")
            else:
                print(f"Failed to process: {error}")
                skipped_count += 1
            
            # Print stderr if available and not empty
            if error and success:
                print(f"--- Warnings/Errors from {dir_name} ---", file=sys.stderr)
                print(error.strip(), file=sys.stderr)
                print(f"--- End Warnings/Errors from {dir_name} ---", file=sys.stderr)
                
            print(f"{'=' * 30}")

    print("\nParallel analysis run complete.")
    print(f"Successfully processed: {processed_count} directories.")
    print(f"Skipped/Failed: {skipped_count} directories.")

if __name__ == "__main__":
    # Get base path from command line arguments, or use current directory
    base_path_arg = sys.argv[1] if len(sys.argv) > 1 else '.'
    
    # Get max workers from command line arguments, or use None (auto-detect)
    max_workers = None
    if len(sys.argv) > 2:
        try:
            max_workers = int(sys.argv[2])
            if max_workers <= 0:
                print("Warning: Invalid worker count. Using auto-detected CPU count instead.")
                max_workers = None
        except ValueError:
            print("Warning: Invalid worker count. Using auto-detected CPU count instead.")
    
    print(f"Starting parallel analysis run in base directory: '{Path(base_path_arg).resolve()}'")
    if max_workers:
        print(f"Using {max_workers} worker processes")
    else:
        print(f"Using {multiprocessing.cpu_count()} worker processes (auto-detected)")
        
    run_analysis_in_parallel(base_path_arg, max_workers)