import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy import sparse
import traceback

def run_ensemble_calculation(data_dir, cleanup=False, compression_method='compressed_npz', 
                           precision='float32', downsample_factor=2):
    """
    Calculate and save ensemble averages for fluorescence data with optimized file sizes.
    
    Parameters:
    -----------
    data_dir : str
        Name of the data directory
    cleanup : bool
        Whether to remove individual seed files after processing
    compression_method : str
        Method to use for compression ('npz', 'compressed_npz', 'sparse', 'hdf5')
    precision : str
        Floating point precision to use ('float64', 'float32', 'float16')
    downsample_factor : int
        Factor by which to downsample data (1 = no downsampling)
        
    Returns:
    --------
    bool : True if calculation was successful, False otherwise
    """
    try:
        print(f"Using data directory: {data_dir}")
        print(f"Compression method: {compression_method}")
        print(f"Data precision: {precision}")
        print(f"Downsample factor: {downsample_factor}")
        
        # Define the exact models to process based on the filenames you provided
        models = [
            'FL_decay_monomer',         # Basic model
            'FL_decay_monomer_ct',     # RP1 model
            'FL_decay_monomer_ct_RP1',   # RP1_CT model
            'FL_decay_monomer_RP1',   # RP1 model
            'FL_decay_monomer_CT'
        ]
        
        # Find all valid unique seeds
        print("Finding valid seed files...")
        seeds = find_unique_valid_seeds(data_dir)
        
        if not seeds:
            print(f"No valid seed files found in {data_dir} directory!")
            return False
        
        print(f"Found {len(seeds)} valid unique seed files")
        if len(seeds) > 10:
            print(f"First 10 seeds: {seeds[:10]}...")
        else:
            print(f"Seeds: {seeds}")
        
        # Create directory for ensemble-averaged data
        ensemble_dir = os.path.join(data_dir, "Ensemble")
        os.makedirs(ensemble_dir, exist_ok=True)
        
        # Process each model
        for model in models:
            print(f"Processing model: {model}")
            
            # Check if enough seed files exist for this model
            seed_files = []
            for seed in seeds:
                file_path = os.path.join(data_dir, f"{model}_{seed}.npz")
                if os.path.exists(file_path):
                    seed_files.append(file_path)
            
            if len(seed_files) == 0:
                print(f"No seed files found for model {model}. Skipping.")
                continue
                
            print(f"Found {len(seed_files)} seed files for model {model}")
            
            # Load data from all seed files
            try:
                # Load and process all seed files for this model
                data_combined = load_seed_data(model, seeds, data_dir)
                
                if not data_combined['list_fl_t_all']:
                    print(f"No valid data found for model {model}. Skipping.")
                    continue
                
                # Calculate ensemble average
                ensemble_avg = calculate_ensemble_average(data_combined)
                
                # Calculate standard deviation
                ensemble_std = calculate_std_deviation(data_combined)
                
                # Base output filename
                output_file = os.path.join(ensemble_dir, f"{model}_ensemble.npz")
                
                # Save data with optimized file size
                save_optimized_data(
                    output_file,
                    data_combined['list_t'],
                    data_combined['lambda_axis'],
                    ensemble_avg,
                    ensemble_std,
                    data_combined['seeds'],
                    compression=compression_method,
                    precision=precision,
                    downsample_factor=downsample_factor
                )
                
                print(f"Saved ensemble-averaged data for {model} to {output_file}")
                
                # Cleanup original files if requested
                if cleanup:
                    for seed in data_combined['seeds']:
                        seed_file = os.path.join(data_dir, f"{model}_{seed}.npz")
                        if os.path.exists(seed_file):
                            try:
                                os.remove(seed_file)
                                print(f"Removed processed file: {seed_file}")
                            except Exception as e:
                                print(f"Could not remove file {seed_file}: {e}")
                
            except Exception as e:
                print(f"Error processing {model}: {e}")
                traceback.print_exc()
                continue
        
        print("Ensemble average calculation complete.")
        return True
        
    except Exception as e:
        print(f"Error in ensemble calculation: {e}")
        traceback.print_exc()
        return False


def find_unique_valid_seeds(data_dir):
    """
    Find all unique, valid seed files for all three model types.
    
    Parameters:
    -----------
    data_dir : str
        Name of the data directory
        
    Returns:
    --------
    list : Sorted list of valid seed numbers
    """
    valid_seeds = set()
    
    # Check all three specific file patterns
    file_patterns = [
        "FL_decay_monomer_*.npz",
        "FL_decay_monomer_ct_*.npz",
        "FL_decay_monomer_ct_RP1_*.npz",
        "FL_decay_monomer_RP1_*.npz",
        "FL_decay_monomer_CT_*.npz"
    ]
    
    for pattern in file_patterns:
        files = glob.glob(os.path.join(data_dir, pattern))
        
        for file in files:
            basename = os.path.basename(file)
            
            # Skip ensemble files
            if "ensemble" in basename:
                continue
                
            try:
                # Extract seed value based on file pattern
                if basename.startswith("FL_decay_monomer_ct_RP1_"):
                    seed = int(basename.split("_")[-1].split(".")[0])
                elif basename.startswith("FL_decay_monomer_RP1_"):
                    # Make sure we don't capture RP1_CT files here
                    if not basename.startswith("FL_decay_monomer_ct_RP1_"):
                        seed = int(basename.split("_")[-1].split(".")[0])
                    else:
                        continue
                elif basename.startswith("FL_decay_monomer_"):
                    # Make sure we don't capture RP1 or RP1_CT files here
                    if not basename.startswith("FL_decay_monomer_RP1_"):
                        seed = int(basename.split("_")[-1].split(".")[0])
                    else:
                        continue
                else:
                    continue
                
                # Verify this is a valid NPZ file
                try:
                    data = np.load(file)
                    # Check if it contains the expected arrays
                    if 'list_t' in data and 'lambda_axis' in data and 'list_fl_t' in data:
                        valid_seeds.add(seed)
                    data.close()
                except Exception as e:
                    print(f"Skipping corrupt file {file}: {e}")
            except ValueError:
                # Not a properly formatted seed file
                pass
    
    return sorted(list(valid_seeds))


def load_seed_data(model, seeds, data_dir):
    """
    Load data from seed files for a specific model.
    
    Parameters:
    -----------
    model : str
        Model name (e.g., 'FL_decay_monomer')
    seeds : list
        List of seed values to load
    data_dir : str
        Name of the data directory
    
    Returns:
    --------
    dict : Dictionary containing the combined data
    """
    data_combined = {
        'list_t': None,
        'lambda_axis': None,
        'list_fl_t_all': [],
        'seeds': []
    }
    
    for seed in seeds:
        filename = os.path.join(data_dir, f"{model}_{seed}.npz")
        
        if os.path.exists(filename):
            try:
                data = np.load(filename)
                
                # Store time and wavelength axis from the first file
                if data_combined['list_t'] is None:
                    data_combined['list_t'] = data['list_t']
                    data_combined['lambda_axis'] = data['lambda_axis']
                
                # Add fluorescence data
                data_combined['list_fl_t_all'].append(data['list_fl_t'])
                data_combined['seeds'].append(seed)
                data.close()
            except Exception as e:
                print(f"Error loading {filename}: {e}")
    
    print(f"Successfully loaded {len(data_combined['seeds'])} files for {model}")
    
    return data_combined


def calculate_ensemble_average(data_combined):
    """
    Calculate the ensemble average of fluorescence data.
    
    Parameters:
    -----------
    data_combined : dict
        Dictionary with combined data from all seeds
        
    Returns:
    --------
    numpy.ndarray : Ensemble averaged fluorescence data
    """
    if not data_combined['list_fl_t_all']:
        raise ValueError("No data files were found to average")
    
    # Convert list to numpy array
    all_data = np.array(data_combined['list_fl_t_all'])
    
    # Calculate mean along the first axis (across different seeds)
    ensemble_average = np.mean(all_data, axis=0)
    
    return ensemble_average


def calculate_std_deviation(data_combined):
    """
    Calculate the standard deviation of fluorescence data across seeds.
    
    Parameters:
    -----------
    data_combined : dict
        Dictionary with combined data from all seeds
        
    Returns:
    --------
    numpy.ndarray : Standard deviation of fluorescence data
    """
    if not data_combined['list_fl_t_all']:
        raise ValueError("No data files were found to calculate standard deviation")
    
    # Convert list to numpy array
    all_data = np.array(data_combined['list_fl_t_all'])
    
    # Calculate standard deviation along the first axis (across different seeds)
    std_dev = np.std(all_data, axis=0)
    
    return std_dev


def save_optimized_data(output_file, list_t, lambda_axis, ensemble_avg, ensemble_std, seeds, 
                      compression='compressed_npz', precision='float32', downsample_factor=1):
    """
    Save data with optimized file size using compression, precision reduction, and downsampling.
    
    Parameters:
    -----------
    output_file : str
        Path to save the output file
    list_t, lambda_axis : numpy.ndarray
        Time and wavelength axes
    ensemble_avg, ensemble_std : numpy.ndarray
        Ensemble average and standard deviation data
    seeds : list
        List of seed values used
    compression : str
        Compression method ('npz', 'compressed_npz', 'sparse', 'hdf5')
    precision : str
        Data precision ('float64', 'float32', 'float16')
    downsample_factor : int
        Factor by which to downsample the data (1 = no downsampling)
    """
    # Apply downsampling if requested
    if downsample_factor > 1:
        list_t = list_t[::downsample_factor]
        lambda_axis = lambda_axis[::downsample_factor]
        ensemble_avg = ensemble_avg[::downsample_factor, ::downsample_factor]
        ensemble_std = ensemble_std[::downsample_factor, ::downsample_factor]
    
    # Convert to the desired precision
    dtype = np.dtype(precision)
    ensemble_avg = ensemble_avg.astype(dtype)
    ensemble_std = ensemble_std.astype(dtype)
    list_t = list_t.astype(dtype)
    lambda_axis = lambda_axis.astype(dtype)
    
    # Save with the requested compression method
    if compression == 'npz':
        np.savez(
            output_file,
            list_t=list_t,
            lambda_axis=lambda_axis,
            ensemble_avg=ensemble_avg,
            ensemble_std=ensemble_std,
            seeds=seeds
        )
    elif compression == 'compressed_npz':
        np.savez_compressed(
            output_file,
            list_t=list_t,
            lambda_axis=lambda_axis,
            ensemble_avg=ensemble_avg,
            ensemble_std=ensemble_std,
            seeds=seeds
        )
    elif compression == 'sparse':
        # Check if data is sparse enough to benefit from sparse storage
        if np.count_nonzero(ensemble_avg) / ensemble_avg.size < 0.5:
            sparse_avg = sparse.csr_matrix(ensemble_avg)
            sparse_std = sparse.csr_matrix(ensemble_std)
            
            sparse.save_npz(f"{os.path.splitext(output_file)[0]}_avg_sparse.npz", sparse_avg)
            sparse.save_npz(f"{os.path.splitext(output_file)[0]}_std_sparse.npz", sparse_std)
            
            # Save axis data separately
            np.savez_compressed(
                f"{os.path.splitext(output_file)[0]}_axes.npz",
                list_t=list_t,
                lambda_axis=lambda_axis,
                seeds=seeds
            )
        else:
            # Fall back to compressed npz if data isn't sparse
            np.savez_compressed(
                output_file,
                list_t=list_t,
                lambda_axis=lambda_axis,
                ensemble_avg=ensemble_avg,
                ensemble_std=ensemble_std,
                seeds=seeds
            )
    elif compression == 'hdf5':
        try:
            import h5py
            with h5py.File(f"{os.path.splitext(output_file)[0]}.h5", 'w') as f:
                f.create_dataset('list_t', data=list_t, compression='gzip', compression_opts=9)
                f.create_dataset('lambda_axis', data=lambda_axis, compression='gzip', compression_opts=9)
                f.create_dataset('ensemble_avg', data=ensemble_avg, compression='gzip', compression_opts=9)
                f.create_dataset('ensemble_std', data=ensemble_std, compression='gzip', compression_opts=9)
                f.create_dataset('seeds', data=seeds)
        except ImportError:
            print("h5py not available. Falling back to compressed NPZ.")
            np.savez_compressed(
                output_file,
                list_t=list_t,
                lambda_axis=lambda_axis,
                ensemble_avg=ensemble_avg,
                ensemble_std=ensemble_std,
                seeds=seeds
            )


# # If this script is run directly (not imported)
# if __name__ == "__main__":
#     # Example usage
#     data_dir = 'Data_directory'  # Replace with your data directory
#     run_ensemble_calculation(data_dir, cleanup=False)