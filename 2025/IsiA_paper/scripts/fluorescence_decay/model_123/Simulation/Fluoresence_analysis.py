import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson
from scipy.optimize import curve_fit
import os
from pathlib import Path
import argparse
import matplotlib as mpl
from matplotlib.gridspec import GridSpec
import pandas as pd
from sklearn.metrics import r2_score

class FluorescenceData:
    """
    Enhanced fluorescence data analysis for publication-quality figures
    and data export capabilities.
    """
    def __init__(self, main_dir_name, components=None, exp_dir=None, start_time=0.1, params=None):
        """
        Initialize the data processor with specified components from the directory.
        
        Parameters:
        -----------
        main_dir_name : str
            Directory containing the ensemble data files
        components : list or None, optional
            List of components to use. Options: ['monomer', 'ct', 'rp1_ct']
            If None, all available components will be used.
        exp_dir : str, optional
            Directory containing experimental data files
        start_time : float, optional
            Starting time point in nanoseconds (default: 0.1 ns)
        params : dict, optional
            Dictionary of simulation parameters to display in plot titles
            Example: {'TIME_FL': [16], 'TIME_ISC': [5], 'TIME_to_CT': [0], 'TIME_to_RP1': [0.2], 'COUPLING_value': [250]}
        """
        self.main_dir = main_dir_name
        self.exp_dir = exp_dir
        self.start_time = start_time
        self.requested_components = components
        self.params = params or {}  # Store parameters or empty dict if None
        
        # Create necessary directories
        self.output_dir = Path(f'Analysis_Figures/Ensemble/{self.main_dir}')
        self.data_export_dir = Path(f'Analysis_Data/{self.main_dir}')
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.data_export_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize fit parameters tracking
        self.fit_results = []
        
        # Set publication quality plot style
        self.setup_publication_style()
        
        # Discover available components and load data
        self.discover_available_components()
        self.load_component_data()
        
        # Load experimental data if path is provided
        self.exp_data = None
        if exp_dir:
            try:
                exp_path = os.path.join(exp_dir, 'Fl_IsiA_monomer_300k_Gabriela_2025.npy')
                if os.path.exists(exp_path):
                    self.exp_data = np.load(exp_path)
                    print(f"Loaded experimental data from {exp_path}")
                else:
                    print(f"Experimental data file not found at {exp_path}")
            except Exception as e:
                print(f"Error loading experimental data: {e}")
        
        # Convert time to ns for better readability
        self.time_ns = np.array(self.list_t) / 1e6
        
        # Find index for time > start_time ns
        if self.start_time > 0:
            self.start_idx = np.where(self.time_ns > self.start_time)[0][0]
            self.time_ns_trimmed = self.time_ns[self.start_idx:]
            print(f"Starting analysis from t = {self.start_time} ns (skipping earlier time points)")
        else:
            # Include all time points if start_time is 0
            self.start_idx = 0
            self.time_ns_trimmed = self.time_ns
            print("Including all time points in analysis (starting from t = 0 ns)")
        
        # Process basic integrated profiles for available components
        self.calculate_basic_profiles()
        
        # Add experimental values to fit results
        self.add_experimental_values_to_fit_results()
        
        # Print summary of loaded components
        self.print_summary()
        
        # Print parameters if provided
        if self.params:
            print("\nSimulation Parameters:")
            for param, value in self.params.items():
                print(f"  {param}: {value}")

    def _generate_param_display_string(self):
        """
        Generate a formatted string to display simulation parameters in plots.
        
        Returns:
        --------
        str
            Formatted parameter string for display
        """
        if not hasattr(self, 'params') or not self.params:
            return ""
        
        param_strings = []
        
        # Common parameter mapping with prettier names
        param_mapping = {
            'TIME_FL': 'FL decay',
            'TIME_ISC': 'ISC decay', 
            'TIME_to_CT': 'CT decay',
            'TIME_to_RP1': 'RP1 decay',
            'COUPLING_value': 'Coupling'
        }
        
        # Units for parameters
        param_units = {
            'TIME_FL': 'ps',
            'TIME_ISC': 'ps',
            'TIME_to_CT': 'ps',
            'TIME_to_RP1': 'ps',
            'COUPLING_value': 'cm⁻¹'
        }
        
        # Format each parameter with its value and unit
        for param, values in self.params.items():
            if param in param_mapping:
                pretty_name = param_mapping[param]
                unit = param_units.get(param, '')
                
                # Handle both single values and lists
                if isinstance(values, list):
                    if len(values) == 1:
                        value_str = f"{values[0]}"
                    else:
                        value_str = f"{values}"
                else:
                    value_str = f"{values}"
                    
                # Create formatted string
                if unit:
                    param_strings.append(f"{pretty_name}: {value_str} {unit}")
                else:
                    param_strings.append(f"{pretty_name}: {value_str}")
        
        # Combine parameters with commas
        return ", ".join(param_strings)

    def setup_publication_style(self):
        """Set up matplotlib style for publication-quality figures (Nature/Science style)."""
        plt.style.use('default')  # Reset to default style first
        
        # Publication quality settings for Nature/Science style
        mpl.rcParams.update({
            # Font settings
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
            'font.size': 9,
            'mathtext.fontset': 'stixsans',
            
            # Axes settings
            'axes.linewidth': 0.8,
            'axes.labelsize': 9,
            'axes.titlesize': 10,
            'axes.spines.right': True,
            'axes.spines.top': True,
            
            # Tick settings
            'xtick.major.width': 0.8,
            'ytick.major.width': 0.8,
            'xtick.minor.width': 0.6,
            'ytick.minor.width': 0.6,
            'xtick.major.size': 3.0,
            'ytick.major.size': 3.0,
            'xtick.minor.size': 1.5,
            'ytick.minor.size': 1.5,
            'xtick.labelsize': 8,
            'ytick.labelsize': 8,
            'xtick.direction': 'in',
            'ytick.direction': 'in',
            'xtick.major.pad': 3.0,
            'ytick.major.pad': 3.0,
            'xtick.top': True,      # Show ticks on top
            'ytick.right': True,    # Show ticks on right
            
            # Legend settings
            'legend.fontsize': 8,
            'legend.frameon': False,
            'legend.edgecolor': '0.8',
            'legend.borderpad': 0.2,
            'legend.labelspacing': 0.2,
            'legend.handlelength': 1.0,
            'legend.handletextpad': 0.5,
            
            # Figure settings
            'figure.figsize': (3.5, 2.625),  # Nature single column width
            'figure.dpi': 300,
            'savefig.dpi': 600,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.02,
            # 'savefig.format': 'tiff',  # High-quality TIFF format for publication
            
            # Grid and other settings
            'grid.alpha': 0.3,
            'lines.linewidth': 1.0,
            'lines.markersize': 3.0
        })
    
    def add_experimental_values_to_fit_results(self):
        """Add the experimental biexponential fit values to the results list"""
        # Experimental values from Gabriela's data
        self.fit_results.append({
            'name': 'Experiment',
            'monomer_pct': None,  # Not applicable for experimental data
            'a1': 23.7,           # 23.7%
            'tau1': 0.21,         # 0.21 ns
            'a2': 76.3,           # 76.3%
            'tau2': 3.89,         # 3.89 ns
            'tau_avg': 3.02,      # 3.02 ns
            'r2': 1.0             # Perfect fit (reference)
        })
        
    def discover_available_components(self):
        """Discover which component data files are available in the directory."""
        ensemble_dir = os.path.join(self.main_dir, 'Ensemble')
        
        # Define file patterns for each component
        component_files = {
            'monomer': 'FL_decay_monomer_ensemble.npz',
            'ct': 'FL_decay_monomer_ct_ensemble.npz',
            'CT' : 'FL_decay_monomer_CT_ensemble.npz',
            'ct_rp1': 'FL_decay_monomer_ct_RP1_ensemble.npz',
            'rp1': 'FL_decay_monomer_RP1_ensemble.npz'
        }
        
        # Check which files exist
        self.available_components = {}
        for component, filename in component_files.items():
            filepath = os.path.join(ensemble_dir, filename)
            self.available_components[component] = os.path.exists(filepath)
        
        # If specific components were requested, check they're available
        if self.requested_components:
            for component in self.requested_components:
                if component not in self.available_components:
                    raise ValueError(f"Requested component '{component}' is not recognized. Valid options: monomer, ct, rp1_ct")
                if not self.available_components[component]:
                    print(f"Warning: Requested component '{component}' data file not found.")
            
            # Use only the requested components that are available
            self.active_components = {comp: self.available_components[comp] 
                                     for comp in self.requested_components 
                                     if self.available_components[comp]}
        else:
            # Use all available components
            self.active_components = self.available_components
        
        # Ensure we have at least one component
        if not any(self.active_components.values()):
            raise FileNotFoundError("No component data files found. Ensure your data directory contains at least one ensemble NPZ file.")
            
        print('='*50)
        print("Available components:", ", ".join([comp for comp, available in self.available_components.items() if available]))
        print("Active components:", ", ".join([comp for comp, active in self.active_components.items() if active]))
        print('='*50)

    def load_component_data(self):
        """Load data for active components from NPZ files."""
        ensemble_dir = os.path.join(self.main_dir, 'Ensemble')
        
        # Initialize component data with placeholders
        self.list_fl_t_monomer = None
        self.list_fl_t_ct = None
        self.list_fl_t_ct_rp1 = None
        self.list_fl_t_rp1 = None
        
        # Data time and wavelength arrays (should be the same for all files)
        self.list_t = None
        self.lambda_axis = None
        
        # Load each component's data if active
        for component, active in self.active_components.items():
            if not active:
                continue
                
            if component == 'monomer':
                file_path = os.path.join(ensemble_dir, 'FL_decay_monomer_ensemble.npz')
                data = np.load(file_path)
                self.list_fl_t_monomer = data['ensemble_avg']
                # Set time and wavelength arrays if not set yet
                if self.list_t is None:
                    self.list_t = data['list_t']
                    self.lambda_axis = data['lambda_axis']
                print(f"Loaded monomer data from {file_path}")
                
            elif component == 'ct' or component == 'CT':
                # Try lowercase version first
                file_path = os.path.join(ensemble_dir, 'FL_decay_monomer_ct_ensemble.npz')
                if not os.path.exists(file_path):
                    # Try uppercase version if lowercase doesn't exist
                    file_path = os.path.join(ensemble_dir, 'FL_decay_monomer_CT_ensemble.npz')
                
                if os.path.exists(file_path):
                    data = np.load(file_path)
                    self.list_fl_t_ct = data['ensemble_avg']
                    # Set time and wavelength arrays if not set yet
                    if self.list_t is None:
                        self.list_t = data['list_t']
                        self.lambda_axis = data['lambda_axis']
                    print(f"Loaded CT data from {file_path}")
                else:
                    print("Warning: Neither FL_decay_monomer_ct_ensemble.npz nor FL_decay_monomer_CT_ensemble.npz found")
                
            elif component == 'ct_rp1':
                file_path = os.path.join(ensemble_dir, 'FL_decay_monomer_ct_RP1_ensemble.npz')
                data = np.load(file_path)
                self.list_fl_t_ct_rp1 = data['ensemble_avg']
                # Set time and wavelength arrays if not set yet
                if self.list_t is None:
                    self.list_t = data['list_t']
                    self.lambda_axis = data['lambda_axis']
                print(f"Loaded RP1+CT data from {file_path}")
            
            elif component == 'rp1':
                file_path = os.path.join(ensemble_dir, 'FL_decay_monomer_RP1_ensemble.npz')
                data = np.load(file_path)
                self.list_fl_t_rp1 = data['ensemble_avg']
                # Set time and wavelength arrays if not set yet
                if self.list_t is None:
                    self.list_t = data['list_t']
                    self.lambda_axis = data['lambda_axis']
                print(f"Loaded RP1 data from {file_path}")
        
        # Create placeholder data for inactive components
        shape = None
        if self.list_fl_t_monomer is not None:
            shape = self.list_fl_t_monomer.shape
        elif self.list_fl_t_ct is not None:
            shape = self.list_fl_t_ct.shape
        elif self.list_fl_t_ct_rp1 is not None:
            shape = self.list_fl_t_ct_rp1.shape
        elif self.list_fl_t_rp1 is not None:
            shape = self.list_fl_t_rp1.shape
        
        if shape is not None:
            if self.list_fl_t_monomer is None:
                self.list_fl_t_monomer = np.zeros(shape)
                print("Created placeholder monomer data (zeros)")
            if self.list_fl_t_ct is None:
                self.list_fl_t_ct = np.zeros(shape)
                print("Created placeholder CT data (zeros)")
            if self.list_fl_t_ct_rp1 is None:
                self.list_fl_t_ct_rp1 = np.zeros(shape)
                print("Created placeholder RP1+CT data (zeros)")
            if self.list_fl_t_rp1 is None:
                self.list_fl_t_rp1 = np.zeros(shape)
                print("Created placeholder RP1 data (zeros)")
        else:
            raise ValueError("No valid data loaded for any component.")

    
    def print_summary(self):
        """Print a summary of the loaded data."""
        print("\n" + "="*50)
        print(f"Data Summary for {self.main_dir}")
        print("="*50)
        print("Available components:")
        for comp, available in self.available_components.items():
            status = "Available" if available else "Not found"
            active = " (Active)" if comp in self.active_components and self.active_components[comp] else ""
            print(f"  - {comp}: {status}{active}")
        
        print(f"Time points: {len(self.list_t)}")
        print(f"Wavelength points: {len(self.lambda_axis)}")
        if self.list_fl_t_monomer is not None:
            print(f"Monomer data shape: {self.list_fl_t_monomer.shape}")
        if self.list_fl_t_ct is not None:
            print(f"CT data shape: {self.list_fl_t_ct.shape}")
        if self.list_fl_t_ct_rp1 is not None:
            print(f"RP1+CT data shape: {self.list_fl_t_ct_rp1.shape}")
        if self.list_fl_t_rp1 is not None:
            print(f"RP1 data shape: {self.list_fl_t_rp1.shape}")
        print(f"Experimental data: {'Loaded' if self.exp_data is not None else 'Not available'}")
        print("="*50 + "\n")
    
    def integrate_over_wavelength(self, fl_data):
        """Integrate fluorescence data over wavelength axis."""
        # Assuming fl_data has shape (time_points, wavelength_points)
        fl_integrated = np.zeros(fl_data.shape[0])
        
        # Integrate each time point's spectrum over wavelength
        for i in range(fl_data.shape[0]):
            fl_integrated[i] = -simpson(fl_data[i, :], self.lambda_axis)
        
        return fl_integrated
    
    def integrate_over_time(self, fl_data):
        """Integrate fluorescence data over time axis."""
        # Assuming fl_data has shape (time_points, wavelength_points)
        fl_integrated = np.zeros(fl_data.shape[1])
        
        # Integrate each wavelength's time trace over time
        for i in range(fl_data.shape[1]):
            fl_integrated[i] = simpson(fl_data[:, i], self.list_t)
        
        return fl_integrated
    
    def normalize_fl_data(self, data):
        """Normalize 2D fluorescence data."""
        # Take absolute values to ensure positivity
        abs_data = np.abs(data)
        # Normalize to the maximum value
        return abs_data / np.max(abs_data) if np.max(abs_data) > 0 else abs_data
    
    def calculate_basic_profiles(self):
        """Calculate basic integrated profiles for active components."""
        # Time-resolved profiles for all components
        self.fl_time_monomer = self.integrate_over_wavelength(self.list_fl_t_monomer)[self.start_idx:]
        self.fl_time_ct = self.integrate_over_wavelength(self.list_fl_t_ct)[self.start_idx:]
        self.fl_time_rp1_ct = self.integrate_over_wavelength(self.list_fl_t_ct_rp1)[self.start_idx:]
        self.fl_time_rp1 = self.integrate_over_wavelength(self.list_fl_t_rp1)[self.start_idx:]
        
        # Wavelength-resolved profiles for all components
        self.fl_wave_monomer = self.integrate_over_time(self.list_fl_t_monomer)
        self.fl_wave_ct = self.integrate_over_time(self.list_fl_t_ct)
        self.fl_wave_rp1_ct = self.integrate_over_time(self.list_fl_t_ct_rp1)
        self.fl_wave_rp1 = self.integrate_over_time(self.list_fl_t_rp1)
        
        # Normalize wavelength profiles
        self.norm_monomer = self.fl_wave_monomer / np.max(self.fl_wave_monomer) if np.max(self.fl_wave_monomer) > 0 else np.zeros_like(self.fl_wave_monomer)
        self.norm_ct = self.fl_wave_ct / np.max(self.fl_wave_ct) if np.max(self.fl_wave_ct) > 0 else np.zeros_like(self.fl_wave_ct)
        self.norm_rp1_ct = self.fl_wave_rp1_ct / np.max(self.fl_wave_rp1_ct) if np.max(self.fl_wave_rp1_ct) > 0 else np.zeros_like(self.fl_wave_rp1_ct)
        self.norm_rp1 = self.fl_wave_rp1 / np.max(self.fl_wave_rp1) if np.max(self.fl_wave_rp1) > 0 else np.zeros_like(self.fl_wave_rp1)

    def get_exp_biexponential_decay(self, time_array):
        """
        Generate experimental bi-exponential decay curve.
        
        Parameters:
        -----------
        time_array : ndarray
            Time points to evaluate the decay at
        
        Returns:
        --------
        ndarray
            Bi-exponential decay curve
        """
        # Parameters from the table
        a1, tau1 = 0.237, 0.21  # First component amplitude and lifetime
        a2, tau2 = 0.763, 3.89  # Second component amplitude and lifetime
        
        # Ensure the time array starts from 0 for the decay calculation
        time_offset = time_array[0]
        adj_time = time_array - time_offset
        
        # Calculate decay
        decay = a1 * np.exp(-adj_time / tau1) + a2 * np.exp(-adj_time / tau2)
        
        # Normalize to start at 1.0
        return decay / decay[0] if len(decay) > 0 else decay

    def biexponential_fit(self, profile, monomer_pct=None, label=None):
        """
        Fit a profile with a bi-exponential model and calculate R².
        
        Parameters:
        -----------
        profile : ndarray
            Normalized fluorescence profile to fit
        monomer_pct : float, optional
            Percentage of monomer in the mixture
        label : str, optional
            Label for the data being fit
        
        Returns:
        --------
        tuple
            (a1, tau1, a2, tau2, tau_avg, r2) - fit parameters and R²
        """
        # Define bi-exponential model function
        def biexponential(t, a1, tau1, a2, tau2):
            return a1 * np.exp(-t / tau1) + a2 * np.exp(-t / tau2)
        
        # Time for fitting (starting from 0)
        time_for_fit = self.time_ns_trimmed - self.time_ns_trimmed[0]
        
        # Initial parameter guess: [a1, tau1, a2, tau2]
        p0 = [0.5, 0.2, 0.5, 2.0]
        bounds = ([0.01, 0.01, 0.01, 0.1], [0.99, 1.0, 0.99, 10.0])
        
        try:
            # Fit the profile with a biexponential model
            params, _ = curve_fit(biexponential, time_for_fit, profile, 
                                p0=p0, bounds=bounds)
            
            # Normalize amplitudes to sum to 1
            a1, tau1, a2, tau2 = params
            sum_a = a1 + a2
            a1 /= sum_a
            a2 /= sum_a
            
            # Calculate average lifetime
            tau_avg = (a1 * tau1 + a2 * tau2) / (a1 + a2)
            
            # Calculate the predicted values for R-squared
            y_pred = biexponential(time_for_fit, a1, tau1, a2, tau2)
            r2 = r2_score(profile, y_pred)
            
            # Store the fit results if monomer_pct is provided
            if monomer_pct is not None:
                self.fit_results.append({
                    'name': label if label else f"Monomer {monomer_pct}%",
                    'monomer_pct': monomer_pct,
                    'a1': a1 * 100,  # Convert to percentage
                    'tau1': tau1,
                    'a2': a2 * 100,  # Convert to percentage
                    'tau2': tau2,
                    'tau_avg': tau_avg,
                    'r2': r2
                })
            
            return a1, tau1, a2, tau2, tau_avg, r2
        
        except Exception as e:
            print(f"Error in biexponential fit: {e}")
            
            # Add failed fit with None values if monomer_pct is provided
            if monomer_pct is not None:
                self.fit_results.append({
                    'name': label if label else f"Monomer {monomer_pct}%",
                    'monomer_pct': monomer_pct,
                    'a1': None,
                    'tau1': None,
                    'a2': None,
                    'tau2': None,
                    'tau_avg': None,
                    'r2': None
                })
            
            return None
    
    def save_fit_results_to_csv(self):
        """
        Save the fit results to a CSV file.
        
        Returns:
        --------
        str
            Path to the saved CSV file
        """
        if not self.fit_results:
            print("No fit results to save.")
            return None
            
        # Convert to DataFrame
        df = pd.DataFrame(self.fit_results)
        
        # Ensure output directory exists
        os.makedirs(self.data_export_dir, exist_ok=True)
        
        # Save to CSV
        csv_path = os.path.join(self.data_export_dir, 'biexponential_fit_results.csv')
        df.to_csv(csv_path, index=False)
        
        print(f"Saved fit results to: {csv_path}")
        
        # Also create a formatted table version
        formatted_path = os.path.join(self.data_export_dir, 'biexponential_fit_table.csv')
        
        # Format the data for publication
        formatted_data = []
        for row in self.fit_results:
            formatted_row = {}
            formatted_row['Name'] = row['name']
            formatted_row['Monomer %'] = row['monomer_pct'] if row['monomer_pct'] is not None else '-'
            
            # Format parameters with appropriate precision
            if row['a1'] is not None:
                formatted_row['a1 (%)'] = f"{row['a1']:.1f}"
                formatted_row['τ1 (ns)'] = f"{row['tau1']:.2f}"
                formatted_row['a2 (%)'] = f"{row['a2']:.1f}"
                formatted_row['τ2 (ns)'] = f"{row['tau2']:.2f}"
                formatted_row['τavg (ns)'] = f"{row['tau_avg']:.2f}"
                formatted_row['R²'] = f"{row['r2']:.3f}"
            else:
                formatted_row['a1 (%)'] = '-'
                formatted_row['τ1 (ns)'] = '-'
                formatted_row['a2 (%)'] = '-'
                formatted_row['τ2 (ns)'] = '-'
                formatted_row['τavg (ns)'] = '-'
                formatted_row['R²'] = '-'
                
            formatted_data.append(formatted_row)
            
        # Convert to DataFrame and save
        formatted_df = pd.DataFrame(formatted_data)
        formatted_df.to_csv(formatted_path, index=False)
        
        print(f"Saved formatted fit results to: {formatted_path}")
        
        return csv_path

    def plot_multiple_ratios(self, ratio_dict, save=True, show=True):
        """
        Plot multiple component mixtures on the same figure for comparison.
        Ensures each mixture is saved with a unique filename.
        
        Parameters:
        -----------
        ratio_dict : dict
            Dictionary with mixture names as keys and tuples of 
            (monomer_pct, ct_pct, rp1_pct, rp1_only_pct) as values.
        save : bool, optional
            Whether to save the figure and data
        show : bool, optional
            Whether to show the figure
                    
        Returns:
        --------
        tuple
            (fig_time, fig_wave) - Time-resolved and wavelength-resolved figures
        """
        # Create separate figures for time-resolved and wavelength-resolved plots
        fig_time = plt.figure(figsize=(7, 5))
        ax_time = fig_time.add_subplot(111)
        
        fig_wave = plt.figure(figsize=(7, 5))
        ax_wave = fig_wave.add_subplot(111)
        
        # Store data for all ratios
        time_data = {}
        wave_data = {}
        
        # Define a nice color gradient for different ratios - using viridis colormap
        cmap = plt.cm.viridis
        colors = cmap(np.linspace(0, 0.9, len(ratio_dict)))
        
        # Add experimental data to both plots
        # For time-resolved plot
        exp_decay = self.get_exp_biexponential_decay(self.time_ns_trimmed)
        ax_time.plot(self.time_ns_trimmed, exp_decay, 'k--', 
                linewidth=1.2, label='Experiment (Bi-exp)')
        
        # For wavelength-resolved plot
        if self.exp_data is not None:
            # Sort experimental data
            exp_x = self.exp_data[:, 0]
            exp_y = self.exp_data[:, 1]
            exp_sort_indices = np.argsort(exp_x)
            
            sorted_exp_x = exp_x[exp_sort_indices]
            sorted_exp_y = exp_y[exp_sort_indices]
            
            # Normalize
            exp_data_norm = sorted_exp_y / np.max(sorted_exp_y)
            
            ax_wave.plot(sorted_exp_x, exp_data_norm, 'k--', 
                    linewidth=1.2, label='Experiment')
        
        # Process each ratio
        for i, (name, ratio_values) in enumerate(ratio_dict.items()):
            # Handle both 3-value and 4-value tuples
            if len(ratio_values) == 4:
                monomer_pct, ct_pct, rp1_pct, rp1_only_pct = ratio_values
            else:
                monomer_pct, ct_pct, rp1_pct = ratio_values
                rp1_only_pct = 0
            
            # Calculate the mixture
            x = monomer_pct / 100     # Monomer fraction
            y = ct_pct / 100          # CT fraction
            z = rp1_pct / 100         # RP1+CT fraction
            w = rp1_only_pct / 100    # RP1 only fraction
            
            # Use the name directly from the dictionary as the label
            # This ensures consistency with the user-specified format
            label = name
            
            print(f"Processing ratio {name}: {monomer_pct}% M, {ct_pct}% CT, {rp1_pct}% RP1+CT, {rp1_only_pct}% RP1")
            
            # Create the mixture
            list_fl_t_mixture = np.zeros_like(self.list_fl_t_monomer)
        
            # Add monomer component if available and requested
            if x > 0 and np.max(self.list_fl_t_monomer) > 0:
                list_fl_t_mixture += self.list_fl_t_monomer * x
                
            # Add CT component if available and requested
            if y > 0 and np.max(self.list_fl_t_ct) > 0:
                list_fl_t_mixture += self.list_fl_t_ct * y
                
            # Add RP1+CT component if available and requested
            if z > 0 and np.max(self.list_fl_t_ct_rp1) > 0:
                list_fl_t_mixture += self.list_fl_t_ct_rp1 * z
                
            # Add RP1 only component if available and requested
            if w > 0 and np.max(self.list_fl_t_rp1) > 0:
                list_fl_t_mixture += self.list_fl_t_rp1 * w
            
            # Calculate time-resolved and wavelength-resolved profiles
            fl_time_mixture = self.integrate_over_wavelength(list_fl_t_mixture)[self.start_idx:]
            fl_wave_mixture = self.integrate_over_time(list_fl_t_mixture)
            
            # Store data
            time_data[name] = {'profile': fl_time_mixture, 'label': label}
            wave_data[name] = {'profile': fl_wave_mixture, 'label': label}
            
            # Plot time-resolved data
            if np.max(fl_time_mixture) > 0:
                # Normalize to start at 1.0
                fl_time_norm = fl_time_mixture / fl_time_mixture[0]
                ax_time.plot(self.time_ns_trimmed, fl_time_norm, '-', 
                        color=colors[i], linewidth=1.5, label=label)
                
                # Fit with biexponential model and store results
                if hasattr(self, 'biexponential_fit'):
                    self.biexponential_fit(fl_time_norm, monomer_pct, name)
            
            # Plot wavelength-resolved data
            if np.max(fl_wave_mixture) > 0:
                # Normalize by max value
                fl_wave_norm = fl_wave_mixture / np.max(fl_wave_mixture)
                
                # Sort the wavelength data
                sort_indices = np.argsort(self.lambda_axis)
                sorted_lambda = self.lambda_axis[sort_indices]
                sorted_fl_wave = fl_wave_norm[sort_indices]
                
                # Create a mask for near-zero values
                threshold = 1e-6
                valid_indices = sorted_fl_wave > threshold
                
                # Plot if there are valid points
                if np.any(valid_indices):
                    ax_wave.plot(sorted_lambda[valid_indices], sorted_fl_wave[valid_indices], '-', 
                            color=colors[i], linewidth=1.5, label=label)
            
            # Save individual data if requested
            if save:
                if hasattr(self, 'save_simulation_data'):
                    # Important: Save each mixture with a unique label to avoid overwriting
                    self.save_simulation_data(fl_time_mixture, fl_wave_mixture, label)
        
        # Configure time-resolved plot
        ax_time.set_xlabel('Time (ns)')
        ax_time.set_ylabel('Normalized Intensity')
        
        # Generate parameter display string for the title
        param_str = ""
        if hasattr(self, 'params') and self.params:
            param_parts = []
            # Format each parameter with its value and unit
            if 'TIME_FL' in self.params:
                param_parts.append(f"FL decay: {self.params['TIME_FL']} ps")
            if 'TIME_ISC' in self.params:
                param_parts.append(f"ISC decay: {self.params['TIME_ISC']} ps")
            if 'TIME_to_CT' in self.params:
                param_parts.append(f"CT decay: {self.params['TIME_to_CT']} ps")
            if 'TIME_to_RP1' in self.params:
                param_parts.append(f"RP1 decay: {self.params['TIME_to_RP1']} ps")
            if 'COUPLING_value' in self.params:
                param_parts.append(f"Coupling: {self.params['COUPLING_value']} cm⁻¹")
            
            param_str = ", ".join(param_parts)
        
        # Set the titles using the parameter string instead of default titles
        if param_str:
            ax_time.set_title("Time-Resolved Fluorescence Comparison")
            fig_time.suptitle(param_str, fontsize=9, y=0.99)
            
            ax_wave.set_title("Wavelength-Resolved Fluorescence Comparison")
            fig_wave.suptitle(param_str, fontsize=9, y=0.99)
        else:
            ax_time.set_title('Time-Resolved Fluorescence Comparison')
            ax_wave.set_title('Wavelength-Resolved Fluorescence Comparison')
        
        ax_time.set_xlim(self.time_ns_trimmed[0], min(8, max(self.time_ns_trimmed)))
        ax_time.set_ylim(0, 1.05)
        ax_time.minorticks_on()
        ax_time.grid(True, alpha=0.3)
        
        # Configure wavelength-resolved plot
        ax_wave.set_xlabel('Wavelength (nm)')
        ax_wave.set_ylabel('Normalized Intensity')
        
        # Focus on the region of interest (660-740 nm)
        ax_wave.set_xlim(660, 740)
        ax_wave.set_ylim(0, 1.05)
        ax_wave.minorticks_on()
        ax_wave.grid(True, alpha=0.3)
        
        # Add legends with better placement
        ax_time.legend(loc='best', fontsize=8, frameon=False)
        ax_wave.legend(loc='best', fontsize=8, frameon=False)
        
        # Adjust layout
        fig_time.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for suptitle
        fig_wave.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for suptitle
        
        # Save the figures
        if save:
            # Generate parameter suffix for filenames
            param_suffix = ""
            if hasattr(self, 'params') and self.params:
                fl_val = str(self.params.get('TIME_FL', ['na']))
                isc_val = str(self.params.get('TIME_ISC', ['na']))
                rp_val = str(self.params.get('TIME_to_RP1', ['na']))
                coupling = str(self.params.get('COUPLING_value', ['na']))
                param_suffix = f"_FL{fl_val}_ISC{isc_val}_RP{rp_val}_Coup{coupling}_start_time_{self.start_time}"
            
            # Time-resolved figure
            time_png_path = os.path.join(self.output_dir, f"multiple_ratios_time_comparison{param_suffix}.png")
            
            # Wavelength-resolved figure
            wave_png_path = os.path.join(self.output_dir, f"multiple_ratios_wave_comparison{param_suffix}.png")
            
            # Save as PNG
            fig_time.savefig(time_png_path, dpi=300, format='png', bbox_inches='tight')
            fig_wave.savefig(wave_png_path, dpi=300, format='png', bbox_inches='tight')
            
            print(f"Saved time-resolved comparison to {time_png_path}")
            print(f"Saved wavelength-resolved comparison to {wave_png_path}")
            
            # Save the fit results to CSV file if the method exists
            if hasattr(self, 'save_fit_results_to_csv'):
                self.save_fit_results_to_csv()
        
        if show:
            plt.show()
        else:
            plt.close(fig_time)
            plt.close(fig_wave)
        
        return fig_time, fig_wave

    def add_experimental_values_to_fit_results(self):
        """Add the experimental biexponential fit values to the results list"""
        # Initialize fit_results list if it doesn't exist
        if not hasattr(self, 'fit_results'):
            self.fit_results = []
        
        # Experimental values from Gabriela's data
        self.fit_results.append({
            'name': 'Experiment',
            'monomer_pct': None,  # Not applicable for experimental data
            'a1': 23.7,           # 23.7%
            'tau1': 0.21,         # 0.21 ns
            'a2': 76.3,           # 76.3%
            'tau2': 3.89,         # 3.89 ns
            'tau_avg': 3.02,      # 3.02 ns
            'r2': 1.0             # Perfect fit (reference)
        })

    def biexponential_fit(self, profile, monomer_pct=None, label=None):
        """
        Fit a profile with a bi-exponential model and calculate R².
        
        Parameters:
        -----------
        profile : ndarray
            Normalized fluorescence profile to fit
        monomer_pct : float, optional
            Percentage of monomer in the mixture
        label : str, optional
            Label for the data being fit
        
        Returns:
        --------
        tuple
            (a1, tau1, a2, tau2, tau_avg, r2) - fit parameters and R²
        """
        # Initialize fit_results list if it doesn't exist
        if not hasattr(self, 'fit_results'):
            self.fit_results = []
            # Add experimental values
            self.add_experimental_values_to_fit_results()
        
        # Import required libraries
        from scipy.optimize import curve_fit
        from sklearn.metrics import r2_score
        import numpy as np
        
        # Define bi-exponential model function
        def biexponential(t, a1, tau1, a2, tau2):
            return a1 * np.exp(-t / tau1) + a2 * np.exp(-t / tau2)
        
        # Time for fitting (starting from 0)
        time_for_fit = self.time_ns_trimmed - self.time_ns_trimmed[0]
        
        # Initial parameter guess: [a1, tau1, a2, tau2]
        p0 = [0.5, 0.2, 0.5, 2.0]
        bounds = ([0.01, 0.01, 0.01, 0.1], [0.99, 1.0, 0.99, 10.0])
        
        try:
            # Fit the profile with a biexponential model
            params, _ = curve_fit(biexponential, time_for_fit, profile, 
                                p0=p0, bounds=bounds)
            
            # Normalize amplitudes to sum to 1
            a1, tau1, a2, tau2 = params
            sum_a = a1 + a2
            a1 /= sum_a
            a2 /= sum_a
            
            # Calculate average lifetime
            tau_avg = (a1 * tau1 + a2 * tau2) / (a1 + a2)
            
            # Calculate the predicted values for R-squared
            y_pred = biexponential(time_for_fit, a1, tau1, a2, tau2)
            r2 = r2_score(profile, y_pred)
            
            # Store the fit results if monomer_pct is provided
            if monomer_pct is not None:
                self.fit_results.append({
                    'name': label if label else f"Monomer {monomer_pct}%",
                    'monomer_pct': monomer_pct,
                    'a1': a1 * 100,  # Convert to percentage
                    'tau1': tau1,
                    'a2': a2 * 100,  # Convert to percentage
                    'tau2': tau2,
                    'tau_avg': tau_avg,
                    'r2': r2
                })
            
            return a1, tau1, a2, tau2, tau_avg, r2
        
        except Exception as e:
            print(f"Error in biexponential fit: {e}")
            
            # Add failed fit with None values if monomer_pct is provided
            if monomer_pct is not None:
                self.fit_results.append({
                    'name': label if label else f"Monomer {monomer_pct}%",
                    'monomer_pct': monomer_pct,
                    'a1': None,
                    'tau1': None,
                    'a2': None,
                    'tau2': None,
                    'tau_avg': None,
                    'r2': None
                })
            
            return None

    def save_fit_results_to_csv(self):
        """
        Save the fit results to a CSV file.
        
        Returns:
        --------
        str
            Path to the saved CSV file
        """
        import pandas as pd
        import os
        
        if not hasattr(self, 'fit_results') or not self.fit_results:
            print("No fit results to save.")
            return None
            
        # Convert to DataFrame
        df = pd.DataFrame(self.fit_results)
        
        # Ensure output directory exists
        os.makedirs(self.data_export_dir, exist_ok=True)
        
        # Save to CSV
        csv_path = os.path.join(self.data_export_dir, 'biexponential_fit_results.csv')
        df.to_csv(csv_path, index=False)
        
        print(f"Saved fit results to: {csv_path}")
        
        # Also create a formatted table version
        formatted_path = os.path.join(self.data_export_dir, 'biexponential_fit_table.csv')
        
        # Format the data for publication
        formatted_data = []
        for row in self.fit_results:
            formatted_row = {}
            formatted_row['Name'] = row['name']
            formatted_row['Monomer %'] = row['monomer_pct'] if row['monomer_pct'] is not None else '-'
            
            # Format parameters with appropriate precision
            if row['a1'] is not None:
                formatted_row['a1 (%)'] = f"{row['a1']:.1f}"
                formatted_row['τ1 (ns)'] = f"{row['tau1']:.2f}"
                formatted_row['a2 (%)'] = f"{row['a2']:.1f}"
                formatted_row['τ2 (ns)'] = f"{row['tau2']:.2f}"
                formatted_row['τavg (ns)'] = f"{row['tau_avg']:.2f}"
                formatted_row['R²'] = f"{row['r2']:.3f}"
            else:
                formatted_row['a1 (%)'] = '-'
                formatted_row['τ1 (ns)'] = '-'
                formatted_row['a2 (%)'] = '-'
                formatted_row['τ2 (ns)'] = '-'
                formatted_row['τavg (ns)'] = '-'
                formatted_row['R²'] = '-'
                
            formatted_data.append(formatted_row)
            
        # Convert to DataFrame and save
        formatted_df = pd.DataFrame(formatted_data)
        formatted_df.to_csv(formatted_path, index=False)
        
        print(f"Saved formatted fit results to: {formatted_path}")
        
        return csv_path

    def save_simulation_data(self, fl_time_mixture, fl_wave_mixture, label):
        """
        Save simulation data to NPY files with unique names for each mixture.
        
        Parameters:
        -----------
        fl_time_mixture : ndarray
            Time-resolved fluorescence data
        fl_wave_mixture : ndarray
            Wavelength-resolved fluorescence data
        label : str
            Descriptive label for the mixture
        """
        import numpy as np
        import os
        import re
        
        # Extract percentages using regex - match "M 0%" format
        monomer_match = re.search(r'M\s+(\d+)%', label)
        ct_match = re.search(r'CT\s+(\d+)%', label)
        rp1_match = re.search(r'RP1\s+(\d+)%', label)
        rp1_ct_match = re.search(r'RP1\+CT\s+(\d+)%', label)
        
        # Set component values (default to 0 if not found)
        monomer_pct = int(monomer_match.group(1)) if monomer_match else 0
        ct_pct = int(ct_match.group(1)) if ct_match else 0
        rp1_pct = int(rp1_match.group(1)) if rp1_match else 0
        rp1_ct_pct = int(rp1_ct_match.group(1)) if rp1_ct_match else 0
        
        # Generate a unique filename that explicitly includes all four component percentages
        # All four components should sum to 100%
        filename = f"Monomer_{monomer_pct}pct_CT_{ct_pct}pct_RP1_{rp1_pct}pct_CTRP1_{rp1_ct_pct}pct"

        # # Add RP components if they're used
        # if rp1_ct_pct > 0:
        #     filename += f"_{rp1_ct_pct}pct_RP1CT"
        # if rp1_pct > 0:
        #     filename += f"_{rp1_pct}pct_RP1"
        
        # Add start time suffix
        start_time_suffix = f"_start_time_{self.start_time}"
        
        # Prepare time-resolved data with time axis
        time_data = np.column_stack((self.time_ns_trimmed, fl_time_mixture))
        
        # Prepare wavelength-resolved data with wavelength axis
        wave_data = np.column_stack((self.lambda_axis, fl_wave_mixture))
        
        # Generate unique filenames
        time_file = os.path.join(self.data_export_dir, f"time_resolved_{filename}{start_time_suffix}.npy")
        wave_file = os.path.join(self.data_export_dir, f"wavelength_resolved_{filename}{start_time_suffix}.npy")
        
        # Check if the files already exist to avoid overwriting
        if os.path.exists(time_file):
            print(f"WARNING: File {time_file} already exists, adding timestamp to avoid overwriting")
            import time
            timestamp = int(time.time())
            time_file = os.path.join(self.data_export_dir, f"time_resolved_{filename}_{timestamp}{start_time_suffix}.npy")
            wave_file = os.path.join(self.data_export_dir, f"wavelength_resolved_{filename}_{timestamp}{start_time_suffix}.npy")
        
        # Save the files
        np.save(time_file, time_data)
        np.save(wave_file, wave_data)
        
        print(f"Saved time-resolved data to {time_file}")
        print(f"Saved wavelength-resolved data to {wave_file}")
        
    def create_specific_mixture(self, monomer_pct, ct_pct, rp1_pct, rp1_only_pct=0, save=True, show=True):
        """
        Create and plot a specific mixture of components with consistent labeling.
        
        Parameters:
        -----------
        monomer_pct : float
            Percentage of monomer (0-100)
        ct_pct : float
            Percentage of CT (0-100)
        rp1_pct : float
            Percentage of RP1+CT (0-100)
        rp1_only_pct : float, optional
            Percentage of RP1 only (0-100), defaults to 0
        save : bool, optional
            Whether to save the figure and data
        show : bool, optional
            Whether to show the figure
        """
        # Verify percentages sum to 100
        total = monomer_pct + ct_pct + rp1_pct + rp1_only_pct
        if abs(total - 100) > 0.001:
            print(f"Warning: Percentages must sum to 100. Current sum: {total}")
            # Normalize to 100%
            factor = 100 / total
            monomer_pct *= factor
            ct_pct *= factor
            rp1_pct *= factor
            rp1_only_pct *= factor
            print(f"Adjusted to: Monomer={monomer_pct:.1f}%, CT={ct_pct:.1f}%, RP1+CT={rp1_pct:.1f}%, RP1={rp1_only_pct:.1f}%")
        
        # Convert to fractions
        x = monomer_pct / 100     # Monomer fraction
        y = ct_pct / 100          # CT fraction
        z = rp1_pct / 100         # RP1+CT fraction
        w = rp1_only_pct / 100    # RP1 only fraction
        
        print(f"Creating mixture: x={x:.2f} (Monomer), y={y:.2f} (CT), z={z:.2f} (RP1+CT), w={w:.2f} (RP1)")
        
        # Create the mixture - ensure components are combined correctly
        list_fl_t_mixture = np.zeros_like(self.list_fl_t_monomer)
        
        # Add monomer component if available and requested
        if x > 0 and np.max(self.list_fl_t_monomer) > 0:
            list_fl_t_mixture += self.list_fl_t_monomer * x
            
        # Add CT component if available and requested
        if y > 0 and np.max(self.list_fl_t_ct) > 0:
            list_fl_t_mixture += self.list_fl_t_ct * y
            
        # Add RP1+CT component if available and requested
        if z > 0 and np.max(self.list_fl_t_ct_rp1) > 0:
            list_fl_t_mixture += self.list_fl_t_ct_rp1 * z
        
        # Add RP1 only component if available and requested
        if w > 0 and np.max(self.list_fl_t_rp1) > 0:
            list_fl_t_mixture += self.list_fl_t_rp1 * w
        
        # Calculate time-resolved and wavelength-resolved profiles
        fl_time_mixture = self.integrate_over_wavelength(list_fl_t_mixture)[self.start_idx:]
        fl_wave_mixture = self.integrate_over_time(list_fl_t_mixture)
        
        # Create a consistent label in the exact format seen in the CSV file
        # This ensures file naming will be consistent with other methods
        monomer_int = int(monomer_pct)
        ct_int = int(ct_pct)
        rp1_ct_int = int(rp1_pct)
        rp1_only_int = int(rp1_only_pct)
        
        # Format matching the CSV file in Image 2
        label = f'Monomer {monomer_int}%: {monomer_int}% M, {ct_int}% CT, {rp1_ct_int}% RP1+CT, {rp1_only_int}% RP1'
        
        # Save data if requested
        if save and hasattr(self, 'save_simulation_data'):
            self.save_simulation_data(fl_time_mixture, fl_wave_mixture, label)
        
        # Create the plot
        self.plot_specific_mixture(fl_time_mixture, fl_wave_mixture, label, save, show)
        
        return fl_time_mixture, fl_wave_mixture

    def plot_specific_mixture(self, fl_time_mixture, fl_wave_mixture, label, save=True, show=True):
        """
        Create publication-quality plots (time-resolved and wavelength-resolved)
        for a specific fluorescence mixture, comparing with experimental data.
        
        Parameters:
        -----------
        fl_time_mixture : ndarray
            1D array of time-resolved fluorescence intensity for the mixture
        fl_wave_mixture : ndarray
            1D array of wavelength-resolved fluorescence intensity for the mixture
        label : str
            Descriptive label for the mixture (used in legends, filenames)
        save : bool, optional
            Whether to save the figure (default: True)
        show : bool, optional
            Whether to display the figure (default: True)

        Returns:
        --------
        matplotlib.figure.Figure or None
            The generated figure object if plotting is successful, otherwise None
        """
        # --- 1. Setup Figure and Axes ---
        fig = plt.figure(figsize=(8, 3.5))  # Wider figure for side-by-side plots
        
        # Use GridSpec with proper spacing
        gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 1], wspace=0.3, 
                    left=0.1, right=0.95, bottom=0.15, top=0.8)
        
        ax1 = fig.add_subplot(gs[0, 0])  # Time-resolved plot (left)
        ax2 = fig.add_subplot(gs[0, 1])  # Wavelength-resolved plot (right)
        
        # Define consistent colors
        mixture_color = '#69b3a2'  # Teal-greenish color for simulation
        exp_color = '#000000'      # Black for experiment
        
        # --- 2. Add parameter title at the top of the figure ---
        param_str = self._generate_param_display_string()
        if param_str:
            fig.suptitle(param_str, fontsize=9, y=0.98)
        
        # --- 3. Time-Resolved Plot (ax1) ---
        plot_time_success = False
        if fl_time_mixture is not None and self.time_ns_trimmed is not None and len(fl_time_mixture) > 0 and fl_time_mixture[0] > 1e-9:
            # Normalize time-resolved data
            fl_time_norm = fl_time_mixture / fl_time_mixture[0]
            
            # Plot simulation line
            ax1.plot(self.time_ns_trimmed, fl_time_norm, '-',
                    color=mixture_color, linewidth=1.5)
            
            # Plot experimental line
            exp_decay = self.get_exp_biexponential_decay(self.time_ns_trimmed)
            if exp_decay is not None and len(exp_decay) == len(self.time_ns_trimmed):
                ax1.plot(self.time_ns_trimmed, exp_decay, '--',
                        color=exp_color, linewidth=1.5)
            
            # Fit result text box
            fit_result = self.biexponential_fit(fl_time_norm)
            if fit_result is not None:
                a1, tau1, a2, tau2, tau_avg, r2 = fit_result  # Fixed unpacking
                fit_text = (f'$a_1$={a1*100:.1f}%, $\\tau_1$={tau1:.2f} ns\n'
                            f'$a_2$={a2*100:.1f}%, $\\tau_2$={tau2:.2f} ns\n'
                            f'$\\tau_{{avg}}$={tau_avg:.2f} ns')
                
                # Position text box in top right
                ax1.text(0.95, 0.95, fit_text, transform=ax1.transAxes,
                        verticalalignment='top', horizontalalignment='right',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, edgecolor='lightgray'),
                        fontsize=8)
            
            # Axes configuration
            ax1.set_xlim(0, min(8, self.time_ns_trimmed[-1]))
            ax1.set_ylim(0, 1.05)
            ax1.set_xlabel('Time (ns)')
            ax1.set_ylabel('Normalized Intensity')
            ax1.minorticks_on()
            ax1.grid(True, alpha=0.3)
            plot_time_success = True
        else:
            ax1.text(0.5, 0.5, "No valid time-resolved data to plot",
                    horizontalalignment='center', verticalalignment='center',
                    transform=ax1.transAxes, fontsize=9, color='grey')
            ax1.set_xlabel('Time (ns)')
            ax1.set_ylabel('Normalized Intensity')
        
        # --- 4. Wavelength-Resolved Plot (ax2) ---
        plot_wave_success = False
        if fl_wave_mixture is not None and self.lambda_axis is not None and len(fl_wave_mixture) > 0 and np.max(fl_wave_mixture) > 1e-9:
            # Normalize and sort wavelength data
            fl_wave_norm = fl_wave_mixture / np.max(fl_wave_mixture)
            sort_indices = np.argsort(self.lambda_axis)
            sorted_lambda = self.lambda_axis[sort_indices]
            sorted_fl_wave_norm = fl_wave_norm[sort_indices]
            
            # Plot simulation line
            ax2.plot(sorted_lambda, sorted_fl_wave_norm, '-',
                    color=mixture_color, linewidth=1.5)
            
            # Plot experimental line if available
            if self.exp_data is not None and self.exp_data.shape[1] >= 2:
                try:
                    exp_x, exp_y = self.exp_data[:, 0], self.exp_data[:, 1]
                    exp_sort_indices = np.argsort(exp_x)
                    sorted_exp_x, sorted_exp_y = exp_x[exp_sort_indices], exp_y[exp_sort_indices]
                    max_exp_y = np.max(sorted_exp_y)
                    if max_exp_y > 1e-9:
                        exp_data_norm = sorted_exp_y / max_exp_y
                        ax2.plot(sorted_exp_x, exp_data_norm, '--',
                                color=exp_color, linewidth=1.5)
                except Exception as e:
                    print(f"Error plotting exp spectrum: {e}")
            
            # Axes configuration
            ax2.set_xlim(660, 750)
            ax2.set_ylim(0, 1.05)
            ax2.set_xlabel('Wavelength (nm)')
            ax2.set_ylabel('Normalized Intensity')
            ax2.minorticks_on()
            ax2.grid(True, alpha=0.3)
            plot_wave_success = True
        else:
            ax2.text(0.5, 0.5, "No valid wavelength-resolved data to plot",
                    horizontalalignment='center', verticalalignment='center',
                    transform=ax2.transAxes, fontsize=9, color='grey')
            ax2.set_xlabel('Wavelength (nm)')
            ax2.set_ylabel('Normalized Intensity')
        
        # --- 5. Add Combined Legend at Bottom ---
        if plot_time_success or plot_wave_success:
            # Create proxy artists for the legend
            sim_line = plt.Line2D([0], [0], color=mixture_color, lw=1.5)
            exp_line = plt.Line2D([0], [0], color=exp_color, lw=1.5, linestyle='--')
            
            # Position legend centered at bottom
            fig.legend(
                [sim_line, exp_line],
                [f'Sim: {label}', 'Exp. Ref.'],
                loc='lower center',
                bbox_to_anchor=(0.5, 0.02),
                ncol=2,
                frameon=False,
                fontsize=8
            )
        
        # --- 6. Adjust layout ---
        plt.tight_layout(rect=[0, 0.08, 1, 0.9])  # Leave space for title and legend
        
        # --- 7. Save and Show ---
        if save and (plot_time_success or plot_wave_success):
            # Create clean filename
            filename_base = label.replace('%', 'pct').replace('+','_').replace(', ', '_').replace(' ', '_').replace(':','_')
            filename_base = "".join(c for c in filename_base if c.isalnum() or c in ('_', '-')).strip()[:100]
            
            # Add parameter suffix to filename
            param_suffix = ""
            if hasattr(self, 'params') and self.params:
                fl_val = str(self.params.get('TIME_FL', ['na']))
                isc_val = str(self.params.get('TIME_ISC', ['na']))
                rp_val = str(self.params.get('TIME_to_RP1', ['na']))
                coupling = str(self.params.get('COUPLING_value', ['na']))
                param_suffix = f"_FL{fl_val}_ISC{isc_val}_RP{rp_val}_Coup{coupling}_start_time_{self.start_time}"
            
            # Save the figure
            png_path = os.path.join(self.output_dir, f"plot_{filename_base}{param_suffix}.png")
            try:
                fig.savefig(png_path, dpi=300, format='png', bbox_inches='tight')
                print(f"Saved figure to {png_path}")
            except Exception as e:
                print(f"Error saving figure: {e}")
        
        if show:
            plt.show()
        else:
            plt.close(fig)
        
        return fig if (plot_time_success or plot_wave_success) else None