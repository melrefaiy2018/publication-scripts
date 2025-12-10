import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit

# --- True Science/Nature Journal Style ---
plt.style.use('default')  # Start fresh
mpl.rcParams.update({
    # Font - exactly what Science/Nature uses
    'font.family': ['Arial', 'Helvetica'],
    'font.size': 7,
    'axes.labelsize': 7,
    'axes.titlesize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'legend.fontsize': 6,
    'font.weight': 'normal',
    
    # Minimal, clean lines
    'axes.linewidth': 0.5,
    'lines.linewidth': 1.0,
    'patch.linewidth': 0.5,
    
    # Clean ticks - minimal and unobtrusive
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.minor.size': 1.2,
    'ytick.minor.size': 1.2,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.minor.width': 0.3,
    'ytick.minor.width': 0.3,
    'xtick.direction': 'out',  # Science/Nature typically use outward ticks
    'ytick.direction': 'out',
    'xtick.minor.visible': False,  # Keep it minimal
    'ytick.minor.visible': False,
    
    # Completely clean - no unnecessary elements
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.spines.left': True,
    'axes.spines.bottom': True,
    'axes.grid': False,
    'axes.axisbelow': True,
    
    # Pure white background
    'axes.facecolor': 'white',
    'figure.facecolor': 'white',
    'savefig.facecolor': 'white',
    
    # High quality output
    'figure.dpi': 150,
    'savefig.dpi': 600,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.05,
    
    # Clean legend
    'legend.frameon': False,  # No frame around legend
    'legend.numpoints': 1,
    'legend.scatterpoints': 1,
    'legend.columnspacing': 1.0,
    'legend.handlelength': 1.5,
})

# Science/Nature color scheme - very clean and minimal
COLORS = {
    'sim': '#000000',        # Black for simulation (primary data)
    'exp': '#808080',        # Gray for experimental
    'fit': '#FF0000',        # Red for fits (when needed)
    'accent': '#0000FF',     # Blue accent
}

# --- Functions (keeping original functionality) ---
def bi_exponential_decay_plot_func(t, a1_pct, tau1, a2_pct, tau2, tau_avg=None):
    a1 = a1_pct / 100.0
    a2 = a2_pct / 100.0
    return a1 * np.exp(-t / tau1) + a2 * np.exp(-t / tau2)

def bi_exponential_decay_fit_func(t, a1_frac, tau1, a2_frac, tau2):
    return a1_frac * np.exp(-t / tau1) + a2_frac * np.exp(-t / tau2)

def extract_params_from_name(run_name):
    params = {}
    direct_matches = re.findall(r'([A-Z]+)(\d+\.?\d*)', run_name)
    for key, value in direct_matches:
        params[key] = float(value) if '.' in value else int(value)
    return params

def compute_mse(y_true, y_pred):
    return np.mean((y_true - y_pred) ** 2)

def get_region_data(wl, inten, region=(710, 740)):
    mask = (wl >= region[0]) & (wl <= region[1])
    return wl[mask], inten[mask]

def calculate_auc(wl, inten):
    return np.trapz(inten, wl)

def area_error(sim_auc, exp_auc):
    return abs(sim_auc - exp_auc) / exp_auc

def clean_axes(ax):
    """Apply Science/Nature journal clean styling"""
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Clean, thin remaining spines
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    
    # Minimal ticks
    ax.tick_params(axis='both', which='major', 
                   length=2.5, width=0.5, direction='out',
                   top=False, right=False, 
                   color='black', labelcolor='black')

def format_panel_label(ax, label, x=-0.15, y=1.05):
    """Add clean panel labels like in Science/Nature"""
    ax.text(x, y, label, transform=ax.transAxes, 
            fontsize=8, fontweight='bold', va='bottom', ha='right')

# --- Main Analysis Function ---
def analyze_simulations(base_dir='Simulation',
                        exp_wavelength_file=None,
                        exp_time_params=None,
                        wavelength_range=(650, 760),
                        N_best=10):

    os.makedirs('outputs', exist_ok=True)
    os.makedirs('plots/best_models', exist_ok=True)

    run_dirs = glob.glob(os.path.join(base_dir, 'Run_*'))
    run_dirs = [d for d in run_dirs if os.path.isdir(d)]

    time_data_map = {}
    wavelength_data_map = {}
    skipped_runs = []
    valid_runs = set()

    print(f"Found {len(run_dirs)} run directories to check...")

    # Data loading (keeping original logic)
    for run_dir in run_dirs:
        run_name = os.path.basename(run_dir)
        data_subdir = os.path.join(run_dir, 'Data')
        
        if not os.path.exists(data_subdir):
            skipped_runs.append(f"{run_name} (no data directory)")
            continue

        time_files = glob.glob(os.path.join(data_subdir, 'time_resolved_*.npy'))
        wavelength_files = glob.glob(os.path.join(data_subdir, 'wavelength_resolved_*.npy'))
        
        if not time_files and not wavelength_files:
            skipped_runs.append(f"{run_name} (no data files)")
            continue
        
        run_has_valid_data = False
            
        if time_files:
            try:
                time_data = np.load(time_files[0])
                if time_data.shape[0] > 0 and time_data.shape[1] >= 2:
                    indices = np.argsort(time_data[:, 0])
                    time_data_map[run_name] = time_data[indices]
                    run_has_valid_data = True
            except Exception as e:
                print(f"Warning: Failed to load time data for {run_name}: {e}")

        if wavelength_files:
            try:
                wavelength_data = np.load(wavelength_files[0])
                if wavelength_data.shape[0] > 0 and wavelength_data.shape[1] >= 2:
                    indices = np.argsort(wavelength_data[:, 0])
                    wavelength_data_map[run_name] = wavelength_data[indices]
                    run_has_valid_data = True
            except Exception as e:
                print(f"Warning: Failed to load wavelength data for {run_name}: {e}")
        
        if run_has_valid_data:
            valid_runs.add(run_name)
        else:
            skipped_runs.append(f"{run_name} (invalid data)")
    
    if skipped_runs:
        print(f"Skipped {len(skipped_runs)} runs due to missing/invalid data")
    
    print(f"Processing {len(valid_runs)} runs with valid data")
    
    if len(valid_runs) == 0:
        print("Error: No valid simulation data found.")
        return
    
    if len(valid_runs) < N_best:
        N_best = len(valid_runs)

    # Load experimental data
    exp_wl_data = None
    if exp_wavelength_file:
        exp_wl_data_raw = np.load(exp_wavelength_file)
        indices = np.argsort(exp_wl_data_raw[:, 0])
        exp_wl_data = exp_wl_data_raw[indices]
        exp_wl_data[:,1] /= np.max(exp_wl_data[:,1])

    # Error calculation (keeping original logic)
    errors = []
    for run_name in valid_runs:
        time_mse = np.nan
        wavelength_mse = np.nan
        region_auc_err = np.nan
        tau_avg_err = np.nan
        
        if run_name in time_data_map:
            t_sim = time_data_map[run_name][:, 0]
            pop_sim = time_data_map[run_name][:, 1] / np.max(time_data_map[run_name][:, 1])
            exp_decay = bi_exponential_decay_plot_func(t_sim, **exp_time_params)
            time_mse = compute_mse(exp_decay, pop_sim)
            
            try:
                popt, _ = curve_fit(bi_exponential_decay_fit_func, t_sim, pop_sim,
                                  p0=[0.5, 1.0, 0.5, 5.0],
                                  bounds=([0, 0, 0, 0], [1, np.inf, 1, np.inf]),
                                  maxfev=5000)
                a1_fit, tau1_fit, a2_fit, tau2_fit = popt
                tau_avg_sim = (a1_fit * tau1_fit + a2_fit * tau2_fit) / (a1_fit + a2_fit)
                tau_avg_exp = exp_time_params['tau_avg']
                tau_avg_err = abs(tau_avg_sim - tau_avg_exp) / tau_avg_exp
            except:
                tau_avg_err = np.nan

        if run_name in wavelength_data_map and exp_wl_data is not None:
            wl_sim = wavelength_data_map[run_name][:, 0]
            int_sim = wavelength_data_map[run_name][:, 1] / np.max(wavelength_data_map[run_name][:, 1])
            exp_interp = np.interp(wl_sim, exp_wl_data[:,0], exp_wl_data[:,1])
            wavelength_mse = compute_mse(exp_interp, int_sim)

            sim_wl_region, sim_inten_region = get_region_data(wl_sim, int_sim)
            exp_wl_region, exp_inten_region = get_region_data(exp_wl_data[:, 0], exp_wl_data[:, 1])
            sim_auc = calculate_auc(sim_wl_region, sim_inten_region)
            exp_auc = calculate_auc(exp_wl_region, exp_inten_region)
            region_auc_err = area_error(sim_auc, exp_auc)

        # Calculate score
        alpha, beta, gamma = 0.1, 0.5, 0.4
        score_components = []
        if not np.isnan(time_mse): score_components.append(alpha * time_mse)
        if not np.isnan(region_auc_err): score_components.append(beta * region_auc_err)
        if not np.isnan(tau_avg_err): score_components.append(gamma * tau_avg_err)
        
        balanced_score = sum(score_components) if score_components else np.inf
        errors.append((run_name, time_mse, wavelength_mse, region_auc_err, tau_avg_err, balanced_score))

    errors_sorted = sorted(errors, key=lambda x: x[5])
    best_run_names = [e[0] for e in errors_sorted[:N_best]]

    # Save results
    rows = []
    for rank, (run_name, time_mse, wavelength_mse, region_auc_err, tau_avg_err, balanced_score) in enumerate(errors_sorted[:N_best], start=1):
        params = extract_params_from_name(run_name)
        row = {'Rank': rank, 'Run Name': run_name, **params,
               'Time MSE': time_mse, 'Wavelength MSE': wavelength_mse,
               'Region AUC Error': region_auc_err, 'Tau Avg Error': tau_avg_err,
               'Score': balanced_score}
        rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv('outputs/best_model_errors.csv', index=False)

    # --- Create Science/Nature Style Plots ---
    for idx, run_name in enumerate(best_run_names):
        has_time_data = run_name in time_data_map
        has_wavelength_data = run_name in wavelength_data_map
        
        if not has_time_data and not has_wavelength_data:
            continue
        
        # Science/Nature figure sizing - exact journal specifications
        if has_time_data and has_wavelength_data:
            # Two-column figure: 183mm = 7.2 inches
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.8))
            plt.subplots_adjust(wspace=0.4)
            
            # Add panel labels
            format_panel_label(ax1, 'a')
            format_panel_label(ax2, 'b')
            
        elif has_time_data:
            # Single column: 89mm = 3.5 inches  
            fig, ax1 = plt.subplots(1, 1, figsize=(3.5, 2.8))
            ax2 = None
        else:
            fig, ax2 = plt.subplots(1, 1, figsize=(3.5, 2.8))
            ax1 = None
        
        # Time-resolved plot
        if has_time_data and ax1 is not None:
            t = time_data_map[run_name][:, 0]
            pop = time_data_map[run_name][:, 1] / np.max(time_data_map[run_name][:, 1])
            
            exp_t = np.linspace(0, np.max(t), 500)
            exp_pop = bi_exponential_decay_plot_func(exp_t, **exp_time_params)
            
            # Clean, minimal plotting - Science/Nature style
            ax1.plot(t, pop, color=COLORS['sim'], linewidth=1.0, label='Simulation')
            ax1.plot(exp_t, exp_pop, color=COLORS['exp'], linewidth=1.0, 
                    linestyle='--', label='Experiment')

            # Minimal parameter display (if fit successful)
            try:
                popt, _ = curve_fit(bi_exponential_decay_fit_func, t, pop,
                                  p0=[0.5, 1.0, 0.5, 5.0],
                                  bounds=([0, 0, 0, 0], [1, np.inf, 1, np.inf]),
                                  maxfev=5000)
                a1_fit, tau1_fit, a2_fit, tau2_fit = popt
                tau_avg = (a1_fit * tau1_fit + a2_fit * tau2_fit) / (a1_fit + a2_fit)
                
                # Very minimal text - just key result
                ax1.text(0.6, 0.8, f'τ = {tau_avg:.1f} ns', 
                        transform=ax1.transAxes, fontsize=6, 
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
                                 edgecolor='none', alpha=0.8))
            except:
                pass

            # Clean axes and labels
            ax1.set_xlabel('Time (ns)')
            ax1.set_ylabel('Normalized fluorescence')
            clean_axes(ax1)
            
            # Minimal legend
            ax1.legend(frameon=False, loc='upper right')

        # Wavelength-resolved plot  
        if has_wavelength_data and ax2 is not None and exp_wl_data is not None:
            wl = wavelength_data_map[run_name][:, 0]
            inten = wavelength_data_map[run_name][:, 1] / np.max(wavelength_data_map[run_name][:, 1])
            
            exp_wl = exp_wl_data[:, 0]
            exp_inten = exp_wl_data[:, 1]

            # Clean plotting
            ax2.plot(wl, inten, color=COLORS['sim'], linewidth=1.0, label='Simulation')
            ax2.plot(exp_wl, exp_inten, color=COLORS['exp'], linewidth=1.0, 
                    linestyle='--', label='Experiment')
            
            ax2.set_xlim(wavelength_range)
            ax2.set_xlabel('Wavelength (nm)')
            ax2.set_ylabel('Normalized intensity')
            clean_axes(ax2)
            
            # Minimal legend
            ax2.legend(frameon=False, loc='upper right')

        # No title - Science/Nature figures typically don't have titles
        # Information goes in the caption
        
        plt.tight_layout()
        
        # Save in publication formats
        base_filename = f'plots/best_models/figure_{idx+1:02d}'
        plt.savefig(f'{base_filename}.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(f'{base_filename}.png', dpi=600, bbox_inches='tight')
        plt.close()

    # Error scatter plot - also in clean style
    time_mses = [e[1] for e in errors_sorted if not np.isnan(e[1])]
    wl_mses = [e[2] for e in errors_sorted if not np.isnan(e[2])]
    
    if time_mses and wl_mses:
        fig, ax = plt.subplots(figsize=(3.5, 3.5))
        
        # Clean scatter plot
        ax.scatter(time_mses, wl_mses, s=12, alpha=0.6, 
                  color='black', edgecolors='none')
        
        # Highlight best models
        best_time = [e[1] for e in errors_sorted[:N_best] if not np.isnan(e[1])]
        best_wl = [e[2] for e in errors_sorted[:N_best] if not np.isnan(e[2])]
        
        if best_time and best_wl:
            ax.scatter(best_time, best_wl, s=16, 
                      color='red', edgecolors='none', alpha=0.8)
        
        ax.set_xlabel('Time MSE')
        ax.set_ylabel('Wavelength MSE')
        clean_axes(ax)
        
        plt.tight_layout()
        plt.savefig('plots/error_comparison.pdf', dpi=600, bbox_inches='tight')
        plt.savefig('plots/error_comparison.png', dpi=600, bbox_inches='tight')
        plt.show()

    print(f"\nGenerated clean, publication-ready figures for {len(best_run_names)} models")
    print("Figures saved in Science/Nature journal style")

# --- Run Analysis ---
if __name__ == "__main__":
    analyze_simulations(
        base_dir='Simulation',
        exp_wavelength_file='/stor/work/Raccah/ActiveProjects/IsiA/CT/Apr25/NeededFiles/Exp_spectra_gabriela/Fl_IsiA_monomer_300k_Gabriela_2025.npy',
        exp_time_params={'a1_pct':23.7, 'tau1':0.21, 'a2_pct':76.3, 'tau2':3.89, 'tau_avg':3.02},
        wavelength_range=(650, 760),
        N_best=10
    )