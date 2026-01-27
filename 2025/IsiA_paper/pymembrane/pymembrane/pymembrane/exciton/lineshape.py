from pathlib import Path

import numba
import numpy as np

from pymembrane.util.physical_constants import kB, hbar


class PigmentLineshape_Kubo:
    """Kubo lineshape function for pigment spectroscopy.
    
    This class calculates lineshape functions using the Kubo formalism based on spectral 
    density parameters that define system-bath coupling. The class supports multiple spectral 
    density models (Drude-Lorentz, Brownian oscillator, Renger) and can include explicit 
    vibronic contributions.
    
    Attributes:
        t_axis (np.ndarray): Time axis for lineshape calculations in femtoseconds.
        e_lambda (float): Total reorganization energy in cm^-1.
        renger_gamma (bool): Flag indicating if Renger spectral density is used.
        explicit_vib (bool): Flag indicating if explicit vibronic modes are included.
    """
    
    t_axis = np.arange(0, 2500, 0.1) # Units: fs
    # t_axis = np.arange(0, 1000, 0.1) # Units: fs

    def __init__(self, list_sd_param, explicit_vib=None, scratch_dir=None, clear_scratch=False):
        """Initialize the PigmentLineshape_Kubo class.
        
        Args:
            list_sd_param (list): List of spectral density parameters. Each element is a list 
                containing [function_gt, function_lambda, dict_params] where:
                - function_gt: Function to calculate g(t), takes (t_axis, temp, **params)
                - function_lambda: Function to calculate reorganization energy, takes (**params)
                - dict_params: Dictionary of parameters for the functions
                
                Example:
                    [[lineshape_drudelorentz,
                      reorganization_drudelorentz,
                      {'w_axis': w_axis,
                       'e_lambda': e_lambda,  # Units: cm^-1
                       'gamma': gamma}],      # Units: cm^-1
                     [lineshape_brownian_oscillator,
                      reorganization_brownian_oscillator,
                      {'w_axis': w_axis,
                       'e_lambda': e_lambda,  # Units: cm^-1
                       'gamma': gamma,        # Units: cm^-1
                       'omega': omega}],      # Units: cm^-1
                     [lineshape_renger,
                      reorganization_renger,
                      {'w_axis': w_axis,
                       'S0': S0,
                       's1': s1,
                       's2': s2,
                       'w1': w1,              # Units: cm^-1
                       'w2': w2}]]            # Units: cm^-1
                       
            explicit_vib (tuple, optional): Tuple of (huang_rhys_factors, frequencies) for 
                explicit vibronic modes. Defaults to None.
            scratch_dir (str or Path, optional): Directory path for caching computed lineshape 
                functions. Defaults to None.
            clear_scratch (bool, optional): If True, archive existing cached files by renaming 
                them with '#' prefix. If False, load precomputed values. Defaults to False.
        """
        self.__list_sd_param = list_sd_param
        self.__gt_by_temp = {}
        self.e_lambda = np.sum([func2(**dict_param) for (func1, func2, dict_param) in list_sd_param])

        self.renger_gamma = False
        for (func1, func2, dict_param) in list_sd_param:
            if func1 is lineshape_renger:
                self.renger_gamma = True
                self.renger_gamma_t = -np.abs(self.t_axis)/250
                self.renger_Jw_param = dict_param

        self.explicit_vib = False
        self._explicit_vib_s = None
        self._explicit_vib_w = None
        self._dict_vib_absorption = {}
        self._dict_vib_fluorescence = {}
        if explicit_vib is not None:
            self.explicit_vib = True
            self._explicit_vib_s, self._explicit_vib_w = explicit_vib
            self._explicit_vib_s = np.array(self._explicit_vib_s)
            self._explicit_vib_w = np.array(self._explicit_vib_w)

        self.dir_scratch = Path(scratch_dir) if scratch_dir else None
        # To save a histroy from the generated lineshape npy files.
        if self.dir_scratch:
            # Load data of pre-computed values
            if self.dir_scratch.is_dir():
                if clear_scratch:
                    # Rename files to start with '#' before deleting files without the '#' prefix
                    for file in self.dir_scratch.iterdir():
                        if not file.name.startswith('#'):
                            new_name = '#' + file.name
                            count = 1
                            while (self.dir_scratch / new_name).exists():
                                new_name = f'#{count}_{file.name}'
                                count += 1
                            file.rename(self.dir_scratch / new_name)
                else:
                    # Load precomputed values
                    list_files = list(self.dir_scratch.glob('T*.npy'))
                    for file in list_files:
                        temp = int(file.stem[1:])
                        gt_temp = np.load(file)
                        if len(gt_temp) == len(self.t_axis):
                            self.__gt_by_temp[temp] = gt_temp
                        else:
                            print(f"Attempted to load lineshape function from {self.dir_scratch / 'T*.npy'}:")
                            print(
                                f'length of gt from file ({len(gt_temp)}) different from current t_axis ({len(self.t_axis)})')
                            print(f'Resetting scratch directory to None!')
                            self.dir_scratch = None
            else:
                print(f'Error: Scratch path {self.dir_scratch} is not a directory. Resetting to None.')
                self.dir_scratch = None

    @property
    def list_sd_param(self):
        """Get the list of spectral density parameters.
        
        Returns:
            list: List of spectral density parameters.
        """
        return self.__list_sd_param

    def g_t(self, temp):
        """Calculate or retrieve the lineshape function g(t) at a given temperature.
        
        Args:
            temp (float): Temperature in Kelvin.
            
        Returns:
            np.ndarray: Complex-valued lineshape function g(t).
        """
        if not (temp in self.__gt_by_temp.keys()):
            self.__gt_by_temp[temp] = self._calculate_gt_by_temp(temp)
        return self.__gt_by_temp[temp]

    def _calculate_gt_by_temp(self, temp):
        """Calculate the lineshape function g(t) for a specific temperature.
        
        Sums contributions from all spectral density models and optionally caches 
        the result to disk.
        
        Args:
            temp (float): Temperature in Kelvin.
            
        Returns:
            np.ndarray: Complex-valued lineshape function g(t).
        """
        g_t = np.zeros_like(self.t_axis) + 0j
        for (func1, func2, dict_param) in self.list_sd_param:
            g_t += func1(self.t_axis, temp, **dict_param)

        if self.dir_scratch:
            np.save(self.dir_scratch / f'T{int(temp)}.npy', g_t)
        return g_t

    def vibronic_absorption(self, temp):
        """Calculate the vibronic absorption lineshape in the time domain.
        
        Computes the absorption lineshape including explicit vibronic contributions 
        from discrete modes.
        
        Args:
            temp (float): Temperature in Kelvin.
            
        Returns:
            np.ndarray: Complex-valued time-domain absorption lineshape with vibronic 
                structure.
        """
        if (not (temp in self._dict_vib_absorption.keys())) and self.explicit_vib:
            t_axis = self.t_axis
            abs_t = 0 * t_axis + 0j
            for (s_vib, omega_vib) in zip(self._explicit_vib_s, self._explicit_vib_w):
                abs_t += s_vib * np.exp(-1j * (omega_vib / hbar) * t_axis) * np.exp(-self.g_t(temp))
            self._dict_vib_absorption[temp] = abs_t

        return self._dict_vib_absorption[temp]

    def vibronic_fluorescence(self, temp):
        """Calculate the vibronic fluorescence lineshape in the time domain.
        
        Computes the fluorescence lineshape including explicit vibronic contributions 
        from discrete modes with proper Stokes shift.
        
        Args:
            temp (float): Temperature in Kelvin.
            
        Returns:
            np.ndarray: Complex-valued time-domain fluorescence lineshape with vibronic 
                structure.
        """
        if (not (temp in self._dict_vib_fluorescence.keys())) and self.explicit_vib:
            t_axis = self.t_axis
            fl_t = 0 * t_axis + 0j
            for (s_vib, omega_vib) in zip(self._explicit_vib_s, self._explicit_vib_w):
                fl_t += s_vib * np.exp(1j * (omega_vib / hbar) * t_axis) * np.exp(
                    1j * self.e_lambda / hbar * t_axis) * np.exp(-np.conj(self.g_t(temp)))

            self._dict_vib_fluorescence[temp] = fl_t

        return self._dict_vib_fluorescence[temp]

@numba.njit()
def numeric_integrate_sd(t_axis, w_axis, sd_over_w2, beta):
    """Numerically integrate the spectral density to compute g(t).
    
    Performs numerical integration of the spectral density using the trapezoidal rule 
    to calculate both real and imaginary parts of the lineshape function g(t).
    
    Args:
        t_axis (np.ndarray): Time axis in femtoseconds.
        w_axis (np.ndarray): Frequency axis in cm^-1.
        sd_over_w2 (np.ndarray): Spectral density divided by frequency squared.
        beta (float): Inverse temperature (1/kB*T) in inverse energy units.
        
    Returns:
        np.ndarray: Complex-valued lineshape function g(t).
    """
    g=[]
    for t in t_axis:
         real_part = sd_over_w2 * ((1 - np.cos(w_axis * t / hbar)) /
                                   np.tanh(0.5 * w_axis * beta))
         imaginary_part = sd_over_w2 * (np.sin(w_axis * t / hbar) - w_axis * t / hbar)
         g.append(np.trapz(y=real_part, x=w_axis) / np.pi + 1j / np.pi *
                  np.trapz(y=imaginary_part, x=w_axis))
    return np.array(g, dtype=np.complex128)

def lineshape_drudelorentz(t_axis, temp, w_axis, e_lambda, gamma, **kwargs):
    """Calculate the lineshape function using the Drude-Lorentz spectral density.
    
    Args:
        t_axis (np.ndarray): Time axis in femtoseconds.
        temp (float): Temperature in Kelvin.
        w_axis (np.ndarray): Frequency axis in cm^-1.
        e_lambda (float or array-like): Reorganization energy in cm^-1. Can be a single 
            value or array for multiple modes.
        gamma (float or array-like): Damping parameter in cm^-1. Can be a single value 
            or array for multiple modes.
        **kwargs: Additional keyword arguments (ignored).
        
    Returns:
        np.ndarray: Complex-valued lineshape function g(t).
    """
    beta = 1 / (kB * temp)
    sd_over_w2 = np.divide(spectraldensity_drudelorentz(w_axis, e_lambda, gamma), w_axis ** 2)
    return numeric_integrate_sd(t_axis, w_axis, sd_over_w2, beta)

def spectraldensity_drudelorentz(w_axis, e_lambda, gamma):
    """Calculate the Drude-Lorentz spectral density.
    
    The Drude-Lorentz spectral density models overdamped oscillator modes and is 
    commonly used for describing solvent relaxation.
    
    Args:
        w_axis (np.ndarray): Frequency axis in cm^-1.
        e_lambda (float or array-like): Reorganization energy in cm^-1. Can be a single 
            value or array for multiple modes.
        gamma (float or array-like): Damping parameter in cm^-1. Can be a single value 
            or array for multiple modes.
            
    Returns:
        np.ndarray: Spectral density J(ω) in appropriate units.
        
    Raises:
        ValueError: If e_lambda and gamma are arrays with different lengths.
    """
    if all([hasattr(func_input, '__iter__') for func_input in [e_lambda, gamma]]):
        e_lambda = np.array(e_lambda)
        gamma = np.array(gamma)
        if len(e_lambda) != len(gamma):
            raise ValueError(
                f"Length mismatch: e_lambda has {len(e_lambda)} elements but "
                f"gamma has {len(gamma)} elements. Both arrays must have the same length."
            )
        w_axis = np.array(w_axis, dtype=np.float64)
        sd = 0.0 * w_axis
        for (e_i, gamma_i) in zip(e_lambda, gamma):
            sd += 2 * e_i * gamma_i * w_axis / (w_axis ** 2 + gamma_i ** 2)
    else:
        sd = 2 * e_lambda * gamma * w_axis / (w_axis ** 2 + gamma ** 2)
    return sd

def reorganization_drudelorentz(e_lambda, gamma, **kwargs):
    """Calculate the reorganization energy for Drude-Lorentz spectral density.
    
    Args:
        e_lambda (float or array-like): Reorganization energy in cm^-1. Can be a single 
            value or array for multiple modes.
        gamma (float or array-like): Damping parameter in cm^-1 (not used in calculation 
            but required for consistency).
        **kwargs: Additional keyword arguments (ignored).
        
    Returns:
        float: Total reorganization energy in cm^-1.
    """
    if all([hasattr(func_input, '__iter__') for func_input in [e_lambda, gamma]]):
        return np.sum(e_lambda)
    else:
        return e_lambda

def lineshape_brownian_oscillator(t_axis, temp, w_axis, e_lambda, gamma, omega, **kwargs):
    """Calculate the lineshape function using the Brownian oscillator spectral density.
    
    Args:
        t_axis (np.ndarray): Time axis in femtoseconds.
        temp (float): Temperature in Kelvin.
        w_axis (np.ndarray): Frequency axis in cm^-1.
        e_lambda (float or array-like): Reorganization energy in cm^-1. Can be a single 
            value or array for multiple modes.
        gamma (float or array-like): Damping parameter in cm^-1. Can be a single value 
            or array for multiple modes.
        omega (float or array-like): Characteristic frequency in cm^-1. Can be a single 
            value or array for multiple modes.
        **kwargs: Additional keyword arguments (ignored).
        
    Returns:
        np.ndarray: Complex-valued lineshape function g(t).
    """
    beta = 1 / (kB * temp)
    sd_over_w2 = np.divide(spectraldensity_brownian_oscillator(w_axis, e_lambda, gamma, omega), w_axis ** 2)
    return numeric_integrate_sd(t_axis, w_axis, sd_over_w2, beta)

def spectraldensity_brownian_oscillator(w_axis, e_lambda, gamma, omega):
    """Calculate the Brownian oscillator spectral density.
    
    The Brownian oscillator spectral density models damped harmonic oscillator modes 
    and is suitable for underdamped vibrational modes.
    
    Args:
        w_axis (np.ndarray): Frequency axis in cm^-1.
        e_lambda (float or array-like): Reorganization energy in cm^-1. Can be a single 
            value or array for multiple modes.
        gamma (float or array-like): Damping parameter in cm^-1. Can be a single value 
            or array for multiple modes.
        omega (float or array-like): Characteristic frequency in cm^-1. Can be a single 
            value or array for multiple modes.
            
    Returns:
        np.ndarray: Spectral density J(ω) in appropriate units.
        
    Raises:
        ValueError: If e_lambda, gamma, and omega are arrays with different lengths.
    """
    if all([hasattr(func_input, '__iter__') for func_input in [e_lambda, gamma, omega]]):
        e_lambda = np.array(e_lambda)
        gamma = np.array(gamma)
        omega = np.array(omega)
        if not (len(e_lambda) == len(gamma) == len(omega)):
            raise ValueError(
                f"Length mismatch: e_lambda has {len(e_lambda)} elements, "
                f"gamma has {len(gamma)} elements, and omega has {len(omega)} elements. "
                f"All arrays must have the same length."
            )
        w_axis = np.array(w_axis, dtype=np.float64)
        sd = 0.0 * w_axis
        for (e_i, gamma_i, omega_i) in zip(e_lambda, gamma, omega):
            omega_sq = omega_i ** 2
            sd += 2 * e_i * omega_sq * w_axis * gamma_i / (omega_sq * gamma_i ** 2 + (w_axis ** 2 - omega_sq) ** 2)
    else:
        omega_sq = omega ** 2
        sd = 2 * e_lambda * omega_sq * w_axis * gamma / (omega_sq * gamma ** 2 + (w_axis ** 2 - omega_sq) ** 2)
    return sd

def reorganization_brownian_oscillator(e_lambda, gamma, omega, **kwargs):
    """Calculate the reorganization energy for Brownian oscillator spectral density.
    
    Args:
        e_lambda (float or array-like): Reorganization energy in cm^-1. Can be a single 
            value or array for multiple modes.
        gamma (float or array-like): Damping parameter in cm^-1 (not used in calculation 
            but required for consistency).
        omega (float or array-like): Characteristic frequency in cm^-1 (not used in 
            calculation but required for consistency).
        **kwargs: Additional keyword arguments (ignored).
        
    Returns:
        float: Total reorganization energy in cm^-1.
    """
    if all([hasattr(func_input, '__iter__') for func_input in [e_lambda, gamma, omega]]):
        return np.sum(e_lambda)
    else:
        return e_lambda

@numba.njit()
def reorganization_renger(w_axis, S0, s1, s2, w1, w2):
    """Calculate the reorganization energy for Renger spectral density.
    
    Computes the reorganization energy by integrating J(ω)/ω over the frequency axis.
    
    Args:
        w_axis (np.ndarray): Frequency axis in cm^-1.
        S0 (float): Overall scaling parameter for the spectral density.
        s1 (float): Scaling parameter for the first component.
        s2 (float): Scaling parameter for the second component.
        w1 (float): Characteristic frequency for the first component in cm^-1.
        w2 (float): Characteristic frequency for the second component in cm^-1.
        
    Returns:
        float: Reorganization energy in cm^-1.
    """
    return 1/np.pi*np.trapz(y=spectraldensity_renger(w_axis, S0, s1, s2, w1, w2)/w_axis, x=w_axis)

def lineshape_renger(t_axis, temp, w_axis, S0, s1, s2, w1, w2, **kwargs):
    """Calculate the lineshape function using the Renger spectral density.
    
    The Renger spectral density is designed to reproduce protein-pigment coupling 
    in photosynthetic systems with a realistic frequency dependence.
    
    Args:
        t_axis (np.ndarray): Time axis in femtoseconds.
        temp (float): Temperature in Kelvin.
        w_axis (np.ndarray): Frequency axis in cm^-1.
        S0 (float): Overall scaling parameter for the spectral density.
        s1 (float): Scaling parameter for the first component.
        s2 (float): Scaling parameter for the second component.
        w1 (float): Characteristic frequency for the first component in cm^-1.
        w2 (float): Characteristic frequency for the second component in cm^-1.
        **kwargs: Additional keyword arguments (ignored).
        
    Returns:
        np.ndarray: Complex-valued lineshape function g(t).
    """
    beta = 1 / (kB * temp)
    sd_over_w2 = np.divide(spectraldensity_renger(w_axis, S0, s1, s2, w1, w2), w_axis ** 2)
    list_gt = numeric_integrate_sd(t_axis, w_axis, sd_over_w2, beta)
    return list_gt

@numba.njit()
def spectraldensity_renger(w_axis, S0, s1, s2, w1, w2):
    """Calculate the Renger spectral density.
    
    The Renger spectral density uses a superohmic form with exponential cutoffs 
    to model protein-pigment interactions in photosynthetic complexes.
    
    Args:
        w_axis (np.ndarray): Frequency axis in cm^-1.
        S0 (float): Overall scaling parameter for the spectral density.
        s1 (float): Scaling parameter for the first component.
        s2 (float): Scaling parameter for the second component.
        w1 (float): Characteristic frequency for the first component in cm^-1.
        w2 (float): Characteristic frequency for the second component in cm^-1.
        
    Returns:
        np.ndarray: Spectral density J(ω) in appropriate units.
    """
    prefactor = np.pi*S0 / (s1 + s2) * np.power(w_axis, 5) / 10080
    sd_1 = s1 / w1 ** 4 * np.exp(-np.sqrt(w_axis / w1))
    sd_2 = s2 / w2 ** 4 * np.exp(-np.sqrt(w_axis / w2))
    return prefactor * (sd_1 + sd_2)

