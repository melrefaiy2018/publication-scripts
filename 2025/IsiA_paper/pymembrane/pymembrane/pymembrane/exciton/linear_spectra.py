"""Linear spectroscopy calculations for excitonic systems.

This module provides methods for calculating linear optical spectra including
absorption, fluorescence, and linear dichroism for excitonic aggregates.
"""

import numba
import numpy as np

from pymembrane.exciton.exciton_aggregate import ExcitonicAggregate
from pymembrane.util.physical_constants import hbar



class LinearSpectra(ExcitonicAggregate):

    def __init__(self, structure):
        '''
        In this initialization we assume that the AtomicProtein object has some specific
        properties defined:
        1. Lineshape for each pigment
        2. Domain for each pigment


        Parameters
        ----------
        structure
        '''
        super().__init__(structure)

    @property
    def w_axis(self):
        """Calculate the frequency axis for spectral calculations.
        
        Returns:
            np.ndarray: Frequency axis in units of energy, centered at the
                central frequency of the Hamiltonian.
        """
        t_axis = self.list_domain[0].list_pigments[0].t_axis
        return (2 * np.pi * hbar * np.fft.fftfreq(len(t_axis), t_axis[1] - t_axis[0])) + self.H2_hamiltonian.central_frequency

    def calc_absorption(self, temp):
        """Calculate the total absorption spectrum.
        
        Computes the absorption spectrum by summing contributions from all
        excitonic states, weighted by their transition dipole moments.
        
        Args:
            temp (float): Temperature in Kelvin.
            
        Returns:
            tuple: A tuple containing:
                - np.ndarray: Frequency axis.
                - np.ndarray: Absorption spectrum intensity.
        """
        for (index_domain, domain) in enumerate(self.list_domain):
            list_abs_by_exciton = np.transpose(domain.abs_and_fl_exc_by_temp(temp)[1])
            list_dipole_by_exciton = domain.dipole_exciton
            for (index_exciton, (abs_exciton, dipole_exciton)) in enumerate(
                    zip(list_abs_by_exciton, list_dipole_by_exciton)):
                if (index_domain, index_exciton) == (0, 0):
                    absorption_t = np.linalg.norm(dipole_exciton) ** 2 * abs_exciton
                else:
                    absorption_t += np.linalg.norm(dipole_exciton) ** 2 * abs_exciton

            if domain._explicit_vib:
                absorption_t += domain.calculate_abs_vib_wdipole(temp)

        return np.fft.ifftshift(self.w_axis), np.real(np.fft.ifftshift(1 / (2 * np.pi) * np.fft.ifft(absorption_t)))

    def calc_linear_dichroism(self, temp, reference_direction):
        """Calculate the linear dichroism (LD) spectrum.
        
        Computes the LD spectrum based on the orientation of transition dipole
        moments relative to a reference direction. The LD signal is proportional
        to (1 - 3cos²θ), where θ is the angle between the dipole and the reference.
        
        Args:
            temp (float): Temperature in Kelvin.
            reference_direction (np.ndarray): Reference direction vector (e.g., 
                membrane normal or light polarization direction).
                
        Returns:
            tuple: A tuple containing:
                - np.ndarray: Frequency axis.
                - np.ndarray: Linear dichroism spectrum intensity.
        """
        LD_spectrum = None

        for (index_domain, domain) in enumerate(self.list_domain):
            list_abs_by_exciton = np.transpose(domain.abs_and_fl_exc_by_temp(temp)[1])
            list_dipole_by_exciton = domain.dipole_exciton

            for (index_exciton, (abs_exciton, dipole_exciton)) in enumerate(
                    zip(list_abs_by_exciton, list_dipole_by_exciton)):
                # Calculate the angle between the dipole and the reference direction
                theta = self.calculate_angle(dipole_exciton, reference_direction)

                # Calculate the orientation-dependent factor
                orientation_factor = (1 - 3 * np.cos(theta) ** 2)

                # Calculate the contribution to the LD spectrum
                contribution = np.linalg.norm(dipole_exciton) ** 2 * orientation_factor * abs_exciton

                # Initialize or accumulate the LD spectrum
                if LD_spectrum is None:
                    LD_spectrum = contribution
                else:
                    LD_spectrum += contribution

                # Include vibrational contributions if explicitly modeled
            if domain._explicit_vib:
                LD_spectrum += domain.calculate_abs_vib_wdipole_for_LD(temp, orientation_factor)

        # Return the frequency axis and the LD spectrum calculated via Fourier transform
        return np.fft.ifftshift(self.w_axis), np.real(np.fft.ifftshift(1 / (2 * np.pi) * np.fft.ifft(LD_spectrum)))
    
    
    def calculate_angle(self, dipole_exciton, reference_direction):
        """Calculate the angle between an exciton dipole and a reference direction.
        
        Args:
            dipole_exciton (np.ndarray): Excitonic transition dipole moment vector.
            reference_direction (np.ndarray): Reference direction vector.
            
        Returns:
            float: Angle θ in radians between the two vectors.
        """
        # Implement the calculation of the angle theta between the dipole_exciton and the reference_direction.
        # One possible implementation if both vectors are numpy arrays:
        cos_theta = np.dot(dipole_exciton, reference_direction) / (
                np.linalg.norm(dipole_exciton) * np.linalg.norm(reference_direction))
        theta = np.arccos(cos_theta)
        return theta

    def calc_fluorescence(self, temp, scale=None):
        """Calculate the fluorescence spectrum.
        
        Computes the thermally-weighted fluorescence spectrum by summing
        contributions from all excitonic states according to their Boltzmann
        populations at the given temperature.
        
        Args:
            temp (float): Temperature in Kelvin.
            scale (list, optional): List of scaling factors for each domain's
                thermal population. If None, raw thermal populations are used.
                
        Returns:
            tuple: A tuple containing:
                - np.ndarray: Frequency axis.
                - np.ndarray: Fluorescence spectrum intensity.
        """
        partition_function = 0
        fluorescence_t = np.zeros(len(self.list_domain[0].list_pigments[0].t_axis),
                                  dtype=np.complex128)
        for (index_domain, domain) in enumerate(self.list_domain):
            list_fl_by_exciton = np.transpose(domain.abs_and_fl_exc_by_temp(temp)[2])
            if domain._explicit_vib:
                list_fl_vib_by_exciton = np.transpose(domain.calculate_fl_vib_wdipole(temp))
            else:
                list_fl_vib_by_exciton = np.zeros_like(list_fl_by_exciton)
            list_dipole_by_exciton = np.array(domain.dipole_exciton, dtype=np.complex128)
            if scale is None:
                list_raw_thermal_pop = np.array(domain._raw_thermal_pop(temp), dtype=np.complex128)
            else:
                list_raw_thermal_pop = np.array(scale[index_domain]*domain.thermal_pop(temp))
            if np.sum(np.abs(list_dipole_by_exciton)) > 0:
                fluorescence_t, partition_function = numba_calc_fluorescence_contributions(list_fl_by_exciton,
                                                                                           list_dipole_by_exciton,
                                                                                           list_raw_thermal_pop,
                                                                                           list_fl_vib_by_exciton,
                                                                                           partition_function,
                                                                                           fluorescence_t,
                                                                                           domain._explicit_vib)

        if scale is not None: partition_function = 1
        return np.fft.ifftshift(self.w_axis), np.real(
            np.fft.ifftshift(1 / (2 * np.pi) * np.fft.ifft(fluorescence_t / partition_function)))
    
@numba.njit()
def numba_calc_fluorescence_contributions(list_fl_by_exciton,
                                          list_dipole_by_exciton,
                                          list_raw_thermal_pop,
                                          list_fl_vib_by_exciton,
                                          partition_function,
                                          fluorescence_t,
                                          flag_explicit_vib):
    """Calculate fluorescence contributions using Just-In-Time compilation (JIT compilation).
    
    This is a Numba-accelerated helper function that computes fluorescence
    contributions from excitonic states, including electronic and vibrational
    components.
    
    Args:
        list_fl_by_exciton (np.ndarray): Fluorescence lineshape for each exciton.
        list_dipole_by_exciton (np.ndarray): Transition dipole moments for each exciton.
        list_raw_thermal_pop (np.ndarray): Thermal populations for each exciton.
        list_fl_vib_by_exciton (np.ndarray): Vibrational fluorescence contributions.
        partition_function (float): Current partition function value.
        fluorescence_t (np.ndarray): Accumulated fluorescence in time domain.
        flag_explicit_vib (bool): Whether to include explicit vibrational contributions.
        
    Returns:
        tuple: A tuple containing:
            - np.ndarray: Updated fluorescence in time domain.
            - float: Updated partition function.
    """

    for (index_exciton, (fl_exciton, dipole_exciton, thermal_pop)) in enumerate(zip(list_fl_by_exciton,
                                                                                    list_dipole_by_exciton,
                                                                                    list_raw_thermal_pop)):
        partition_function += thermal_pop
        intermediate = thermal_pop * np.linalg.norm(dipole_exciton) ** 2 * fl_exciton
        fluorescence_t += intermediate
        if flag_explicit_vib:
            fluorescence_t = fluorescence_t + thermal_pop * list_fl_vib_by_exciton[:, index_exciton]

    return fluorescence_t, partition_function