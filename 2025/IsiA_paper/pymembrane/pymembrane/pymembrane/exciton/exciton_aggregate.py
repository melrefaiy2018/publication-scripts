import copy
from functools import cached_property

import numba
import numpy as np

from pymembrane.exciton.lineshape import spectraldensity_renger
from pymembrane.util.physical_constants import hbar, kB


@numba.njit()
def n_w(w_axis, temp):
    return 1 / (np.exp(w_axis / (kB * temp)) - 1)


class ExcitonDomain:

    def __init__(self, list_pigments, H2_domain):
        """
        This is a class that represents a domain of pigments that are electronically coupled.

        Parameters
        ----------
        list_pigments : list
            A list of PigmentAtomic objects that are in this domain.
        H2_domain : numpy.ndarray
            The Hamiltonian of the domain.
        """
        self.list_pigments = list_pigments
        self.list_names = [pigment.name for pigment in self.list_pigments]
        self.__rc = np.any(['PHO' in name for name in self.list_names])
        self._explicit_vib = False
        if any([pigment._explicit_vib_s is not None for pigment in self.list_pigments]):
            self._explicit_vib = True

        # Construct the exciton states for the domain
        # * NOte: self._U2_site_exc[:,i] is the ith exciton.
        self.__H2_domain = H2_domain
        self._list_e_exc, self._U2_site_exc = np.linalg.eigh(H2_domain)

        # Calculate Exciton Reorganization Energy
        self._e_lambda_exciton = np.dot(self.U2_exc_site ** 4, [spectra_pigment.e_lambda
                                                                for spectra_pigment
                                                                in self.list_pigments])

        self._dipole_exciton = None
        self._dict_abs_fl_by_temp = {}

        self._dict_therm_pop_by_temp = {}
        self._dict_tau_by_temp = {}

        self._n_cla = np.sum([1 for name in self.list_names if 'CLA' in name])
        self._n_clb = np.sum([1 for name in self.list_names if 'CHL' in name])

    def __repr__(self):
        return f'Kubo Domain ({self.name}): \n {[pigment.name for pigment in self.list_pigments]}'

    def __len__(self):
        return len(self.list_pigments)

    @property
    def domain_type(self):
        return 'Exciton'

    @property
    def location(self):
        return np.mean([pig.location for pig in self.list_pigments], axis=0)

    @cached_property
    def name(self):
        return f'{list(set([pig.name.split("_")[0] for pig in self.list_pigments]))}_{self.list_pigments[0]}'

    def thermal_pop(self, temp):
        if temp not in self._dict_therm_pop_by_temp.keys():
            therm_pop = self._raw_thermal_pop(temp)
            therm_pop = np.array(therm_pop)/np.sum(therm_pop)
            self._dict_therm_pop_by_temp[temp] = therm_pop

        return self._dict_therm_pop_by_temp[temp]

    def thermal_mu2(self, temp):
        return np.sum([p_therm*np.linalg.norm(dipole)**2 for (p_therm, dipole) in zip(self.thermal_pop(temp),
                                                                                  self.dipole_exciton)])

    def _raw_thermal_pop(self, temp):
        return [np.exp(-(energy-e_lambda)/(kB*temp)) for (energy, e_lambda) in zip(self._list_e_exc, self.e_lambda_exciton)]

    def G2_gt_site(self, temp):
        # This will construct a [N_pig, N_time] array
        return np.array([spectra_pigment.g_t(temp) for spectra_pigment in self.list_pigments])

    @property
    def e_lambda_exciton(self):
        return self._e_lambda_exciton

    def calculate_lifetime_by_temp(self, temp):
        # Calculate lifetime (assumes all spectral densities are the same)
        # Calculate Lifetime Broadening
        # -----------------------------
        # I have set the correlation distance to 0 in order to simplify the calculation
        # of lifetimes.
        renger_jw_param = copy.deepcopy(self.list_pigments[0]._lineshape.renger_Jw_param)
        J_w = spectraldensity_renger
        list_tau = []
        for index_exciton in np.arange(len(self._list_e_exc)):
            tau_exciton = 0
            for index_exciton_N in np.arange(len(self._list_e_exc)):
                if index_exciton_N != index_exciton:
                    omega_MN = self._list_e_exc[index_exciton] - \
                               self._list_e_exc[index_exciton_N]
                    gamma_MN = np.dot(
                        np.abs(self._U2_site_exc[:, index_exciton]) ** 2,
                        np.abs(self._U2_site_exc[:, index_exciton_N]) ** 2)
                    if omega_MN > 0:
                        renger_jw_param['w_axis'] = np.array([omega_MN])
                        tau_exciton += gamma_MN * (
                                1 + n_w(omega_MN, temp)) * J_w(**renger_jw_param)[0]
                    else:
                        renger_jw_param['w_axis'] = np.array([-omega_MN])
                        tau_exciton += gamma_MN * (
                            n_w(-omega_MN, temp)) * J_w(**renger_jw_param)[0]
            list_tau.append(np.max([hbar/1000, tau_exciton]))
        return list_tau

    def abs_and_fl_exc_by_temp(self, temp):
        if not (temp in self._dict_abs_fl_by_temp.keys()):
            self._dict_abs_fl_by_temp[temp] = self.calculate_abs_fl_exc(temp)

        return self._dict_abs_fl_by_temp[temp]

    def calculate_abs_fl_exc(self, temp):
        t_axis = self.list_pigments[0].t_axis
        T2_array = np.transpose(np.tile(t_axis, [self.n_exc, 1]))
        coherent_phase = self._list_e_exc/hbar*T2_array
        reorganization_correction = 2*self.e_lambda_exciton/hbar*T2_array
        G2_gt_exciton = self.G2_gt_exciton(temp)

        # NOTE: The code below assumes that if any pigment has renger's lineshape in this
        # domain, all pigments have Renger's lineshape in this domain.
        renger_gamma_t = np.zeros_like(G2_gt_exciton)
        if self.list_pigments[0]._lineshape.renger_gamma:
            t_axis = self.list_pigments[0].t_axis
            list_tau = self.calculate_lifetime_by_temp(temp)
            for (index, tau) in enumerate(list_tau):
                #TODO: Fix naming convention
                #      (list_tau is really gamma since it is in units of energy)
                renger_gamma_t[:,index] = (t_axis*tau/hbar)

        spectrum_abs_t = np.exp(-1j*coherent_phase - G2_gt_exciton - renger_gamma_t)
        spectrum_fl_t  = np.exp(-1j*(coherent_phase - reorganization_correction) - np.conj(G2_gt_exciton)- renger_gamma_t)
        return (t_axis, spectrum_abs_t, spectrum_fl_t)

    def calculate_abs_vib_wdipole(self, temp):
        t_axis = self.list_pigments[0].t_axis
        absorption_t = np.zeros(len(t_axis), dtype=np.complex128)
        for (E_n, pigment) in zip(np.diag(self.H2_domain), self.list_pigments):
            dipole_pig = pigment.dipole
            if pigment.explicit_vib:
                abs_vib_t = pigment.vibronic_absorption(temp)
                absorption_t += np.linalg.norm(dipole_pig) ** 2 * np.exp(-1j*E_n/hbar*t_axis)*abs_vib_t

        return absorption_t

    def calculate_fl_vib_wdipole(self, temp):
        t_axis = self.list_pigments[0].t_axis
        n_exc = len(self._list_e_exc)
        fluorescence_t = np.zeros([n_exc, len(t_axis)], dtype=np.complex128)
        for index_exc in np.arange(n_exc):
            for (index_site, pigment) in enumerate(self.list_pigments):
                dipole_pig = pigment.dipole

                if pigment.explicit_vib:
                    E_exc = self._list_e_exc[index_exc]
                    lambda_exc = self.e_lambda_exciton[index_exc]
                    c_site_exc = self._U2_site_exc[index_site, index_exc]

                    coherent_phase = E_exc/hbar * t_axis
                    reorganization_correction = (lambda_exc) / hbar * t_axis
                    fl_t = pigment.vibronic_fluorescence(temp)
                    fluorescence_t[index_exc, :] += (np.linalg.norm(dipole_pig) ** 2 * (c_site_exc) ** 2
                                                    * np.exp(-1j * (coherent_phase - reorganization_correction))
                                                    * fl_t)

        return fluorescence_t

    def G2_gt_exciton(self, temp):
        # Lineshape Function: G_IIII(t) = \sum_n |U_In|^4 g_n(t)
        # Returns an object that is [N_exc, N_time] array
        G2_iiii_t = self.U2_exc_site ** 4 @ self.G2_gt_site(temp)
        return np.transpose(G2_iiii_t)

    @property
    def n_exc(self):
        return len(self._list_e_exc)

    @property
    def U2_exc_site(self):
        return np.transpose(self._U2_site_exc)

    @property
    def U2_site_exc(self):
        return self._U2_site_exc

    @property
    def dipole_exciton(self):
        if self._dipole_exciton is None:
            self._dipole_exciton = np.dot(self.U2_exc_site, self.dipole_site)

        return self._dipole_exciton

    @property
    def dipole_site(self):
        return np.array([pigment.dict_data['dipole_mag'] * pigment.get_dipole_dir()
                         for pigment in self.list_pigments])

    @property
    def H2_domain(self):
        return self.__H2_domain

    @property
    def is_rc(self):
        return self.__rc


class ExcitonicAggregate:

    def __init__(self, structure):
        self.structure = structure
        self.dict_pigments = structure.dict_pigments

        # Construct lists of pigments by domain
        self.dict_pigments_by_domain_index = {}
        for (key, pigment) in self.dict_pigments.items():
            if pigment.domain_index in self.dict_pigments_by_domain_index.keys():
                self.dict_pigments_by_domain_index[pigment.domain_index].append(pigment)
            else:
                self.dict_pigments_by_domain_index[pigment.domain_index] = [pigment]

        # Construct Spectral Domain objects
        self.dict_domain_by_index = {}
        for (index_domain, list_pig) in self.dict_pigments_by_domain_index.items():
            list_names = [pigment.name for pigment in list_pig]
            H2_domain = self.H2_hamiltonian.subset_by_name(list_names, list_names)
            self.dict_domain_by_index[index_domain] = ExcitonDomain(list_pig, H2_domain)

        # Construct Domain List
        # NOTE: For this to work correctly the domain indices need to be a contiguous set of
        # integers starting from 0.
        # TODO: This should either be checked or enforced.
        list_domain_index = list(self.dict_domain_by_index.keys())
        list_domain_index.sort()
        self.list_domain = [self.dict_domain_by_index[index] for index in list_domain_index]

    @property
    def H2_hamiltonian(self):
        return self.structure.H2_hamiltonian

    @property
    def n_dom(self):
        return len(self.list_domain)
