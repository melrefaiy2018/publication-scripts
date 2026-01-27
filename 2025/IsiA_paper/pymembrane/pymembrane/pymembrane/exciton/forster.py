"""Exciton transfer routines based on generalized Forster theory.

This module provides routines for computing exciton transfer rates
between donor and acceptor domains using a generalized Forster
formalism. The functions here expect domain-like objects from the
package and a site-basis Hamiltonian; they return either a matrix of
exciton-to-exciton rates or a thermalized scalar rate.

Notes:
    - Domain-like objects are expected to expose attributes and methods
      used throughout the package (for example: ``list_pigments``,
      ``U2_exc_site``, ``U2_site_exc``, ``n_exc``,
      ``abs_and_fl_exc_by_temp``, ``thermal_pop`` and a boolean
      ``_explicit_vib``).
"""

import copy

import numpy as np
from typing import Any

from pymembrane.util.physical_constants import hbar


def generalized_forster_exciton(domain_acceptor: Any, domain_donor: Any, temp: float, H2_hamiltonian: Any) -> np.ndarray:
    """Compute exciton-to-exciton generalized Forster rates.

    Calculates the matrix of rates ``k_{a<-d}`` from each donor exciton
    ``d`` to each acceptor exciton ``a`` using the time-domain
    generalized Forster expression:

        k_{a<-d} = |V_{a,d}|^2 / hbar^2 * \\int dt A_a(t) F_d^*(t)

    Vibronic correction terms are added when either domain declares
    ``_explicit_vib``; these corrections are accumulated in the
    ``J2_ad`` (acceptor vibronic) and ``L2_ad`` (donor vibronic)
    contributions.

    Args:
        domain_acceptor: Acceptor-domain object. Expected attributes/methods
            include: ``list_pigments``, ``U2_exc_site``, ``U2_site_exc``,
            ``n_exc``, ``abs_and_fl_exc_by_temp``, ``thermal_pop``, and
            ``_explicit_vib`` when vibronic states exist.
        domain_donor: Donor-domain object with the same interface as the
            acceptor.
        temp: Temperature in Kelvin used to evaluate lineshapes.
        H2_hamiltonian: Site-basis Hamiltonian object supporting
            ``subset_by_name`` to extract donor-acceptor couplings.

    Returns:
        numpy.ndarray: Array of shape (n_acceptor_excitons, n_donor_excitons)
            containing exciton-to-exciton rates ``k_{a<-d}``.

    Raises:
        AttributeError: If the provided domain objects do not implement the
            expected attributes or methods.

    Notes:
        - Internal contributions: ``I2_ad`` (coherent overlap),
          ``J2_ad`` (acceptor vibronic corrections), ``L2_ad`` (donor
          vibronic corrections).
        - Time integrals are computed using numpy's ``trapz`` over the
          time axis returned by the domains' ``abs_and_fl_exc_by_temp``.

    Example:
        >>> K = generalized_forster_exciton(acc_domain, don_domain, 300.0, H2)
        >>> K.shape
        (acc_domain.n_exc, don_domain.n_exc)
    """
    # Construct Site Basis Hamiltonian
    # --------------------------------
    H2_site_ad = H2_hamiltonian.subset_by_name(domain_acceptor.list_names,
                                               domain_donor.list_names)

    # Construct Exciton Basis Hamiltonian
    # -----------------------------------
    H2_site_ad_fc = copy.copy(H2_site_ad)
    for (index_donor, site_donor) in enumerate(domain_donor.list_pigments):
        for (index_acc, site_acceptor) in enumerate(domain_acceptor.list_pigments):
            H2_site_ad_fc[index_acc, index_donor] = H2_site_ad[index_acc, index_donor]*site_donor.fc_00*site_acceptor.fc_00
    H2_exc_ad = domain_acceptor.U2_exc_site@H2_site_ad_fc@domain_donor.U2_site_exc

    if domain_acceptor._explicit_vib:
        # If we are including explicit vibrations, then can precalculate some of the vibronic
        # Hamiltonian terms.
        H2_avib_dexc = np.zeros_like(H2_exc_ad)
        for alpha_exc in np.arange(domain_donor.n_exc):
            for n_site in np.arange(domain_acceptor.n_exc):
                # H2_avib_dexc[n_site, alpha_exc] += np.sum(domain_donor.U2_site_exc[:, alpha_exc]
                #                                            * H2_site_ad_fc[n_site, :])
                H2_avib_dexc[n_site, alpha_exc] += np.sum(domain_donor.U2_site_exc[:, alpha_exc]
                                                           * H2_site_ad[n_site, :]
                                                           * np.array([domain_donor.list_pigments[m_site].fc_00
                                                                         for m_site in np.arange(domain_donor.n_exc)]))

    if domain_donor._explicit_vib:
        H2_aexc_dvib = np.zeros_like(H2_exc_ad)
        for m_site in np.arange(domain_donor.n_exc):
            for beta_exc in np.arange(domain_acceptor.n_exc):
                # H2_aexc_dvib[beta_exc, m_site] += np.sum(domain_acceptor.U2_site_exc[:, beta_exc]
                #                                          * H2_site_ad_fc[:,m_site])
                H2_aexc_dvib[beta_exc, m_site] += np.sum(domain_acceptor.U2_site_exc[:, beta_exc]
                                                         * H2_site_ad[:,m_site]
                                                         * np.array([domain_acceptor.list_pigments[m_site].fc_00
                                                                         for m_site in np.arange(domain_acceptor.n_exc)]))

    # Construct the overlap integrals
    # -------------------------------
    I2_ad = np.zeros([domain_acceptor.n_exc, domain_donor.n_exc])
    J2_ad = np.zeros([domain_acceptor.n_exc, domain_donor.n_exc])
    L2_ad = np.zeros([domain_acceptor.n_exc, domain_donor.n_exc])
    t_axis, A2_donor, F2_donor = domain_donor.abs_and_fl_exc_by_temp(temp)
    A2_acceptor = domain_acceptor.abs_and_fl_exc_by_temp(temp)[1]

    # The block of code below might be replaced by a numba function
    for index_donor in np.arange(domain_donor.n_exc):
        for index_acceptor in np.arange(domain_acceptor.n_exc):
            I2_ad[index_acceptor, index_donor] = 2*np.real(np.trapz(np.conj(F2_donor[:, index_donor])
                                                                    *A2_acceptor[:, index_acceptor],
                                                                    x=t_axis))
            #print(f'exciton {index_acceptor}<--exciton {index_donor} overlap: {I2_ad[index_acceptor, index_donor]*np.abs(H2_exc_ad[index_acceptor, index_donor])**2/hbar**2}')

            # If the acceptor domain has explicit vibrations, then we will add the correction term
            # that accounts for vibronic <-- exciton transitions.
            # 2/hbar**2 * SUM_n |c_n^a|^2 * |H2_exc_vib[exc, n_site]|^2 * SUM_\nu FC_n(0,\nu) * [integral over lineshapes]
            if domain_acceptor._explicit_vib:
                for (n_site, site_energy) in enumerate(np.diag(domain_acceptor.H2_domain)):
                    contribution = 2 * np.abs(domain_acceptor.U2_site_exc[n_site,index_acceptor])**2
                    contribution *= (H2_avib_dexc[n_site, index_donor])**2/hbar**2
                    if contribution == 0:
                        contribution = 0
                    else:
                        contribution *= np.real(np.trapz(np.conj(F2_donor[:, index_donor])
                                                        * domain_acceptor.list_pigments[n_site].vibronic_absorption(temp)
                                                        * np.exp(-1j*site_energy/hbar*t_axis),
                                                        x=t_axis))

                        #print(f'exciton {index_acceptor} vibron {n_site} <--exciton {index_donor} overlap: {contribution}')
                    J2_ad[index_acceptor, index_donor] += contribution

            if domain_donor._explicit_vib:
                for m_site in np.arange(domain_donor.n_exc):
                    energy_exc = domain_donor._list_e_exc[index_donor]
                    lambda_exc = domain_donor.e_lambda_exciton[index_donor]
                    coherent_phase = energy_exc/hbar * t_axis
                    reorganization_correction = (lambda_exc) / hbar * t_axis

                    contribution = 2 * np.abs(domain_donor.U2_site_exc[m_site,index_donor])**2
                    contribution *= (H2_aexc_dvib[index_acceptor, m_site])**2/hbar**2
                    if contribution == 0:
                        contribution = 0
                    else:
                        contribution *= np.real(np.trapz(np.conj(domain_donor.list_pigments[m_site].vibronic_fluorescence(temp)
                                                             * np.exp(-1j * (coherent_phase - reorganization_correction)))
                                                     * A2_acceptor[:,index_acceptor],
                                                     x=t_axis))
                    L2_ad[index_acceptor, index_donor] += contribution

                    #print(f'exciton {index_acceptor} <--exciton {index_donor} vibron {m_site} overlap: {contribution}')

    # Calculate rates of transport (Aexc<--Dexc)
    # ------------------------------------------
    K2_ad = I2_ad * np.abs(H2_exc_ad)**2/hbar**2 + J2_ad + L2_ad

    return K2_ad


def generalized_forster_thermal(domain_acceptor: Any, domain_donor: Any, temp: float, H2_hamiltonian: Any) -> float:
    """Compute the donor->acceptor thermalized generalized Forster rate.

    This function first computes the exciton-to-exciton rate matrix via
    :func:`generalized_forster_exciton` and then constructs the thermalized
    scalar rate by weighting donor exciton rates with the donor thermal
    population P(d):

        K_{A<-D} = \\sum_{d in D} P(d) \\sum_{a in A} k_{a<-d}

    Args:
        domain_acceptor: Acceptor-domain object (see
            ``generalized_forster_exciton``).
        domain_donor: Donor-domain object (see ``generalized_forster_exciton``).
        temp: Temperature in Kelvin used for lineshapes and thermal populations.
        H2_hamiltonian: Site-basis Hamiltonian used to extract donor-acceptor
            couplings.

    Returns:
        float: Thermalized donor->acceptor rate (scalar). Units match the
            project's convention (typically s^-1 when SI units are used).

    Notes:
        - This routine will call ``domain_donor.thermal_pop(temp)``; the
          returned population is expected to sum to 1. If not, the output will
          scale accordingly.

    Example:
        >>> k_dom = generalized_forster_thermal(acc_domain, don_domain, 300.0, H2)
    """
    # Calculate rates of transport (Aexc<--Dexc)
    # ------------------------------------------
    K2_ad = generalized_forster_exciton(domain_acceptor, domain_donor, temp, H2_hamiltonian)

    # Construct thermalized rates (Adom<--Ddom)
    # -----------------------------------------
    p_therm = domain_donor.thermal_pop(temp)
    kdom_ad = 0

    # The block of code below should be replaced by a numba function
    for index_donor in np.arange(domain_donor.n_exc):
        for index_acceptor in np.arange(domain_acceptor.n_exc):
            kdom_ad += K2_ad[index_acceptor, index_donor]*p_therm[index_donor]

    return kdom_ad

