import numpy as np

from pymembrane.exciton.rate_matrix import RateMatrix

class NonlinearSpectra(RateMatrix):
    def __init__(self, structure, temp, r_neighbor=70, k_ij=None):
        super().__init__(structure,
                         temp,
                         r_neighbor,
                         k_ij)

    def calculate_time_resolved_fluorescence(self, P1_0, N_t, dt, N_save, N_resolve=None):
        """
        This method calculates the time-resolved fluorescence profile by weighting the
        thermalized domain fluorescence by their population.
        """
        if N_resolve is None:
            N_resolve = N_save
        t_axis, P2_pop_t = self.time_evolve(P1_0, N_t+N_resolve, dt)

        list_fl_t = []
        list_t = []
        w_axis = []
        for index_save in np.arange(0,N_t,N_save):
            list_t.append(t_axis[index_save])
            P1_pop_t = np.sum(P2_pop_t[index_save:(index_save+N_resolve)], axis=0)
            w_axis, fl_t = self.calc_fluorescence(self.temp, scale=P1_pop_t)
            list_fl_t.append(fl_t)

        return list_t, w_axis, list_fl_t

    def calculate_tcspc(self, P1_0, N_t, dt, N_save, N_resolve=None, method_fast=True):
        """Calculate time-correlated single photon counting (TCSPC) signal.
        
        Computes the TCSPC decay trace by integrating the fluorescence intensity
        over time bins. Two methods are available: a full spectral integration
        or a fast method that uses thermally-averaged transition dipole moments.
        
        Args:
            P1_0 (np.ndarray): Initial population distribution across domains.
            N_t (int): Total number of time steps for population evolution.
            dt (float): Time step size for integration.
            N_save (int): Interval between saved time bins (in number of steps).
            N_resolve (int, optional): Number of time steps to integrate over for 
                each bin to calculate the thermalized population. If None, defaults 
                to N_save.
            method_fast (bool, optional): If True, uses the fast method based on 
                thermal transition dipole moments. If False, uses full spectral 
                integration. Defaults to True.
        
        Returns:
            tuple: A tuple containing:
                - list_t (list): List of time values at the start of each bin.
                - tcspc (np.ndarray): Array of integrated photon counts for each time bin.
        
        Notes:
            The fast method computes TCSPC using the weighted sum of thermal transition
            dipole moment squared (thermal_mu2) for each domain, which is computationally
            more efficient than integrating the full fluorescence spectrum.
        """
        if not method_fast:
            list_t, w_axis, list_fl_t = self.calculate_time_resolved_fluorescence(P1_0, N_t, dt, N_save, N_resolve)
            tcspc = np.sum(list_fl_t, axis=1)
        else:
            if N_resolve is None:
                N_resolve = N_save
            t_axis, P2_pop_t = self.time_evolve(P1_0, N_t + N_resolve, dt)
            M1_mu_them = [domain.thermal_mu2(self.temp) for domain in self.list_domain]
            tcspc = []
            list_t = []
            for index_save in np.arange(0,N_t,N_save):
                list_t.append(t_axis[index_save])
                P1_pop_t = np.sum(P2_pop_t[index_save:(index_save + N_resolve)], axis=0)
                tcspc.append(np.sum(P1_pop_t*M1_mu_them))

        return list_t, np.array(tcspc)