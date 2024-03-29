"""
Fitter for Carrier and Sideband Rabi Flopping
"""
from motional_distribution import motional_distribution as md
import numpy as np
from scipy.special.orthogonal import eval_genlaguerre as laguerre
import warnings
class rabi_flop_time_evolution(object):
    def __init__(self, sideband_order, eta, nmax = 1000):
        self.sideband_order = sideband_order # 0 for carrier, 1 for 1st sideband etc
        self.eta = eta
        self.nmax = nmax
        self.rabi_coupling = self.compute_rabi_coupling()
    def compute_rabi_coupling(self):
        '''
        Rabi couplings, see Leibfried (2003), eq:70
        '''
        eta = self.eta
        if self.sideband_order == 0:
            coupling_func = lambda n: np.exp(-1./2*eta**2) * laguerre(n, 0, eta**2)
        elif self.sideband_order == 1:
            coupling_func = lambda n: np.exp(-1./2*eta**2) * eta**(1)*(1./(n+1.))**0.5 * laguerre(n, 1, eta**2)
        elif self.sideband_order == 2:
            coupling_func = lambda n: np.exp(-1./2*eta**2) * eta**(2)*(1./((n+1.)*(n+2)))**0.5 * laguerre(n, 2, eta**2)
        elif self.sideband_order == 3:
            coupling_func = lambda n: np.exp(-1./2*eta**2) * eta**(3)*(1./((n+1)*(n+2)*(n+3)))**0.5 * laguerre(n, 3 , eta**2) 
        elif self.sideband_order == 4:
            coupling_func = lambda n: np.exp(-1./2*eta**2) * eta**(4)*(1./((n+1)*(n+2)*(n+3)*(n+4)))**0.5 * laguerre(n, 4 , eta**2)
        elif self.sideband_order == 5:
            coupling_func = lambda n: np.exp(-1./2*eta**2) * eta**(5)*(1./((n+1)*(n+2)*(n+3)*(n+4)*(n+5)))**0.5 * laguerre(n, 5 , eta**2)      
        elif self.sideband_order == -1:
            coupling_func = lambda n: 0 if n == 0 else np.exp(-1./2*eta**2) * eta**(1)*(1./(n))**0.5 * laguerre(n - 1, 1, eta**2)
        elif self.sideband_order == -2:
            coupling_func = lambda n: 0 if n <= 1 else np.exp(-1./2*eta**2) * eta**(2)*(1./((n)*(n-1.)))**0.5 * laguerre(n - 2, 2, eta**2)
        elif self.sideband_order == -3:
            coupling_func = lambda n: 0 if n <= 2 else np.exp(-1./2*eta**2) * eta**(3)*(1./((n)*(n-1.)*(n-2)))**0.5 * laguerre(n -3, 3, eta**2)
        elif self.sideband_order == -4:
            coupling_func = lambda n: 0 if n <= 3 else np.exp(-1./2*eta**2) * eta**(4)*(1./((n)*(n-1.)*(n-2)*(n-3)))**0.5 * laguerre(n -4, 4, eta**2)
        elif self.sideband_order == -5:
            coupling_func = lambda n: 0 if n <= 4 else np.exp(-1./2*eta**2) * eta**(5)*(1./((n)*(n-1.)*(n-2)*(n-3)*(n-4)))**0.5 * laguerre(n -5, 5, eta**2)
        else:
            raise NotImplementedError("Can't calculate rabi couplings sideband order {}".format(self.sideband_order))
        return np.array([coupling_func(n) for n in range(self.nmax)])
        
    def compute_evolution_thermal(self, nbar, delta, time_2pi, t, excitation_scaling = 1.):
        '''returns the state evolution for temperature nbar, detuning delta, rabi frequency T_Rabi for times t'''
        print 'Compute Evolution Thermal'
        omega = self.rabi_coupling
        ones = np.ones_like(t)
        p_n = md.thermal(nbar, self.nmax)
        if 1 - p_n.sum() > 1e-6:
            raise Exception ('Hilbert space too small, missing population')
        if delta == 0:
            #prevents division by zero if delta == 0, omega == 0
            effective_omega = 1
        else:
            effective_omega = omega/np.sqrt(omega**2 + delta**2)
        result = np.outer(p_n * effective_omega, ones) * (np.sin( np.outer( np.sqrt(omega**2+delta**2)*np.pi/time_2pi, t ))**2)
        result = np.sum(result, axis = 0)
        result = excitation_scaling * result
        return result
    
    def compute_evolution_coherent(self, nbar, alpha, delta, time_2pi, t, excitation_scaling = 1.):
        '''returns the state evolution for temperature nbar, detuning delta, rabi frequency T_Rabi for times t'''
        print 'Compute Evolution Coherent'
        omega = self.rabi_coupling
        ones = np.ones_like(t)
        p_n = md.displaced_thermal(alpha, nbar, self.nmax)
        print 'got array'
        if 1 - p_n.sum() > 1e-6:
            raise Exception ('Hilbert space too small, missing population')
        if delta == 0:
            #prevents division by zero if delta == 0, omega == 0
            effective_omega = 1
        else:
            effective_omega = np.abs(omega)/np.sqrt(omega**2+delta**2)
        print 'looking for result'
        # This equation was fixed from the previous version 
        result = np.outer(p_n * effective_omega, ones) * (np.sin( np.outer( np.sqrt(omega**2+delta**2)*np.pi/time_2pi, t ))**2)
        print 'summing'
        result = np.sum(result, axis = 0)
        print 'scaling'
        result = excitation_scaling * result
        print 'done'
	return result
    
if __name__ == '__main__':
    from matplotlib import pyplot
    def thermal_example():
        eta = 0.05
        nbar = 3
        times = np.linspace(0, 300, 100)
        te = rabi_flop_time_evolution(-1 ,eta)
        prob_red = te.compute_evolution_thermal(nbar = nbar, delta = 0, time_2pi = 15, t = times)
        te = rabi_flop_time_evolution(1 ,eta)
        prob_blue = te.compute_evolution_thermal(nbar = nbar, delta = 0, time_2pi = 15, t = times)
        te = rabi_flop_time_evolution(0 ,eta)
        prob_car = te.compute_evolution_thermal(nbar = nbar, delta = 0, time_2pi = 15, t = times)
        pyplot.plot(times, prob_red, label = 'red sideband')
        pyplot.plot(times, prob_blue, label = 'carrier')
        pyplot.plot(times, prob_car, label = 'blue sideband')
        pyplot.legend()
        pyplot.show()
        #extracting back the temperature from the ratio of sidebands
        ratio_nbar = prob_red[1:] /(prob_blue[1:] - prob_red[1:])
        print 'nbar from sideband ratio', ratio_nbar[0]
    
    def displaced_thermal_example():
        eta = 0.05
        for alpha in np.linspace(0, 25, 6):
            times = np.linspace(0, 50, 300)
            pyplot.subplot(211)
            te = rabi_flop_time_evolution(-1 ,eta)
            prob_blue = te.compute_evolution_coherent(nbar = 3, alpha = alpha, delta = 0, time_2pi = 15, t = times)
            pyplot.plot(times, prob_blue, label = 'alpha a = {}'.format(alpha))
            pyplot.subplot(212)
            te = rabi_flop_time_evolution(-2 ,eta)
            prob_blue = te.compute_evolution_coherent(nbar = 3, alpha = alpha, delta = 0, time_2pi = 15, t = times)
            pyplot.plot(times, prob_blue, label = 'alpha a = {}'.format(alpha))
        pyplot.subplot(211)
        pyplot.title('First order red sideband', fontsize = 24)
        pyplot.subplot(212)
        pyplot.title('Second order red sideband', fontsize = 24)
        # Will error if run on Mac OSX
        # pyplot.tight_layout()
        pyplot.legend(prop={'size':16})
        pyplot.xlabel('Excitation Time (arb)', fontsize = 20)
        pyplot.ylabel('Excitation', fontsize = 20)
        pyplot.tick_params(axis='both', which='major', labelsize=14)
        pyplot.show()
    
    def display_rabi_coupling():
        eta = 0.05
        for sideband in [-1,-2,-3,-4,-5]:
            te = rabi_flop_time_evolution(sideband ,eta)
            couplings = te.compute_rabi_coupling()
            couplings = np.abs(couplings)
            pyplot.plot(couplings, label = 'sideband {}'.format(sideband))
        pyplot.title('Rabi Coupling Strength', fontsize = 24)
        pyplot.xlabel('n', fontsize = 24)
        pyplot.ylabel(u'$\Omega_n$ / $\Omega_0$', fontsize = 24)
        pyplot.legend()
        pyplot.show()

        
    thermal_example()  
    displaced_thermal_example()
    display_rabi_coupling()
    

