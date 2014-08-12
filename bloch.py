# from __future__ import division
from scipy.optimize import *
from get_data import *
import lmfit as lmfit 
from datafit import DataFit
import numpy as np
from scipy.integrate import odeint
import math
class Bloch(DataFit):
    name = 'Bloch'
    def __init__(self):
        DataFit.__init__(self)
        self.all_parameters=['t_excitation', 'height','gamma','ohm','shift']
        # t_excitation can be retrieved from [Parameter 10] in $PATH/LabRAD/cct/data/Experiments.dir/Spectrum729.dir/2014Jun17.dir/1903_14.dir
        # label = Spectrum.manual_excitation_time
        # data = Value(40.0, 'us')
        self.excitation_time=0
        self.guess_dict={ 't_excitation': lambda : 40e-6 , # in sec 
			'height': self.guess_height	,
            'shift': lambda : 0.0,
            'gamma': lambda :  1e-4 ,
            'ohm': lambda: 2*np.pi*1000 # in Hz
		}
        # Don't need to set data, if given data using the new datafit class (raw is passed in)
        #self.parameters = self.guess_dict
        # self.raw= raw
        # self.setData(raw)
    
    def guess(self, param):
        return self.guess_dict[param]()

    def guess_height(self):
        # max at zero corrsponds to multiplicative coeficient
        # height is a scale factor
        # Since the max is not 1 . Use a amplitude scale factor of (max in data)/ (max in model)
        print ("guess_height: {}".format(max(self.dataY)))
        return max(self.dataY)

    def guess_shift(self):
        #algorithm that finds the time (x value) corresponding to max height (ymax) 
        shift =0.0
        for i in range(len(self.dataX)):
            if (self.dataY[i]==max(self.dataY)):
                shift= self.dataX[i]
        return (shift)
        
    def model (self,params,x):
        print ("inside model")
        t_excitation = 40.0e-6
        t_excitation = params['t_excitation'].value
        print ("t_excitation: ".format(t_excitation))
        height = params['height'].value
        shift = params['shift'].value
        gamma = params['gamma'].value
        ohm = params['ohm'].value
        exc = []
        #step = math.floor((50*ohm)*2/len(x)-1) #Need to truncate this value
        step=(50*ohm)*2/(len(x))
        # Scan delta over the whole range but also want to keep the delta_list (which is the size of Y model data)
        delta_list = np.arange(-50*ohm, 50*ohm, step)
        print ("delta_list: {}".format(len(delta_list)))
        for delta_i in delta_list:
            print ("delta_i: {}".format(delta_i))
            # Solve for the atomic inversion R3, normalize it to get excitation probability to fit data 
            p_exc = height*((blochSolver(delta_i, gamma = gamma,T_0=t_excitation)+1)/2)
            exc.append(p_exc)
        # plt.plot(delta_list, exc)
        return exc

ohm = 2*np.pi*1000 #Hz
#T_1 is now called gamma. It is the parameter that determines damping on r3.
gamma = 1e-5 #Set T as infinite (large)
# We don't care about T_2, it is taken to be very small
#Initial Condition (In z direction, initially R_3= -1)
R0 = [0,0,-1]
T_2pi = 10. #excitation time 
def blochSolver(delta, gamma = 1, T_0 = 0.6/ohm):
    def f(R, t):
        r1 = R[0]
        r2 = R[1]
        r3 = R[2]
        # Coupled ODE, as 3 different functions 
        # assume T2-->0
        f1 = - delta*r2
        f2 = delta*r1 + ohm*r3
        f3 = - r3/gamma -ohm*r2
        return [f1,f2,f3]
    t = np.linspace(0, T_0,50) # Time evolution
    # Manually minimize mxstep, if too large then Integration unsucessful.
    sol = odeint(f, R0, t,mxstep=300000)#,full_output=1) # use full output for debugging
    # print ("sol:"+str(sol))
    # r3 = sol[0][:,2]
    r3 = sol[:,2]
    print ("r3: ".format(r3))
    return r3[-1]
