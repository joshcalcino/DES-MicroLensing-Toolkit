import MicroLensingGenerator
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from scipy.interpolate import interp1d
import sys
import pandas as pd
import matplotlib as mpl

""" 

    """


"""  Argument inputs:
        t_0, p, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type"""
def nike(MJD_list, tMax) :
    if len(MJD_list) == 1:
        raise Exception ("length of MJD_list is ", len(MJD_list))
    
    x = MicroLensingGenerator.GenerateMicrolensingEvent(tMax, 0.1, 220, 30, 4.5, 0.5, MJD_list, 12, 1, 1)
    # x = MicroLensingGenerator.GenerateMicrolensingEvent(1, 0.1, 220, 30, 4.5, 0.5, np.arange(1, 1000, 2), 12, 1, 1)
    t = x.times
    u = x.get_u(t)
    A = x.A
    delta_mag = x.get_delta_mag(t)  # change in the magnitude of the star due to the lensing
    
    # show the light curve
    mpl.rcParams['figure.figsize'] = (8, 5)
    sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')


    #interp = interp1d(u, A, bounds_error =  False, kind = 'linear')
    #utime = np.arange(u.min(),u.max(),0.001)
    #dm = 2.5*np.log10(interp(utime))
    print "u: ", u
    print "delta_mag: ", delta_mag
    
    interp = interp1d(u, delta_mag, bounds_error =  False, kind = 'linear')
    t_time = np.arange(t.min(),t.max(),0.01)
    dm = interp(t_time)
    
    print("dm: ", dm)
    
   #plot a vector of MJDs instead of a list of MJDs 

    plt.clf()
    plt.plot(t_time, dm, c = "r")
    plt.scatter(t, u, s = 20, c = "b")
    #plt.plot(t_time, dm, c="r")
    #plt.scatter(u, delta_mag, s=20,c="b")
    #plt.ylim(0, 0.5)
    plt.xlabel("MJD")
    plt.ylabel("Magnitude Difference")
    plt.grid()
    #plt.xlim(0, 5)
    plt.title("Fake plots")
#    return u, delta_mag

"""
plt.plot(u, 2.5*np.log10(interp(u)))
#plt.plot(u, interp2(delta_mag))
plt.ylim(0, 0.5)
plt.xlabel("Einstein Radii")
plt.ylabel("Magnitude Difference")
plt.grid()
plt.title("Fake graph with chosen inputs")
plt.show()
"""
