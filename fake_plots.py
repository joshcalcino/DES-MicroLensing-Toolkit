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


"""  Argument inputs:
        t_0, p, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type"""
def nike(event) :
    if len(event.times) == 1:
        raise Exception ("Length of MJD_list is: ", len(event.times))
    
    #event = MicroLensingGenerator.GenerateMLEvent(tMax, 0.1, 220, 30, 4.5, 0.5, MJD_list, 12)
    t = event.times
    u = event.u
    A = event.A
    delta_mag = event.delta_mag  # change in the magnitude of the star due to the lensing
    
    # show the light curve
    mpl.rcParams['figure.figsize'] = (8, 5)
    sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')


    #interp = interp1d(u, A, bounds_error =  False, kind = 'linear')
    #utime = np.arange(u.min(),u.max(),0.001)
    #dm = 2.5*np.log10(interp(utime))
    
    t_time = np.arange(t.min(),t.max(),0.01)
    interp = interp1d(t, delta_mag, bounds_error =  False, kind = 'linear')
    dm = interp(t_time)
    
    
   #plot a vector of MJDs instead of a list of MJDs 

    plt.clf()
    plt.plot(t_time, dm, c = "r")
    plt.scatter(t, delta_mag, s = 20, c = "b")
    #plt.plot(t_time, dm, c="r")
    #plt.scatter(u, delta_mag, s=20,c="b")
   # plt.ylim(0, 0.1)
    plt.xlabel("Time in MJD")
    plt.ylabel("Magnitude Difference")
    plt.grid()
   # plt.xlim(0, 5)
    plt.title("Fake plots")
   # return event
#    return u, delta_mag

def get_dm(MJD_list, t_max):
    x = MicroLensingGenerator.GenerateMLEvent(t_max, 0.1, 220, 30, 4.5, 0.5, MJD_list, 12, 1, 1)
    A = x.A
    t = x.times
    u = x.get_u()
    return A

def quickTest():
    x = MicroLensingGenerator.GenerateMLEvent(56870, 0.1, 220, 30, 4.5, 0.5, [56877,56234,56589,56900,56824], 12,'z', 1, 1)
    return x

def clear():
    plt.clf()

def plot_many(event):
    if len(event.times) == 1:
        raise Exception ("Length of MJD_list is: ", len(event.times))
    
    t = event.times
    u = event.u
    A = event.A
    delta_mag = event.delta_mag  # change in the magnitude of the star due to the lensing
    
    #mpl.rcParams['figure.figsize'] = (8, 5)
    #sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')

    t_time = np.arange(t.min(),t.max(),0.01)
    interp = interp1d(t, delta_mag, bounds_error =  False, kind = 'linear')
    dm = interp(t_time)
    
    plt.plot(t_time, dm, c = "r")
    plt.scatter(t, delta_mag, s = 20, c = "b")
   # plt.ylim(0, 0.1)
    plt.xlabel("Time in MJD")
    plt.ylabel("Magnitude Difference")
    plt.grid()
   # plt.xlim(0, 5)
    plt.title("Fake plots")
#    return u, delta_mag
    return 0
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
