import MicroLensingGenerator
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from scipy.interpolate import interp1d

def nike() :
    x = MicroLensingGenerator.GenerateMicrolensingEvent(1, 0.1, 220, 30, 4.5, 0.5, np.arange(1, 100, 2), 12, 1, 1)
    t = x.times
    u = x.get_u(t)
    A = x.A
    delta_mag = x.get_delta_mag(t)  # change in the magnitude of the star due to the lensing

    # show the light curve
    interp = interp1d(u, A, bounds_error =  False, kind = 'linear')
    utime = np.arange(u.min(),u.max(),0.01)
    dm = 2.5*np.log10(interp(utime))
    
    
    plt.clf()
    plt.plot(utime, dm, c="r")
    plt.scatter(u, delta_mag, s=20,c="b")
    plt.ylim(0, 0.5)
    plt.xlabel("Einstein Radii")
    plt.ylabel("Magnitude Difference")
    plt.grid()
    plt.title("Fake plots")
    return u, delta_mag




