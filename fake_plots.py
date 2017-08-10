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


def nike(event) :
    plt.clf()

def clear():
    plt.clf()

def plot_many(event, color = "rosybrown"):
    if len(event.times) == 1:
        raise Exception ("Length of MJD_list is: ", len(event.times))
    plt.ion()
 
    t = event.times
    b = event.bandpass
    t_0 = event.t_0    
    m_0 = event.m_0
    mag = event.light_curve
    error = event.light_curve_error
    delta_mag = event.delta_mag  # change in the magnitude of the star due to the lensing
    
    mpl.rcParams['figure.figsize'] = (8, 5)
    t_time = np.arange(t.min(),t.max(),0.01)

    interp = interp1d(t, delta_mag, bounds_error =  False, kind = 'linear')
    dm = interp(t_time) 

    interp2 = interp1d(t, mag, bounds_error = False, kind = 'linear') 
    cm = interp2(t_time)    

    interp3 = interp1d(t, m_0, bounds_error = False, kind = 'linear')
    mm = interp3(t_time)

    plt.axvline(x=t_0, label="t_0", color ='orange') 
    #plt.plot(t_time, dm, c = "rosybrown", zorder=1)
    plt.plot(t_time, cm, c = "g", zorder = 1)
    plt.plot(t_time, mm, c = 'r', zorder = 1)

    #Mag
    for i in range(0, len(b), 1):
        if b[i] == 'r':
            plt.scatter(t[i], mag[i], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
            plt.errorbar(t[i], mag[i], error[i], capsize = 4, c='tomato', elinewidth=1 )
        if b[i] == 'Y':
            plt.scatter(t[i], mag[i], s= 50, c = "gold",marker='8', edgecolor='black', zorder=6)
            plt.errorbar(t[i], mag[i], error[i], capsize = 4, c='goldenrod', elinewidth=1)
        if b[i] == 'z':
            plt.scatter(t[i], mag[i], s= 50, c = "m",marker='p', edgecolor='black', zorder=4)
            plt.errorbar(t[i], mag[i], error[i], capsize =4, c='mediumorchid', elinewidth =1)
        if b[i] == 'g':
            plt.scatter(t[i], mag[i], s= 50, c = "g",marker='D', edgecolor='black', zorder=3)
            plt.errorbar(t[i], mag[i], error[i], capsize =4, c='lime', elinewidth=1)
        if b[i] == 'i':
            plt.scatter(t[i], mag[i], s= 50, c = "mediumblue",marker='>', edgecolor='black', zorder=5)
            plt.errorbar(t[i], mag[i], error[i], capsize =4, c='royalblue', elinewidth=1)
    #MJD
    for i in range(0, len(b), 1):
        if b[i] == 'r':
            plt.scatter(t[i],  m_0[i], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
            plt.errorbar(t[i], m_0[i],error[i], capsize =4, c='tomato', elinewidth=1 )
        if b[i] == 'Y':
            plt.scatter(t[i],  m_0[i], s= 50, c = "gold",marker='8', edgecolor='black', zorder=6)
            plt.errorbar(t[i], m_0[i],error[i], capsize =4, c='goldenrod', elinewidth=1)
        if b[i] == 'z':
            plt.scatter(t[i],  m_0[i], s= 50, c = "m",marker='p', edgecolor='black', zorder=4)
            plt.errorbar(t[i], m_0[i],error[i], capsize =4, c='mediumorchid', elinewidth =1)
        if b[i] == 'g':
            plt.scatter(t[i],  m_0[i], s= 50, c = "g",marker='D', edgecolor='black', zorder=3)
            plt.errorbar(t[i], m_0[i],error[i], capsize =4, c='lime', elinewidth=1)
        if b[i] == 'i':
            plt.scatter(t[i],  m_0[i], s= 50, c = "mediumblue",marker='>', edgecolor='black', zorder=5)
            plt.errorbar(t[i], m_0[i],error[i], capsize =4, c='royalblue', elinewidth=1)
    """
    #DelMag
    for i in range(0, len(b), 1):
        if b[i] == 'r':
            plt.scatter(t[i],  delta_mag[i], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
            plt.errorbar(t[i], delta_mag[i],error[i], capsize =4, c='tomato', elinewidth=1 )
        if b[i] == 'Y':
            plt.scatter(t[i],  delta_mag[i], s= 50, c = "gold",marker='8', edgecolor='black', zorder=6)
            plt.errorbar(t[i], delta_mag[i],error[i], capsize =4, c='goldenrod', elinewidth=1)
        if b[i] == 'z':
            plt.scatter(t[i],  delta_mag[i], s= 50, c = "m",marker='p', edgecolor='black', zorder=4)
            plt.errorbar(t[i], delta_mag[i],error[i], capsize =4, c='mediumorchid', elinewidth =1)
        if b[i] == 'g':
            plt.scatter(t[i],  delta_mag[i], s= 50, c = "g",marker='D', edgecolor='black', zorder=3)
            plt.errorbar(t[i], delta_mag[i],error[i], capsize =4, c='lime', elinewidth=1)
        if b[i] == 'i':
            plt.scatter(t[i],  delta_mag[i], s= 50, c = "mediumblue",marker='>', edgecolor='black', zorder=5)
            plt.errorbar(t[i], delta_mag[i],error[i], capsize =4, c='royalblue', elinewidth=1)
        """
    plt.xlabel("MJD")
    plt.ylabel("Magnitude")
    plt.title("Microlensing Curves")

    #plt.xlim(56900, 57000)
    #plt.ylim(0, 0.1)

    plt.grid()
    plt.show()
    return 0
