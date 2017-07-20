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
    if len(event.times) == 1:
        raise Exception ("Length of MJD_list is: ", len(event.times))
    
    t = event.times
    u = event.u
    A = event.A
    b = event.bandpass
    delta_mag = event.delta_mag  
    
    # show the light curve
    mpl.rcParams['figure.figsize'] = (8, 5)
    
    t_time = np.arange(t.min(),t.max(),0.01)
    interp = interp1d(t, delta_mag, bounds_error =  False, kind = 'linear')
    dm = interp(t_time)
    
    
   #plot a vector of MJDs instead of a list of MJDs 

    for i in range(0, len(b), 1):
        if b[i] == 'r':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
        if b[i] == 'Y':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "gold",marker='8', edgecolor='black', zorder=6)
        if b[i] == 'g':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "g",marker='D', edgecolor='black', zorder=3)
        if b[i] == 'z':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "m",marker='p', edgecolor='black', zorder=4)
        if b[i] == 'i':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "mediumblue",marker='>', edgecolor='black', zorder=5)
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

def quickTest(m):
    x = MicroLensingGenerator.GenerateMLEvent(56870, 0.1, 220, 30, 4.5, 0.5, [56877,56234,56589,56900,56824], m,'r', 1, 1)
    return x

def clear():
    plt.clf()

def plot_many(event):
    if len(event.times) == 1:
        raise Exception ("Length of MJD_list is: ", len(event.times))
    plt.ion()
 
    t = event.times
    b = event.bandpass
    t_0 = event.t_0     

    u = event.u
    A = event.A
    delta_mag = event.delta_mag  # change in the magnitude of the star due to the lensing
    
    mpl.rcParams['figure.figsize'] = (8, 5)

    t_time = np.arange(t.min(),t.max(),0.01)
    interp = interp1d(t, delta_mag, bounds_error =  False, kind = 'linear')
    dm = interp(t_time)
    
    plt.axvline(x=t_0, label="t_0", color ='orange') 
    plt.plot(t_time, dm, c = "rosybrown", zorder=1)
    for i in range(0, len(b), 1):
        if b[i] == 'r':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
        if b[i] == 'Y':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "gold",marker='8', edgecolor='black', zorder=6)
        if b[i] == 'g':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "g",marker='D', edgecolor='black', zorder=3)
        if b[i] == 'z':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "m",marker='p', edgecolor='black', zorder=4)
        if b[i] == 'i':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "mediumblue",marker='>', edgecolor='black', zorder=5)
    plt.xlabel("MJD")
    plt.ylabel("Magnitude Difference")
    plt.title("Microlensing Curves")

    #plt.xlim(56900, 57000)
    #plt.ylim(0, 0.1)

    plt.grid()
    plt.show()
    return 0
