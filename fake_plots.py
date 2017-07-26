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

def plot_many(event, color = "rosybrown"):
    if len(event.times) == 1:
        raise Exception ("Length of MJD_list is: ", len(event.times))
    plt.ion()
 
    t = event.times
    b = event.bandpass
    t_0 = event.t_0    
    curves = event.all_curves
    error_bars_r = event.generate_noise[0] 
    error_bars_Y = event.generate_noise[1] 
    error_bars_g = event.generate_noise[2] 
    error_bars_z = event.generate_noise[3] 
    error_bars_i = event.generate_noise[4] 

    print "lens of r,y,g,z,i:", len(error_bars_r), len(error_bars_Y), len(error_bars_g), len(error_bars_z), len(error_bars_i)
    #error_bars = event.generate_noise
    u = event.u
    A = event.A
    delta_mag = event.delta_mag  # change in the magnitude of the star due to the lensing
    
    mpl.rcParams['figure.figsize'] = (8, 5)

    t_time = np.arange(t.min(),t.max(),0.01)
    interp = interp1d(t, delta_mag, bounds_error =  False, kind = 'linear')
    dm = interp(t_time)

    interp2 = interp1d(t, curves, bounds_error = False, kind = 'linear') 
    cm = interp2(t_time)    

    r_index = 0
    Y_index = 0
    g_index = 0
    z_index = 0
    i_index = 0

    r_avg = 0
    g_avg = 0
    i_avg = 0

    plt.axvline(x=t_0, label="t_0", color ='orange') 
    plt.plot(t_time, dm, c = color, zorder=1)
    #plt.plot(t_time, cm, c = "rosybrown", zorder=1)
    for i in range(0, len(b), 1):
        if b[i] == 'r':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
            plt.errorbar(t[i],delta_mag[i],error_bars_r[r_index], capsize =4, c='tomato', elinewidth=1 )
            #plt.scatter(t[i],event.light_curve_r[r_index], s= 50, c = "r",marker='s', edgecolor='black', zorder=2)
            #plt.errorbar(t[i],event.light_curve_r[r_index],error_bars_r[r_index], capsize =4, c='tomato', elinewidth=1 )
            #plt.scatter(t[i], event.m_0[i] - 20.168972, c='r')
            #print "r-i:", event.m_0[i] - 20.168972
            r_avg += event.m_0[i] - 20.168972
            r_index += 1
        if b[i] == 'Y':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "gold",marker='8', edgecolor='black', zorder=6)
            plt.errorbar(t[i], delta_mag[i],error_bars_Y[Y_index], capsize =4, c='goldenrod', elinewidth=1)
            #plt.scatter(t[i],event.light_curve_Y[Y_index], s= 50, c = "gold",marker='s', edgecolor='black', zorder=2)
            #plt.errorbar(t[i],event.light_curve_Y[Y_index],error_bars_Y[Y_index], capsize =4, c='goldenrod', elinewidth=1 )
            #print "error bars Y:", error_bars_Y[Y_index]
            Y_index += 1
        if b[i] == 'z':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "m",marker='p', edgecolor='black', zorder=4)
            plt.errorbar(t[i], delta_mag[i],error_bars_z[z_index], capsize =4, c='mediumorchid', elinewidth =1)
            #plt.scatter(t[i],event.light_curve_z[z_index], s= 50, c = "m",marker='s', edgecolor='black', zorder=2)
            #plt.errorbar(t[i],event.light_curve_z[z_index],error_bars_z[z_index], capsize =4, c='mediumorchid', elinewidth=1 )
            z_index += 1
        if b[i] == 'g':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "g",marker='D', edgecolor='black', zorder=3)
            plt.errorbar(t[i], delta_mag[i],error_bars_g[g_index], capsize =4, c='lime', elinewidth=1)
            #plt.scatter(t[i],event.light_curve_g[g_index], s= 50, c = "g",marker='s', edgecolor='black', zorder=2)
            #plt.errorbar(t[i],event.light_curve_g[g_index],error_bars_g[g_index], capsize =4, c='lime', elinewidth=1 )
            #plt.scatter(t[i], event.m_0[i] - 20.561211, c='g')
            g_avg += event.m_0[i] - 20.561211
            #print "g-r:", event.m_0[i] - 20.561211
            g_index += 1
        if b[i] == 'i':
            plt.scatter(t[i],delta_mag[i], s= 50, c = "mediumblue",marker='>', edgecolor='black', zorder=5)
            plt.errorbar(t[i], delta_mag[i],error_bars_i[i_index], capsize =4, c='royalblue', elinewidth=1)
            #plt.scatter(t[i],event.light_curve_i[i_index], s= 50, c = "mediumblue",marker='s', edgecolor='black', zorder=2)
            #plt.errorbar(t[i],event.light_curve_i[i_index],error_bars_i[i_index], capsize =4, c='royalblue', elinewidth=1 )
            #plt.scatter(t[i], event.m_0[i] - 20.003294, c='mediumblue')
            i_avg += event.m_0[i] - 20.003294
            #print "i-z:", event.m_0[i] - 20.003294
            i_index += 1

    #print "r-i", r_avg/r_index
    #print "g-r", g_avg/g_index
    #print "i-z", i_avg/i_index
    plt.xlabel("MJD")
    plt.ylabel("Magnitude")
    plt.title("Microlensing Curves")

    #plt.xlim(56900, 57000)
    #plt.ylim(0, 0.1)

    plt.grid()
    plt.show()
    return 0
