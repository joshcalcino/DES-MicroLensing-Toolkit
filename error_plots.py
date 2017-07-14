import numpy as np
import matplotlib.pyplot as plt
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


def error_plots(magnitude, magerr, t_eff, magnitude2, magerr2, t_eff2):

    mag_plot = []
    magerr_plot = []
    t_eff_plot = []

    mag_plot2 = []
    magerr_plot2 = []
    t_eff_plot2 = []


    for line in magnitude:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        mag_plot.append(current_line)

    for line in magerr:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        magerr_plot.append(current_line)
    
    total1 = 0
    for line in t_eff:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        total1 = total1 + current_line
        t_eff_plot.append(current_line)

    for line in magnitude2:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        mag_plot2.append(current_line)
        
    for line in magerr2:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        magerr_plot2.append(current_line)
       
    total2 = 0
    for line in t_eff2:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        total2 = total2 + current_line
        t_eff_plot2.append(current_line)
    
    total1 = total1/len(t_eff_plot)
    total2 = total2/len(t_eff_plot2)
    

    print("Total1: " + str(total1))
    print("Total2: " + str(total2))

    mag_plot = np.array(mag_plot)
    magerr_plot = np.array(magerr_plot)
    t_eff_plot = np.array(t_eff_plot)

    ix = np.argsort(mag_plot)
    mag_plot = mag_plot[ix]
    magerr_plot = magerr_plot[ix]

    #eventually feed in the magnitudes instead of the min and max below
    #x_axis = np.arange(mag_plot.min(), mag_plot.max(), 0.001)
    #interp = interp1d(mag_plot, magerr_plot, bounds_error = False)
    #interp_mag = interp(x_axis)

    plt.clf()

    x_axis = mag_plot
    interp = interp1d(mag_plot, magerr_plot, bounds_error = False)
    interp_mag = interp(x_axis)

    error = interp_mag + (-2.5 * np.log10(t_eff_plot / total1 )*.5)
    
    print error
    stopity

    #ratio of t_eff: denominator is average t_eff for plot, numerator is from each individual magnitude measurement
    #-2.5*log(that)*.5
    #interp(mag_list)+(-2.5*log(t_eff_list/t_eff)*.5 is the error?
    #plt.scatter(mag_plot2, magerr_plot2, color = 'red')
    plt.scatter(mag_plot, magerr_plot, color = 'blue', s = 5) #interp1d
    plt.plot(x_axis, interp_mag, c = "r")
    #plt.scatter(np.array(mag_plot)+(-.86*.5), magerr_plot, color = 'blue', alpha = .1)
    plt.title("Magnitude and Error of selection of stars")
    plt.xlabel("Magnitude")
    plt.ylabel("Error")
    #plt.axhline(y = 0.0085, color = 'r')
    plt.savefig("t_eff_overlap.png")
    plt.show()
