from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from scipy.interpolate import interp1d
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

def get_errors(mag_plot, t_eff, magerr_plot, mjd_list):

    print t_eff

    total = 0
    for line in range(0, len(t_eff)):
        total = total + t_eff[line]

    avg_t_eff = total/len(t_eff)
    print(avg_t_eff)
    avg_t_eff  = np.mean(t_eff)
    print(avg_t_eff)

    mag_plot = np.array(mag_plot)
    magerr_plot = np.array(magerr_plot)
    t_eff_plot = np.array(t_eff)

    ix = np.argsort(mag_plot)
    mag_plot = mag_plot[ix]
    magerr_plot = magerr_plot[ix]

    plt.clf()

    x_axis = mag_plot
    interp = interp1d(mag_plot, magerr_plot, bounds_error = False)
    interp_mag = interp(x_axis)

    error = interp_mag + (-2.5 * np.log10(t_eff / total )*.5)

    plt.errorbar(mjd_list, mag_plot, error, fmt = 'o', color = 'green', ecolor = 'r', capsize = 20, elinewidth = 1.5)
    plt.scatter(mjd_list, mag_plot, c = 'r')
    plt.savefig("errorbarsyay.png")
    plt.show()
    #return error
