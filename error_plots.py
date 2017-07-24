import numpy as np
import matplotlib.pyplot as plt
import MicroLensingGenerator
import matplotlib.pyplot as plt
import matplotlib
import pickle
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from scipy.interpolate import interp1d
import sys
import pandas as pd
import matplotlib as mpl


def error_plots(teff_low, teff_hi,input_file = "Y_error_data.txt", bandpass = 'Y'):

    mag,magerr,teff,quick_id = np.genfromtxt(input_file,unpack=True)


    mag_plot = []
    magerr_plot = []
    t_eff_plot = []

    mag_plot2 = []
    magerr_plot2 = []
    t_eff_plot2 = []

    """
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

    print(len(t_eff_plot))

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
    """ 
    #total1 = total1/len(t_eff_plot)
    #total2 = total2/len(t_eff_plot2)
    

    #print("Total1: " + str(total1))
    #print("Total2: " + str(total2))

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

    # training set the high teff set
    teff_low_ref = 0.7
    teff_hi_ref = 0.9
    teff_ref = (teff_low_ref+teff_hi_ref)/2.
    ix, = np.where((teff >= teff_low_ref ) & ( teff < teff_hi_ref))
    mag_training = mag[ix]
    magerr_training=magerr[ix]
    teff_training = teff[ix]

    interp = interp1d(mag_training, magerr_training, bounds_error = False)

    ix, = np.where((teff >= teff_low ) & ( teff < teff_hi))
    #ix, = np.where((mag > 21) & (mag < 22))
    mag_testing = mag[ix]
    magerr_testing=magerr[ix]
    teff_testing = teff[ix]

    predicted_mag_error = interp(mag_testing)
    fd = open("magerr_model_{}.pickle".format(bandpass), 'wb')
    pickle.dump([mag_training, magerr_training], fd)
    fd.close()

    error = predicted_mag_error + (-2.5 * np.log10(teff_testing / teff_ref )*.5)
    print (teff_testing/teff_ref).min()
    print (teff_testing/teff_ref).max()
    print (teff_testing/teff_ref).mean()


    #predicted_mag_error = np.log10(teff_testing/teff_ref)
    #magerr_testing = (-2.5 * np.log10(teff_testing / teff_ref )*.5)

    #predicted_mag_error = teff_testing
    #magerr_testing =  magerr_testingi
    error_file = pickle.load(open("magerr_model_{}.pickle".format(bandpass), 'rb'))
    interp_woo = interp1d(error_file[0], error_file[1], bounds_error = False)
    timey = np.arange(error_file[0].min(),error_file[0].max(),0.05)
    plot_me = interp_woo(np.sort(mag_testing))
    wimey = interp_woo(timey)
    #print error

    print(np.sort(mag_testing))

    denominator = (np.square(magerr_plot)*np.log(10)**2)*90
    numerator = 6.25*t_eff_plot*5
    calculated = abs(numerator/denominator)
    #calculated = (6.25)/(np.square(magerr_plot)*np.log(10)**2)*t_eff_plot
    #print calculated
    plt.scatter(mag_testing, magerr_testing, color = 'blue', s = 15)
    plt.plot(np.sort(mag_testing), plot_me, c = 'crimson')
    #plt.plot(timey, wimey, c = 'r')
    #plt.scatter(mag_training, magerr_training, color = 'blue')
    plt.xlabel("Stellar Magnitude")
    plt.ylabel("Measured Magnitude Error")
    plt.suptitle("Interpolated and Measured Error for ~500 stars")
    #plt.plot([0, 1], [0, 1], c = 'r')
    #plt.xlim(0, .3)
    #plt.ylim(0, .3)
    
    #plt.xlim(0, predicted_mag_error.max())
    #plt.ylim(0, magerr_testing.max())
    #plt.scatter(mag_plot, calculated, color = 'red')
    #plt.scatter(mag_plot, error, color = 'green')
    plt.show()
    plt.savefig("error_for_paper_06-1.5.png")
    
    """

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

    """
