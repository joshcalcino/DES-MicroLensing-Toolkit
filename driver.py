import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator
import fake_plots
import getData
import star
import getHPIX
import fitsio
import os

def load_data():
    hpix = getHPIX.pix() #list of all pixels in survey
    index = 0
    for pix in hpix:
        if pix != "11737": continue
        mjd_list, mag_list, magerr_list, error_list, teff_list, ra_list, dec_list, band_list, objID_list = \
            np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        data = getData.getData(pix) #160,000 objects with seperate obs
        objID = data.uniqueIDs #list of all objects
        for i in range(4,30,1): #for i in range(1, len(objID)).
            mjd = data.get_timesByIDs(objID[i])i
            mag = #get mag for each bandpass?
            magerr = data.get_magerr(objID[i])
            error = #get error bars
            t_eff = data.get_t_eff(objID[i])
            ra = data.get_RA(objID[i])
            dec = data.get_DEC(objID[i])
            bandpass = data.get_bandpass(objID[i])
            is_star = data.isStar(objID[i])
            """Ds = data.get_Ds(objID[i])
            curve_type = data.get_curve_type(objID[i])"""

             if is_star:
                for j in mjd:
                    MJD = np.asarray(mjd)
                    MAG = np.asarray(mag)
                    MAGERR = np.asarray(magerr)
                    ERROR = np.asarray(error)
                    TEFF = np.asarray(t_eff)
                    RA = np.asarray(ra)
                    DEC = np.asarray(dec)
                    BAND = np.asarray(bandpass)
                    obj_id = np.asarray(objID[i]).astype(int)
                    obj_id =obj_id*np.ones(m0.size).astype(int)

                    mjd_list = np.append(mjd_list, MJD)
                    mag_list = np.append(mag_list, MAG)
                    magerr_list = np.append(magerr_list, MAGERR)
                    error_list = np.append(error_list, ERROR)
                    teff_list = np.append(teff_list, TEFF)
`                   ra_list = np.append(ra_list, RA)
                    dec_list = np.append(dec_list, DEC)
                    band_list = np.append(band_list, BAND)
                    objID_list = np.append(objID_list, obj_id)
        
        filters = np.zeros(m0_list.size)
        ix = band_list == "u"; filters[ix] = 0
        ix = band_list == "g"; filters[ix] = 1
        ix = band_list == "r"; filters[ix] = 2
        ix = band_list == "i"; filters[ix] = 3
        ix = band_list == "z"; filters[ix] = 4
        ix = band_list == "Y"; filters[ix] = 5
        ix = band_list == "y"; filters[ix] = 5

        save_data(mjd_list, mag_list, magerr_list, error_list, teff_list, ra_list, dec_list, filters, objID_list, pix)

def nike(mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, index):
    curves = []
    size =  len(ra_list)
    for i in range(size): #for i in data:
        print "i:", i
        #ppend(mjd_list, MJD)
        #pppend(mjd_list, MJD)pend(mjd_list, MJD)

        if i != index: continue
        star_set = star.star()
        #calculated variables
        t_eff= teff_list[i]
        m_0 = m0_list[i]
        ra = ra_list[i]
        dec = dec_list[i]
        objID = objID_list[i]
        mjd = mjd_list[i]
        band = band_list[i]
        
        curves = star_set.get_curves(mjd,band , t_eff, m_0)
    #plots(curves)
    return curves

def plots(event, start =0, stop = 2, step=1):
    print "event len:", len(event)
    fake_plots.clear()
    for i in range(start, stop, step):
        print "i", i
        fake_plots.plot_many(event[i])    
    return 0

def save_data(mjd_list, mag_list, magerr_list, error_list, teff_list, ra_list, dec_list, filters, objID_list, pix):
    file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/test/ml_curves" + str(pix) + ".fits"
    if os.path.exists(file_name):
        os.remove(file_name)
        print "removed the file!"
    print "cool!"
        
    fits = fitsio.FITS(file_name,'rw')
    array_list = [mjd_list, mag_list, magerr_list, error_list, teff_list, ra_list, dec_list, filters, objID_list]
    names = ['mjd', 'mag', 'magerr', 'error', 'teff', 'ra', 'dec', 'filters', 'objID']
    fits.write(array_list, names=names, overwrite = True)
    print "saved!"        
