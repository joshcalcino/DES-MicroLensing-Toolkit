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
    mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list  = \
        [], [], [], [], [], [], []
    for pix in hpix:
        if pix != "11737": continue
        data = getData.getData(pix) #160,000 objects with seperate obs
        objID = data.uniqueIDs #list of all objects
        for i in range(4,5,1): #for i in data:
            #variables from data
            mjd = data.get_timesByIDs(objID[i])
            t_eff = data.get_t_eff(objID[i])
            m_0 = data.get_mag(objID[i]) 
            is_star = data.isStar(objID[i])
            RA = data.get_RA(objID[i])
            DEC = data.get_DEC(objID[i])
            bandpass = data.get_bandpass(objID[i])
            print "mjd", len(mjd)
            print "RA", len(RA)
            print "bandpass", len(bandpass)

            if is_star:
                mjd_list.append( mjd )
                teff_list.append( t_eff )
                m0_list.append( m_0 )
                ra_list.append( RA )
                dec_list.append( DEC )
                objID_list.append( objID[i] )
                band_list.append( bandpass )
            """
            Ds = data.get_Ds(objID[i])
            curve_type = data.get_curve_type(objID[i])
            """
        #save_data(mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, pix)    
    return mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list

def nike(mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, index):
    curves = []
    size =  len(ra_list)
    for i in range(size): #for i in data:
        print "i:", i
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

def plot_many( start, stop, step=1):
    index = 0
    while start < stop:
         fake_plots.plot_many(event_list[start])
         start += step
         index += 1
    return index

def plot1(event):
    fake_plots.clear()
    fake_plots.plot_many(event)    
    return 0

def save_data(self, mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, pix):
    mjd_array = np.asarray(mjd_list)
    teff_array = np.asarray(teff_list)
    m0_array = np.asarray(m0_list) 
    ra_array = np.asarray(ra_list)
    dec_array = np.asarray(dec_list)
    objID_array = np.asarray(objID_list)
        
    self.file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/test/ml_curves" + str(pix) + ".fits"
    if os.path.exists(self.file_name):
        os.remove(self.file_name)
        print "removed the file!"
    print "cool!"
        
    fits = fitsio.FITS(self.file_name,'rw')
    array_list = [mjd_array, teff_array, m0_array, ra_array, dec_array, objID_array]
    names = ['mjd_array', 'teff_array', 'm0_array', 'ra_array', 'dec_array', 'objID_array'] 
    fits.write(array_list, names=names, overwrite = True)
    print "saved!"        
