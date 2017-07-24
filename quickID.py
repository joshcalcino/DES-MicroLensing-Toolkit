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
import time

start_time = time.time()
        
def collectData():
    hpix = getHPIX.pix() #list of all pixels in survey
    index = 0
    for pix in hpix: #for pix in hpix
        if pix != "11737": continue
        m0_list, magerr_list, teff_list, band_list, objID_list, mjd_list = \
            np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
            #[], [], [], [], []
        data = getData.getData(pix) #160,000 objects with seperate obs
        objID = data.uniqueIDs #list of all objects
        #print len(objID)
        for i in range(1,len(objID)): #for i in range(4,30,1):
            mjd = data.get_timesByIDs(objID[i])
            m_0 = data.get_mag(objID[i])
            magerr = data.get_magerr(objID[i])
            t_eff = data.get_t_eff(objID[i])
            bandpass = data.get_bandpass(objID[i])
            is_star = data.isStar(objID[i])

            if is_star:
                for j in mjd:
                    #m0_list.append( m_0 )
                    #magerr_list.append( magerr )
                    #teff_list.append( t_eff )
                    #band_list.append( bandpass )
                    #objID_list.append( objID[i] )
                    m0 = np.asarray(m_0)
                    magerr = np.asarray(magerr)
                    teff = np.asarray(t_eff)
                    band = np.asarray(bandpass)
                    MJD = np.asarray(mjd)    
                    obj_id = np.asarray(objID[i]).astype(int)
                    obj_id =obj_id*np.ones(m0.size).astype(int)
                    m0_list = np.append(m0_list, m0)
                    magerr_list = np.append(magerr_list, magerr)
                    teff_list = np.append(teff_list, teff)
                    band_list = np.append(band_list, band)
                    objID_list = np.append(objID_list, obj_id)
                    mjd_list = np.append(mjd_list, MJD)
                """
                m0_list = np.asarray(m_0)
                magerr_list = np.asarray(magerr)
                teff_list = np.asarray(t_eff)
                band_list = np.asarray(bandpass)
                objID_list = np.asarray(objID[i])
                print m0_list
                print teff_list
                print objID_list

                np.append(m0_list, m_0)
                np.append(magerr_list, magerr)
                np.append(teff_list, t_eff)
                np.append(band_list, bandpass)
                np.append(objID_list, objID[i])
                print m0_list
                print magerr_list

                m0_list.append( m_0 )
                magerr_list.append( magerr )
                teff_list.append( t_eff )
                band_list.append( bandpass )
                objID_list.append( objID[i] )"""
                

        filters = np.zeros(m0_list.size)
        ix = band_list == "u"; filters[ix] = 0
        ix = band_list == "g"; filters[ix] = 1
        ix = band_list == "r"; filters[ix] = 2
        ix = band_list == "i"; filters[ix] = 3
        ix = band_list == "z"; filters[ix] = 4
        ix = band_list == "Y"; filters[ix] = 5
        ix = band_list == "y"; filters[ix] = 5

        #print m0_list, magerr_list, teff_list, filters, objID_list, mjd_list, pix
        save_data(m0_list, magerr_list, teff_list, filters, objID_list, mjd_list, pix)

    # this is the end of the loop of all pixels

def save_data(m0_list, magerr_list, teff_list, filters, objID_list, mjd_list, pix):
    #m0_array = np.asarray(m0_list)
    #magerr_array = np.asarray(magerr_list)
    #teff_array = np.asarray(teff_list)
    #band_array = np.asarray(band_list)
    #objID_array = np.asarray(objID_list)
    #try:
        #pix = int(pix[0])
    #except:
        #pix = int(pix)
    
    file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/teff_practice/pixel" + str(pix) + ".fits"
    if os.path.exists(file_name):
        os.remove(file_name)
        print "removed the file!"
    print "cool!"

    fits = fitsio.FITS(file_name,'rw')
    #array_list = [m0_array, magerr_array, teff_array, objID_array]
    array_list = [m0_list, magerr_list, teff_list, filters, objID_list, mjd_list]
    #names = ['m0_array', 'magerr_array', 'teff_array', 'objID_array']
    names = ['m0', 'magerr', 'teff', 'filters', 'objID', 'mjd']
    fits.write(array_list, names=names, overwrite = True)
    print "saved!"
    print "My program took",(time.time() - start_time)/60, "minutes to run" 
                

