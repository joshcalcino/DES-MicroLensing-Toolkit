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
        
def collectData():
    hpix = getHPIX.pix() #list of all pixels in survey
    index = 0
    m0_list, magerr_list, teff_list, band_list, objID_list = \
        [], [], [], [], []
        #np.array([]), np.array([]), np.array([]), np.array([]), np.array([])  
    for pix in range(11737, 11741): #for pix in hpix
        #if pix != "11737": continue
        data = getData.getData(pix) #160,000 objects with seperate obs
        objID = data.uniqueIDs #list of all objects
        for i in range(4,30,1): #for i in data:
            mjd = data.get_timesByIDs(objID[i])
            m_0 = data.get_mag(objID[i])
            magerr = data.get_magerr(objID[i])
            t_eff = data.get_t_eff(objID[i])
            bandpass = data.get_bandpass(objID[i])
            is_star = data.isStar(objID[i])

            if is_star:
                for i in mjd:
                    m0_list.append( m_0 )
                    magerr_list.append( magerr )
                    teff_list.append( t_eff )
                    band_list.append( bandpass )
                    objID_list.append( objID[i] )

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
                objID_list.append( objID[i] )
                """

        save_data(m0_list, magerr_list, teff_list, objID_list, pix)

def save_data(m0_list, magerr_list, teff_list, objID_list, pix):
    m0_array = np.asarray(m0_list)
    magerr_array = np.asarray(magerr_list)
    teff_array = np.asarray(teff_list)
    #band_array = np.asarray(band_list)
    objID_array = np.asarray(objID_list)

    file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/teff_practice/pixel" + str(pix) + ".fits"
    if os.path.exists(file_name):
        os.remove(file_name)
        print "removed the file!"
    print "cool!"

    fits = fitsio.FITS(file_name,'rw')
    array_list = [m0_array, magerr_array, teff_array, objID_array]
    #array_list = [m0_list, magerr_list, teff_list, band_list, objID_list]
    names = ['m0_array', 'magerr_array', 'teff_array', 'objID_array']
    #names = ['m0_list', 'magerr_list', 'teff_list', 'band_list', 'objID_list']
    fits.write(array_list, names=names, overwrite = True)
    print "saved!" 
                

