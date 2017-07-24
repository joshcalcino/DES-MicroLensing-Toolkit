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

class driver(object):

    def __init__(self):
        self.hpix = getHPIX.pix() #list of all pixels in survey 
        self.data = self.all_data()
        print "in driver- teff:", self.data[0].get_t_eff(self.data[0].uniqueIDs[0])
    
    def mishelleCode():
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
`                       ra_list = np.append(ra_list, RA)
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

    def nike(self, test):
        for pix in self.hpix:
            data = getData.getData(pix) #160,000 objects with seperate obs
            objID = data.uniqueIDs #list of all objects
            for i in range(0,10,1): #for i in data:
                #variables from data
                mjd = data.get_timesByIDs(objID[i])
                t_eff = data.get_t_eff(objID[i])
                m_0 = data.get_m_0(objID[i]) 
                is_star = data.isStar(objID[i], mjd)
                """
                Ds = data.get_Ds(objID[i])
                curve_type = data.get_curve_type(objID[i])
                """

                star = star.star()
                
                #calculated variables
                if is_star = True:
                    star_event = star.get_curves(self.mjd, t_eff, m_0)
                    self.event_list.append(star_event)
                #call the save to fits file method
                #self.plot_many(0,100)
            index =+ 1
        print "index:", index
        return 0


    def all_MJDs(self):
        for i in range(0,len(data),1):
            for index in range(0,len(data),1):
                for bandpass in ['y','r','g']:
                    mjd = data[i].get_MJD(index, bandpass)


   # @profile
    def plot_many(self, start, stop, step=1):
        fake_plots.clear()
        index = 0
        while start < stop:
            fake_plots.plot_many(self.event_list[start])
            start += step
            index += 1
        return index

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
