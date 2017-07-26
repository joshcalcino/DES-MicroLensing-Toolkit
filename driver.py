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
            mjd = data.get_timesByIDs(objID[i])
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

                    #star_event = star.get_curves(mjd, teff, m_0)
                    #for k in len(star_event):
                        #mag_r = star_event[k].light_curve_ + "bandpass"
                        #mag_Y = star_event[k].light_curve_Y

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
                    ra_list = np.append(ra_list, RA)
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

def load_data(pixel="11200", test_ID = 11173700000150):
    hpix = getHPIX.pix() #list of all pixels in survey
    index = 0
    mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list  = \
        [], [], [], [], [], [], []
    for pix in hpix:
        if pix != pixel: continue
        data = getData.getData(pix) #160,000 objects with seperate obs
        objID = data.uniqueIDs #list of all objects
        for i in range(0,1,1): #for i in data:
            #variables from data
            objID = [test_ID]
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
    tmp = []; import pickle
    curves = []
    size =  len(ra_list)
    # i is a given star_id or quick_object_id

    # this is the approach when you hae giant single vectors of quanitties
    #idlist = np.unique(objID_list)
    #for qid in idlist:
        #ix = qid == objID_list
        #teff = teff_list[ix]   # this isnt' a lsit of vectors but a giant flat vector
        # bandpasses = bandpass_list[ix]
        #for band in np.unique(bandpasses):  # or for filter in range(0,5):  or band in ["g","r","i","z","Y"] :
            #iy = bandpass_list[ix] == band
            #rband = mag_list[ix][iy]

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
   
        curves = star_set.get_curves(mjd, band, objID, t_eff, m_0)
        tmp2 = dict()
        for l in range(0,len(curves),1):
            tmp2["mag"]  = curves[l].light_curve
            tmp2["magerr"]  = curves[l].light_curve_error
            tmp2["mjd"]  = curves[l].times
            tmp2["delmag"]  = curves[l].delta_mag
            tmp2["objID"] = curves[l].quickid
        tmp.append(tmp2)
        #self.save_data(curves)
    pickle.dump(tmp, open("data_july_25.pickle","wb"))
    return curves

def nikex(mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, index):
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
   
        curves = star_set.get_curvesx(mjd,band , t_eff, m_0)
    #plots(curves)
    return curves

def nikeu(mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, index):
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
   
        curves = star_set.get_curvesu(mjd,band , t_eff, m_0)
    #plots(curves)
    return curves

def nikeM(mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, index):
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
   
        curves = star_set.get_curvesM(mjd,band , t_eff, m_0)
    #plots(curves)
    return curves

def plots(event, start =0, stop = 2, step=1):
    print "event len:", len(event)
    fake_plots.clear()
    for i in range(start, stop, step):
        print "i", i
        fake_plots.plot_many(event[i])    
    return 0

def plots2(event1, event2, event3, start, stop, step):
    fake_plots.clear()
    e1 = len(event1)
    e2 = len(event2)
    e3 = len(event3)
    for i in range(start, stop, step):
        if i < e1:
            fake_plots.plot_many(event1[i], "blue")
        if i < e2:
            fake_plots.plot_many(event2[i], "orange")
        if i < e3:    
            fake_plots.plot_many(event3[i], "red")
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
