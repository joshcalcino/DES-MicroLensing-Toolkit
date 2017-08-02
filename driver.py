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

def load_data(pixel="11200", test_ID = 11120000000150):
    hpix = getHPIX.pix() #list of all pixels in survey
    index = 0
    mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list  = \
        [], [], [], [], [], [], []
    for pix in hpix:
        if pix != pixel: continue
        data = getData.getData(pix) #160,000 objects with seperate obs
        objID = data.uniqueIDs #list of all objects
        for i in range(0,2,1): #for i in objID:
            #variables from data
            #objID = [test_ID]
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
    return mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, pix

def marika(low, high):
    t1 = time.time()
    numStars = []
    hpix = getHPIX.pix()
    if high > len(hpix): raise Exception
    for i in range(0,len(hpix),1):
        if (i >= low) and (i < high):
            data = getData.getData(hpix[i])
            stars = data.count_stars()
            objs = data.uniqueIDs
            numStars.append([hpix[i], stars, len(objs)])
    t2 = time.time()
    seconds = t2 - t1
    print "Marika took", seconds, " seconds to run. Or", seconds/60, "minutes." 
    return numStars

def mishelle(index = 11200):
    data = getData.getData(index)
    uniq = data.uniqueIDs
    return uniq

def nike(data, index=0):
    mjd_list = data[0]     
    teff_list = data[1]
    m0_list = data[2]
    ra_list = data[3]
    dec_list = data[4]
    objID_list = data[5]
    band_list = data[6]
    pix = data[7]

    curves = []
    size =  len(ra_list)

    for i in range(size): #for i in data:
        print "i:", i
        #if i != index: continue
        star_set = star.star()
        #calculated variables
        t_eff= teff_list[i]
        m_0 = m0_list[i]
        ra = ra_list[i]
        dec = dec_list[i]
        objID = objID_list[i]
        mjd = mjd_list[i]
        band = band_list[i]
   
        tmp = star_set.get_curves(mjd, band, objID, ra, dec, t_eff, m_0)
        curves.append(tmp)
    save_data(curves, pix)
    #pickle.dump(tmp, open("data_july_25.pickle","wb"))
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

def save_data(data, pix):
    objID_list = [] 
    ra_list = []
    dec_list = []
    mag_list = [] #final magnitude, aka lightcurve
    magerr_list = [] #error bar
    mjd_list = []
    bandpass_list = []

    for list_of_curves in data:
        for curve in list_of_curves:
            objID_list = np.append(objID_list, curve.quickid)
            ra_list = np.append(ra_list, curve.ra)
            dec_list = np.append(dec_list,curve.dec)
            mag_list = np.append(mag_list, curve.light_curve)
            magerr_list = np.append(magerr_list, curve.light_curve_error)
            mjd_list = np.append(mjd_list, curve.times)
            bandpass_list = np.append(bandpass_list, curve.bandpass)

    file_name = "/home/s1/marika/data/DES-MicroLensing-Toolkit/fitsData/test/ml_curves" + str(pix) + ".fits"
    if os.path.exists(file_name):
        os.remove(file_name)
        print "removed the file!"
    print "cool!"
        
    fits = fitsio.FITS(file_name,'rw')
    print 1
    array_list = [objID_list, ra_list,dec_list,mag_list,magerr_list,mjd_list,bandpass_list]
    print 2
    names = ['OBJID','RA','DEC','MAG','MAGERR','MJD','BAND']  
    print 3
    fits.write(array_list, names=names, overwrite = True)
    print "saved!"

        
