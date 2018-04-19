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
import random
import pickle
from scipy.interpolate import interp1d
from matplotlib.ticker import FormatStrFormatter
import data_practice as dp

start_time = time.time()
"""
p = number of pixels to sample
o = number of objects to sample
nike():
    load_data() # loads one pixel of data, all stars
    get_curves()
    save_data()
"""
#runs one hpix, and loads data, then ets curves and adds to list, then saves data
def nike(p = 2, o = 5): #when ready to run, remove 'p=2, o=5'
    hpix = getHPIX.pix()
    hpix = ['9280', '9281']
    pixies = random.sample(hpix, p) #when ready to run, remove this line
    plotable = []
    print "Pixels:", pixies
    for pix in pixies: #when ready to run over ALL pixels, replace 'pixies' with 'hpix'
        f_name = pix + 'test' #when ready to run remove this line
        t0 = time.time()
        data = load_data(pix, o) #when ready to run, change 'o' to 100
        t1 = time.time()
        curves = get_curves(data)
        plotable.append( curves )
        t2 = time.time()
        files = save_data(curves, f_name) #when ready to run replace 'f_name' with 'pix'
        t3 = time.time()
        data_time = round(t1-t0, 2)
        curves_time = round(t2-t1, 2)
        save_time = round(t3-t2, 2)
        print "Pixel", pix, "Complete:"
        print "Load Data:", data_time, "seconds"
        print "Curves:", curves_time, "seconds"
        print "Save File:", save_time, "seconds"
    return plotable

""" Currently load_data only works for 1 pixel at a time. """
def load_data(pix, sample_size = 2):
    mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list, pix_list  = \
        [], [], [], [], [], [], [], []
    data = getData.getData(pix) #13,000-250,000 objects with seperate obs
    objID = data.sIDs #list of all stars
    objs = random.sample(objID, sample_size)
    print "Object IDs:", objs
    
    for ID in objs: 
        #variables from data
        data.loadStar( ID )
        mjd = data.getMJD( ID ) 
        t_eff = data.get_t_eff( ID )
        m_0 = data.get_mag( ID )
        RA = data.get_RA( ID )
        DEC = data.get_DEC( ID )
        band = data.get_band( ID )

        # appends data information to list
        mjd_list.append( mjd ) 
        teff_list.append( t_eff )
        m0_list.append( m_0 )
        ra_list.append( RA )
        dec_list.append( DEC )
        objID_list.append( ID )
        band_list.append( band )

    return mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list

""" Takes data from all objects within one pixel."""
def get_curves(data): 
    mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, band_list = data
    objects =  len(mjd_list) 

    curves = []
    print "Number of Objects:", objects

    for i in range(0, objects, 1):
        print "Object", i, ":", objID_list[i] 
        t_eff= np.asarray( teff_list[i] )
        m_0 = np.asarray( m0_list[i] )
        ra = np.asarray( ra_list[i] )
        dec = np.asarray( dec_list[i] )
        objID = np.asarray( objID_list[i] )
        mjd = np.asarray( mjd_list[i] )
        band = np.asarray( band_list[i] )
 
        lcs = create_lightcurves(mjd, band, objID, ra, dec, t_eff, m_0) #45,000 MLEs
        curves.append( lcs )
    return curves

""" Takes DES data and creates Microlensing events. Outputs a list of lightcurves per star."""
def create_lightcurves(MJD_list, band, objID, ra, dec, t_eff = 1, m_0=30, Ds=5,curve_type =1):
    lightcurve = []
    v_t = 220 #velocity = 220 m/s
    lcID = 0 #lightcurve ID
    trange = [56535, 56747, 56992, 57234, 57430] # Dates: 8.31.13, 3.31.14, 12.1.14, 7.31.15, 2.12.16   
    for t_0 in trange:
        #if t_0 != 56535: continue
        print "t_0 index: ", t_0, lcID
        for u_0 in np.arange(0,2,.2):
            for x in np.arange(0.1,1,0.1):
                for M in range(10,101,10):
                    lightcurve.append(MicroLensingGenerator.GenerateMLEvent(
                        t_0, u_0, v_t, M, Ds, x, MJD_list, m_0, band, objID, ra, dec, lcID, t_eff, curve_type))
                    lcID += 1
    print "total index:", lcID
    return lightcurve

""" Saves data from all curves, of all objects in one pixel."""
def save_data(pixel_data, pix):
    objID_list, ra_list, dec_list, mag_list, magerr_list, mjd_list, bandpass_list, mlens_list, \
        t0_list, u0_list, x_list, lcid_list = [], [], [], [], [], [], [], [], [], [], [], []

    #This for loop appends data from each curve to a single row.
    for obj in pixel_data: #pixel_data needs to be one pixel worth of data
        for curve in obj: #obj needs to be one object per pixel
            lcid_list = np.append( lcid_list, curve.lcid )
            objID_list = np.append( objID_list, curve.quickid )
            ra_list = np.append( ra_list, curve.ra )
            dec_list = np.append( dec_list, curve.dec )
            mag_list = np.append( mag_list, curve.light_curve )
            magerr_list = np.append( magerr_list, curve.light_curve_error )
            mjd_list = np.append( mjd_list, curve.times )
            bandpass_list = np.append( bandpass_list, curve.bandpass )
            mlens_list = np.append( mlens_list, curve.mlens )
            t0_list = np.append( t0_list, curve.t0 )
            u0_list = np.append( u0_list, curve.u0 )
            x_list = np.append( x_list, curve.xx )
            
    file_name = "/home/s1/marika/data/DES-MicroLensing-Toolkit/fitsData/test/ml_curves" + str(pix) + ".fits"
    if os.path.exists(file_name):
        os.remove(file_name)
        print "removed the file!"
    print "cool!"
        
    fits = fitsio.FITS(file_name,'rw')
    print 1
    array_list = [lcid_list, objID_list, ra_list, dec_list, mag_list, magerr_list, mjd_list, \
         bandpass_list, mlens_list, t0_list, u0_list, x_list]
    print 2
    column_headers = ['LCID', 'OBJID','RA','DEC','MAG','MAGERR','MJD','BAND', 'M_LENS', 'T_0', 'U_0', 'X']  
    print 3
    fits.write(array_list, names=column_headers)
    print "saved!"

def plots(curves, start =0, stop = 2, step=1):
    print "event len:", len(curves)
    fake_plots.clear()
    for i in range(start, stop, step):
        print "i", i
        fake_plots.plot_many(curves[i])    

def plots2(lc1, lc2, lc3, start, stop, step):
    fake_plots.clear()
    l1 = len(lc1)
    l2 = len(lc2)
    l3 = len(lc3)
    for i in range(start, stop, step):
        if i < l1:
            fake_plots.plot_many(lc1[i], "blue")
        if i < l2:
            fake_plots.plot_many(lc2[i], "chartreuse2")
        if i < l3:    
            fake_plots.plot_many(lc3[i], "red")




