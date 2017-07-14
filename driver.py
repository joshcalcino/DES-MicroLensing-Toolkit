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
from astropy.io import fits
import pyfits

class driver(object):

    def __init__(self):

	counts = np.array([312, 334, 308, 317])					
      	names = np.array(['NGC1', 'NGC2', 'NGC3', 'NGC4'])
	c1 = pyfits.Column(name='target', format='10A', array=names) 			
	c2 = pyfits.Column(name='counts', format='J', unit='DN', array=counts)	
	c3 = pyfits.Column(name='notes', format='A10')					
	c4 = pyfits.Column(name='spectrum', format='1000E')
	c5 = pyfits.Column(name='flag', format='L', array=[True, False, True, True])	#these 5 lines ^^^ make columns
 	cols = pyfits.ColDefs([c1, c2, c3, c4, c5])					#column definitions
	tbhdu = pyfits.BinTableHDU.from_columns(cols)					#create a new table
	hdu = pyfits.PrimaryHDU() 							#create a primary HDU object to hold all lists / tables
	hdulist = pyfits.HDUList([hdu, tbhdu])  					#create a HDUList containing the primary HDU and newly created table extension
	self.hpix = getHPIX.pix()							#list of all pixels in survey 

        for i in self.hpix:
            self.file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/test/lc_curves" + str(i) + ".fits"
            hdulist.writeto(self.file_name, clobber=True)				#writes the primary HDU and newly created table extension into a file for each pixel

    def nike(self):
        for pix in self.hpix:
            data = getData.getData(pix) #160,000 objects with seperate obs
            final_IDs = data.isStar() #list of objects that are stars
            for i in range(0, 10, 1): #goes through the star IDs, just 10 for now
		
		#variables from data
                mjd = data.get_timesByIDs(final_ID[i]) 
                t_eff = data.get_t_eff(final_ID[i]) 
                mag_list = data.get_m_0(final_ID[i])  
                #Ds = data.get_Ds(objID[i]) 
                #curve_type = data.get_curve_type(objID[i]

                #calculated variables
		star = star.star()
 		star_events = star.get_curves(mjd, t_eff, mag_list) #returns 600,000 light curves for the object  
                self.save_data(pix, star_events)

                #call the save to fits file method
                #self.plot_many(0,100)

            index =+ 1
        print "index:", index 
        return 0
    
    def save_data(self, pix, events):
        
        #get the different columns that we want to return
        #for i in events:
            #ID = events[i].objID
            #RA = ....
            #DEC = ...
            #MJD_LIST = ...
            #MAG_LIST = event[i].light_curve
            #           filename, data, header=None, checksum=False, verify=True, **kwargs
            #fits.append(pix, ID, OBJECT_ID)
            #fits.append(pix, RA, RA)
            #fits.append(pix, DEC, DEC)
            #fits.append(pix, MAG_LIST[j], MAG)
            #fits.append(pix, MJD_LIST[j], MDJ_LIST)

        #save all data per pixel, that way, there are only 1800 total files
        return "saved!!!"
