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

class driver(object):

    def __init__(self):
        self.hpix = getHPIX.pix() #list of all pixels in survey 
       # print "in driver- teff:", data.get_t_eff(data.uniqueIDs[0])

    def nike(self):
        for pix in self.hpix:
            data = getData.getData(pix) #160,000 objects with seperate obs
            final_IDs = data.star_list() #list of objects that are stars
            for i in range(0,10,1): #final_IDs:

                #variables from data
                mjd = data.get_timesByIDs(final_ID[i])
                t_eff = data.get_t_eff(final_ID[i])
                mag_list = data.get_m_0(final_ID[i]) 
                """
                Ds = data.get_Ds(objID[i])
                curve_type = data.get_curve_type(objID[i])
                """

                star = star.star()

                #calculated variables 
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
