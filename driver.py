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

    def save_data():
        stuff = 0
        self.fits_file.append(stuff)
        return "saved!!!"
