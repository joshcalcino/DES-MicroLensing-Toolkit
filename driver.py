import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator
import fake_plots
import getData
import parameters


class driver(object):

    def __init__(self):
        self.hpix = self.all_hpix() 
        self.data = self.all_data()
        teff = self.data[0].get_t_eff(self.data[0].uniqueIDs[0])
        print teff

    def nike(self, test):
        index = 0
            # for i in data.uniqueIDs:
        for j in self.data:
            #self.mjd = data.get_timesByIDs(i)
            data = self.data[j]
            for i in range(0,10,1):
                self.mjd = data.get_timesByIDs(data.uniqueIDs[i])   
                star = parameters.star() 
                t_eff = data.get_t_eff(data.uniqueIDs[i])
                star_event = star.get_curves(self.mjd, t_eff)  
                self.event_list.append(star_event)
            #self.plot_many(0,100)
            index =+ 1
        print "index:", index 
        return 0

    def all_MJDs(self):
        for i in range(0,len(data),1):
            for index in range(0,len(data),1): 
                for bandpass in ['y','r','g']:
                    mjd = data[i].get_MJD(index, bandpass)

            
    def all_hpix(self):
        hpix_list = []
       # for hpix in range(11737, 11738, 1): #05441-11878
       #     hpix_list.append(hpix)
        hpix_list.append(11737)
        return hpix_list

    def all_data(self):
        files = []
        for i in self.hpix:
            data = getData.getData(i)
            files.append(data)
        return files 

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
