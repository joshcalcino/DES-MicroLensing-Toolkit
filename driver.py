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

    """
    def __init__(self, events):
        self.event_list = events
    """
    def nike(self):
        index = 0
        for hpix in range(11737,11738,1):
            data = getData.getData(hpix)
            # for i in data.uniqueIDs:
            for i in range(0,10,1):
                #self.mjd = data.get_timesByIDs(i)
                self.mjd = data.get_timesByIDs(data.uniqueIDs[i])   
                self.star = parameters.star() 
                self.event_list = self.star.get_curves(self.mjd)  
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
        return 0


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
