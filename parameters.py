import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator
import fake_plots
import getData
import itertools

class loop(object):
 
    #need to loop over Parameter Limits:
    #Mlens   [10,100,1]             -Mass of the deflector
    #t_0     (MJD_min-1424, MJD_max+1424, 1)  -time of maximum light amplification
    #u_0     [0,2, 0.01]            -impact parameter
    #v_t     [20,220km/s, 100]      -transverse velocity of the lensing object
    #Ds      20 kpc                 -distance to the source, based on actually calculations, not varied
    #x       (0,1,0.05)             -percentage of distance from Dl to Ds (Dl/Ds)
    #m_0                            -average magnitude of source
    #t_E     (438, 1424,1)          -converted from (1.2-3.9) yrs to MJD    

    """ star_object(MJD_list): takes  a list of integers and returns a list of Mircolensing events. """
    def star_object(self, MJD_list, m_0=30, Ds=20, t_eff=0,curve_type =1):
        self.lightcurve = []             

#        t_eff = 0                   #need to pick a variable
#        curve_type = 1              #PAC curve
#        m_0 = 30                    #eventually will need to call data.get_m_0()
#        Ds = 20                     #eventually will need to call data.get_Ds()

        """
        "The loop below takes approximately 30 minutes per mass producing about 12 trillion  light curves. "
        urange = self.get_drange(0,2,.01)
        x_range = self.get_drange(0,1,0.05)
        for m in range(10, 100,1):
            for u_0 in urange:
                for v_t in range(20,220, 100):
                    for x in x_range:
                        for t_E in range(438,1424,1):
                            for t_0 in range(int(min(MJD_list)-1424), int(max(MJD_list)+1424), 1):        
                                lightcurve.append(MicroLensingGenerator.GenerateMLEvent(t_0, u_0, v_t, m, Ds, x, MJD_list, m_0, t_eff, curve_type))
        return lightcurve 
        """

        """The loop below takes approximately 3 minutes per mass producing 640,200 light curves. The complete loop would produce 6.4 million light curves in 30 minutes."""
        urange = self.get_drange(0,2,.5) #.1
        x_range = self.get_drange(0,1,.5) #.1
        index = 0
        v_t = 220
        for m in range(10, 101, 20): #10
            print "m, number of events: ", m,  index
            for u_0 in urange:
                for x in x_range:
                    for t_E in range(438,1424,50): #30
                        for t_0 in range(int(min(MJD_list)-1424), int(max(MJD_list)+1424), 50): #30        
                            self.lightcurve.append(MicroLensingGenerator.GenerateMLEvent(t_0, u_0, v_t, m, Ds, x, MJD_list, m_0, t_eff, curve_type))
                            index += 1
        return self.lightcurve 

    """ get_drange(start, stop, step): takes 3 float parameters and returns a list from start to stop with the step interval. """
    def get_drange(self, start, stop, step): #function that returns a list of decimals in a given range
        r = []
        while start < stop:
            r.append(start)
            start += step
            start = round(start,2)
        return r
"""
    def plot_compare(start, length):
        #compares index lightcurves
        for i in length:
            fake_plots.nike(self.lightcurve[start].times, self.lightcurves[start].t_0)
        return 0
"""
 
