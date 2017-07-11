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
class star(object):

 
    #need to loop over Parameter Limits:
    #Mlens   [10,100,1]             -Mass of the deflector
    #t_0     (MJD_min-1424, MJD_max+1424, 1)  -time of maximum light amplification == should range based on life of survey
    #u_0     [0,2, 0.01]            -impact parameter
    #v_t     [20,220km/s, 100]      -transverse velocity of the lensing object
    #Ds      20 kpc                 -distance to the source, based on actually calculations, not varied
    #x       (0,1,0.05)             -percentage of distance from Dl to Ds (Dl/Ds)
    #m_0                            -average magnitude of source

    """ star_object(MJD_list): takes  a list of integers and returns a list of Mircolensing events. """
    def get_curves(self, MJD_list, m_0=30, Ds=20, t_eff=1,curve_type =1):
        self.lightcurve = []             

        m_0 = 15                    #eventually will need to call data.get_m_0()
#        Ds = 20                     #eventually will need to call data.get_Ds()

        """The loop below takes approximately 3 minutes per mass producing 640,200 light curves. The complete loop would produce 6.4 million light curves in 30 minutes."""
        urange = self.get_drange(0,2,.1) #.1
        x_range = self.get_drange(0.1,1,.1) #.1
        mrange = range(10,101,10)
        M_step = len(mrange)
        x_step = len(x_range)
        u_step = len(urange)
        v_t = 220
        index = 0
        print "u step:", u_step
        print "x step:", x_step
        print "M step:", M_step
       # for t_0 in range(int(min(MJD_list)-365), int(max(MJD_list)+365), 20): #30        
        for t_0 in range(56900, 57000, 5): #30        
            print "t_0 index: ", t_0, index
            for u_0 in urange:
                for x in x_range:
                    for M in range(10, 101, 10): #10
                        self.lightcurve.append(MicroLensingGenerator.GenerateMLEvent(t_0, u_0, v_t, M, Ds, x, MJD_list, m_0, t_eff, curve_type))
                        index += 1
        print "total index:", index
        return self.lightcurve 

    def get_drange(self, start, stop, step): #function that returns a list of decimals in a given range
    """ get_drange(start, stop, step): takes 3 float parameters and returns a list from start to stop with the step interval. """
        r = []
        while start < stop:
            r.append(start)
            start += step
            start = round(start,2)
        return r
 
