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

    def get_drange(self, start, stop, step): #function that returns a list of decimals in a given range
        print "start: ", start, "stop: ", stop, "step: ", step
        r = []
        while start < stop:
            r.append(start)
            start += step
            start = round(start,2)
        return r
 
    def star_object(self, MJD_list):
        #assuming getData has already been called
        #get the MJD_list
        #loop over all of the variables above
        #in each loop, call an instance of object "MicroLensingGenerator"
        t_eff = 0 #need to pick a variable
        curve_type = 1 #PAC curve
        lightcurve = []
        m_0 = 30
        Ds = 20
       # MJD_list = [56915, 56945]
        urange = self.get_drange(0,2,.01)
        x_range = self.get_drange(0,1,0.05)
        """     
        for i in range(0,5,1):
            lightcurve.append(MicroLensingGenerator.GenerateMicrolensingEvent(i, i, i, i, i, i, MJD_list, i, i, i))
        print lightcurve """
        for m in range(10, 100,1):
            for u_0 in urange:
                for v_t in range(20,220, 100):
                    for x in x_range:
                        for t_E in range(438,1424,1):
                            for t_0 in range(int(min(MJD_list)-1424), int(max(MJD_list)+1424), 1):        
                                lightcurve.append(MicroLensingGenerator.GenerateMicrolensingEvent(t_0, u_0, v_t, m, Ds, x, MJD_list, m_0, t_eff, curve_type))
       #                         print "curve t_0: ", t_0
        return lightcurve 
