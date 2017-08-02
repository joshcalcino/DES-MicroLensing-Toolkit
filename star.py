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

    def __init__(self):
        self.lightcurve = []
        #self.final_mag_list = [] 
        """
        Need to loop over Parameter Limits:

        Mlens   [10,100,1]             -Mass of the deflector
        t_0     (55000, 59000, 50)      -time of maximum light amplification
                                             should range based on life of survey
        u_0     [0, 2, 0.01]            -impact parameter
        x       (0,1,0.05)             -percentage of distance from Dl to Ds (Dl/Ds)
        bandpass   {g,r,i,z,Y}
        v_t     [220km/s]              -transverse velocity of the lensing object
        Ds      20 kpc                 -distance to the source, based on actual calculations, not varied
        m_0                            -average magnitude of source
        """

    """ Takes a various arguments and returns a list of Mircolensing events. """

    def get_curves(self, MJD_list, bandpass, objID, ra, dec, t_eff = 1, m_0=30, Ds=5,curve_type =1):
        v_t = 220
        index = 0
        t_0 = 56961        
        # for t_0 in range(int(min(MJD_list)-365), int(max(MJD_list)+365), 20): #30        
        """
        for t_0 in range(int(min(MJD_list))-100, int(max(MJD_list))+100, 200):        
            print "t_0 index: ", t_0, index
            for u_0 in np.arange(0,2,.1):
        for number in range(0, 1, 1):
            for u_0 in np.linspace(0, 2, 10):
                for x in np.arange(0.1,1,.1):
                    for M in range(10,100,10): 
        """
        #for t_0 in range(int(min(MJD_list))-100, int(max(MJD_list))+100, 200): #oct 31, 2014 (MJD = 56961)        
        for marika in range(0,1,1):
            t_0 = 56961
            print "t_0 index: ", t_0, index
            for u_0 in np.arange(0,2,.2):
                for x in np.arange(0.1,1,0.1):
                    for M in range(10,101,10): 
                        self.lightcurve.append(MicroLensingGenerator.GenerateMLEvent(
                            t_0, u_0, v_t, M, Ds, x, MJD_list, m_0, bandpass, objID, ra, dec, t_eff, curve_type))
                        #self.final_mag_list.append(self.lightcurve)
                        index += 1
        print "total index:", index
        """
        for marika in range(0,1,1):
            t_0 = 56961
            M = 50
            u_0 = .5
            #x = .5
            print "t_0 index: ", t_0, index
            for x in np.arange(0.1,1,.1):
                self.lightcurve.append(MicroLensingGenerator.GenerateMLEvent(
                            t_0, u_0, v_t, M, Ds, x, MJD_list, m_0, bandpass, t_eff, curve_type))
                index += 1
        """
        return self.lightcurve

