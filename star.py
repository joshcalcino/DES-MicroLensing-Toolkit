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
        lcID = 0
        trange = [56535, 56747, 56992, 57234, 57430] # Dates: 8.31.13, 3.31.14, 12.1.14, 7.31.15, 2.12.16   
        for t_0 in trange:
            if t_0 != 56535: continue
            print "t_0 index: ", t_0, lcID
            for u_0 in np.arange(0,2,.2):
                for x in np.arange(0.1,1,0.1):
                    for M in range(10,101,10): 
                        self.lightcurve.append(MicroLensingGenerator.GenerateMLEvent(
                            t_0, u_0, v_t, M, Ds, x, MJD_list, m_0, bandpass, objID, ra, dec, lcID, t_eff, curve_type))
                        lcID += 1
        print "total index:", lcID
        return self.lightcurve

