
import numpy as np
    import sys
    import pandas as pd
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    import MicroLensingGenerator
    import fake_plots
    import getData

class loop (object):
 
    #need to loop over Parameter Limits:
    #Mlens   [10,100,1]             -Mass of the deflector
    #t_0     (t_E min, t_E max, 1)  -time of maximum light amplification
    #u_0     [0,2, 0.01]            -impact parameter
    #v_t     [20,220km/s, 100]      -transverse velocity of the lensing object
    #Ds      30 kpc                 -distance to the source, based on actually calculations, not varied
    #x       (0,1,0.1)              -percentage of distance from Dl to Ds (Dl/Ds)
    #m_0     (0,30,1)               -average magnitude of source
    #t_E     (438, 1424,1)          -converted from (1.2-3.9) yrs to MJD    
 
    def star_object(hpix=11373)
        #assuming getData has already been called
        #get the MJD_list
        #loop over all of the variables above
        #in each loop, call an instance of object "MicroLensingGenerator"
        Ds = 30 #kpc
        MJD_list = getData.getMJD_list();
        t_eff = 0 #need to pick a variable
        curve_type = 1 #PAC curve
        lightcurve = []
        index = 1
        for m in range(10, 100,1):
            index = index + 1
            for t_0 in range(t1, t2, 0.01):
                index = index +  1
                for u_0 in range(0,2,0.01):
                    index = index + 1
                    for v_t in range(20,220, 100):
                        index = index + 1
                        for x in range(1,1,0.1):
                            index = index + 1
                            for m_0 in range(0,30,1):
                                index = index + 1
                                for t_E in range(438,1424,1):
                                    index = index + 1
                                    lightcurve[index] = MicroLensingGenerator.GenerateMicrolensingEvent(t_0, u_0, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type)

        """
        for all mass from 10 to 100
            for all t_0 from min to max
                for all 
                    for all
                        lightcurve1 = microlesningenerator.generatoonioanoxa( m1, t0, )
        """
        return 0 


