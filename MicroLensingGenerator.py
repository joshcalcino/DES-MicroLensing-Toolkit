import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

#from desqcat import load_hpx, load_cat, load_cat_epochs
##Celeste 11:04 36 June
class GenerateMicrolensingEvent(object):
    """This class will generate a random microlensing event. All of the parameters are randomly generated."""

    def __init__(self, t_0, p, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type):
        print "init 1a"
        self.M_lens = M_lens                        #Mass of the lens, solar masses 
        self.Ds = Ds                                #Distance between observer and source, kpc
        self.ImpactParameter = p                    #Minimum distance between object and undeflected line of sight at time, t_0
        self.t_0 = t_0                              #Time of maximum light distortion
        self.m_0 = m_0                              #Average magnitude of individual source
        self.t_eff = t_eff                          
        self.times = self.get_MJD_list(MJD_list) ##formerly []         #list of times source was observed
        self.x = x                                  #ratio between Ds and Dl; percentage of Dl compared to Ds
        self.V_t = V_t*(1.496e-8)*86400             #transverse velocity, possibly determined via Gaia, in AU/day, input in km/s
        self.t_E = self.get_t_E()                   #lensing timescale
        self.curve_type = curve_type
        self.A = self.get_A()
        print "init 1b"

    def get_MJD_list(self, identity):
        sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
        from desqcat import load_hpx, load_cat, load_cat_epochs
        mpl.rcParams['figure.figsize'] = (8, 5)
        hpix = identity
        cat_wide = load_cat(hpix)
        cat = load_cat(hpix, long=True)
        epochs = {}
        epochs['g'] = load_hpx(hpix, 'g') #change later to loop through all bands
        cat_cols = ['QUICK_OBJECT_ID', 'RA', 'DEC', 'HPX2048', 'BAND',
            'NEPOCHS', 'FLAGS', 'WAVG_FLAGS', 'EXPNUM',
            'WAVG_MAG_PSF', 'WAVG_MAGERR_PSF',
            'WAVG_MAG_AUTO', 'WAVG_MAGERR_AUTO']
        epoch_cols = ['QUICK_OBJECT_ID', 'EXPNUM', 'MJD_OBS', 'EXPTIME', 'BAND', 'T_EFF',
              'MAG_PSF', 'MAGERR_PSF', 'MAG_AUTO', 'MAGERR_AUTO']
        ecat = load_cat_epochs(hpix, cat_cols, epoch_cols)
        ecat = ecat.query('MAG_PSF < 30')
        list_times = ecat['QUICK_OBJECT_ID']
        obj_expnum_counts = ecat[['QUICK_OBJECT_ID', 'EXPNUM', 'BAND']].groupby(['QUICK_OBJECT_ID', 'EXPNUM'], as_index=False).count()
        obj_expnum_counts.columns = ['QUICK_OBJECT_ID', 'EXPNUM', 'COUNTS']
        duplicated_objects = obj_expnum_counts.QUICK_OBJECT_ID[obj_expnum_counts.COUNTS>1]
        ecat = ecat[np.in1d(ecat.QUICK_OBJECT_ID.values, duplicated_objects.values, invert=True)]
        
        quick_id = list_times[0]
        myobj_df = ecat.loc[quick_id]
        myobj_r = ecat.query("QUICK_OBJECT_ID==" + str(quick_id) + " & BAND=='r'")[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        index = 0
        quick_id = list_times[index]
        while len(myobj_r['MJD_OBS']) == 1:
            index = index + 1
            quick_id = list_times[index]
            myobj_df = ecat.loc[quick_id]
            myobj_r = ecat.query("QUICK_OBJECT_ID==" + str(quick_id) + " & BAND=='r'")[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]

        len_id = len(myobj_r['MJD_OBS'])

        quick_id_array = np.zeros(len_id)

        for time in range(0, len_id):
            quick_id_array[time] = myobj_r['MJD_OBS'][time]
        print("times are: ")
        print(quick_id_array)
        return quick_id_array        

    def get_A(self):
        print "get_A 2"
        t = self.times
        u = self.get_u(t)
        A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
        print "u is: " + str(u)
        print "A is: " + str(A)
        print "A is " + str(len(A)) + " long"
        return A        

    def get_curve_type(self):
        print "get_curve_type 3"
        curve_type = self.curve_type
        #Put curve type loop here, return value of curve type
        #1. Paci 2. Ellipse 3. Parallax 4. Cluster
        return 0

    def get_r_E(self):  # r_E is the Einstein ring radius in units of.. not sure yet
        print "get_r_E 4"
        m_denominator=1.0
        d_denominator=10.
        M = self.M_lens
        Ds = self.Ds
        x = self.x
        r_E = 4.54 * np.sqrt(M / m_denominator)*np.sqrt(Ds/d_denominator)*(np.sqrt((x*(1-x)))/(0.5)) #in AU
        """
        r_E = 4.848e-9*Ds*( 0.902 * np.sqrt(M / const.M_sun.value) * np.sqrt(10000 / (x*Ds)) * np.sqrt(
            1 - x))  # in milli arcseconds, now in whatever units Ds is in (km)
        """
        print "r_E is: " + str(r_E)
        return r_E

    def get_t_E(self):  # time it takes the source to move a distance equal to the Einstein ring radius
        print "get_t_E 5"
        t_E = self.get_r_E() / self.V_t #needs the same units as r_E to get out seconds
        print 't_E is'
        print t_E
        return t_E

    def get_u(self, t):
        print "get U 6"
        p = self.ImpactParameter
        t_E = self.get_t_E()
        print 't_E is: ' + str(t_E)
        t_0 = self.t_0
        u = np.sqrt(p ** 2 + ((t - t_0) / t_E) ** 2)
        print "u is " + str(len(u)) + " long"
        return u

    def get_delta_mag(self, t):  # change in the magnitude of the star due to the lensing
        print "get delta mag 7"
        A = self.A
        delta_mag = 2.5 * np.log10(A)
        print "delta_mag is: ",delta_mag
        return delta_mag

    def generate_data(self):
        print "generate data 8"
        t = self.times
        delta_mag_list = self.get_delta_mag(t)
        mag_list = self.m_0 + delta_mag_list
        final_mag_list = mag_list + self.generate_noise(t)
        return final_mag_list

    def save_data(self, data):  # save the data as a text file
        print "save data 9"
        delta_mag = np.reshape(data['delta_mag'], (len(data['delta_mag']), 1))
        time = np.reshape(data['time'], (len(data['time']), 1))
        data = np.concatenate((delta_mag, time), axis=1)
        np.savetxt('sample_microevent.txt', data)

    def generate_noise(self, t):
        print "generate noise 10"
        noise = 0 
        return noise
