import numpy as np
import sys
import pickle
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import data_practice as dp

class GenerateMLEvent(object):
    """This class will generate a microlensing event.
            It takes 10 parameters: 
            t_0, p, v_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type
    """

    def __init__(self, t_0, u_0, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff = 0, curve_type = 1, bandpass):
        self.M_lens = M_lens                #Lens mass,             solar masses 
        self.Ds = Ds                        #Dist to source         kpc
        self.x = x                          #% of Dl compared to Ds  (0, 1)
        self.r_E = self.get_r_E()
        self.V_t = V_t * 5.775e-4             #transverse velocity, in AU/day, input in km/s
        self.t_E = self.get_t_E()           #lensing timescale
        self.u_0 = u_0                      #min dist betwn obj and line of sight at t_0, unitless
        self.t_0 = t_0                      #Time of maximum light distortion
        self.times = MJD_list               #list of times source was observed, in days
        self.u = self.get_u()
        self.A = self.get_A()               #amplification of magnification
        self.m_0 = m_0                      #Avg magnitude of src   
        self.t_eff = t_eff
        self.curve_type = curve_type
        self.delta_mag = self.get_delta_mag()   #calculates change in magnitude
        self.error_file = pickle.load(open("magerr_model_{}.pickle".format(bandpass), 'rb'))
        self.interp = interp1d(self.error_file[0], self.error_file[1])
        self.light_curve = self.generate_data() #list of mag at times accounting for noise, delta and initial magnitudes
        self.generate_noise = generate_noise()

    """ get_r_E(): Calculates radius of the Einstein ring given M, Ds, and x."""
    def get_r_E(self):  # r_E is the Einstein ring radius in units of.. not sure yet
        m_denominator = 1.0
        d_denominator = 10.
        M = self.M_lens
        Ds = self.Ds
        x = self.x
        r_E = 4.54 * np.sqrt(M / m_denominator)*np.sqrt(Ds/d_denominator)*(np.sqrt((x*(1-x)))/(0.5)) #in AU
        return r_E

    """ get_t_E(): Calculates time it take source to travel the radius of the Einstein ring given r_E and V_t. """
    def get_t_E(self):  # time it takes the source to move a distance equal to the Einstein ring radius
        t_E = self.r_E / self.V_t #needs the same units as r_E to get out seconds
        return t_E

    """ get_u(): Calculates distance from source to lens' line of sight given u_0, t_E, and t_0. """
    def get_u(self):
        t = self.times
        p = self.u_0
        t_E = self.t_E              #convert to days for MJD division
        t_0 = self.t_0
        u = np.zeros(len(t))
        for i in range(0, len(t)):
            sign = 1
            if t[i] - t_0 < 0:
                sign = 1
            u[i] = sign * np.sqrt(p ** 2 + ((t[i] - t_0) / t_E) ** 2)
        return u

    """ get_A(): Calculates the amplitude of the light curve due to a lensing event given u. """
    def get_A(self):
        u = self.u
        A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
        return A        

    """ get_delta_mag(): Calculates the change in magnitude of the light source due to lensing, given A.  """
    def get_delta_mag(self):  # change in the magnitude of the star due to the lensing
        A = self.A
        delta_mag = 2.5 * np.log10(np.absolute(A))
        return delta_mag

    """ generate_noise(t): Calculates noise due to interference given t. """
    def generate_noise(self):
        self.interp(self.light_curve)
        return mag_sigma_list 

    """ generate_data(): Calculates the resulting change in magnitude of the source (including compensation for noise) given initial mag and change in mag. """
    def generate_data(self):
        final_mag_list = self.m_0 + self.delta_mag # + self.generate_noise
        return final_mag_list 

    def get_curve_type(self):
        #Put curve type loop here, return value of curve type
        #1. Paci 2. Ellipse 3. Parallax 4. Cluster
        return 1

    def save_data(self, data):  # save the data as a text file
        delta_mag = np.reshape(data['delta_mag'], (len(data['delta_mag']), 1))
        time = np.reshape(data['time'], (len(data['time']), 1))
        data = np.concatenate((delta_mag, time), axis=1)
        np.savetxt('sample_microevent.txt', data)

