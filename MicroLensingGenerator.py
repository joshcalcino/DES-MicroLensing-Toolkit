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
            It takes 11 parameters: 
            t_0, u_0, v_t, M_lens, Ds, x, MJD_list, m_0, bandpass, t_eff, curve_type
    """

    def __init__(self, t_0, u_0, V_t, M_lens, Ds, x, MJD_list, m_0, bandpass, objID, ra, dec, t_eff = 0, curve_type = 1):
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
        self.m_0 = m_0                      #magnitude per observation, list   
        self.t_eff = t_eff
        self.bandpass = bandpass
        self.curve_type = curve_type
        self.delta_mag = self.get_delta_mag()   #calculates change in magnitude
        self.interp_r, self.interp_Y, self.interp_g, self.interp_z, self.interp_i = self.get_error_files()
        self.light_curve = self.generate_data(bandpass) #list of mag at times, accounting for delta and initial magnitudes
        self.light_curve_error = self.generate_error(bandpass)
        self.quickid = np.ones(bandpass.size)*objID
        self.ra = np.ones(bandpass.size)*ra
        self.dec = np.ones(bandpass.size)*dec

    def get_error_files(self):
        fr = pickle.load(open("magerr_model_r.pickle", 'rb'))
        fY = pickle.load(open("magerr_model_Y.pickle", 'rb'))
        fg = pickle.load(open("magerr_model_g.pickle", 'rb'))
        fz= pickle.load(open("magerr_model_z.pickle", 'rb')) 
        fi = pickle.load(open("magerr_model_i.pickle", 'rb'))
        er, eY, eg, ez, ei = self.get_interps(fr, fY, fg, fz, fi)
        return er, eY, eg, ez, ei

    def get_interps(self, er, eY, eg, ez, ei): 
        ir = interp1d(er[0], er[1], bounds_error = False)    
        iY = interp1d(eY[0], eY[1], bounds_error = False)    
        ig = interp1d(eg[0], eg[1], bounds_error = False)    
        iz = interp1d(ez[0], ez[1], bounds_error = False)    
        ii = interp1d(ei[0], ei[1], bounds_error = False)    
        return ir, iY, ig, iz, ii 

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
    def generate_error(self, bandpass):
        hack_index = 0
        ir, = np.where(bandpass=='r')
        iY, = np.where(bandpass=='Y')
        ig, = np.where(bandpass=='g')
        iz, = np.where(bandpass=='z')
        ii, = np.where(bandpass=='i')
        nr = self.interp_r(self.light_curve[ir])
        nY = self.interp_Y(self.light_curve[iY])
        ng = self.interp_g(self.light_curve[ig])
        nz = self.interp_z(self.light_curve[iz])
        ni = self.interp_i(self.light_curve[ii])
        #Below is a HACK. FIX THIS!!!
        ix = ((ng < 0.005) | np.isnan(ng)) 
        if np.any(ix): ng[ix] = 0.005
        ix = ((nr < 0.005) | np.isnan(nr)) 
        if np.any(ix): nr[ix] = 0.005
        ix = ((ni < 0.005) | np.isnan(ni)) 
        if np.any(ix): ni[ix] = 0.005
        ix = ((nz < 0.005) | np.isnan(nz)) 
        if np.any(ix): nz[ix] = 0.005
        ix = ((nY < 0.005) | np.isnan(nY)) 
        if np.any(ix): nY[ix] = 0.005
        #Above is a HACK. FIX THIS!!
        error = np.zeros(bandpass.size)
        error[ig] = ng
        error[ir] = nr
        error[ii] = ni
        error[iz] = nz
        error[iY] = nY
        return error

    """ generate_data(): Calculates the resulting change in magnitude of the source (including compensation for noise) given initial mag and change in mag. """
    def generate_data(self, bandpass):
        ir, = np.where(bandpass=='r')
        iY, = np.where(bandpass=='Y')
        ig, = np.where(bandpass=='g')
        iz, = np.where(bandpass=='z')
        ii, = np.where(bandpass=='i')
        final_mag_list_r = self.m_0[ir] + self.delta_mag[ir] # + self.generate_noise
        final_mag_list_Y = self.m_0[iY] + self.delta_mag[iY] # + self.generate_noise
        final_mag_list_g = self.m_0[ig] + self.delta_mag[ig] # + self.generate_noise
        final_mag_list_z = self.m_0[iz] + self.delta_mag[iz] # + self.generate_noise
        final_mag_list_i = self.m_0[ii] + self.delta_mag[ii] # + self.generate_noise
        #return final_mag_list_r, final_mag_list_Y,final_mag_list_g,final_mag_list_z,final_mag_list_i 
        final_mag_list = np.zeros(bandpass.size)
        final_mag_list[ig] = self.m_0[ig]+self.delta_mag[ig]
        final_mag_list[ir] = self.m_0[ir]+self.delta_mag[ir]
        final_mag_list[ii] = self.m_0[ii]+self.delta_mag[ii]
        final_mag_list[iz] = self.m_0[iz]+self.delta_mag[iz]
        final_mag_list[iY] = self.m_0[iY]+self.delta_mag[iY]
        return final_mag_list
    """ 
    was in the structure as self.all_curves 
    def get_curves(self):
        lc = self.m_0 + self.delta_mag
        return lc
    """
    def get_curve_type(self):
        #Put curve type loop here, return value of curve type
        #1. Paci 2. Ellipse 3. Parallax 4. Cluster
        return 1

    def save_data(self, data):  # save the data as a text file
        delta_mag = np.reshape(data['delta_mag'], (len(data['delta_mag']), 1))
        time = np.reshape(data['time'], (len(data['time']), 1))
        data = np.concatenate((delta_mag, time), axis=1)
        np.savetxt('sample_microevent.txt', data)

