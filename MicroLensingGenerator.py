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
            self.t_0 = t_0                              #Time of maximum light distortion
            self.ImpactParameter = p                    #min dist betwn obj and line of sight at t_0, unitless
            print(V_t)
            self.M_lens = M_lens                        #Lens mass,             solar masses 
            #self.V_t = V_t*(1.496e-8)*86400             #transverse velocity, in AU/day, input in km/s
            self.V_t = V_t*5.775e-4             #transverse velocity, in AU/day, input in km/s
            self.Ds = Ds                                #Dist to source         kpc
            self.x = x                                  #% of Dl compared to Ds  (0, 1)
            #self.times = self.get_MJD_list(MJD_list)    #list of times source was observed
            self.times = MJD_list                      #list of times source was observed, in days
            #self.hpix = MJD_list 
            self.m_0 = m_0                              #Avg magnitude of src   
            self.t_eff = t_eff                          
            self.t_E = self.get_t_E()                   #lensing timescale
            self.A = self.get_A()
            print("velocity: ", self.V_t)                       # Pac curve
    

        def get_A(self):
            t = self.times
            u = self.get_u(t)
            A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
            print("A: ", A)
            return A        
    

        def get_curve_type(self):
            curve_type = self.curve_type
            #Put curve type loop here, return value of curve type
            #1. Paci 2. Ellipse 3. Parallax 4. Cluster
            return 0
    

        def get_r_E(self):  # r_E is the Einstein ring radius in units of.. not sure yet
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
            print("r_E: ", r_E)
            return r_E
    

        def get_t_E(self):  # time it takes the source to move a distance equal to the Einstein ring radius
            t_E = self.get_r_E() / self.V_t #needs the same units as r_E to get out seconds
            print("t_E: ", t_E)
            return t_E
    

        def get_u(self, t):
            p = self.ImpactParameter
            t_E = self.get_t_E() #convert to days for MJD division
            t_0 = self.t_0
            u = np.zeros(len(t))
            for i in range(0, len(t)):
                sign = 1
                if t[i] - t_0 < 0:
                    sign = 1
                u[i] = sign * np.sqrt(p ** 2 + ((t[i] - t_0) / t_E) ** 2)
            print("u in deep class: ", u)
            print("p :", p)
            print("t_0: ", t_0)
            return u
    

        def get_delta_mag(self, t):  # change in the magnitude of the star due to the lensing
            A = self.A
            delta_mag = 2.5 * np.log10(np.absolute(A))
            return delta_mag
    

        def generate_data(self):
            t = self.times
            delta_mag_list = self.get_delta_mag(t)
            mag_list = self.m_0 + delta_mag_list
            final_mag_list = mag_list + self.generate_noise(t)
            return final_mag_list
    

        def save_data(self, data):  # save the data as a text file
            delta_mag = np.reshape(data['delta_mag'], (len(data['delta_mag']), 1))
            time = np.reshape(data['time'], (len(data['time']), 1))
            data = np.concatenate((delta_mag, time), axis=1)
            np.savetxt('sample_microevent.txt', data)
    

        def generate_noise(self, t):
            noise = 0 
            return noise

