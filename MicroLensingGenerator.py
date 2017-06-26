import numpy as np
import astropy.constants as const

class GenerateMicrolensingEvent(object):
    """This class will generate a random microlensing event. All of the parameters are randomly generated.
    def get_pi_rel(self):
        pi_rel = np.random.uniform(low=0.01, high=0.1) # pi_rel (\pi_{\rm rel}) is the relative parallax between the
        # lens and the source, measured in milli-arcseconds.
        return pi_rel
    """

    def __init__(self, t_0, p, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type):
        #self.Ds, self.Dl, self.Drel = self.get_Dist()
        self.curve_type = curve_type
        self.M_lens = M_lens
        self.Ds = Ds
        self.ImpactParameter = self.get_ImpactParameter(p)
        self.t_0 = self.get_t_0(t_0)
        self.m_0 = m_0
        self.t_eff = t_eff
        self.times = self.get_times(MJD_list) ##formerly []
        self.x = x
        self.V_t = V_t
        self.t_E = self.get_t_E() 
   
    def get_curve_type(self):
        curve_type = self.curve_type
        #Put curve type loop here, return value of curve type
        #1. Paci 2. Ellipse 3. Parallax 4. Cluster
        return 0

    def get_Mass(self):
        M = np.random.uniform(low=0.1,
                              high=100) * const.M_sun.value  # the mass of the lens, relative to solar mass.
        return M

    def get_t_0(self, t0):  # time of maximum approach, in years since the start of DES.. for now
        t_0 = t0
        print 't_0 is'
        print t_0
        return t_0

    def get_ImpactParameter(self, impact):
        p = impact
        print 'p is'
        print p
        return p

    def get_r_E(self):  # r_E is the Einstein ring radius in units of.. not sure yet
        M = self.M_lens
        Ds = self.Ds
        x = self.x
        # r_E = np.sqrt( (4*const.G * M)/const.c ** 2 * ( (Ds - Dl)/(Ds*Dl) ) )
        r_E = 0.902 * np.sqrt(M / const.M_sun.value) * np.sqrt(10000 / (x*Ds)) * np.sqrt(
            1 - x)  # in milli arcseconds
        return r_E


    def get_t_E(self):  # time it takes the source to move a distance equal to the Einstein ring radius
        t_E = self.get_r_E() / self.V_t #needs the same units as r_E to get out seconds
        print 't_E is'
        print t_E
        return t_E

    def get_u(self, t):
        p = self.get_ImpactParameter(self.ImpactParameter)
        t_E = self.get_t_E()
        t_0 = self.get_t_0(self.t_0)
        u = np.sqrt(p ** 2 + ((t - t_0) / t_E) ** 2)
        return u

    def get_final_mag(self, t):
        final_mag = self.get_delta_mag(t)+self.m_0
        return final_mag

    def get_delta_mag(self, t):  # change in the magnitude of the star due to the lensing
        u = self.get_u(t)
        A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
        delta_mag = 2.5 * np.log10(A)
        return delta_mag

    def get_times(self, MJD_list):
        self.times = MJD_list
        return self.times

    def generate_data(self):
        t = self.get_times(self.times)
        delta_mag_list = self.get_delta_mag(t)
        mag_list = self.m_0 + delta_mag_list
        final_mag_list = mag_list + self.generate_noise(t)
        return final_mag_list

    def simple_fake(self):
        t = self.times()
        p = 0.7
        t_0 = 2.
        t_E = 0.5
        u = np.sqrt(p ** 2 + ((t - t_0) / t_E) ** 2)
        A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
        delta_mag = 2.5 * np.log10(A)
        data = {"delta_mag": delta_mag, "time": t}
        return data

    def save_data(self, data):  # save the data as a text file
        delta_mag = np.reshape(data['delta_mag'], (len(data['delta_mag']), 1))
        time = np.reshape(data['time'], (len(data['time']), 1))
        data = np.concatenate((delta_mag, time), axis=1)
        np.savetxt('sample_microevent.txt', data)


    def generate_noise(self, t):
        noise = 0 
        return noise
