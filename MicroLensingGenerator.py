import numpy as np
##Celeste 11:04 36 June
class GenerateMicrolensingEvent(object):
    """This class will generate a random microlensing event. All of the parameters are randomly generated."""

    def __init__(self, t_0, p, V_t, M_lens, Ds, x, MJD_list, m_0, t_eff, curve_type):
        print "init 1a"
        self.t_0 = t_0                              #Time of maximum light distortion
        self.ImpactParameter = p                    #min dist betwn obj and line of sight at t_0, unitless
        self.V_t = V_t*(1.496e-8)*86400             #transverse velocity, in AU/day, input in km/s
        self.M_lens = M_lens                        #Lens mass,             solar masses 
        self.Ds = Ds                                #Dist to source         kpc
        self.x = x                                  #% of Dl compared to Ds  (0, 1)
        self.times = MJD_list ##formerly []         #list of times source was observed, in days
        self.m_0 = m_0                              #Avg magnitude of src   
        self.t_eff = t_eff                          
        self.curve_type = curve_type
        self.t_E = self.get_t_E()                   #lensing timescale
        self.A = self.get_A()
        print "init 1b"

    def get_A(self):
        print "get_A 2"
        t = self.times
        u = self.get_u(t)
        A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
        print "A is " + str(len(A)) + " long"
        print "u is " + str(len(u)) + " long"
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
        return t_E

    def get_u(self, t):
        print "get U 6"
        p = self.ImpactParameter
        t_E = self.get_t_E()
        print 't_E is: ' + str(t_E)
        t_0 = self.t_0
        u = np.sqrt(p ** 2 + ((t - t_0) / t_E) ** 2)
        return u

    def get_delta_mag(self, t):  # change in the magnitude of the star due to the lensing
        print "get delta mag 7"
        A = self.A
        delta_mag = 2.5 * np.log10(A)
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
