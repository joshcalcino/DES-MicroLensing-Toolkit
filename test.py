import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

# g = np.loadtxt('desgriz/des_g.dat')
# r = np.loadtxt('desgriz/des_r.dat')
# i = np.loadtxt('desgriz/des_i.dat')
# z = np.loadtxt('desgriz/des_z.dat')


class GenerateObservations(object):

    """Given the data from a star, this class will then generate the expected observations
    over the DES time period."""

    '''Define some constants that will be used throughout the class'''

    def __init__(self):
        """Initialise the sensitivity function for each band."""
        print 'Initialising GenerateObservations'
        self.g_band, self.g_band_min, self.g_band_max = self.load_g_band_sf()
        self.r_band, self.r_band_min, self.r_band_max = self.load_r_band_sf()
        self.i_band, self.i_band_min, self.i_band_max = self.load_i_band_sf()
        self.z_band, self.z_band_min, self.z_band_max = self.load_z_band_sf()

    def load_g_band_sf(self):
        """Loads the sensitivity function for the g_band and returns an interpolated version of it."""
        g = np.loadtxt('desgriz/des_g.dat') # the wavelength is in Angstroms so convert to meters below
        g_band = interp1d(g[:, 0]*1e-10, g[:, 1])
        return g_band, np.min(g[:, 0])*1e-10, np.max(g[:, 0])*1e-10

    def load_r_band_sf(self):
        """Loads the sensitivity function for the r_band and returns an interpolated version of it."""
        r = np.loadtxt('desgriz/des_r.dat') # the wavelength is in Angstroms so convert to meters below
        r_band = interp1d(r[:, 0]*1e-10, r[:, 1])
        return r_band, np.min(r[:, 0])*1e-10, np.max(r[:, 0])*1e-10

    def load_i_band_sf(self):
        """Loads the sensitivity function for the i_band and returns an interpolated version of it."""
        i = np.loadtxt('desgriz/des_i.dat') # the wavelength is in Angstroms so convert to meters below
        i_band = interp1d(i[:, 0]*1e-10, i[:, 1])
        return i_band, np.min(i[:, 0])*1e-10, np.max(i[:, 0])*1e-10

    def load_z_band_sf(self):
        """Loads the sensitivity function for the g_band and returns an interpolated version of it."""
        z = np.loadtxt('desgriz/des_z.dat') # the wavelength is in Angstroms so convert to meters below
        z_band = interp1d(z[:, 0]*1e-10, z[:, 1])
        return z_band, np.min(z[:, 0])*1e-10, np.max(z[:, 0])*1e-10

    '''Now generate the actual sensitivity functions that can be called with scipy.integrate.quad'''

    def get_g_band_sf(self, wavelength):
        """The actual sensitivity function for the g band flux."""
        g_band_sf = np.zeros(len(wavelength))
        for i in range(0, len(wavelength)):
            if wavelength[i] > self.g_band_max or wavelength[i] < self.g_band_min:
                g_band_sf[i] = 0.
            else:
                g_band_sf[i] = self.g_band(wavelength[i])
        return g_band_sf

    def get_r_band_sf(self, wavelength):
        """The actual sensitivity function for the g band flux."""
        r_band_sf = np.zeros(len(wavelength))
        for i in range(0, len(wavelength)):
            if wavelength[i] > self.r_band_max or wavelength[i] < self.r_band_min:
                r_band_sf[i] = 0.
            else:
                r_band_sf[i] = self.r_band(wavelength[i])
        return r_band_sf

    def get_i_band_sf(self, wavelength):
        """The actual sensitivity function for the g band flux."""
        i_band_sf = np.zeros(len(wavelength))
        for i in range(0, len(wavelength)):
            if wavelength[i] > self.i_band_max or wavelength[i] < self.i_band_min:
                i_band_sf[i] = 0.
            else:
                i_band_sf[i] = self.i_band(wavelength[i])
        return i_band_sf

    def get_z_band_sf(self, wavelength):
        """The actual sensitivity function for the g band flux."""
        z_band_sf = np.zeros(len(wavelength))
        for i in range(0, len(wavelength)):
            if wavelength[i] > self.z_band_max or wavelength[i] < self.z_band_min:
                z_band_sf[i] = 0.
            else:
                z_band_sf[i] = self.z_band(wavelength[i])
        return z_band_sf


ob = GenerateObservations()

xnew = np.linspace(3.e-7, 1.1e-6, 1000)
g_band_sf = ob.get_g_band_sf(xnew)
r_band_sf = ob.get_r_band_sf(xnew)
i_band_sf = ob.get_i_band_sf(xnew)
z_band_sf = ob.get_z_band_sf(xnew)

plt.plot(xnew, g_band_sf, 'g', xnew, r_band_sf, 'r', xnew, i_band_sf, 'b', xnew, z_band_sf, 'k')
plt.show()

'''
x = np.linspace(0, 1500, 500)

g = g_sensitivity_function(x)
r = r_sensitivity_function(x)
i = i_sensitivity_function(x)
z = z_sensitivity_function(x)

plt.plot(x, g, 'g', x, r, 'r', x, i, 'b', x, z, 'k')
plt.show()
'''