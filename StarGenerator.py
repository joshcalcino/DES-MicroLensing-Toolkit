import numpy as np
import astropy.constants as const
from scipy.integrate import quad
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


class GenerateStandardStar(object):

    """Generate a catalogue of normal stars based off of values from:
    http://adsabs.harvard.edu/abs/1991Ap%26SS.181..313D
    Assuming that a star is a black body to determine the luminosity in different bands..

    Stars are assumed to be on the main sequence, or else this would get too complicated."""

    '''Define constants that are being used throughout the code.'''

    h = const.h.value
    c = const.c.value
    k = const.k_B.value
    R_sun = const.R_sun.value
    L_sun = const.L_sun.value
    M_sun = const.M_sun.value

    def __init__(self):
        self.M = np.random.uniform(low=0.5, high=17.5) * const.M_sun.value
        self.R = self.get_radius(self.M)
        self.L = self.get_luminosity(self.M)
        self.T = self.get_temperature()
        self.r = self.get_distance()

    def get_distance(self):
        """For now this will just generate a uniform random number, will eventually replace
        with something to resemble the density of stars in th DES field"""
        dist = np.random.uniform(low=10, high=10) * const.pc.value
        return dist

    def get_radius(self, M):
        M = M + M * self.generate_fluctuations()
        if M < 1.66:
            R = 1.06 * (M/self.M_sun) ** 0.945 * self.R_sun
        else:
            R = 1.33 * (M/self.M_sun) ** 0.555 * self.R_sun
        return R

    def get_luminosity(self, M):
        if M < 0.7:
            M = M + M * self.generate_fluctuations()
            L = 0.35 * (M/const.M_sun.value) ** 2.62 * const.L_sun.value
        else:
            M = M + M * self.generate_fluctuations()
            L = 1.02 * (M/const.M_sun.value) ** 3.92 * const.L_sun.value
        return L

    def get_temperature(self):
        T = (self.L/(4*np.pi*const.sigma_sb.value*self.R**2)) ** (1./4)
        return T

    def get_flux(self, wavelength):
        flux = (2 * np.pi * self.h * self.c ** 2 / wavelength ** 5) / (np.exp(self.h * self.c / (wavelength * self.k
            * self.T)) - 1) * (self.R / self.r) ** 2
        return flux

    def generate_fluctuations(self):
        fluct = np.random.normal(0, 0.02)
        return fluct


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
        if wavelength > self.g_band_max or wavelength < self.g_band_min:
            return 0.
        else:
            return self.g_band(wavelength)

    def get_r_band_sf(self, wavelength):
        if wavelength > self.r_band_max or wavelength < self.r_band_min:
            return 0.
        else:
            return self.r_band(wavelength)

    def get_i_band_sf(self, wavelength):
        if wavelength > self.i_band_max or wavelength < self.i_band_min:
            return 0.
        else:
            return self.i_band(wavelength)

    def get_z_band_sf(self, wavelength):
        if wavelength > self.z_band_max or wavelength < self.z_band_min:
            return 0.
        else:
            return self.z_band(wavelength)

    def get_g_monochromatic_flux(self, wavelength, flux):
        g_band_sf = self.get_g_band_sf(wavelength)
        return flux(wavelength) * g_band_sf

    def get_r_monochromatic_flux(self, wavelength, flux):
        r_band_sf = self.get_r_band_sf(wavelength)
        return flux(wavelength) * r_band_sf

    def get_i_monochromatic_flux(self, wavelength, flux):
        i_band_sf = self.get_i_band_sf(wavelength)
        return flux(wavelength) * i_band_sf

    def get_z_monochromatic_flux(self, wavelength, flux):
        z_band_sf = self.get_z_band_sf(wavelength)
        return flux(wavelength) * z_band_sf

    def get_apparent_magnitudes(self, flux):
        """This function will generate the apparent magnitudes that we measure, with some added noise and error bars.
        The Cx constants represent the arbitrary 0 point for the magnitude scale, traditionally taken as the flux
        coming from the star Vega. Since the only variable needed to simulate a star in GenerateStandardStar is the
        mass of the star, the Cx constants are determined by the flux coming from a star with the same mass as Vega."""
        Cg = -23.914949483915933
        Cr = -25.122907959259308
        Ci = -24.111869728236627
        Cz = -25.367460437973467
        A1 = quad(self.get_i_monochromatic_flux, 0, 1., args=flux, points=[800 * 1e-9], epsrel=1.e-5)
        A2 = quad(self.get_i_monochromatic_flux, 0, 1., args=flux, points=[800 * 1e-9], epsrel=1.e-5)
        g_band, tmp = -2.5 * np.log10(quad(self.get_g_monochromatic_flux, 0, 1., args=flux, points=[500*1e-9])) + Cg
        r_band, tmp = -2.5 * np.log10(quad(self.get_r_monochromatic_flux, 0, 1., args=flux, points=[600*1e-9])) + Cr
        i_band, tmp = -2.5 * np.log10(quad(self.get_i_monochromatic_flux, 0, 1., args=flux, points=[800*1e-9])) + Ci
        z_band, tmp = -2.5 * np.log10(quad(self.get_z_monochromatic_flux, 0, 1., args=flux, points=[900*1e-9])) + Cz
        print A1
        print A2
        # t = self.generate_observation_times()
        # g_band, r_band = g_band + self.generate_noise(t), r_band + self.generate_noise(t)
        # i_band, z_band = i_band + self.generate_noise(t), z_band + self.generate_noise(t)
        return [g_band, r_band, i_band, z_band]

    def generate_noise(self, t):
        """Need to decide whether or not the noise is introduced here or when the star is generated.."""
        noise = np.random.normal(0.0, 0.02, len(t))
        return noise

    def generate_errors(self, t):
        noise = np.random.normal(0.0, 0.02, len(t))
        return noise

    def generate_observation_times(self):
        t = np.linspace(0, 5, 10)
        return t




Star = GenerateStandardStar()
print 'The mass is'
print Star.M / const.M_sun.value
print 'The Luminsoity is'
print Star.L / const.L_sun.value
print 'The temperature is'
print Star.T

Observations = GenerateObservations()
print Observations.get_apparent_magnitudes(Star.get_flux)


# xdata = np.linspace(0, 1e-6, 1000)
# g_band, r_band = Observations.get_g_band_sf(xdata), Observations.get_r_band_sf(xdata)
# i_band, z_band = Observations.get_i_band_sf(xdata), Observations.get_z_band_sf(xdata)
#
# plt.plot(xdata, g_band, 'g', xdata, r_band, 'r', xdata, i_band, 'b', xdata, z_band, 'k')
# plt.show()

"""
    Assuming that the sensitivity functions are step functions, so pretty much a 0th order approximation.
    Values are obtained from here:
    http://www.stsci.edu/hst/wfc3/documents/handbooks/currentIHB/c06_uvis06.html

    Should find out what these sensitivity functions actually are.

    def get_g_monochromatic_flux(self, wavelength):
        print wavelength
        h = const.h.value
        c = const.c.value
        k = const.k_B.value
        peak = 477.3 * 10**(-9)
        width = 134.4 * 10**(-9)
        g = np.exp(- (wavelength - peak) ** 2 / (2 * width))
        flux = (2 * np.pi * h * c ** 2 / wavelength ** 5) / (np.exp(h * c / (wavelength * k * self.T)) - 1) * (
            self.R / self.r) ** 2
        return flux * g

    def get_r_monochromatic_flux(self, wavelength):
        h = const.h.value
        c = const.c.value
        k = const.k_B.value
        peak = 624.2 * 10**(-9)
        width = 146.3  * 10**(-9)
        r = np.exp(- (wavelength - peak) ** 2 / (2 * width))
        flux = (2 * np.pi * h * c ** 2 / wavelength ** 5) / (np.exp(h * c / (wavelength * k * self.T)) - 1) * (
            self.R / self.r) ** 2
        return flux * r

    def get_i_monochromatic_flux(self, wavelength):
        h = const.h.value
        c = const.c.value
        k = const.k_B.value
        peak = 764.7 * 10**(-9)
        width = 117.1 * 10**(-9)
        i = np.exp(- (wavelength - peak) ** 2 / (2 * width))
        flux = (2 * np.pi * h * c ** 2 / wavelength ** 5) / (np.exp(h * c / (wavelength * k * self.T)) - 1) * (
            self.R / self.r) ** 2
        return flux * i

    def get_z_monochromatic_flux(self, wavelength):
        h = const.h.value
        c = const.c.value
        k = const.k_B.value
        peak = 916.6 * 10**(-9)
        width = 118.2 * 10**(-9)
        z = np.exp(- (wavelength - peak) ** 2 / (2 * width))
        flux = (2 * np.pi * h * c ** 2 / wavelength ** 5) / (np.exp(h * c / (wavelength * k * self.T)) - 1) * (
            self.R / self.r) ** 2
        return flux * z

    def get_apparent_magnitudes(self):
        Cg = 0.
        Cr = 0.
        Ci = 0.
        Cz = 0.
        g_band, tmp = -2.5*np.log10(quad(self.get_g_monochromatic_flux, 0, 1., points=[600*10**(-9)])) + Cg
        r_band, tmp = -2.5 * np.log10(quad(self.get_r_monochromatic_flux, 0, 1., points=[600*10**(-9)])) + Cr
        i_band, tmp = -2.5 * np.log10(quad(self.get_i_monochromatic_flux, 0, 1., points=[600*10**(-9)])) + Ci
        z_band, tmp = -2.5 * np.log10(quad(self.get_z_monochromatic_flux, 0, 1., points=[600*10**(-9)])) + Cz
        return [g_band, r_band, i_band, z_band]"""