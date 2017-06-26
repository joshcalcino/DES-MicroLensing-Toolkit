import MicroLensingGenerator
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from scipy.interpolate import interp1d


x = MicroLensingGenerator.GenerateMicrolensingEvent(50, 0.1, 1e-15, 30, 4.5, 0.5, np.arange(1, 100, 0.01), 12, 1, 1)

t = x.times
u = x.get_u(t)
A = x.A
interp = interp1d(u, A, bounds_error =  False, kind = 'linear')

plt.plot(u, 2.5*np.log10(interp(u)))
plt.ylim(0, 0.5)
plt.xlabel("Einstein Radii")
plt.ylabel("Magnitude Difference")
plt.grid()
plt.show()


