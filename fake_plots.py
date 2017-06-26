import MicroLensingGenerator
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
import numpy as np
import astropy.table as t
import matplotlib.image as img
from scipy.optimize import newton
from scipy.interpolate import interp1d

<<<<<<< HEAD

x = MicroLensingGenerator.GenerateMicrolensingEvent(1, 0.1, 220, 30, 4.5, 0.5, np.arange(1, 100, 0.01), 12, 1, 1)

t = x.times
u = x.get_u(t)
A = x.A
delta_mag = x.get_delta_mag(t)
=======
x = MicroLensingGenerator.GenerateMicrolensingEvent(50, 0.1, 1e-15, 30, 4.5, 0.5, np.arange(1, 100, 0.01), 12, 1, 1)
                                                   #t_0,u_0, v_t, M_lens,D_s, x,  t[list],     m_0, t_eff, curve_type  
t = x.times
u = x.get_u(t)
A = x.A
delta_mag = x.get_delta_mag(t)  # change in the magnitude of the star due to the lensing
>>>>>>> 6e85e7e20e88e8363d22d3f0fc63bb9e9a8b9856
interp = interp1d(u, A, bounds_error =  False, kind = 'linear')

plt.plot(u, 2.5*np.log10(interp(u)))
#plt.plot(u, interp2(delta_mag))
plt.ylim(0, 0.5)
plt.xlabel("Einstein Radii")
plt.ylabel("Magnitude Difference")
plt.grid()
plt.title("Fake graph with chosen inputs")
plt.show()

