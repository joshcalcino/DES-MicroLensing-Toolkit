import numpy as np

sigma_noise = 0.03


''' First we test to see if the light-curve has a significant bump over
what we would expect to be generated by the noise alone. '''


def std_dev_above_noise(data):
    mag = data['mag']
    avg_mag = np.average(mag)
    variance = np.sum((mag - avg_mag)**2)/len(mag)
    sigma = np.sqrt(variance)
    return sigma > 3*sigma_noise

''' '''

