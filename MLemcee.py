from testfitter import *
import logging
import numpy as np
import emcee
from chainconsumer import ChainConsumer
from MicroLensingGenerator import *
import matplotlib.pyplot as plt

MicroEvent = GenerateMicrolensingEvent()
#data = MicroEvent.simple_fake()
data = MicroEvent.generate_data()
delta_mag_data = data['delta_mag']
t = data['time']
errors = data['errors']


def get_delta_mag(t, p, t_0, t_max):
    #t_max = 2.
    #t_0 = 0.5
    u = np.sqrt(p ** 2 + ((t - t_max) / t_0) ** 2)
    A = (u ** 2 + 2) / (u * np.sqrt(u ** 2 + 4))
    delta_mag = 2.5 * np.log10(A)
    return delta_mag


def log_posterior(params):
    p = params[0]
    t_0 = params[1]
    t_max = params[2]
    if p < 0 or t_0 < 0: #  or t_max < 0:
        return -np.inf
    delta_mag = get_delta_mag(t, p, t_0, t_max)
    chi2 = ((delta_mag_data - delta_mag)**2) / errors**2
    L = -0.5 * chi2.sum()
    if isinstance(L, complex) or not isinstance(L, float):
        return -np.inf
    return L


def start(num_walkers=1):
    x = np.abs(np.random.normal(loc=1., scale=2., size=(num_walkers, 3)))
    print "Starting positions are:"
    print x
    return x

# logging.basicConfig(level=logging.DEBUG)
# sampler = EnsembleSampler(save_interval=60, temp_dir="here", num_steps=10000)
# output = sampler.fit({"log_posterior": log_posterior, "start": start, "uid": "hello"})

p_real = MicroEvent.ImpactParameter
t_0_real = MicroEvent.t_0
t_max_real = MicroEvent.t_max

nwalkers = 100
ndim = 3

start_pos = start(nwalkers)

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[])
sampler.run_mcmc(start_pos, 1000)

output = sampler.flatchain[100:]

c = ChainConsumer().add_chain(output, parameters=["$p$", "$t_0$", "$t_{\\rm max}$"])
truth = [p_real, t_0_real, t_max_real]
c.plot(filename="data_constraints.png", truth=truth)
print "printing c.get_summary()"
print(c.get_summary())


fit_data = c.get_summary()
# print(fit_data['$p$'])
# print 'ayy lmao'
# print(fit_data['$p$'][1])

p_fit, t_0_fit, t_max_fit = fit_data['$p$'][1], fit_data['$t_0$'][1], fit_data['$t_{\\rm max}$'][1]
tt = np.linspace(np.min(t), np.max(t), 200)
delta_mag_real = get_delta_mag(tt, p_real, t_0_real, t_max_real)
delta_mag_fit = get_delta_mag(tt, p_fit, t_0_fit, t_max_fit)

fig = plt.figure(2)
under, = plt.plot(tt, delta_mag_real, label='Underlying Light-curve')
fitted, = plt.plot(tt, delta_mag_fit, label='Fitted Light-curve')
data = plt.errorbar(t, delta_mag_data, yerr=errors, fmt='o', label='Simulated Light-curve')
plt.legend(handles=[under, fitted, data], loc=1)

plt.show()
