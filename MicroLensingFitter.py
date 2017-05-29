from mcmc import *
from testfitter import *
import numpy as np
import MicroLensingGenerator as mlg
import matplotlib.pyplot as plt
import astropy.constants as const

class Fitted(object): # Setting up the mcmc stuff
    def __init__(self, debug=False):
        self._debug = debug
    def getIndex(self, param):
        for i, x in enumerate(self.getParams()):
            if x[0] == param:
                return i + 3
        return None
    def debug(self, string):
        if self._debug:
            print(string)
            sys.stdout.flush()
    def getNumParams(self):
        raise NotImplementedError
    def getChi2(self, params):
        raise NotImplementedError
    def getParams(self):
        raise NotImplementedError

class SimpleMicroLensingFitter(Fitted):
    def __init__(self, debug=True):
        super(SimpleMicroLensingFitter, self).__init__(debug=debug)
        self.params = [('p', 0.0, 1.4, r'$p$')]  # Telling the mcmc code what parameters we want to run over
                        # ('name', lower limit, upper limit, TeX string)

    def getNumParams(self):  # For the mcmc.py
        return len(self.params)

    def getParams(self):  # Getting the parameters
        return self.params

    def setData(self, res):
        self.res = res
        self.nsn = len(res['time'])  # the number of data points being fitted

    def get_delta_mag(self, params, t):
        p = params[0]
        t_max = 2.
        t_0 = 0.5
        u = np.sqrt(p ** 2 + ((t - t_max) / t_0) ** 2)
        A = (u**2 + 2)/( u*np.sqrt( u**2 + 4 ) )
        delta_mag = 2.5*np.log10(A)
        return delta_mag

    def getChi2(self, params):  # Calculate the chi2 without errors for now
        delta_mag = self.res['delta_mag']
        time = self.res['time']
        delta_mag_fit = self.get_delta_mag(params, time)
        chi2 = (np.exp(delta_mag_fit/2.5) - np.exp(delta_mag/2.5)/0.01)**2
        print chi2
        print chi2.sum()
        plt.plot(data['time'], data['delta_mag'], data['time'], delta_mag_fit)
        plt.show()
        return chi2.sum()
'''
class MicroLensingFitter(Fitted):
    def __init__(self, debug=True):
        super(MicroLensingFitter, self).__init__(debug=debug)
        self.params = [('M', 100, 100, r'$M$'), ('Ds', 1000, 1000, r'$D_{\rm s}$'),
                       ('Dl', 500, 500, r'$D_{\rm l}$'), ('V', 150, 150, r'$V$'),
                       ('p', 0.0, 1.4, r'$p$'), ('t_max', 2.5, 2.5, r'$t_{\rm max}$')]  # Telling the mcmc code what parameters we want to run over
                        # ('name', lower limit, upper limit, TeX string)

    def getNumParams(self):  # For the mcmc.py
        return len(self.params)

    def getParams(self):  # Getting the parameters
        return self.params

    def setData(self, res):
        self.res = res
        self.nsn = len(res['time'])  # the number of data points being fitted

    def get_Dist(self, params):
        Ds = params[1] # Ds is the distance from the observer to the source.
        Dl = params[2] # Dl is the distance from the observer to the lens.
        print Ds
        print Dl
        Drel = Ds - Dl # Drel is the relative distance between the lens and the source.
        return [Ds, Dl, Drel]

    def get_Mass(self, params):
        M = params[0] # the mass of the lens, relative to solar mass.
        return M

    def get_ImpactParameter(self):
        p = np.random.uniform(low=0.0, high=0.8) # p is the impact paramater, lower p indicates a smaller lens to source
        # relative parallax as seen from the observer. Dimensionless unit measured in units of Einstein ring radiues.
        return p

    def get_r_E(self, params): # r_E is the Einstein ring radius in units of.. not sure yet
        M = self.get_Mass(params)
        Ds, Dl, Drel = self.get_Dist(params)
        #r_E = np.sqrt( (4*const.G * M)/const.c ** 2 * ( (Ds - Dl)/(Ds*Dl) ) )
        r_E = 0.902 * np.sqrt(M/1.0) * np.sqrt(10000/Dl) * np.sqrt(1 - Dl/Ds) # in milli arcseconds
        return r_E

    def get_r_dot(self, params): # relative proper motion between the lens and the source
        Ds, Dl, Drel = self.get_Dist(params)
        V = params[3] # relative transverse velocity of the lens with respect to the source
        r_dot = 4.22 * (V/200.0) * (10000/Ds)
        return r_dot

    def get_t_0(self, params): # time it takes the source to move a distance equal to the Einstein ring radius
        t_0 = self.get_r_E(params)/self.get_r_dot(params)
        return t_0

    def get_u(self, params, t):
        p = params[4]
        t_0 = self.get_t_0(params)
        t_max = params[5]
        u = np.sqrt( p**2 + ( (t - t_max)/t_0 )**2 )
        return u

    def get_delta_mag(self, params, t): # change in the magnitude of the star due to the lensing
        u = self.get_u(params, t)
        A = (u**2 + 2)/( u*np.sqrt( u**2 + 4 ) )
        delta_mag = 2.5*np.log10(A)
        return delta_mag

    def do_plot(self):
        plt.plot(data['time'], data['delta_mag'], data['time'], delta_mag_fit)
        plt.show()

    def getChi2(self, params):  # Calculate the chi2 without errors for now
        delta_mag = self.res['delta_mag']
        time = self.res['time']
        delta_mag_fit = self.get_delta_mag(params, time)
        chi2 = ((delta_mag_fit - delta_mag)**2).sum()
        self.do_plot()
        return chi2

MicroEventGenerator = mlg.GenerateMicrolensingEvent()
MicroEvent = MicroEventGenerator.simple_fake()
MicroEventGenerator.save_data(MicroEvent)

'''
# End of SupernovaFitter class
# Begin the commands to run and analyse the chains

fitter = SimpleMicroLensingFitter()
MicroEvent = GenerateMicrolensingEvent()
data = MicroEvent.generate_data()
data = MicroEvent.simple_fake()
fitter.setData(data)  # Loading in the SN Ia data


#plt.scatter(data['time'], data['delta_mag'])
#plt.plot(time, delta_mag)
#plt.show()

cambMCMCManager = CambMCMCManager('SimpleMicroEventFit', fitter, debug=True)
cambMCMCManager.configureMCMC(numCalibrations=10, calibrationLength=1000, thinning=1, maxSteps=20000)
cambMCMCManager.configureSaving(stepsPerSave=1000)

for walk in range(0, 4):  # uncomment this when all your chains are finished and you want to produce a plot and get bounds
    cambMCMCManager.doWalk(walk)  # This command runs the chain in the brackets, they take a while so I do each chain in a different terminal

cambMCMCManager.consolidateData(uid='SimpleMicroEventFit')  # The name of the saved chains
print(cambMCMCManager.getParameterBounds())  # prints the bounds obtained by the mcmc chains
cambMCMCManager.testConvergence(uid='SimpleMicroEventFit')  # doing tests on the chains
cambMCMCManager.plotResults(plotLine=True, filename='SimpleMicroEventFit')  # plot the results and save them as pdf and png

#MicroEvent = GenerateMicrolensingEvent()
#data = MicroEvent.generate_data()
#print data
#delta_mag = MicroEvent.simple_fake()
#time = MicroEvent.generate_times()
# plt.scatter(time, delta_mag)
#plt.plot(time, delta_mag)
#plt.show()
#print MicroEvent.get_t_0()
#print delta_mag
'''