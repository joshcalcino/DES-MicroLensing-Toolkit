import numpy as np
import counting_stars as stars
#import decam2hp
import numpy as np
import healpy as hp
import healpy as hp
import matplotlib.pyplot as plt;

def Aug3(gs=100):
    nsides=32
    ra = 0
    dec = 0
    phi_test = ra*2*np.pi/360.
    theta_test = (90-dec)*2*np.pi/360
    whatever = hp.ang2pix(nsides, theta_test, phi_test)
    #print(whatever)
    data = stars.counting_stars()
    pixel, counts = data.counting_plot()
    theta, phi = hp.pix2ang(nsides, pixel)
    #theta, phi = hp.pix2ang(nsides, whatever)
    #print(str(theta) + " " + str(theta_test))
    #print(str(phi) + " " + str(phi_test))
    s_ra = phi*360./(2*np.pi)
    #print("theta: ", theta)
    s_dec = 90 - (theta*360./(2.*np.pi))
    #print(s_ra)
    #print(s_dec)
    plt.clf()
    ix = s_ra > 180; s_ra[ix]=s_ra[ix]-360
    ix, = np.where((s_ra>-100)&(s_ra<100))
    plotmap(s_ra[ix], s_dec[ix], counts[ix],gs=gs)
    return
    plt.hexbin(s_ra, s_dec, counts, bins = "log")
    plt.colorbar()
    plt.ylabel("Dec")
    plt.xlabel("RA")
    #plt.savefig("counts_on_sky.png")

def plotmap(ra,dec,vals,vmin=0,vmax=500000,save=False,file="",gs=100) :
    from equalArea import mcplot
    from equalArea import mcbryde
    plt.clf();
    plt.axes().set_aspect('equal')
    x,y=mcbryde.mcbryde(ra,dec)
    #Add number of stars to the plt eventually
    #plt.hexbin(x,y,vals,gridsize=gs)
    mcplot.plot(x,y,vals,vmin=vmin,vmax=vmax,cmap="jet",save=save,file=file,gridsize=gs, sum = True)
    plt.ylabel("Dec")
    plt.xlabel("RA")
