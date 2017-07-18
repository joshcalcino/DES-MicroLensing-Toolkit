import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator
import fake_plots
import getData
import star
import getHPIX
import fitsio

class driver(object):

    def __init__(self):
        self.star = star.star()
        self.file_name = ""
        self.hpix = getHPIX.pix() #list of all pixels in survey 
    """
        for i in self.hpix:
            self.file_name = "/home/s1/marika/data/DES-MicroLensing-Toolkit/fitsData/test/lc_curves" + str(i) + ".fits"
            pyfits.writeto(self.file_name, a, clobber=True)
    """
       # print "in driver- teff:", data.get_t_eff(data.uniqueIDs[0])

    def nike(self):
        for pix in range(11737, 11738, 1): #in aelf.hpix:
            data = getData.getData(pix) #160,000 objects with seperate obs
            self.file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/test/lc_curves" + str(pix) + ".fits"
            fits = fitsio.FITS(self.file_name,'rw')
            #empty arrays
            mjd_total = np.array([])
            mag_total = np.array([])
            magerrs = np.array([])
            qid = np.array([])
            ra = np.array([])
            dec = np.array([])
            for i in range(0,10,1): #final_IDs:
                ID = data.uniqueIDs[i]
                if data.isStar(ID):
                
                    #variables from data
                    mjd = data.get_timesByIDs(ID)
                    t_eff = data.get_t_eff(ID)
                    mag_list = data.get_m_0(ID)
                    RA = data.get_RA(ID)
                    DEC = data.get_DEC(ID)
                    #ra.append(RA)
                    #dec.append(DEC) 
                    """
                    Ds = data.get_Ds(objID[i])
                    curve_type = data.get_curve_type(objID[i])
                    """
		    
                    #calculated variables 
                    star_events, final_mag_list = self.star.get_curves(mjd, t_eff, mag_list) #returns 600,000 light curves for the object 
                     
                    self.save_data(pix, star_events, ID,  mjd, RA, DEC, final_mag_list)

                    #call the save to fits file method
                    #self.plot_many(0,100)

	    index =+ 1
        print "index:", index 
        return 0

    
    def save_data(self, pix, star_events, ID,  mjd, RA, DEC, final_mag_list):
        mjd_total = np.array([mjd])
        mag_total = np.array([final_mag_list])
        magerrs = np.array([])
        qid = np.array([ID])
        ra = np.array([RA])
        dec = np.array([DEC])

        array_list = [mjd_total, final_mag_list, magerrs, qid, ra, dec]
        names = ['mjd_total', 'final_mag_list', 'magerrs', 'qid', 'ra', 'dec'] 
        fits.write(array_list, names=names)


         """mjd_total = np.array([])
            loop over stars:
                mjd, mags, magerrs, qid, ra, dec = marikas_code()
                mjd_total = np.array(mjd_total, mjd)
                mag_total = np.array(mag_total, mjd)
            write)fits(mag_total, mag_err_total....)
            mjd, mag, magerr, qid, ra, dec
            mjd, mag, magerr, qid, ra, dec
            mjd, mag, magerr, qid, ra, dec
            mjd, mag, magerr, qid, ra, dec
            50000., 17., 0.1, 10, 0, 0.
            50001., 17.1, 0.1, 10, 0, 0.
            50002., 17.2, 0.1, 10, 0, 0.
            501002., 19.2, 0.3, 1010, 10, 10.
            501002., 19.2, 0.3, 1010, 10, 10.
            501002., 19.2, 0.3, 1010, 10, 10.
            ix, = np.where(qid == 10)
            nrows = nstars* nlightcurves/star * number_of_mags/star
                10000    600,000              len(mjd_list)"""


 
        #data_table = pyfits.open(file_name.fits)

        #for i in final_mag_list: #len(final_mag_list) ~ 600,000
            #data_event = [ID, RA, DEC, mjd, final_mag_list[i]]
            #data_table.append(file_name, data_event)
            
        """
        #Example:
        #$col1 = pyfits.Column(name = 'wave', format = 'D', array = wv_full_test)
        #$col2 = pyfits.Column(name = 'fx', format = 'D', array = fx_full_test)
        #$col3 = pyfits.Column(name = 'fx_error', format = 'D', array = (wv_full_test*.1))

        #$cols = pyfits.ColDefs([col1, col2, col3])
        #$tbhdu = pyfits.new_table(cols)
        #$hdu = pyfits.PrimaryHDU()
        #$thdulist = pyfits.HDUList([hdu, tbhdu])
        #$tbhdu.writeto('filename.fits')
        """

        #get the different columns that we want to return
        #for i in events:
            #RA = RA
            #DEC = DEC
            #MJD_LIST = mjd
            #MAG_LIST = final_mag_list
            #           filename, data, header=None, checksum=False, verify=True, **kwargs
            #fits.append(pix, ID, OBJECT_ID)
            #fits.append(pix, RA, RA)
            #fits.append(pix, DEC, DEC)
            #fits.append(pix, MAG_LIST[j], MAG)
            #fits.append(pix, MJD_LIST[j], MDJ_LIST)

        #save all data per pixel, that way, there are only 1800 total files
        return "saved!!!"
