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
from astropy.io import fits
import pyfits

class driver(object):

    def __init__(self):
        self.star = star.star()
        #self.clearDir()
        self.file_name = ""
        #a = np.array([1,2,3,4,5])
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
            hdulist = fits.HDUList() 
            #final_IDs = data.star_list() #list of objects that are stars
            #print "individual stars in pixel:", len(final_IDs)
            for i in range(0,10,1): #final_IDs:
                ID = data.uniqueIDs[i]
                if data.isStar(ID):
                
                    #variables from data
                    mjd = data.get_timesByIDs(ID)
                    t_eff = data.get_t_eff(ID)
                    mag_list = data.get_m_0(ID)
                    RA = data.get_RA(ID)
                    DEC = data.get_DEC(ID) 
                    """
                    Ds = data.get_Ds(objID[i])
                    curve_type = data.get_curve_type(objID[i])
                    """
		    
                    #calculated variables 
                    star_events, final_mag_list = self.star.get_curves(mjd, t_eff, mag_list) #returns 600,000 light curves for the object  
                    
                    #self.save_data(pix, star_events, ID,  mjd, RA, DEC, final_mag_list)

                    #call the save to fits file method
                    #self.plot_many(0,100)

                    """
                    each row of the table is a lightcurve
                    
                    one column specify what pixel each light curves are in, each row = the pixel number
                    one column specify what object ID the light curves are for, each row = the object ID number
                    one column RA
                    one column DEC
                    one column with a lightcurve ID, each row = new lightcurve ID
                    one column where each row is the mjd list (numpy array?) for a lightcurve
                    one column where each row is the final_mag_list for a lightcurve
                    """
			mjd_total = np.array([])
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
				10000    600,000              len(mjd_list)
		    a = np.array([1,2,3,4,5]) 
                    b = np.array([3.2,5.6,7.8])
                    c = np.array([2,2,2,2,2,2,2,2,22,33])
                    d = np.array([10,2.3])
		    e = np.array([1000,20000,440000,5678,223345])
		    f = np.array([1,2,3])
		    g = np.array([4,5,6])
		
		    c1 = fits.Column(name = "pix", format = "D", array = a)
                    c2 = fits.Column(name = "star_events", format = "D", array = b)
                    c3 = fits.Column(name = "ID", format = "D", array =  c)
		    c4 = fits.Column(name = "mjd", format = "D", array = d)
		    c5 = fits.Column(name = "RA", format = "D", array = e)
                    c6 = fits.Column(name = "DEC", format = "D", array = f)
		    c7 = fits.Column(name = "final_mag_list", format = "D", array = g)
		    cols = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7])
		    table = fits.BinTableHDU.from_columns(cols)  
                    hdulist.append(table)
            self.file_name = "/home/s1/mmironov/DES-MicroLensing-Toolkit/fitsData/test/lc_curves" + str(pix) + ".fits"
	    hdulist.writeto(self.file_name, clobber=True)	
	    index =+ 1
        print "index:", index 
        return 0

    
								


    def from_columns(cls, columns, header=None, nrows=0, fill=False, 
		     **kwargs):
	coldefs = cls._columns_type(columns)
        data = FITS_rec.from_columns(coldefs, nrows=nrows, fill=fill)
        hdu = cls(data=data, header=header, **kwargs)
        coldefs._add_listener(hdu)
        return hdu
 
    def save_data(self, pix, events, ID, mjd, RA, DEC, final_mag_list):
        print "save_data" 
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
