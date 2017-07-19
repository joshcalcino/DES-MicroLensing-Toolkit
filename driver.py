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
        self.hpix = getHPIX.pix() #list of all pixels in survey 
        """
        for i in self.hpix:
            self.file_name = "/home/s1/marika/data/DES-MicroLensing-Toolkit/fitsData/test/lc_curves" + str(i) + ".fits"
            pyfits.writeto(self.file_name, a, clobber=True)
        """

       # print "in driver- teff:", data.get_t_eff(data.uniqueIDs[0])
        """
        #Example:
        #$col1 = pyfits.Column(name = 'ID', format = 'D', array = wv_full_test)
        #$col2 = pyfits.Column(name = 'RA', format = 'D', array = fx_full_test)
        #$col3 = pyfits.Column(name = 'DEC', format = 'D', array = (wv_full_test*.1))
        #$col4 = pyfits.Column(name = 'MAG_LIST', format = 'D', array = wv_full_test)
        #$col5 = pyfits.Column(name = 'MJD_LIST', format = 'D', array = fx_full_test)
        """

    def nike(self):
        index = 0
        for pix in self.hpix:
            if pix != "11737": continue
            data = getData.getData(pix) #160,000 objects with seperate obs
            objID = data.uniqueIDs #list of all objects
            for i in range(4,5,1): #for i in data:
                #variables from data
                mjd = data.get_timesByIDs(objID[i])
                t_eff = data.get_t_eff(objID[i])
                m_0 = data.get_m_0(objID[i]) 
                is_star = data.isStar(objID[i])
                RA = data.get_RA(objID[i])
                DEC = data.get_DEC(objID[i])
                """
                Ds = data.get_Ds(objID[i])
                curve_type = data.get_curve_type(objID[i])
                """

                item = star.star()

                #calculated variables
                if is_star == True:
                    star_event = item.get_curves(mjd, t_eff, m_0) #returns about 36000 light curves 
                    self.save_data(pix, star_event, objID[i], mjd, RA, DEC)
                    self.plots( star_event[0] )
                #call the save to fits file method
                #self.plot_many(0,100)
            index += 1
        print "index:", index 
        return 0

    def nike2(self):
        index = 0
        mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list  = \
            np.array([]), np.array([]), np.array([]), np.array([]), np.array([]), np.array([])
        for pix in self.hpix:
            if pix != "11737": continue
            data = getData.getData(pix) #160,000 objects with seperate obs
            objID = data.uniqueIDs #list of all objects
            for i in range(4,30,1): #for i in data:
                #variables from data
                mjd = data.get_timesByIDs(objID[i])
                t_eff = data.get_t_eff(objID[i])
                m_0 = data.get_m_0(objID[i]) 
                is_star = data.isStar(objID[i])
                RA = data.get_RA(objID[i])
                DEC = data.get_DEC(objID[i])
                if is_star:
                    mjd_list = np.vstack([mjd_list, mjd])
                    teff_list = np.vstack([teff_list, t_eff])
                    m0_list = np.vstack([m0_list, m_0])
                    ra_list = np.vstack([ra_list, RA])
                    dec_list = np.vstack([dec_list, DEC])
                    objID_list = np.vstack([objID_list, objID[i]])
                """
                Ds = data.get_Ds(objID[i])
                curve_type = data.get_curve_type(objID[i])
                """
        return mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list

    def plot_curves(self, mjd_list, teff_list, m0_list, ra_list, dec_list, objID_list, index):
        size =  len(ra_list)
        for i in range(size): #for i in data:
            if i != index: continue
            item = star.star()
            #calculated variables
            mjd = mjd_list[i]
            t_eff= teff_list[i]
            m_0 = m0_list[i]
            ra = ra_list[i]
            dec = dec_list[i]
            objID = objID_list[i]
            print mjd_list
            print type(mjd)
            print np.shape(mjd), np.shape(t_eff), np.shape(m_0),np.shape(ra), np.shape(dec)
            star_event = item.get_curves(mjd, t_eff, m_0) #returns about 36000 light curves 
            print np.shape(star_event)
            self.plots(star_event, 100, 110, 1)
        return 0
    
    def clearDir():
        return 0

    def plot1(self,event):
        fake_plots.clear()
        fake_plots.plot_many(event)    
        return 0

    def plots(self,event, start =0, stop = 2, step=1):
        fake_plots.clear()
        for i in range(start, stop, step):
            fake_plots.plot_many(event[i])    
        return 0

    def save_data(self, pix, events, ID, mjd, RA, DEC):
        print "save_data" 
        
        #file_name = ""
        #file_name = "/home/s1/marika/data/DES-MicroLensing-Toolkit/fitsData/test/lc_curves" + pix + ".fits"
        #data_table = pyfits.open(file_name)

        #for i in range(0, len(final_mag_list),1): #len(final_mag_list) ~ 600,000
            #data_event = [ID, RA, DEC, mjd[i], final_mag_list[i]]

        
        #data_table.append(file_name, data_event)
            
        """
        #$cols = pyfits.ColDefs([col1, col2, col3, col4, col5])
        #$tbhdu = pyfits.new_table(cols)
        #$hdu = pyfits.PrimaryHDU()
        #$thdulist + 'ID'  = pyfits.HDUList([hdu, tbhdu])
        #$tbhdu.writeto('file_name.fits')
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
        return "saved!!!"
