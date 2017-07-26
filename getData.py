import numpy as np
import get_errors
import sys
import pandas as pd
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator

class getData(object):

    def __init__(self, hpix=11737):
        self.list_times, self.uniqueIDs, self.ecat = self.pull_data(hpix)
        self.index = 0
        #self.star_list = self.star_list()

    def get_MJD(self, index =1000):
        quick_id = self.list_times[index]
        MJD_list = self.ecat.query("QUICK_OBJECT_ID== {}".format(quick_id))['MJD_OBS']
        return MJD_list
    
    def get_timesByIDs(self, ID):
        MJD_list = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['MJD_OBS']
        return MJD_list   

    def get_RA(self, IDs):
        ra = self.ecat.query("QUICK_OBJECT_ID== {}".format(IDs))['RA']
        return ra   

    def get_DEC(self, IDs):
        dec = self.ecat.query("QUICK_OBJECT_ID== {}".format(IDs))['DEC']
        return dec   

    def get_t_eff(self, ID):
        t_eff = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['T_EFF']
        return t_eff

    def get_magerr(self, ID):
        magerr  = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['MAGERR_PSF']
        return magerr

    def get_mag(self, ID):
        mag = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['MAG_PSF']
        return mag

    def get_bandpass(self, ID):
        band = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['BAND']
        return band

    def get_spread(self, ID):
        wavg  = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['WAVG_SPREAD_MODEL'] 
        return wavg
       

    def get_spread_err(self, ID):
        spreaderr = self.ecat.query("QUICK_OBJECT_ID== {}".format(ID))['SPREADERR_MODEL']
        return spreaderr
       
    def find20mag(self):
        for x in range(0, len(self.uniqueIDs)):
            Id = self.uniqueIDs[self.index]
            mag = self.get_mag(Id)
            mag = mag[0]
            if (mag < 20.05 and mag > 19.95):
                print("ID is: " + str(Id))
                print("********************")
            self.index += 1

    def grab_details_for_error(self, ID):
        t_eff = self.get_t_eff(ID) 
        magerr = self.get_magerr(ID)
        mag_psf = self.get_mag(ID)
        mjd = self.get_timesByIDs(ID)
        band = self.get_bandpass(ID)
        return t_eff, magerr, mag_psf, mjd, band
    
    def pull_data(self,hpix):
        hpix = int(hpix)
        sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
        from desqcat import load_hpx, load_cat, load_cat_epochs
        cat_wide = load_cat(hpix)
        cat = load_cat(hpix, long=True)
        cat_cols = ['QUICK_OBJECT_ID', 'RA', 'DEC','BAND', 'EXPNUM', 'WAVG_SPREAD_MODEL','SPREADERR_MODEL'] #other options: HPX2048, NEPOCHS, FLAGS, WAVG_FLAGS, EXPNUM, WAVG_MAG_PSF, WAVG_MAGERR_PSF, WAVG_MAG_AUTO, WAVG_MAGERR_AUTO 
        epoch_cols = ['QUICK_OBJECT_ID', 'EXPNUM', 'MJD_OBS', 'BAND', 'T_EFF',
              'MAG_PSF', 'MAGERR_PSF', 'MAG_AUTO', 'MAGERR_AUTO', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL'] #other options: EXPTIME, 
        ecat = load_cat_epochs(hpix, cat_cols, epoch_cols)
        ecat = ecat.query('MAG_PSF < 30')
        obj_expnum_counts = ecat[['QUICK_OBJECT_ID', 'EXPNUM', 'BAND']].groupby(['QUICK_OBJECT_ID', 'EXPNUM'], as_index=False).count()
        obj_expnum_counts.columns = ['QUICK_OBJECT_ID', 'EXPNUM', 'COUNTS']
        duplicated_objects = obj_expnum_counts.QUICK_OBJECT_ID[obj_expnum_counts.COUNTS>1]
        ecat = ecat[np.in1d(ecat.QUICK_OBJECT_ID.values, duplicated_objects.values, invert=True)]
        list_times = ecat['QUICK_OBJECT_ID']
        uniqueIDs = ecat['QUICK_OBJECT_ID'].unique()
        return list_times, uniqueIDs, ecat

    def get_error_details(self, quick_id):
        t_eff, magerr , mag_psf, mjd, bandpass = self.grab_details_for_error(quick_id)
        mag_psf = np.asarray(mag_psf, dtype = float)
        t_eff = np.asarray(t_eff, dtype = float)
        magerr = np.asarray(magerr, dtype = float)
        mjd = np.asarray(mjd, dtype = float)
        bandpass = np.asarray(bandpass, dtype = float)
        testing = get_errors.return_error(mag_psf, t_eff, magerr, mjd, bandpass)
        print("stop")
        #np.savetxt( "newdata.txt", np.array([mag_psf, magerr, t_eff, quick_id]).T, "%5.3f %5.2f %5.2f %d") 
        #print "magerr:", magerr
        #print "Nlist:", N_list
        #N_list = (6.25/((np.square(magerr)*np.log(10)**2)))*5/90*t_eff
       # N_list = 1
        return t_eff, magerr, mag_psf, mjd, bandpass

    def unit_test(self, mjd_list, x):
        u = x.get_u(mjd_list)
        mjd = np.zeros(5)
        expected = np.array([5.22711, 8.27541, 8.27500, 7.98218, 6.88552])
        print expected
        for item in range(0, 5):
            mjd_list[item] = mjd[item]
        print("Expected: ", expected)
        print("Real: ", mjd_list)
        for i in range(0, len(mjd_list)):
            if mjd_list[i] == expected[i]:
                    print("True")
            else:
                return False
        return True 


    def isStar(self, ID): #takes individual object and returns true if is a star, false if galaxy
        test = False
        wavg = self.get_spread(ID)
        spreaderr = self.get_spread_err(ID)
        for n in range(0, len(wavg)):
            if abs(wavg[n]) < (0.003 + spreaderr[n]): #is this for a particular bandpass?
                test = True
        if test == True:
            print "ObjID", ID, "is a star!!"
        else:
            print "ObjID", ID, "is NOT at star."
        return test

    """
    def star_list(self, bandpass = 'g'): #takes ALL ids from data and returns a list of only stars from data
        stars = []
        for ID in self.uniqueIDs:
            #if isStar(ID, bandpass):
                #stars.append(ID)
            isStar_truth = False
            wavg = self.get_spread(ID, bandpass)
            spreaderr = self.get_spread_err(ID, bandpass)
            for n in range(0, len(wavg)):
                if abs(wavg[n])<(0.003 +  spreaderr[n]):
                    isStar_truth = True
            if isStar_truth == True:
                print "found a star!!"
                stars.append(ID)
            else:
                print "no star"
        return stars
    """
