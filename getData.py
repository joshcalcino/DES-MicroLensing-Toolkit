import numpy as np
import get_errors
import sys
import pandas as pd
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator

class getData(object):

    def __init__(self, hpix=9280):
        print "Pixel:", hpix
        self.sIDs, self.star_cat, self.ecat = self.pull_data(hpix)
        self.index = 0
        self.get_starData(self.sIDs[0])

    def pull_data(self,hpix):
        hpix = int(hpix)
        sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
        from desqcat import load_hpx, load_cat, load_cat_epochs
        cat_wide = load_cat(hpix)
        cat_cols = ['QUICK_OBJECT_ID', 'RA', 'DEC','BAND', 'EXPNUM', 'WAVG_MAG_PSF', 'WAVG_MAGERR_PSF', 'WAVG_SPREAD_MODEL','SPREADERR_MODEL'] 
        epoch_cols = ['QUICK_OBJECT_ID','EXPTIME', 'EXPNUM', 'MJD_OBS', 'BAND', 'T_EFF','MAG_PSF', 'MAGERR_PSF']
        ecat = load_cat_epochs(hpix, cat_cols, epoch_cols)

        test_cat = load_cat(hpix, long=False, cols=['QUICK_OBJECT_ID', 'WAVG_MAG_PSF_R','WAVG_SPREAD_MODEL_R','WAVG_SPREADERR_MODEL_R'])
        
        good_qids = test_cat.query('WAVG_MAG_PSF_R < 21.5 & (WAVG_SPREAD_MODEL_R<(0.003+WAVG_SPREADERR_MODEL_R))').QUICK_OBJECT_ID
        star_cat = ecat[ecat.QUICK_OBJECT_ID.isin(good_qids)]

        obj_expnum_counts = star_cat[['QUICK_OBJECT_ID', 'EXPNUM', 'BAND']].groupby(['QUICK_OBJECT_ID', 'EXPNUM'], as_index=False).count()
        obj_expnum_counts.columns = ['QUICK_OBJECT_ID', 'EXPNUM', 'COUNTS']
        duplicated_objects = obj_expnum_counts.QUICK_OBJECT_ID[obj_expnum_counts.COUNTS>1]
        star_cat = star_cat[np.in1d(star_cat.QUICK_OBJECT_ID.values, duplicated_objects.values, invert = True)]
        sIDs = star_cat['QUICK_OBJECT_ID'].unique()
        return sIDs, star_cat, ecat

 def loadStar(self, ID):
        data = self.star_cat.query("QUICK_OBJECT_ID== {}".format(ID))[['MJD_OBS', 'RA','DEC','T_EFF','MAGERR_PSF','MAG_PSF', 'BAND', "WAVG_SPREAD_MODEL", "SPREADERR_MODEL"]]
        self.edata = data
        self.edata_id = ID
        #print "Star ID", ID, "data loaded." 

    def getMJD(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["MJD_OBS"]

    def get_RA(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["RA"]

    def get_DEC(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
       return self.edata["DEC"]

    def get_t_eff(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["T_EFF"]

    def get_magerr(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["MAGERR_PSF"]

    def get_mag(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["MAG_PSF"]

    def get_band(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["BAND"]

    def get_spread(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["WAVG_SPREAD_MODEL"]

    def get_spread_err(self, ID):
        if self.edata_id != ID:
            self.loadStar(ID)
        return self.edata["SPREADERR_MODEL"]

       
    def find20mag(self):
        maglist = []
        sIDs = self.sIDs
        print "sIDs:", len(sIDs)
        for ID in sIDs:
            mag = self.get_mag(ID)
            for i in mag:
                if (i < 20.05 and i > 19.95):
                    maglist.append(ID)
                    break
        print "maglist len:", len(maglist)
        return maglist
 
    def oldFind20Mag(self):
        for x in range(0, len(self.uniqueIDs)):
            Id = self.uniqueIDs[self.index]
            mag = self.get_mag(Id)
            mag = mag[0]
            if (mag < 20.05 and mag > 19.95):
                print("ID is: " + str(Id))
                print("********************")
            self.index += 1

    def get_wavgs(self, ID):
        wg = self.cat_wide.query("QUICK_OBJECT_ID== {}".format(ID))['WAVG_MAG_PSF_G']
        wr = self.cat_wide.query("QUICK_OBJECT_ID== {}".format(ID))['WAVG_MAG_PSF_R']
        wi = self.cat_wide.query("QUICK_OBJECT_ID== {}".format(ID))['WAVG_MAG_PSF_I']
        wz = self.cat_wide.query("QUICK_OBJECT_ID== {}".format(ID))['WAVG_MAG_PSF_Z']
        wy = self.cat_wide.query("QUICK_OBJECT_ID== {}".format(ID))['WAVG_MAG_PSF_Y']
        return wg, wr, wi, wz, wy
        
    def grab_details_for_error(self, ID):
        t_eff = self.get_t_eff(ID) 
        magerr = self.get_magerr(ID)
        mag_psf = self.get_mag(ID)
        mjd = self.get_timesByIDs(ID)
        band = self.get_bandpass(ID)
        return t_eff, magerr, mag_psf, mjd, band
   
    def count_stars(self):
        return self.sIDs.size


    def get_error_details(self, quick_id):
        t_eff, magerr , mag_psf, mjd, bandpass = self.grab_details_for_error(quick_id)
        mag_psf = np.asarray(mag_psf, dtype = float)
        t_eff = np.asarray(t_eff, dtype = float)
        magerr = np.asarray(magerr, dtype = float)
        mjd = np.asarray(mjd, dtype = float)
        bandpass = np.asarray(bandpass, dtype = float)
        testing = get_errors.return_error(mag_psf, t_eff, magerr, mjd, bandpass)
        print("stop")
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
