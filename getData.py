import numpy as np
import get_error
import sys
import pandas as pd
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator

class getData(object):

    def __init__(self, hpix=11737):
        sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
        from desqcat import load_hpx, load_cat, load_cat_epochs
      #  mpl.rcParams['figure.figsize'] = (8, 5)
        cat_wide = load_cat(hpix)
        cat = load_cat(hpix, long=True)
        cat_cols = ['QUICK_OBJECT_ID', 'BAND', 'EXPNUM']
        epoch_cols = ['QUICK_OBJECT_ID', 'EXPNUM', 'MJD_OBS', 'BAND', 'T_EFF',
              'MAG_PSF', 'MAGERR_PSF', 'MAG_AUTO', 'MAGERR_AUTO', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL']
       # cat_cols = ['QUICK_OBJECT_ID', 'RA', 'DEC', 'HPX2048', 'BAND',
       #     'NEPOCHS', 'FLAGS', 'WAVG_FLAGS', 'EXPNUM',
       #     'WAVG_MAG_PSF', 'WAVG_MAGERR_PSF',
       #     'WAVG_MAG_AUTO', 'WAVG_MAGERR_AUTO']
       # epoch_cols = ['QUICK_OBJECT_ID', 'EXPNUM', 'MJD_OBS', 'EXPTIME', 'BAND', 'T_EFF',
       #       'MAG_PSF', 'MAGERR_PSF', 'MAG_AUTO', 'MAGERR_AUTO']
        ecat = load_cat_epochs(hpix, cat_cols, epoch_cols)
        ecat = ecat.query('MAG_PSF < 30')
        self.list_times = ecat['QUICK_OBJECT_ID']
        self.uniqueIDs = self.list_times.unique()
        obj_expnum_counts = ecat[['QUICK_OBJECT_ID', 'EXPNUM', 'BAND']].groupby(['QUICK_OBJECT_ID', 'EXPNUM'], as_index=False).count()
        obj_expnum_counts.columns = ['QUICK_OBJECT_ID', 'EXPNUM', 'COUNTS']
        duplicated_objects = obj_expnum_counts.QUICK_OBJECT_ID[obj_expnum_counts.COUNTS>1]
        self.ecat = ecat[np.in1d(ecat.QUICK_OBJECT_ID.values, duplicated_objects.values, invert=True)]

    def isStar(self, ID):
        #for i in self.uniqueIDs:
            #if star == True:
            #    stars.append(ID)
            #else:
            #    nothing
        return stars

    def get_MJD(self, index =1000, bandpass='g'):
        quick_id = self.list_times[index]
        myobj_df = self.ecat.loc[quick_id]
      #  myobj_r = self.ecat.query("QUICK_OBJECT_ID==" + str(quick_ID) + " & BAND==", bandpass)[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        MJD_list = myobj_r['MJD_OBS']
       # print "MJD_list: ", MJD_list
        return MJD_list
    
    def get_timesByIDs(self, IDs, bandpass='g'):
        myobj_df = self.ecat.loc[IDs]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(IDs, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        MJD_list = myobj_r['MJD_OBS']
        return MJD_list   

    def get_RA(self, IDs, bandpass='g'):
        myobj_df = self.ecat.loc[IDs]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(IDs, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        MJD_list = myobj_r['MJD_OBS']
        return MJD_list   

    def get_DEC(self, IDs, bandpass='g'):
        myobj_df = self.ecat.loc[IDs]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(IDs, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        MJD_list = myobj_r['MJD_OBS']
        return MJD_list   

    def grab_details_for_error(self, quick_id, bandpass = 'r'):
        myobj_df = self.ecat.loc[quick_id]
        #  myobj_r = self.ecat.query("QUICK_OBJECT_ID==" + str(quick_ID) + " & BAND==", bandpass)[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'T_EFF']]
        t_eff = myobj_r['T_EFF']
        magerr = myobj_r['MAGERR_PSF']
        mag_psf = myobj_r['MAG_PSF']
        mjd = myobj_r['MJD_OBS']
        return t_eff, magerr, mag_psf, mjd

    def get_m_0(self, ID):
        maths = 0
        myobj_df = self.ecat.loc[IDs]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(IDs, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        m_0s = myobj_r['MAG_PSF']
        for i in m_0s:
            maths += i
        m_0 = float(maths/len(m_0s))
        return m_0   
        

    def get_error_details(self, quick_id):
        t_eff, magerr , mag_psf, mjd= self.grab_details_for_error(quick_id)
        print "teff:", t_eff
        error = get_error.get_error(mag_psf, t_eff, magerr, mjd)
        #print "magerr:", magerr
        #print "Nlist:", N_list
        #N_list = (6.25/((np.square(magerr)*np.log(10)**2)))*5/90*t_eff
       # N_list = 1
        return t_eff, magerr, mag_psf, mjd

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

    def isStar(self, quick_id_list, mjd):

        for quick_id in range(0, len(mjd)):
            stars = []
            myobj_df = self.ecat.loc[quick_id_list[quick_id]]
            myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL', 'T_EFF']]
            wavg = myobj_r['WAVG_SPREAD_MODEL']
            spreaderr = myobj_r['SPREADERR_MODEL']

            isStar_truth = False

            for n in range(0, len(mjd)):
                if abs(wavg[n])<(0.003 +  spreaderr[n]):
                    isStar_truth = True
                else:
                    isStar_truth = False
                    break

            if isStar_array == True:
                stars.append(quick_id_list[quick_id])
        
        return stars
