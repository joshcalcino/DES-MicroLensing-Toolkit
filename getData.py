import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

class getData(object):

    def __init__(self, hpix=11737, bandpass = 'g'):
        hpix = int(raw_input("Enter hpix value: "))
        bandpass = raw_input("Enter bandpass: ")
        index = int(raw_input("Enter inded: "))
        print hpix, bandpass, index
        sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
        from desqcat import load_hpx, load_cat, load_cat_epochs
        mpl.rcParams['figure.figsize'] = (8, 5)
        cat_wide = load_cat(hpix)
        cat = load_cat(hpix, long=True)
      #  epochs = {}
       # epochs[bandpass] = load_hpx(hpix, bandpass) #change later to loop through all bands
        cat_cols = ['QUICK_OBJECT_ID', 'RA', 'DEC', 'HPX2048', 'BAND',
            'NEPOCHS', 'FLAGS', 'WAVG_FLAGS', 'EXPNUM',
            'WAVG_MAG_PSF', 'WAVG_MAGERR_PSF',
            'WAVG_MAG_AUTO', 'WAVG_MAGERR_AUTO']
        epoch_cols = ['QUICK_OBJECT_ID', 'EXPNUM', 'MJD_OBS', 'EXPTIME', 'BAND', 'T_EFF',
              'MAG_PSF', 'MAGERR_PSF', 'MAG_AUTO', 'MAGERR_AUTO']
        ecat = load_cat_epochs(hpix, cat_cols, epoch_cols)
        ecat = ecat.query('MAG_PSF < 30')
        self.list_times = ecat['QUICK_OBJECT_ID']
        obj_expnum_counts = ecat[['QUICK_OBJECT_ID', 'EXPNUM', 'BAND']].groupby(['QUICK_OBJECT_ID', 'EXPNUM'], as_index=False).count()
        obj_expnum_counts.columns = ['QUICK_OBJECT_ID', 'EXPNUM', 'COUNTS']
        duplicated_objects = obj_expnum_counts.QUICK_OBJECT_ID[obj_expnum_counts.COUNTS>1]
        self.ecat = ecat[np.in1d(ecat.QUICK_OBJECT_ID.values, duplicated_objects.values, invert=True)]

    def get_MJD(self, index =6, bandpass='g'):
        quick_id = self.list_times[index]
        myobj_df = self.ecat.loc[quick_id]
      #  myobj_r = self.ecat.query("QUICK_OBJECT_ID==" + str(quick_ID) + " & BAND==", bandpass)[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        MJD_list = myobj_r['MJD_OBS']
      #  print "MJD_list: ", MJD_list
        return MJD_list

