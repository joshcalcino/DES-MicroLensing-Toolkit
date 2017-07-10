import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator

class data_practice(object):

    def __init__(self, hpix=11737):
        print("starting")
        sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
        from desqcat import load_hpx, load_cat, load_cat_epochs
        mpl.rcParams['figure.figsize'] = (8, 5)
        cat_wide = load_cat(hpix)
        cat = load_cat(hpix, long=True)
        cat_cols = ['QUICK_OBJECT_ID', 'RA', 'DEC', 'HPX2048', 'BAND',
            'NEPOCHS', 'FLAGS', 'WAVG_FLAGS', 'EXPNUM',
            'WAVG_MAG_PSF', 'WAVG_MAGERR_PSF',
            'WAVG_MAG_AUTO', 'WAVG_MAGERR_AUTO', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL']
        #epoch_cols = ['QUICK_OBJECT_ID', 'MAG_PSF', 'MAG_PSF_R']
        epoch_cols = ['QUICK_OBJECT_ID', 'EXPNUM', 'MJD_OBS', 'EXPTIME', 'BAND', 'T_EFF',
              'MAG_PSF', 'MAGERR_PSF', 'MAG_AUTO', 'MAGERR_AUTO', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL']
        ecat = load_cat_epochs(hpix, cat_cols, epoch_cols)
        #ecat = ecat.query('QUICK_OBJECT_ID == 11173700000000')
        ecat = ecat.query('MAG_PSF < 30')
        self.list_times = ecat['QUICK_OBJECT_ID']
        #self.list_wavg = cat['WAVG_SPREAD_MODEL_R']
        #self.list_wavgerr = cat['SPREADERR_MODEL_R']
        obj_expnum_counts = ecat[['QUICK_OBJECT_ID', 'EXPNUM', 'BAND']].groupby(['QUICK_OBJECT_ID', 'EXPNUM'], as_index=False).count()
        obj_expnum_counts.columns = ['QUICK_OBJECT_ID', 'EXPNUM', 'COUNTS']
        duplicated_objects = obj_expnum_counts.QUICK_OBJECT_ID[obj_expnum_counts.COUNTS>1]
        self.ecat = ecat[np.in1d(ecat.QUICK_OBJECT_ID.values, duplicated_objects.values, invert=True)]

    def grab_details(self, quick_id, bandpass='r'):
        myobj_df = self.ecat.loc[quick_id]
      #  myobj_r = self.ecat.query("QUICK_OBJECT_ID==" + str(quick_ID) + " & BAND==", bandpass)[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL']]
        new_list = myobj_r['MAG_PSF']
        wavg = myobj_r['WAVG_SPREAD_MODEL']
        spreaderr = myobj_r['SPREADERR_MODEL']
        #print("wavg: ", wavg)
        #print("spreaderr: ", spreaderr)
        #print "list: ", new_list
        return new_list, wavg, spreaderr


    def avg_mag(self):

        mag_file = open('mag_averages.txt', 'w')
        mag_list = []

        myobj_rband = self.ecat.query("BAND=='{}'".format('r'))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL']]
        self.quick_id_list = myobj_rband['QUICK_OBJECT_ID']     
        #need a list of the times only for the r band!
        quick_id = self.quick_id_list[0]
        current_quick_id = self.quick_id_list[0]
        index = 0
        while index in range(0, len(self.quick_id_list)-1):
            print("index at top: ", index)    
            obj_mag, wavg, spreaderr = self.grab_details(quick_id, bandpass = 'r')
            acceptable_mag = []
            n = 0
            #print("quick id: ", self.quick_id_list[0]) 
            #print("quick id: ", self.quick_id_list[1]) 
            #print("quick id: ", self.quick_id_list[2]) 
            while quick_id == current_quick_id:
                if abs(wavg[n])<(0.003 +  spreaderr[n]) and (index < len(self.quick_id_list)-1):
                    print("Found a star")
                    acceptable_mag.append(obj_mag[n])
                    current_quick_id = self.quick_id_list[index+n+1]
                    n = n+1
                    print("current quick id: ", current_quick_id)
                    print("quick id: ", quick_id)
                else:
                    if index < len(self.quick_id_list)-1:
                        print("Not a star")
                        current_quick_id = self.quick_id_list[index+n+1]
                        n = n+1
                        print("current quick id: ", current_quick_id)
                        print("quick id: ", quick_id)
            index = index + n
            quick_id = current_quick_id
            
            
            total = 0
            if (len(acceptable_mag) > 0):
                for i in range(0, len(acceptable_mag)):
                    total = total + acceptable_mag[i]
                total = total/(len(acceptable_mag))
                print("Index: ", index)
                print("Total: ", total)
                print("**********")
                mag_list.append(total)
                mag_file.write(str(total)+"\n")
            #print("index: ", index)        
            acceptable_mag = []
            #if(index >=2000):
                #print("mag_list", mag_list)
                #return mag_list
        return mag_list





"""

        1. loop through grabdetails on each index
              
            
            1a. store all the magnitudes for the index
            1b. average all magnitudes in the index
            1c. store avg magnitude in the array we return.
         


        Last Step. Return star count in each bin.


"""            
