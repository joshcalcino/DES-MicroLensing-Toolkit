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
        self.cat_wide = load_cat(hpix)
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

    #@profile
    def grab_details(self, quick_id, bandpass):
        myobj_df = self.ecat.loc[quick_id]
      #  myobj_r = self.ecat.query("QUICK_OBJECT_ID==" + str(quick_ID) + " & BAND==", bandpass)[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL', 'T_EFF']]
        new_list = myobj_r['MAG_PSF']
        wavg = myobj_r['WAVG_SPREAD_MODEL']
        spreaderr = myobj_r['SPREADERR_MODEL']
        magerr = myobj_r['MAGERR_PSF']
        t_eff = myobj_r['T_EFF']
        #print("wavg: ", wavg)
        #print("spreaderr: ", spreaderr)
        #print "list: ", new_list
        return new_list, wavg, spreaderr, magerr, t_eff

    def grab_detalis_for_error(self, quick_id, bandpass):
        myobj_df = self.ecat.loc[quick_id]
        #  myobj_r = self.ecat.query("QUICK_OBJECT_ID==" + str(quick_ID) + " & BAND==", bandpass)[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'BAND']]
        myobj_r = self.ecat.query("QUICK_OBJECT_ID== {} & BAND=='{}'".format(quick_id, bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL', 'T_EFF']]
        t_eff = myobj_r['T_EFF']
        magerr = myobj_r['MAGERR_PSF']
        return t_eff, magerr

    def avg_mag_redux(self):
        mag_file = open('mag_avgs_redux.txt', 'w')
        mag_list = []
        test_list = self.cat_wide['MAG_PSF_R']
        print(test_list[0])
        spreaderr = self.cat_wide['SPREADERR_MODEL_R']
        wavg = self.cat_wide['WAVG_SPREAD_MODEL_R']
        #print self.cat_wide.head()
        for n in range(0, len(test_list)):
            if abs(wavg[n])<(0.003 +  spreaderr[n]) and (test_list[n] < 30):
                mag_list.append(test_list[n])
                mag_file.write(str(test_list[n])+"\n")
        #print mag_list
        return mag_list

    #@profile
    def avg_mag(self, outfile = "{}_error_data.txt", bandpass = 'r'):
        mag_file = open('delete_me.txt', 'w')
        avg_file = open('mag_averages_short_2.txt', 'w')
        err_file = open('magerr_short_2.txt', 'w')
        eff_file = open('t_eff_short_2.txt', 'w')
        mag_list = []
        outfile = outfile.format(bandpass)

        myobj_rband = self.ecat.query("BAND=='{}'".format(bandpass))[['MJD_OBS','MAG_PSF', 'MAGERR_PSF', 'QUICK_OBJECT_ID', 'BAND', 'WAVG_SPREAD_MODEL', 'SPREADERR_MODEL']]
        self.quick_id_list = myobj_rband['QUICK_OBJECT_ID']     
        #need a list of the times only for the r band!
        quick_id = self.quick_id_list[0]
        current_quick_id = self.quick_id_list[0]
        index = 0
        error_mag = []
        magerr_final = []
        t_eff_final = []
        quick_id_final  = []

        while index in range(0, len(self.quick_id_list)-1):
            print("index at top: ", index)    
            obj_mag, wavg, spreaderr, magerr, t_eff = self.grab_details(quick_id, bandpass)
            acceptable_mag = []
            n = 0
            print(self.quick_id_list)
            while quick_id == current_quick_id:
                if (len(wavg)>0) and (abs(wavg[n])<(0.003 +  spreaderr[n])) and (index < len(self.quick_id_list)-1):
                    print("Found a star")
                    print("obj mag: " + str(obj_mag[n]))
                    acceptable_mag.append(obj_mag[n])
                    magerr_final.append(magerr[n])
                    error_mag.append(obj_mag[n])
                    t_eff_final.append(t_eff[n])
                    avg_file.write(str(obj_mag[n])+"\n")
                    err_file.write(str(magerr[n])+"\n")
                    eff_file.write(str(t_eff[n])+"\n")
                    quick_id_final.append(self.quick_id_list[index+n])
                    current_quick_id = self.quick_id_list[index+n+1]
                    n = n+1
                    print(len(wavg))
                    print(len(spreaderr))
                    print("current quick id: ", current_quick_id)
                    print("quick id: ", quick_id)
                else:
                    if index < len(self.quick_id_list)-1:
                        #print("Not a star")
                        #print("t_eff at n: : " +  str(t_eff[n]))
                        current_quick_id = self.quick_id_list[index+n+1]
                        n = n+1
                        #print("current quick id: ", current_quick_id)
                        #print("quick id: ", quick_id)
            index = index + n
            quick_id = current_quick_id
            #error_mag = acceptable_mag
            
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
            if(index >=2000):
                #print("mag_list", mag_list)
                np.savetxt(outfile, np.array([error_mag, magerr_final, t_eff_final, quick_id_final]).T,
                        "%7.4f %7.4f %7.4f %d")
                return error_mag, magerr_final, t_eff_final
        return mag_list

    def get_t_eff(self, quick_id):
        t_eff, magerr = self.grab_details_for_error(quick_id)
        N_list = ((-2.5)**2/((magerr)**2*np.log(10)))*5/90*t_eff
        return N_list

"""

qid = self.qiuick_id_list
unique_qid = np.unique(qid)
for uqid in unique_qid:
    index = qid==uqid
    index = np.where(qid==uqid)
    mags = obj_mag[index]



        1. loop through grabdetails on each index
              
            
            1a. store all the magnitudes for the index
            1b. average all magnitudes in the index
            1c. store avg magnitude in the array we return.
         


        Last Step. Return star count in each bin.


"""            
