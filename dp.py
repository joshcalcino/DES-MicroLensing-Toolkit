import numpy as np
import sys
import MicroLensingGenerator
sys.path.append('/data/des51.b/data/neilsen/wide_cadence/python')
from desqcat import load_hpx, load_cat, load_cat_epochs

def avg_mag_redux(hpix):

    cat_wide = load_cat(hpix)
    mag_list = []
    test_list =cat_wide['MAG_PSF_R']
    #print(test_list[0])
    spreaderr = cat_wide['SPREADERR_MODEL_R']
    wavg = cat_wide['WAVG_SPREAD_MODEL_R']
    #print self.cat_wide.head()
    count = 0
    for n in range(0, len(test_list)):
        if abs(wavg[n])<(0.003 +  spreaderr[n]) and (test_list[n] <= 21.5):
            #mag_list.append(test_list[n])
            count += 1
    #print mag_list
    #counts_file.write(str(count) + "\n")
    return count

