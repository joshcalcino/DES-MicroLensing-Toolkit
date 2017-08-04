import dp
import numpy as np
import sys
import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator
import dp
import numpy as np

class counting_stars(object):

    def __init__(self):
        self.count = 0

        de = []
        temp1 = []
        temp2 = []
        final = []

        # imports files lists with structure cat_hpx_#####.fits
        files = os.listdir("/data/des51.b/data/kadrlica/projects/y3q2/v7/cat/")

        # Removes cat_hpx_ from file name
        for i in range(0,len(files),1):
            temp1.append(files[i].replace("cat_hpx_", ""))

        # Removes .fits from file name
        for j in range(0,len(temp1),1):
            temp2.append(temp1[j].replace(".fits", ""))

        # Removes any remaining files that do not have numbers. 
        for k in range(0,len(temp2),1):
            if temp2[k].isdigit() == True:
                final.append(temp2[k])

        self.file_list = final
        #print("Length of list: ", len(self.file_list))

    def one_republic(self):
        for pixel in range (0, len(self.file_list)):
            count_file = open('counts_and_id.txt', 'a')
            temp_count = 0
            current_pixel = int(self.file_list[pixel])
            data = dp.avg_mag_redux(hpix = current_pixel)
            temp_count = data
            self.count += temp_count
            print("***************")
            print("At pixel ", current_pixel)
            print("Count in pixel: ", temp_count)
            print("Current star count: ", self.count)
            print("***************")
            print(str(current_pixel) + "-" + str(temp_count) + "\n")
            count_file.write(str(current_pixel) + "-" + str(temp_count) + "\n")
            count_file.close()
        return self.count

    def progress(self, id_num):
        for i in range(0, len(self.file_list)):
            if self.file_list[i] == id_num:
                print(i)
                return(i)

    def counting_plot(self):
        input_file = 'counts_and_id.txt'
        lately = np.genfromtxt(input_file, unpack = True)
        id_array2 = []
        count2 = []
        with open(input_file) as f:
            #lines = [line.rstrip('\n') for line in open(input_file)]
            lines = f.readlines()


        id_array = np.zeros(0)
        count = np.zeros(0)
        for x in range(0, len(lines)):
            lines[x] = lines[x].rstrip('\n')
            string = str(lines[x])
            strip = string.split("-")
            id_array = np.append(id_array, int(strip[0]))
            count = np.append(count, int(strip[1]))
        #print(id_array)
        #print(count)
        count = count.astype(int)
        id_array = id_array.astype(int)
        """
        plt.ion()
        plt.clf()
        plt.hist(count, bins = 500)
        plt.xlim(0, 50000)
        plt.ylabel("Counts")
        plt.xlabel("Stars in a pixel")
        plt.show()
        plt.close()
        """
        
        return id_array, count

