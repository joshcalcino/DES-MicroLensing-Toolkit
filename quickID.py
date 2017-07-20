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
import os

class collect():
        
    def __init__(self):
        self.hpix = getHPIX.pix() #list of all pixels in survey

    def collectData():
        index = 0
        for pix in self.hpix:
            if pix != "11737": continue
            data = getData.getData(pix) #160,000 objects with seperate obs
            objID = data.uniqueIDs #list of all objects
            for i in range(4,30,1): #for i in data:

