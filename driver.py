import numpy as np
import sys
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import MicroLensingGenerator
import fake_plots
import getData
import parameters


class driver(object):


    def __init__(self, events):
        self.event_list = events

   # @profile
    def plot_many(self, start, stop, step=1):
        fake_plots.clear()
        index = 0
        while start < stop:
            fake_plots.plot_many(self.event_list[start])
            start += step
            index += 1
        return index
