import numpy as np
import matplotlib.pyplot as plt

#appease the mighty Jim Annis and find out how many stars have mag
#above and below 20 in the r band.

#@profile
def cutter(mag_list_file):
    underGrads = []
    grads = []
    hist = []
    
    for line in mag_list_file:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        if current_line <= 20:
            underGrads.append(current_line)
        if current_line > 20:
            grads.append(current_line)
        if current_line <= 21.5:
            hist.append(current_line)
    u20 = len(underGrads)
    o20 = len(grads)

    #plt.hist(hist, bins = range(0, 25 + 1, 1), alpha = 0.5, histtype = 'bar', ec='black')
    plt.hist(hist, bins = 'auto', alpha = 0.5, histtype = 'bar', ec='black')
    plt.title("Number of Counts of each Magnitude")
    plt.yscale('log', nonposy = 'clip')
    plt.xlabel("Magnitude")
    plt.ylabel("Counts")
    #plt.savefig("r_band_all_data.png")
    plt.show()
    return u20, o20
