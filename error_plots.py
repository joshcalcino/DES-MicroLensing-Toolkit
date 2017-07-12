import numpy as np
import matplotlib.pyplot as plt


def error_plots(magnitude, magerr, t_eff):

    mag_plot = []
    magerr_plot = []
    t_eff_plot = []

    for line in magnitude:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        mag_plot.append(current_line)

    for line in magerr:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        magerr_plot.append(current_line)

    for line in t_eff:
        current_line = line
        current_line.rstrip("\n")
        current_line = float(current_line)
        t_eff_plot.append(current_line)

    plt.scatter(mag_plot, magerr_plot)
    plt.title("Magnitude and Error of selection of stars")
    plt.xlabel("Magnitude")
    plt.ylabel("Error")
    plt.axhline(y = 0.0085, color = 'r')
    plt.show()
