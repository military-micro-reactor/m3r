#! /usr/bin/python3
from mcnptools import Meshtal
meshtally = Meshtal("k_eff/meshtal")

from mpl_toolkits.mplot3d import axes3d

import matplotlib.pyplot as plot
from matplotlib import cm
import numpy as np

tally14 = meshtally.GetTally(104)

xbins = tally14.GetXRBins()
ybins = tally14.GetYZBins()
ebins = tally14.GetEBins()

X, Y = np.meshgrid(xbins, ybins)
flux = np.zeros((len(xbins), len(ybins)))
for i in range(len(xbins)):
    for j in range(len(ybins)):
        flux[i][j] = tally14.GetValue(i, j, 0, 0)


figure = plot.figure()
plot.title("M3R Thermal Flux")
axes = figure.add_subplot(111, projection='3d')
surf = axes.plot_surface(X, Y, flux, cmap=cm.coolwarm)
plot.show()
