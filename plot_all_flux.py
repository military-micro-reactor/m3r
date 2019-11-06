#! /usr/bin/python3
from mcnptools import Meshtal
meshtally = Meshtal("meshtal")

from mpl_toolkits.mplot3d import axes3d

import matplotlib.pyplot as plot
from matplotlib import cm
import numpy as np

tally104 = meshtally.GetTally(104)

xbins = tally104.GetXRBins()
ybins = tally104.GetYZBins()
ebins = tally104.GetEBins()
fluxes = np.zeros([len(ebins), len(xbins), len(ybins)])

X, Y = np.meshgrid(xbins, ybins)
for i in range(len(xbins)):
    for j in range(len(ybins)):
        for e in range(len(ebins)):
            fluxes[e][i][j] = tally104.GetValue(i, j, 0, e)


figure = plot.figure()
plot.title("M3R Flux")
axes = figure.add_subplot(2, 1, 1, projection='3d')
axes.plot_surface(X, Y, fluxes[0], cmap=cm.coolwarm)

axes = figure.add_subplot(2, 1, 2, projection='3d')
axes.plot_surface(X, Y, fluxes[4], cmap=cm.coolwarm)
plot.show()
