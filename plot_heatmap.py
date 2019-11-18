#! /usr/bin/python3
from mcnptools import Meshtal
meshtally = Meshtal("k_eff/meshtal")

import matplotlib.pyplot as plot
from matplotlib import cm
from matplotlib.ticker import LogFormatter, LogLocator
import numpy as np

_J = 6.242e12
recoverable_energy = 200/_J  # recoverable energy per fission in joules
nu = 2.44
reactor_power = 10e6       # reactor at 200kWt
I = nu*reactor_power/recoverable_energy
bin_width = 500/200.0

print("MCNP flux to actual flux factor is I={}".format(I))

tally14 = meshtally.GetTally(104)

xbins = tally14.GetXRBins()
ybins = tally14.GetYZBins()
ebins = tally14.GetEBins()

X, Y = np.meshgrid(xbins, ybins)
thermal_flux = np.zeros((len(xbins)-101, len(ybins)-200))
for i in range(100,len(xbins)-101):
    for j in range(100, len(ybins)-100):
        thermal_flux[i][j] = I*tally14.GetValue(i, j, 0, 0)

figure, axes = plot.subplots()
axes.set_title("Flux Tally at 10MW$_t$")
heatmap = axes.imshow(np.log10(thermal_flux+1), cmap=cm.coolwarm,
                      origin='lower', interpolation='bilinear')
axes.set_xlabel("x-axis bin")
axes.set_ylabel("y-axis bin")
colorbar = figure.colorbar(heatmap, ax=axes)
colorbar.ax.set_ylabel("$log_{10}$ of flux")
plot.show()
