#! /usr/bin/python3
from mcnptools import Meshtal
meshtally = Meshtal("k_eff/meshtal")

import matplotlib.pyplot as plot
from matplotlib import cm
from matplotlib.ticker import LogFormatter, LogLocator
import numpy as np

_J = 6.242e12
recoverable_energy = 200/_J #recoverable energy per fission in joules
nu = 2.4
reactor_power = 1.0e7       #reactor at 10MWt
I = nu*reactor_power/recoverable_energy
bin_width=500/200.0

print("MCNP flux to actual flux factor is I={}".format(I))

tally14 = meshtally.GetTally(104)

xbins = tally14.GetXRBins()
ybins = tally14.GetYZBins()
ebins = tally14.GetEBins()

X, Y = np.meshgrid(xbins, ybins)
thermal_flux = np.zeros((len(xbins), len(ybins)))
for i in range(len(xbins)):
    for j in range(len(ybins)):
        thermal_flux[i][j] = I*tally14.GetValue(i, j, 0, len(ebins)-1)

figure, axes = plot.subplots()
axes.set_title("M$^3$R Fast Flux Tally at 10MW$_t$")
heatmap = axes.imshow(np.log10(thermal_flux+1), cmap=cm.coolwarm, origin='lower', interpolation='bilinear')
axes.set_xlabel("x-axis bin")
axes.set_ylabel("y-axis bin")
colorbar = figure.colorbar(heatmap, ax=axes)
colorbar.ax.set_ylabel("$log_{10}$ of flux")
plot.show()
