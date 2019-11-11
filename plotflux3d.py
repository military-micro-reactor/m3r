#! /usr/bin/python3
from mcnptools import Meshtal
meshtally = Meshtal("k_eff/meshtal")

from mpl_toolkits.mplot3d import axes3d

import matplotlib.pyplot as plot
from matplotlib import cm
import numpy as np

_J = 6.242e12
recoverable_energy = 200/_J  # recoverable energy per fission in joules
nu = 2.4
reactor_power = 1.0e7       # reactor at 10MWt
I = nu*reactor_power/recoverable_energy
bin_width = 500/200.0

print("MCNP flux to actual flux factor is I={}".format(I))

tally14 = meshtally.GetTally(104)

xbins = tally14.GetXRBins()
ybins = tally14.GetYZBins()
ebins = tally14.GetEBins()

X, Y = np.meshgrid(xbins, ybins)
flux = np.zeros((len(xbins), len(ybins)))
for i in range(len(xbins)):
    for j in range(len(ybins)):
        flux[i][j] = tally14.GetValue(i, j, 0, 0)


figure = plot.figure(figsize=plot.figaspect(0.5))
figure.suptitle("M$^3$R Thermal Flux at 10MW$_t$")
axes = figure.add_subplot(1, 2, 2, projection='3d')
surface = axes.plot_surface(X, Y, flux, cmap=cm.coolwarm)

axes = figure.add_subplot(1, 2, 1)
heatmap = axes.imshow(np.log10(I*flux+1), cmap=cm.coolwarm,
                      origin='lower', interpolation='bilinear')
axes.set_xlabel("x-axis bin")
axes.set_ylabel("y-axis bin")
colorbar = figure.colorbar(heatmap, ax=axes)
colorbar.ax.set_ylabel("$log_{10}$ of flux")

plot.show()
