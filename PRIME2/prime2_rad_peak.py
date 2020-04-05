#! /usr/bin/python3
from mcnptools import Mctal
tallyfile = Mctal("mctam")

import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
from statistics import mean
from math import sin, cos, pi

i_hat = [0, 6.75]
j_hat = [-5.8457, 3.375]

prime2_basis = np.column_stack((i_hat, j_hat))

core_ij = []

imin = 0
imax = 4
for j in range(-4, 1):
    row_ij = np.concatenate(
        (np.arange(imin, imax + 1),
         j * np.ones(imax - imin + 1))).reshape(2,
                                                imax - imin + 1).transpose()
    core_ij.extend(row_ij)
    imin -= 1

imax -= 1
imin = -4
for j in range(1, 5):
    row_ij = np.concatenate(
        (np.arange(imin, imax + 1),
         j * np.ones(imax - imin + 1))).reshape(2,
                                                imax - imin + 1).transpose()
    core_ij.extend(row_ij)
    imax -= 1

    
coords = np.dot(prime2_basis, np.array(core_ij).transpose())
tally = tallyfile.GetTally(14)
fbins = tally.GetFBins()
tdata = [tally.GetValue(f, 0, 0, 0, 0, 0, 0, 0) for f in range(len(fbins))]

mean_power = mean(tdata)
adjusted_tdata = [k/mean_power for k in tdata]
print("Peaking factor is {:4f}".format(max(adjusted_tdata)))
fig, ax = plot.subplots()
scattered = ax.scatter(coords[0], coords[1], c=adjusted_tdata, cmap=cm.jet, s=1024, marker='H')
colorbar = fig.colorbar(scattered, ax=ax)

ax.set_xlabel("x-position in core [cm]", fontsize=16)
ax.set_ylabel("y-position in core [cm]", fontsize=16)
ax.set_aspect('equal')
ax.set_xlim([-30, 30])
ax.set_ylim([-30, 30])
colorbar.ax.set_ylabel("Peaking Factor", fontsize=16)
colorbar.ax.tick_params(axis='both', which='major', labelsize=14)
# ax.set_title("Radial Peaking Factor Plot for the PRIME Core", fontsize=20)

plot.show()
