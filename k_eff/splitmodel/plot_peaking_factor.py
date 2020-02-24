#! /usr/bin/python3
from mcnptools import Mctal
tallyfile = Mctal("tally_19.2.20")

import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
from statistics import mean
from math import sin, cos, pi

j_hat = [1, 0]
i_hat = [0.5, -0.86603]

core_i_hat = [0.86603, 0.5]
core_j_hat = [0, 1]

assm_basis = np.concatenate((i_hat, j_hat)).reshape(2, 2).transpose()
core_basis = np.dot(8, np.concatenate((core_i_hat,
                                       core_j_hat)).reshape(2, 2).transpose())

core_ij = np.array([[0, 1, 1, 2],
                    [0, 0, 1, 0]])
core_disp = np.dot(core_basis, core_ij)

assm_ij = []

imin = 0
imax = 4
for j in range(-4, 1):
    row_ij = np.concatenate((np.arange(imin, imax+1),
                             j*np.ones(imax-imin+1))).reshape(2,imax-imin+1).transpose()
    assm_ij.extend(row_ij)
    imin -= 1

imax -= 1
imin = -4
for j in range(1, 5):
    row_ij = np.concatenate((np.arange(imin, imax+1),
                             j*np.ones(imax-imin+1))).reshape(2,imax-imin+1).transpose()
    assm_ij.extend(row_ij)
    imax -= 1

x = []
y = []
assm_xy = np.dot(assm_basis, np.array(assm_ij).transpose())
assm_xy = np.add([[core_disp[0][0]], [core_disp[1][0]]], assm_xy)
x.extend(assm_xy[0])
y.extend(assm_xy[1])
for assm in range(1, 4):
    for k in range(0, 6):
        rotate_60 = np.array([[cos(k*pi/3), -sin(k*pi/3)], [sin(k*pi/3), cos(k*pi/3)]])
        assm_xy = np.dot(assm_basis, np.array(assm_ij).transpose())
        assm_xy = np.add([[core_disp[0][assm]],[core_disp[1][assm]]], assm_xy)
        assm_xy = np.dot(rotate_60, assm_xy)
        x.extend(assm_xy[0])
        y.extend(assm_xy[1])
        
power_tally = []
tallies = (114, 124, 134, 144)
for t in range(len(tallies)):
    tally = tallyfile.GetTally(tallies[t])
    num_fbins = len(tally.GetFBins())
    tdata = [tally.GetValue(i,0,0,0,0,0,0,0)
             for i in range(num_fbins)]
    power_tally.extend(tdata)
    if t > 0:
        for i in range(5):
            power_tally.extend(tdata)

weighted_power_tally = power_tally
nonzero = [k for k in weighted_power_tally if k > 0]
mean_power = mean(nonzero)
adjusted_power = [k/mean_power for k in power_tally]
print(max(adjusted_power))

xy_list = np.dot(assm_basis, np.array(assm_ij).transpose())
fig, ax = plot.subplots()
scattered = ax.scatter(x, y, c=adjusted_power, cmap=cm.jet, s=128, marker="h")
colorbar = fig.colorbar(scattered, ax=ax)

ax.set_xlabel("x-position in core", fontsize=16)
ax.set_ylabel("y-position in core", fontsize=16)
ax.set_aspect('equal')
colorbar.ax.set_ylabel("Peaking Factor", fontsize=14)
ax.set_title("Peaking Factor Plot for the PRIME Core", fontsize=20)
fig.tight_layout()
plot.show()
