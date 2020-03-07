#! /usr/bin/python3
from mcnptools import Mctal
tallyfile = Mctal("axial_peaking_tally")

import numpy as np
import matplotlib.pyplot as plot
from scipy import optimize
from matplotlib import cm
from statistics import mean

tally = tallyfile.GetTally(204)
data = [
    tally.GetValue(f, 0, 0, 0, 0, 0, 0, 0, 0)
    for f in range(len(tally.GetFBins()))
]
data = [k/mean(data) for k in data]
x = np.linspace(-118.4, 118.4, len(data))


def cosine_fit_func(x, a, b):
    return a*np.cos(b*x)


fit_params, fit_cov = optimize.curve_fit(cosine_fit_func,
                                         x, data,
                                         p0=[1, 0.001])
print(max(data))
print(fit_params)
fig, ax = plot.subplots()
ax.plot(x, cosine_fit_func(x, fit_params[0], fit_params[1]),
        color='tab:red', label="Axial Peaking Factor, cosine fit", linewidth=4)
ax.scatter(x, data, marker='+')
ax.set_xlabel("z-position in center assembly [cm]", fontsize=16)
ax.set_ylabel("Axial peaking factor", fontsize=16)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.legend()
plot.show()
