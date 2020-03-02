#! /usr/bin/python3
from mcnptools import Mctal
tallyfile = Mctal("axial_peaking_tally_1.3.20")

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
print(fit_params)
fig, ax = plot.subplots()
ax.scatter(x, data, marker='+')
ax.plot(x, cosine_fit_func(x, fit_params[0], fit_params[1]),
        color='tab:red')
plot.show()
