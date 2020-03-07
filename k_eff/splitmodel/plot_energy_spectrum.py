from mcnptools import Mctal
tallyfile = Mctal("mctan")

import numpy as np
import matplotlib.pyplot as plot
from matplotlib import cm
from statistics import mean

tally204 = tallyfile.GetTally(34)

ebins = tally204.GetEBins()
# F D U S M C E T
vals = [tally204.GetValue(0, 0, 0, 0, 0, 0, i, 0) for i in range(len(ebins))]

fig, ax = plot.subplots()
ax.plot(ebins, vals)
ax.set_xlabel("Energy, MeV", fontsize=14)
ax.set_ylabel("Flux per source particle rate", fontsize=14)
plot.xscale('log')
fig.tight_layout()
plot.show()
