#! /usr/bin/python3

from sys import argv
from mcnptools import Mctal

mct = Mctal(argv[1])
tally = mct.GetTally(int(argv[2]))

with open(argv[3], "w+") as outfile:
    ebins = tally.GetEBins()
    for i in range(len(ebins)):
        val = tally.GetValue(0, 0, 0, 0, 0, 0, i, 0)
        outfile.write("{} , {}\n".format(ebins[i], val))
