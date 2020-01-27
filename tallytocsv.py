#! /usr/bin/python3

from sys import argv
from mcnptools import Mctal

mct = Mctal(argv[1])
tally = mct.GetTally(int(argv[2]))

with open(argv[3], "w+") as outfile:
    ebins = tally.GetEBins()
    fbins = tally.GetFBins()
    for i in range(len(ebins)):
        outfile.write("{}, ".format(ebins[i]))
        for f in range(len(fbins)):
            val = tally.GetValue(f, 0, 0, 0, 0, 0, i, 0)
            outfile.write("{}, ".format(val))
        outfile.write("0\n")
