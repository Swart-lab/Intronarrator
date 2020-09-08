#!/usr/bin/env python
from sys import argv
from matplotlib import *
import matplotlib.pyplot as plt
from numpy import *
from matplotlib.ticker import PercentFormatter

intron_lens = []
for line in open(argv[1]):
    if line[0] != '#':
        atoms = line.split()
        if atoms[2] == "intron":
            alen = int(atoms[4]) - int(atoms[3]) +1
            intron_lens.append(alen)

            #if alen == 12:
            #    print(alen, line[:-1])
                

plt.hist(intron_lens, bins=range(10, 35), ec='white', align='left', 
    weights=ones(len(intron_lens))/len(intron_lens))
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

plt.savefig(argv[2])
