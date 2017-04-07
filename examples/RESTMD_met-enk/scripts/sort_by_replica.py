#!/usr/bin/env python2.7

import os,sys
from numpy import *

# Sorts files written as walkers, to new directory ordered
# as replicas. NOTE: Assumes that all files were written on
# the last step of index.

nrep = int(sys.argv[1])
start = int(sys.argv[2])
stop = int(sys.argv[3])

for j in range(stop - start):
    # Read exchange data
    exchanges = genfromtxt('../log.lammps-' + str(j+start),skip_header=1003)
    print exchanges
    for walk in range(nrep):
        # Loop over walkers
        basename = 'WT.' + str(walk) + '.d-' + str(j+start)
        path = '../' + str(walk) + '/'
        fname = path + basename
        for rep in range(nrep):
            # Loop over replicas
            #if (rep == int(exchanges[rep+1])):
            os.system("cp %s ../replica_data/%d/WT-%d" % (fname,exchanges[walk+1],j+start))

sys.exit()
