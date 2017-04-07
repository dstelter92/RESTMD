#!/usr/bin/env python2.7

import os, sys
from numpy import *

##################
### For RESTMD ###
##################

### PARAMETERS ###
NumReplica = 2
start = int(sys.argv[1])
stop = int(sys.argv[2])

steps = 1000000
exchange = 1000
thermo = 50

header = 98 # match to log.lammps.0 files
footer = 57
##################

length = steps/thermo
length_walk = steps/exchange
ratio = exchange/thermo

print "Gather enthalpies from log files, follow exchanges, save data for each replica instead of per walker"

size = ((length)*(stop-start),3*NumReplica+1)
size_walk = ((length)*(stop-start),NumReplica+1)
full_data = empty(size)
full_data[:] = NAN
# Format: step, indx1, indx2,... indxN, walk1, walk2,...walkN, rep1, rep2,...repN

# Save Exchanges
for i in range(stop-start):
    wname = '../log.lammps' + '-' + str(start+i)
    Walk = genfromtxt(wname, skip_header=3, skip_footer=1)
    for j in range(NumReplica):
        for k in range(length):
            indx = k+(i*length)
            if (exchange == steps):
                full_data[:,j+1][indx] = int(Walk[j+1])
            else:
                full_data[:,j+1][indx] = int(Walk[:,j+1][k/ratio])

# Save Walkers
for i in range(NumReplica):
    print "\n===--- Walker", i, "---==="
    for j in range(stop-start):
        hname = '../' + str(i) + '/log.lammps.' + str(i) + '-' + str(start+j)
        Hist = genfromtxt(hname, skip_header=header, skip_footer=footer, usecols = (0,2))
        print "Walker: ", i, "Data index: ", j+start
        print "Step start: ", Hist[0][0], "Step end: ", Hist[-1][0]
        print "Energy: min: ", amin(Hist[:,1]), "max: ", amax(Hist[:,1])

        for k in range(length):
        # Load enthalpy data in global array
            indx = k+(j*length)
            #if (first_flag == 0):
                # Put time in full_data
            full_data[:,0][indx] = Hist[:,0][k]
            full_data[:,i+NumReplica+1][indx] = Hist[:,1][k]
            #print k, k+(i*length), k/ratio

# Save Replicas
for rep in range(NumReplica):
    for indx in range(len(full_data[:,0])):
        for j in range(NumReplica):
            if (full_data[:,j+1][indx] == rep):
                full_data[:,rep+(2*NumReplica)+1][indx] = full_data[:,j+NumReplica+1][indx]


print "\nSaving outputs..."
savetxt('full_replica_data.dat', full_data)
walkers = empty(size_walk, dtype=int)
for rep in range(NumReplica):
    Rname = 'replica-' + str(rep) + '.dat'
    Wname = 'walker-' + str(rep) + '.dat'
    savetxt(Rname, full_data[:, (rep+(2*NumReplica)+1)])
    savetxt(Wname, full_data[:, (rep+NumReplica+1)])
    walkers[:,rep+1] = full_data[:,rep+1].astype(int)

walkers[:,0] = full_data[:,0].astype(int)
savetxt('exchange_list.dat', walkers, fmt='%d')

sys.exit()
