#!/usr/bin/env python2.7

import os, sys
from numpy import *

##################
### For RESTMD ###
##################

### PARAMETERS ###
NumReplica = 4
start = 150
stop = 200

Elo = -300000
Ehi = -100000
binsize = 512
steps = 500000
exchange = 10000
thermo = 100
##################
header = 92
footer = 26

length = steps/thermo
length_walk = steps/exchange
ratio = exchange/thermo
nbins = int(-(Elo - Ehi) / binsize) + 1

print "Gather enthalpies from log files, follow exchanges, save data for each replica instead of per walker"
print "# bins:", nbins, "Elo:", Elo, "Ehi:", Ehi

size = ((length)*(stop-start),3*NumReplica+1)
size_walk = ((length)*(stop-start),NumReplica+1)
full_data = empty(size)
full_data[:] = NAN
# Format: step, walk1, walk2,...walkN, rep1, rep2,...repN

first_flag = 0

for i in range(NumReplica):
    print "\n===--- Walker", i, "---==="
    for j in range(stop-start):
        # sum over data
        hname = 'log/log.lammps.' + str(i) + '-' + str(start+j)
        wname = 'log/log.lammps-' + str(start+j)

        Hist = genfromtxt(hname, skip_header=header, skip_footer=footer, usecols = (0,2))

        print "Walker: ", i, "Data index: ", j
        print "  Step start: ", Hist[0][0], "Step end: ", Hist[-1][0]
        print "  Energy: min: ", amin(Hist[:,1]), "max: ", amax(Hist[:,1])

        Walk = genfromtxt(wname, skip_header=3, skip_footer=1)

        for k in range(length):
        # Load enthalpy data in global array
            indx = k+(j*length)
            if (first_flag == 0):
                full_data[:,0][indx] = Hist[:,0][k]
            full_data[:,i+1][indx] = Hist[:,1][k]
            full_data[:,i+NumReplica+1][indx] = int(Walk[:,i+1][k/ratio])
            #print k, k+(j*length), k/ratio

            for rep in range(NumReplica):
                if (full_data[:,i+NumReplica+1][indx] == rep):
                    full_data[:,rep+(2*NumReplica)+1][indx] = full_data[:,i+1][k]
                if (full_data[:,i+NumReplica+1][indx] == NAN):
                    sys.exit()
    first_flag = 1


print "Saving outputs..."
#savetxt('full_replica_data.dat', full_data)
walkers = empty(size_walk)
for rep in range(NumReplica):
    Rname = 'replica-' + str(rep) + '.dat'
    Wname = 'walker-' + str(rep) + '.dat'
    savetxt(Rname, full_data[:, (rep+(2*NumReplica)+1)])
    savetxt(Wname, full_data[:, (rep+1)])
    walkers[:,rep+1] = full_data[:,rep+NumReplica+1]

walkers[:,0] = full_data[:,0]
savetxt('exchange_list.dat', walkers)

sys.exit()
