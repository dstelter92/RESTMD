#!/usr/bin/env python2.7

import os, sys
from numpy import *

##################
### For RESTMD ###
##################

### PARAMETERS ###
NumReplica = 2
start = 55
stop = 101

Elo = -7000
Ehi = -5500
binsize = 4
steps = 1000
exchange = 1000
thermo = 50
##################
header = 88
footer = 53

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
        hname = str(i) + '/log.lammps.' + str(i) + '-' + str(start+j)
        wname = 'log.lammps' + '-' + str(start+j)

        Hist = genfromtxt(hname, skip_header=header, skip_footer=footer, usecols = (0,2))
        print "Walker: ", i, "Data index: ", j+start
        print "Step start: ", Hist[0][0], "Step end: ", Hist[-1][0]
        print "Energy: min: ", amin(Hist[:,1]), "max: ", amax(Hist[:,1])
        Walk = genfromtxt(wname, skip_header=3, skip_footer=1)


        for k in range(length):
        # Load enthalpy data in global array
            indx = k+(j*length)
            if (first_flag == 0):
                full_data[:,0][indx] = Hist[:,0][k]
            full_data[:,i+1][indx] = Hist[:,1][k]
            if (exchange == steps):
                full_data[:,i+NumReplica+1][indx] = int(Walk[i+1])
            else:
                full_data[:,i+NumReplica+1][indx] = int(Walk[:,i+1][k/ratio])
            #print k, k+(j*length), k/ratio

            for rep in range(NumReplica):
                if (full_data[:,i+NumReplica+1][indx] == rep):
                    full_data[:,rep+(2*NumReplica)+1][indx] = full_data[:,i+1][k]
                if (full_data[:,i+NumReplica+1][indx] == NAN):
                    sys.exit()
    first_flag = 1


print "\nSaving outputs..."
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
