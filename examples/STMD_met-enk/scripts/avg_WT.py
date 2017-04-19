#!/usr/bin/env python2.7

import os,sys
from numpy import *

# Import Ts from WT files, average

start = int(sys.argv[1])
stop = int(sys.argv[2])
rep = int(sys.argv[3])

path = '../0/'
basename = 'WT.0.d'

fname = path + basename + '-' + str(start)
num = stop - start
size = len(loadtxt(fname,))
data = empty((size,num))
final_data = empty((size,2))

final_data[:,0] = loadtxt(fname,usecols=(1,2))[:,0]

for i in range(stop - start):
    fname = path + basename + '-' + str(i+start)
    data[:,i] = loadtxt(fname,usecols=(1,2))[:,1]

for i in range(size):
    final_data[i][1] = average(data[i])

savetxt('Ts_avg.%d.dat'%rep,final_data)
#print final_data

sys.exit()
