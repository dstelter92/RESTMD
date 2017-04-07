#!/bin/bash

###
# Script to run and get all data from Met-Enk RESTMD run.
# Usage:
# ./run_analysis.bash #replicas start_indx stop_indx
# ./run_analysis.bash 2 10 $(cat ../last)
# 
# Typically, you would want to start collecting production data
# after STG4 has started, AND the system has time to equilibrate
# Ts(E). This system at this temperature range is easy to
# equilibrate, but a good way to check is by the flatness of the
# energy histogram, or ensure that the full energy range is
# sampled. In this example, I start data collection after 10 M
# steps. 
###

NREP=$1
START=$2
STOP=$3

# Get enthalpies per replica
python get-hist-enthalpies.py $START $STOP
sleep 2

# Sort WT files by replica
python sort_by_replica.py $NREP $START $STOP
sleep 2

# Average production WT files, sorted above, for each replica
for ((i=0; i<$NREP; i++))
do
  python avg_WT.py $START $STOP $i
done
sleep 2

# Run ST-WHAM using averaged, production WT file to calculate sampling weight
python st-wham_RESTMD.py
sleep 1

exit 0
