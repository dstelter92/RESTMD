#!/bin/bash
#
#$ -l h_rt=7:00:00
#$ -j y
#$ -m e
#$ -N metenk
#$ -pe mpi_16_tasks_per_node 16

# Script to perform a long RESTMD run using PBS submission.
# Script will restart itself, runs for total of 20 M steps
# if the resubmission is uncommented. 

NSLOTS=8
NREPLICA=2

LAST=$(cat last)
PROCPERPART=$(($NSLOTS/$NREPLICA))
PROC=${NREPLICA}x$PROCPERPART

LMP=~/LAMMPS/lammps-dev/src/lmp_ubuntu.openmp

# Run LAMMPS, be sure to edit the path!
if [ $LAST -eq 1 ]
then
  INP=run.in.restmd
else
  INP=restart.in.restmd
fi

mpirun -np $NSLOTS $LMP -p $PROC -i $INP
if [ $? -ne 0 ]; then
  echo "Simulation did not complete."
  exit 1
fi
sleep 2
 
./rename.sh $NREPLICA $LAST

# Update Walkers
tail -1 log.lammps | awk '{for (i=2; i<=NF; i++) print $i}' | tr '\n' ' ' > last_walkers
sed -e "s/ZZZ/$(sed 's:/:\\/:g' last_walkers)/" restart.in.restmd.template > restart.in.restmd
wait ${!}


# Advance index counter
LAST=$(($LAST + 1))
echo $LAST > last
sleep 2

exit 0
