#!/bin/bash
#
#$ -l h_rt=7:00:00
#$ -j y
#$ -m e
#$ -N metenk
#$ -pe mpi_16_tasks_per_node 16

NSLOTS=16
NREPLICA=1

LAST=$(cat last)
PROCPERPART=$(($NSLOTS/$NREPLICA))
PROC=${NREPLICA}x$PROCPERPART

if [ $LAST -eq 1 ]
then
  mpirun -np $NSLOTS /path/to/lmp_exe -in run.in.restmd > output
  ./rename.sh $NREPLICA $LAST
else
  mpirun -np $NSLOTS /path/to/lmp_exe -in restart.in.restmd > output
  ./rename.sh $NREPLICA $LAST
fi
sleep 2

LAST=$(($LAST + 1))
echo $LAST > last
sleep 2

# (uncomment to) resubmit
#if [ $LAST -lt 21 ]
#then
#  qsub submit.sh
#fi

exit 0
