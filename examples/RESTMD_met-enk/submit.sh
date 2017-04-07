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

NSLOTS=16
NREPLICA=2

LAST=$(cat last)
PROCPERPART=$(($NSLOTS/$NREPLICA))
PROC=${NREPLICA}x$PROCPERPART


# Run LAMMPS, be sure to edit the path!
if [ $LAST -eq 1 ]
then
  mpirun -np $NSLOTS /path/to/lmp_exe -p $PROC -in run.in.restmd
  ./rename.sh $NREPLICA $LAST
else
  mpirun -np $NSLOTS /path/to/lmp_exe -p $PROC -in restart.in.restmd
  ./rename.sh $NREPLICA $LAST
fi
sleep 2

# Update Walkers
tail -1 log.lammps | awk '{for (i=2; i<=NF; i++) print $i}' | tr '\n' ' ' > last_walkers
sed -e "s/ZZZ/$(sed 's:/:\\/:g' last_walkers)/" restart.in.restmd.template > restart.in.restmd

# Advance index counter
LAST=$(($LAST + 1))
echo $LAST > last
sleep 2

# (Uncomment to) Resubmit
#if [ $LAST -lt 20 ]
#then
#  qsub submit.sh
#  #./submit.sh
#fi

exit 0
