#!/bin/bash

# Rename files

NREPLICA=$1
LAST=$2

cp log.lammps log.lammps-$LAST
for ((i=0; i<$NREPLICA; i++));
do
  cp WH.$i.d $i/WH.$i.d-$LAST
  cp WT.$i.d $i/WT.$i.d-$LAST
  cp WHP.$i.d $i/WHP.$i.d-$LAST
  cp log.lammps.$i $i/log.lammps.$i-$LAST
  cp oREST.$i.d $i/oREST.$i.d-$LAST
  cp $i.dcd $i/$i-$LAST.dcd
done

sleep 2
exit 0
