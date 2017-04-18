#!/bin/bash

start=$1
stop=$2
init=1          # if init=1, initialize files, otherwise just append

if [ $init -eq 1 ]
then
    > energies.dat
fi

for ((i=$start; i<$stop; i++))
do
    echo -ne "File index: $i"\\r
    grep -v '[a-zA-Z]' ../log.lammps-$i | awk '{if(NF==7 && $1%100==0) print $3/5100}' | head -n -5 | tail -n +9 >> energies.dat
done

exit 0
