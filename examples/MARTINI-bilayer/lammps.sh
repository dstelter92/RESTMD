#!/bin/bash
#
#$ -l h_rt=1:00:00
#$ -j y
#$ -N restmd-test
#$ -pe mpi_16_tasks_per_node 16

export MPI_IMPLEMENTATION=openmpi

mpirun -np $NSLOTS /project/normoliq/dstelter/LAMMPS+STMD/test-build/src/lmp_openmpi.scc -p 2x8 -in run.in.restmd > output
#mpirun -np 1 /project/normoliq/dstelter/LAMMPS+STMD/LAMMPS-r12125/src/lmp_openmpi.scc < in.sw 
