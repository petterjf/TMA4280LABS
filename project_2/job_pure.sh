#!/bin/bash
#PBS -N n4096P1
#PBS -A imf_lille-tma4280
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=20:mpiprocs=16
#PBS -o stdout_n4096P1
#PBS -e stderr_n4096P1
#PBS -M petterjf@stud.ntnu.no
#PBS -q training

cd $PBS_O_WORKDIR
module load gcc/6.3.0
module load openmpi/2.0.1
mpirun pure.out 4096