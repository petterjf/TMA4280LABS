#!/bin/bash
#PBS -N poisson_hybrid_14
#PBS -A imf_lille-tma4280
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:01:00
#PBS -l select=1:ncpus=20:mpiprocs=16
#PBS -o stdout
#PBS -e stderr
#PBS -M petterjf@stud.ntnu.no
#PBS -q training

cd $PBS_O_WORKDIR
module load gcc/6.3.0
module load openmpi/2.0.1
mpirun pure.out 16384