#!/bin/bash
#PBS -l select=1:ncpus=40:ngpus=4:mpiprocs=4
#PBS -l walltime=48:00:00
#PBS -N FARGO3D
#PBS -q gp4
#PBS -P MST107457
#PBS -j oe

module purge
module load intel/2018_u1 
module load mpi/openmpi-3.0.0/intel2018u1
module load cuda/10.0.130

cd $PBS_O_WORKDIR

mpiexec -mca pml cm -mca mtl psm2 ./fargo3d_mpi -t mkl.par > output
