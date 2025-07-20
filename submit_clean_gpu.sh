#!/bin/bash
###### Job name ######
#PBS -N fargo3d_gpu
###### Output files ######
#PBS -o fargo3d_gpu.out
#PBS -e fargo3d_gpu.err
###### Number of nodes and cores ######
#PBS -l select=1:ncpus=yy:mpiprocs=yy:ngpus=yy
###### Specify how many hours do you need ######
#PBS -l walltime=24:00:00
###### Queue name ######
#PBS -q gp1d
###### Sandbox ######
#PBS -W sandbox=PRIVATE
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh

module purge
#module load gcc/12.3.0 openmpi/4.1.4 cuda/11.8
#module load nvhpc/23.7
module load cuda/11.8 intel/2020u4 openmpi/4.1.1
module load python/3.10.2



###### Run your jobs with parameters ######
if [ -n "$PBS_NODEFILE" ]; then
  if [ -f "$PBS_NODEFILE" ]; then
    NPROCS=$(wc -l < "$PBS_NODEFILE")
  fi
fi

# Actual run command will be appended by the function
