#!/bin/bash
###### Job name ######
#PBS -N fargo3d
###### Output files ######
#PBS -o fargo3d.out
#PBS -e fargo3d.err
###### Number of nodes and cores (PBS Pro syntax) ######
#PBS -l select=xx:ncpus=yy:mpiprocs=yy
#PBS -l walltime=48:00:00
###### Queue name ######
#PBS -q xxx
###### Sends mail to yourself when the job begins and ends ######
##PBS -M mlehmann@asiaa.sinica.edu.tw
#PBS -m be
###### Specific the shell types ######
#PBS -S /bin/bash

###### Enter this job's working directory ######
cd $PBS_O_WORKDIR

###### Load modules to setup environment ######
. /etc/profile.d/modules.sh

module purge
module load cuda/11.8 intel/2020u4 openmpi/4.1.1
module load python/3.10.2

###### Run your jobs with parameters ######
if [ -n "$PBS_NODEFILE" ]; then
  if [ -f "$PBS_NODEFILE" ]; then
    NPROCS=$(wc -l < "$PBS_NODEFILE")
  fi
fi

# Actual run command will be appended by the function
