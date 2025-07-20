#!/bin/bash
###### Job name ######
#PBS -N fargo3d
###### Output files ######
#PBS -o fargo3d.out
#PBS -e fargo3d.err
###### Number of nodes and cores ######
#PBS -l nodes=2:ppn=16
#PBS -l walltime=48:00:00
###### Queue name ######
#PBS -q small
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
module add "openmpi/4.0.4_ic19.0"


###### Run your jobs with parameters ######
if [ -n "$PBS_NODEFILE" ]; then
  if [ -f $PBS_NODEFILE ]; then
    NPROCS=`wc -l < $PBS_NODEFILE`
  fi
fi
#
#
#
#
#
#

$OPENMPI_HOME/bin/mpirun -v -machinefile $PBS_NODEFILE -np $NPROCS ./fargo3d_cos_nudz1dm7_bet1dm3_Z1dm1_St1dm2_r2H_z08H_LR_EQ_xl setups/mkl/parfiles/cos_nudz1dm7_bet1dm3_Z1dm1_St1dm2_r2H_z08H_LR_EQ_xl.par  -> cos_nudz1dm7_bet1dm3_Z1dm1_St1dm2_r2H_z08H_LR_EQ_xl.out
