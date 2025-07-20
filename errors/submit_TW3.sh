#!/bin/bash
#SBATCH -A MST110258       # Account name/project number
#SBATCH -J fargo3d          # Job name
#SBATCH -p ct224             # Partition name
#SBATCH --ntasks=224 # (-n) Number of MPI tasks (i.e. processes)
#SBATCH --cpus-per-task=1 # (-c) Number of cores per MPI task
#SBATCH --nodes=4 # (-N) Maximum number of nodes to be allocated
#SBATCH --ntasks-per-node=56 # Maximum number of tasks on each node
#SBATCH --ntasks-per-socket=28 # Maximum number of tasks on each socket
#SBATCH -o %j.out           # Path to the standard output file
#SBATCH -e %j.err           # Path to the standard error ouput file

#module load compiler/intel/2020u4 IntelMPI/2020
#mpiexec.hydra -bootstrap slurm -n 24 /home/user/bin/intel-hello


mpirun ./fargo3d_cos_b1dm3_St1dm2_Z1dm1_R2H_Z08H_HR1k_pres setups/mkl/parfiles/cos_b1dm3_St1dm2_Z1dm1_R2H_Z08H_HR1k_pres.par -> cos_b1dm3_St1dm2_Z1dm1_R2H_Z08H_HR1k_pres.out
