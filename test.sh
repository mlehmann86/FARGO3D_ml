#!/bin/bash

#SBATCH --cpus-per-task=4 

#SBATCH --gres=gpu:8             

#SBATCH --job-name=fargo3d

#SBATCH --output=%x-%j.txt

#SBATCH --error=%x-%j.err

#SBATCH --ntasks-per-node=8

#SBATCH --nodes=1

#SBATCH --time=120:00:00

#SBATCH --partition v100     

#SBATCH --mail-type=BEGIN,END




# Load necessary modules

module load compiler/2022.1.0

module load cuda/11.8.0

module load intel_mpi/2021.6.0


# Set MPI environment variables

export I_MPI_FABRICS=shm:ofi

export I_MPI_DEBUG=5


# Check CUDA installation

nvcc --version

nvidia-smi


# Check MPI installation

mpirun --version


# Ensure correct paths
export PATH=$MPI_HOME/bin:$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$MPI_HOME/lib:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
/ceph/sharedfs/software/compiler/2022.2.0.262/intel/oneapi/mpi/2021.6.0/bin/mpirun ./fargo3d_cos_b1d0_us_St1dm2_Z1dm2_r6H_z08H_fim053_ss203_3D_2PI_stnew_LR150 setups/mkl/parfiles/cos_b1d0_us_St1dm2_Z1dm2_r6H_z08H_fim053_ss203_3D_2PI_stnew_LR150.par -> cos_b1d0_us_St1dm2_Z1dm2_r6H_z08H_fim053_ss203_3D_2PI_stnew_LR150.out
