#!/bin/bash
#SBATCH --cpus-per-task=4 
#SBATCH --gres=gpu:8             
#SBATCH --job-name=fargo3d
#SBATCH --output=%x-%j.txt
#SBATCH --error=%x-%j.err
#SBATCH --ntasks-per-node=8
#SBATCH --nodes=1
#SBATCH --time=120:00:00
#SBATCH --partition xxxx     
#SBATCH --mail-type=BEGIN,END


# Source the provided setup script
#source /ceph/sharedfs/pkg/mpich/4.2.2/gpu/setup.sh

# Purge any existing modules for a clean environment
module purge

# Load a fully compatible GCC, OpenMPI, and CUDA toolchain.
# The `mpi/openmpi-3.1.6/cuda/gcc930` module was explicitly built
# with CUDA awareness using gcc/9.3.0, preventing version conflicts.
# GCC 9.3.0 is fully supported by the CUDA 11.8 toolkit.
module load cuda/12.2.0 mpi/openmpi-4.1.5/cuda-12.2/gcc-11.5.0 gcc/11.5.0

# Load python if needed
module load python3/3.9.23

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
