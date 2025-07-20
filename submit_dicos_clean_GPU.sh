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

# Load the appropriate modules for the selected GPU architecture.
###MODULE_LOAD_COMMANDS###

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
