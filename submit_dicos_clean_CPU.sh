#!/bin/bash
#SBATCH --job-name=fargo3d
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=xxxx
#SBATCH --nodes=1
#SBATCH --ntasks=yy
#SBATCH --time=2-00:00:00

# Load modules to ensure a consistent environment for the job
echo "Loading modules..."
module purge
# Ensure the runtime environment matches the compilation environment
module load gcc/13.1.0
module load mpi/openmpi-5.0.3/gcc13.1.0
echo "Modules loaded."

# The specific 'mpirun' command will be appended to this script
# by the submit_fargo function.
echo "Starting FARGO3D simulation..."

