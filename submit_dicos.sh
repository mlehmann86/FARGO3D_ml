#!/bin/bash
#SBATCH --job-name=fargo3d
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=edr1-al9_large
#SBATCH --nodes=1
#SBATCH --ntasks=192
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

mpirun -np 192 ./fargo3d_cos_bet1d1_gam53_ss15_q1_r0516_nu1dm11_COR_HR150_nu1dm7_2D setups/mkl/mkl.par > cos_bet1d1_gam53_ss15_q1_r0516_nu1dm11_COR_HR150_nu1dm7_2D.out
