#!/bin/bash

cd /home/u1264970/FARGO3D 

makegpu='make SETUP=mkl PARALLEL=1 GPU=1 MPICUDA=1'

#---------------------------------
echo "cleaning up" &
make mrproper
process_id=$!
wait $process_id

#---------------------------------
echo "loading necessary modules" 
module purge &

process_id=$!
wait $process_id

module add "nvidia/cuda/9.2" 

process_id=$!
wait $process_id

module add "openmpi/3.1.4"

process_id=$!
wait $process_id

#---------------------------------
echo "compiling fargo3d with current setup and GPU-MPI CUDA" &

$makegpu
process_id=$!
wait $process_id
echo "compilation finished" &
echo "copying executable and batch file into work directory"

cp /home/u1264970/FARGO3D/fargo3d /work/u1264970/FARGO3D/. & cp /home/u1264970/FARGO3D/submit.sh /work/u1264970/FARGO3D/. 

process_id=$!
wait $process_id

cd /work/u1264970/FARGO3D/ &

process_id=$!
wait $process_id

sbatch submit.sh


echo "Job submitted--------------------."

