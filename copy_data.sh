#!/bin/bash

# --- Define array of simulation names ---
simulations=(
    cos_bet1d4_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D
    cos_bet1d3_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D
    cos_bet1d2_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D
)

# --- Paths ---
source_base="/tiara/home/mlehmann/data/FARGO3D/outputs"
dest_base="/theory/lts/mlehmann/FARGO3D/outputs"
parfile_base="setups/mkl/parfiles"
optfile_base="setups/mkl/optfiles"

# --- Loop over simulations ---
for sim in "${simulations[@]}"; do
    echo "Copying simulation folder for $sim..."
    cp -r "$source_base/$sim" "$dest_base/"

    echo "Copying .par file for $sim..."
    cp "$parfile_base/${sim}.par" "$dest_base/$sim/mkl.par"

    echo "Copying .opt file for $sim..."
    cp "$optfile_base/${sim}.opt" "$dest_base/$sim/mkl.opt"
done

# --- Final large print statement ---
echo "=================================================="
echo "========   ALL SIMULATIONS SUCCESSFULLY COPIED   ========"
echo "=================================================="
