


FILE=/tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare
if [ -f "$FILE" ]; then
    mv /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.par
#else 
#    echo "$FILE does not exist."
fi




gedit /tiara/home/mlehmann/data/FARGO3D/setups/mkl/{mkl.par,mkl3D_us.opt,mkl3D_NODUST_us.opt,condinit3D_us.c,condinit3D_EQ_us.c,mkl.bound3D_us.0,mkl.bound3D_us.1,stockholm.c,algogas.c,substep2_b.c,edamp.c}  &
