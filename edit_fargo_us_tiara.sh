


FILE=/tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare
if [ -f "$FILE" ]; then
    mv /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.par
#else 
#    echo "$FILE does not exist."
fi




gedit /tiara/home/mlehmann/data/FARGO3D/setups/mkl/{mkl.par,mkl2D_us.opt,mkl2D_NODUST_us.opt,condinit2D_EQ_us.c,mkl.bound2D_us_unp.0,mkl.bound2D_us_unp.1,stockholm.c,algogas.c,substep2_b.c,boundaries.txt}  &
