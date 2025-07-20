


FILE=/tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare
if [ -f "$FILE" ]; then
    mv /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.par
#else 
#    echo "$FILE does not exist."
fi




gedit /tiara/home/mlehmann/data/FARGO3D/setups/mkl//{mkl.par,mkl3D.opt,condinit3D.c,condinit3D_EQ.c,mkl.bound3D.0,mkl.bound3D.1,edamp.c,stockholm.c,substep2_b.c,algogas.c}  &
