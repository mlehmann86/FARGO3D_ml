


FILE=/tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare
if [ -f "$FILE" ]; then
    mv /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.pare /tiara/home/mlehmann/data/FARGO3D/setups/mkl/mkl.par
#else 
#    echo "$FILE does not exist."
fi




gedit /tiara/home/mlehmann/data/FARGO3D/setups/mkl/{mkl.par,mkl2D.opt,condinit2D.c,condinit2D_EQ.c}  &
