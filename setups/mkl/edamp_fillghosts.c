//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void Edamp_fillghosts_cpu(real dt) {


//<USER_DEFINED>
  // These are INPUT and OUTPUT because we read the physical zone
  // values and write new values into the ghost zones of the same array.
  INPUT(T_guess);
  INPUT(Div_old);
  OUTPUT(T_guess);
  OUTPUT(Div_old);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
//<\INTERNAL>

//<EXTERNAL>
  real* T_predict = T_guess->field_cpu;
  real* D_old     = Div_old->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny + 2 * NGHY;
  int size_z = Nz + 2 * NGHZ;
//<\EXTERNAL>

//<CONSTANT>
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

  for (k = 0; k < size_z; k++) {
    for (j = 0; j < size_y; j++) {
      for (i = 0; i < size_x; i++) {
//<#>
        ll = i + j * pitch + k * stride;

        int Ny = size_y - 2 * NGHY;
        int Nz = size_z - 2 * NGHZ;

        // This kernel executes on all cells. We check if the current cell
        // is a ghost cell. If so, we fill it. Otherwise, we do nothing.
        // This assumes X is periodic and does not have ghost zones to fill here.
        if (j < NGHY || j >= Ny + NGHY || k < NGHZ || k >= Nz + NGHZ) {

            int j_src = j;
            int k_src = k;

            // Clamp the j-index to the boundary of the physical domain
            if (j < NGHY) j_src = NGHY;
            if (j >= Ny + NGHY) j_src = Ny + NGHY - 1;

            // Clamp the k-index to the boundary of the physical domain
            #ifdef Z
              if (k < NGHZ) k_src = NGHZ;
              if (k >= Nz + NGHZ) k_src = Nz + NGHZ - 1;
            #endif

            // Get the linear index of the source physical cell
            int ll_src = i + j_src * pitch + k_src * stride;

            // Copy the data from the physical source cell to this ghost cell
            T_predict[ll] = T_predict[ll_src];
            D_old[ll]     = D_old[ll_src];
        }
//<\#>
      }
    }
  }
//<\MAIN_LOOP>

}
