//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void dgfloor_cpu() {
  
//<USER_DEFINED>
  int ii;
  for (ii=0; ii<NFLUIDS; ii++) {
    INPUT(Fluids[ii]->Density);
    OUTPUT(Fluids[ii]->Density);
}
//<\USER_DEFINED>

//<EXTERNAL>
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real* rhog  = Fluids[0]->Density->field_cpu;
  real* rhod1 = Fluids[1]->Density->field_cpu;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  int  o;
//<\INTERNAL>
  
//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=1; k<size_z; k++) {
#endif
#ifdef Y
    for(j=1; j<size_y; j++) {
#endif
#ifdef X
      for(i=XIM; i<size_x; i++) {
#endif
//<#>
	

//enforce minimum d/g
//for(o=1; o<NFLUIDS; o++){
//dens[o][l] = dens[o][l]/dens[0][l] > 1.0e-16 ? dens[o][l]:1.0e-16*dens[0][l];
//rho[o][l] = rho[o][l]/rho[0][l] > 1.0e-16 ? rho[o][l]:1.0e-16*rho[0][l];
//}

rhod1[l] = rhod1[l]/rhog[l] > 1.0e-16 ? rhod1[l]:1.0e-16*rhog[l];

//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}

