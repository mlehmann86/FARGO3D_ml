//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void PTorque_cpu(real dt) {

//<USER_DEFINED>
  INPUT(Vx_temp);
  OUTPUT(Vx_temp);

  real orbits, planetmass_taper;
  
  orbits = PhysicalTime/(2.0*M_PI); 

  if (orbits < PTORQUEON) {
  planetmass_taper = 0.0; 
  } else if( (orbits >= PTORQUEON) && (orbits < (PTORQUEON + PTORQUETAPER)) ){
  planetmass_taper = pow(sin(0.5*M_PI*(orbits-PTORQUEON)/PTORQUETAPER), 2.0);
  } else {
  planetmass_taper = 1.0; 
  }

//<\USER_DEFINED>

//<EXTERNAL>
  real* vx = Vx->field_cpu;
  real* vx_temp = Vx_temp->field_cpu;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP;
  int size_y = Ny+2*NGHY;
  int size_z = Nz+2*NGHZ;
  real taper = planetmass_taper; 
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  real x;
  real bigrad;
  real curlyf;
  real lambda;
  real ex1;
  real ex2;
  real h0;
  real omega0; 
  real p1 = 0.0297597;
  real p2 = 1.09770;
  real p3 = 0.938567;
  real p4 = 0.0421186;
  real p5 = 0.902328;
  real p6 = 1.03579;
  real p7 = 0.0981183;
  real p8 = 4.68108;
  real delta=1.3;
  real cconst=0.798; 
//<\INTERNAL>

//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real ASPECTRATIO(1);
// real PTORQUEMASS(1);
//<\CONSTANT>


//<MAIN_LOOP>
  for (k=1; k<size_z; k++) {
    for (j=1; j<size_y; j++) {
      for (i=XIM; i<size_x; i++) {
//<#>
	h0 = ASPECTRATIO*R0;
	omega0 = sqrt(G*MSTAR/R0/R0/R0);
	bigrad = ymed(j)*sin(zmed(k));
	x = (bigrad - R0)/h0;


	ex1 = (x+p2)/p3;
	ex1*= ex1;
	ex2 = (x-p5)/p6;
	ex2*= ex2;
	curlyf = p1*exp(-ex1) + p4*exp(-ex2);
	curlyf*= tanh(p7 - p8*x);
	lambda =-curlyf*pow(omega0*R0, 2.0);
	lambda*= pow(PTORQUEMASS*taper, 2.0);
	lambda*= pow(R0/h0, 4.0)/bigrad;

/*
        lambda = cconst*R0*pow(PTORQUEMASS*taper,2.0);
        lambda*= pow(R0*omega0, 2.0);
        lambda/= 2.0*bigrad*bigrad;
        if (x < -delta){
        lambda *= -pow(R0/(bigrad-R0),4.0);
        } else if (x > delta){
        lambda *=  pow(R0/(bigrad-R0),4.0);
        } else {
        lambda *=0.0; 
        }
*/
        vx_temp[l] = vx_temp[l] + dt*lambda; 
//<\#>
      }
    }
  }
//<\MAIN_LOOP>
}
