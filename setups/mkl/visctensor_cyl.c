//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>

void visctensor_cyl_cpu(){

//<USER_DEFINED>
  INPUT2D(Viscosity_profile);
  INPUT(Density);
#ifdef ALPHAVISCOSITY
  INPUT(Energy);
#endif
#ifdef X
#ifdef COLLISIONPREDICTOR
  INPUT(Vx_half);
#else
  INPUT(Vx);
#endif
  OUTPUT(Mmx);
  OUTPUT(Mpx);
#endif
#ifdef Y
#ifdef COLLISIONPREDICTOR
  INPUT(Vy_half);
#else
  INPUT(Vy);
#endif
  OUTPUT(Mmy);
  OUTPUT(Mpy);
#endif
#ifdef Z
#ifdef COLLISIONPREDICTOR
  INPUT(Vz_half);
#else
  INPUT(Vz);
#endif
  OUTPUT(Mmz);
  OUTPUT(Mpz);
#endif
//<\USER_DEFINED>

//<EXTERNAL>
  real* rho = Density->field_cpu;
#ifdef ALPHAVISCOSITY
  real* energy = Energy->field_cpu;
#endif
#ifdef X
#ifdef COLLISIONPREDICTOR
  real* vx = Vx_half->field_cpu;
#else
  real* vx = Vx->field_cpu;
#endif
#endif
#ifdef Y
#ifdef COLLISIONPREDICTOR
  real* vy = Vy_half->field_cpu;
#else
  real* vy = Vy->field_cpu;
#endif
#endif
#ifdef Z
#ifdef COLLISIONPREDICTOR
  real* vz = Vz_half->field_cpu;
#else
  real* vz = Vz->field_cpu;
#endif
#endif
#ifdef X
  real* tauxx = Mmx->field_cpu;
#endif
#ifdef Y
  real* tauyy = Mmy->field_cpu;
#endif
#ifdef Z
  real* tauzz = Mmz->field_cpu;
#endif
#if defined(X) && defined(Z)
  real* tauxz = Mpx->field_cpu;
#endif
#if defined(Y) && defined(X)
  real* tauyx = Mpy->field_cpu;
#endif
#if defined(Z) && defined(Y)
  real* tauzy = Mpz->field_cpu;
#endif
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = XIP; 
  int size_y = Ny+2*NGHY-1;
  int size_z = Nz+2*NGHZ-1;
  real dx = Dx;
  real ss = SIGMASLOPE;
  real r0 = R0;
  real* viscosity_profile = Viscosity_profile->field_cpu;
  int pitch2d = Pitch2D;
//<\EXTERNAL>

//<INTERNAL>
  int i;
  int j;
  int k;
  real div_v;
  real viscosity;
  real viscositym;
  real viscosityzm;
  real viscosityzmym;
  real bigrad;
//<\INTERNAL>

//<CONSTANT>
// real NU(1);
// real GAMMA(1);
// real ALPHA(1);
// real Sxj(Ny+2*NGHY);
// real Syj(Ny+2*NGHY);
// real Szj(Ny+2*NGHY);
// real Sxk(Nz+2*NGHZ);
// real Syk(Nz+2*NGHZ);
// real Szk(Nz+2*NGHZ);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
// real InvVj(Ny+2*NGHY);
//<\CONSTANT>

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
#ifdef ALPHAVISCOSITY
#ifdef ISOTHERMAL
	viscosity     = ALPHA*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j)/(G*MSTAR));
	viscositym    = ALPHA*0.0625*(energy[l]+energy[lxm]+energy[lym]+energy[lxm-pitch])*(energy[l]+energy[lxm]+energy[lym]+energy[lxm-pitch])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
	viscosityzmym = ALPHA*0.0625*(energy[l]+energy[lym]+energy[lzm]+energy[lym-stride])*(energy[l]+energy[lym]+energy[lzm]+energy[lym-stride])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
	viscosityzm   = ALPHA*0.0625*( energy[l]+energy[lzm]+energy[lxm]+energy[lxm-stride] )*( energy[l]+energy[lzm]+energy[lxm]+energy[lxm-stride] )*sqrt(ymed(j)*ymed(j)*ymed(j)/(G*MSTAR));
#else
	viscosity     = ALPHA*GAMMA*(GAMMA-1.0)*energy[l]/rho[l]*sqrt(ymed(j)*ymed(j)*ymed(j)/(G*MSTAR));
	viscositym    = ALPHA*GAMMA*(GAMMA-1.0)*(energy[l]+energy[lxm]+energy[lym]+energy[lxm-pitch])/(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
	viscosityzmym = ALPHA*GAMMA*(GAMMA-1.0)*(energy[l]+energy[lym]+energy[lzm]+energy[lym-stride])/(rho[l]+rho[lym]+rho[lzm]+rho[lym-stride])*sqrt(ymin(j)*ymin(j)*ymin(j)/(G*MSTAR));
	viscosityzm   = ALPHA*GAMMA*(GAMMA-1.0)*(energy[l]+energy[lzm]+energy[lxm]+energy[lxm-stride])/(rho[l]+rho[lzm]+rho[lxm]+rho[lxm-stride])*sqrt(ymed(j)*ymed(j)*ymed(j)/(G*MSTAR));
#endif
#else

	int l2D_local = j + pitch2d * k;

	viscosity     = viscosity_profile[l2D_local];
	viscositym    = viscosity_profile[l2D_local - 1];
	viscosityzm   = viscosity_profile[l2D_local - pitch2d];
	viscosityzmym = viscosity_profile[l2D_local - 1 - pitch2d];





#endif


//Evaluate centered divergence.
	div_v = 0.0;
#ifdef X
	div_v += (vx[lxp]-vx[l])*SurfX(j,k);
#endif
#ifdef Y
	div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
#endif
#ifdef Z
	div_v += (vz[lzp]*SurfZ(j,k+1)-vz[l]*SurfZ(j,k));
#endif
	div_v *= 2.0/3.0*InvVol(j,k);

#if defined(X)
	tauxx[l] = viscosity*rho[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
#endif
#if defined(Y) && defined(X)
	tauxx[l] += viscosity*rho[l]*(vy[lyp]+vy[l])/ymed(j);
#endif
#ifdef Y
	tauyy[l] = viscosity*rho[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
#endif
#ifdef Z
	tauzz[l] = viscosity*rho[l]*(2.0*(vz[lzp]-vz[l])/(zmin(k+1)-zmin(k)) - div_v);
#endif

#if defined(X) && defined(Z)
	tauxz[l] = viscosityzm*.25*(rho[l]+rho[lzm]+rho[lxm]+rho[lxm-stride])*((vx[l]-vx[lzm])/(zmed(k)-zmed(k-1)) + (vz[l]-vz[lxm])/zone_size_x(j,k)); //centered on lower, left "radial" edge in y
#endif

#if defined(Y) && defined(X)
	tauyx[l] = viscositym*.25*(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z
#endif

#if defined(Z) && defined(Y)
	tauzy[l] = viscosityzmym*.25*(rho[l]+rho[lym]+rho[lzm]+rho[lym-stride])*((vz[l]-vz[lym])/(ymed(j)-ymed(j-1)) + (vy[l]-vy[lzm])/(zmed(k)-zmed(k-1))); //centered on lower, inner edge in x ("azimuthal")
#endif
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
