#include "fargo3d.h"

void Init() {
 int i,j,k;
  int gj,gk;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  real H, cs2, smallp, smallp_eff, smallq;

  real rdnum,rho_tmp,v_tmp, buffer, cutoff;
  int idum;


  real rhog, rhod, bump, dbump, dgmid, FR, bcoeff,HD, HDHG, fg, fd;
  real omega, omega_kep, omega_gas, vphi_gas, vphi_dust, vrad_gas, vrad_dust, vz_gas, vz_dust, vR_gas, vR_dust;
  real r, R3, R, z, DeltaR, DeltaRin, DeltaRout, TaperRin, TaperRout, taper, sharp, vert, Epsilon;
  real eta, St, StPrime, tstop, cs0, cs, rhog0, H0, omega_kep0, dgratio,dgratio_p, gas_energy, bet,delta,alpha,wd;
  real  vR_gas_0,St0,eps0,eta0,a,b,c;


  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
#ifdef Z
  v3  = Vz->field_cpu;
#else
  // v3 is unused in 2D, but good practice to nullify
  v3  = NULL;
#endif


if(Fluidtype == GAS) {
masterprint("--------------------------------------------------------------\n");

masterprint("VERTICALLY UNSTRATIFIED SETUP!\n");

if ((NFLUIDS >= 2)) {

masterprint("DUST + GAS! \n");

#ifdef DUSTVG
masterprint("DUST FEELS VERTICAL GRAVITY!\n");
#endif

#ifndef PRESET
#ifndef EQUILIBRIUM
masterprint("INITIALIZING UNSETTLED DUST LAYER!\n");
#endif
#endif


#ifdef PRESET
masterprint("INITIALIZING PRE-SETTLED DUST LAYER!\n");
#endif

#ifdef EQUILIBRIUM
masterprint("INITIALIZING (PRE-SETTLED) DUST LAYER in EQUILIBRIUM!\n");
#endif

#ifdef FEEDBACK
masterprint("FEEDBACK INCLUDED!\n");
#endif

#ifndef FEEDBACK
masterprint("FEEDBACK NOT INCLUDED!\n");
#endif

}
if ((NFLUIDS == 1)) {
masterprint("GAS ONLY! \n");
}

masterprint("DEBUG: NX = %d, NZ = %d\n", NX, NZ);

if ((NX >= 2) && (NZ >= 2)) {
    masterprint("SIMULATION IS 3D! \n");
} else if ((NX >= 2) && (NZ == 1)) {
    masterprint("SIMULATION IS 2D NON-AXISYMMETRIC! \n");
} else if ((NZ >= 2) && (NX == 1)) {
    masterprint("SIMULATION IS 2D AXISYMMETRIC! \n");
}
masterprint("--------------------------------------------------------------\n");
}

#ifdef Z
masterprint("VERTICAL MOTIONS INCLUDED!\n");
#endif

  sharp     = 16.0;
  DeltaRin  = YMIN - R0;
  DeltaRout = YMAX - R0;
  TaperRin  = DUSTTAPER*DeltaRin;
  TaperRout = DUSTTAPER*DeltaRout;


//this setup assumes p=0 and q=0 and a constant stokes number!
  smallp    = SIGMASLOPE + FLARINGINDEX + 1.0;
  smallq    = 1.0 - 2.0*FLARINGINDEX;

  H0         = ASPECTRATIO*R0;
  omega_kep0 = sqrt(G*MSTAR/R0/R0/R0);
  cs0        = omega_kep0*H0;
  rhog0      = SIGMA0*(1.0 + BUMP_AMP)/sqrt(2.0*M_PI)/H0;
#ifdef Z
masterprint("PERFORMING Z-LOOP!\n");
    for (k=0; k<Nz+2*NGHZ; k++) {
#endif
    for (j=0; j<Ny+2*NGHY; j++) {

#ifndef Z
      // In 2D runs, explicitly set k=0 for the 'l' index macro.
      k = 0;
#endif

//setup below is for spherical coords only (taken from Nelson et al 2013, with p -> -p and q -> -q from their notation)

      //first do densities and vphi
      R = Ymed(j);
      H = ASPECTRATIO*pow(R/R0,FLARINGINDEX)*R;

      R3 = R*R*R;
      omega_kep = sqrt(G*MSTAR/R3);   //  cs2 = pow((sqrt(G*MSTAR/pow(Ymed(j)*sin(Zmed(k)),3.0))*ASPECTRATIO*pow(Ymed(j)*sin(Zmed(k))/R0,FLARINGINDEX)*Ymed(j)*sin(Zmed(k))),2.0)
      cs    = omega_kep*H;
      cs2   = pow(cs, 2.0);           //  cs2 = pow((sqrt(G*MSTAR/pow(Ymed(j)*sin(Zmed(k)),3.0))*ASPECTRATIO*pow(Ymed(j)*sin(Zmed(k))/R0,FLARINGINDEX)*Ymed(j)*sin(Zmed(k))),2.0)


      bump  = 1.0 + BUMP_AMP*exp(-0.5*pow((R-BUMP_RAD)/(BUMP_WIDTH*H0),2.0));
      dbump = -(R-BUMP_RAD)/pow(BUMP_WIDTH*H0,2.0)*(bump - 1.0);
      dbump*=  R/bump;

      smallp_eff = smallp  - dbump; //smallp  \equiv -dln rhog_mid/dln r

//NEW
Epsilon=METALLICITY; //For a fully unsettled layer

#ifdef DUSTVG
#ifdef EQUILIBRIUM
      //delta=NUDZ*pow(R/H,2.0);
      delta=NUDZ/(pow(H,2.0)*omega_kep);

    //  Epsilon=METALLICITY*(ZMAX-ZMIN)/H*sqrt(STOKES1/delta)/sqrt(2.0*M_PI);
      Epsilon=METALLICITY*sqrt(STOKES1/delta)*pow(R/R0,-SIGMASLOPE);
#endif
#ifdef PRESET
        eta       = 0.5*(smallp_eff+smallq)*pow(H/R,2.0);


	HD=0.5*eta*R; //Same as Li & Youdin (2021)
	HDHG=HD/H;

	HDHG*=1.0; //Check influence of initial layer thickness (Stockholm)
if ((HDHGFAC < 1.0)) {
	HDHG=HDHGFAC; //OVERRIDE
}

      Epsilon=METALLICITY*(ZMAX-ZMIN)/H/HDHG/sqrt(2.0*M_PI);
#endif
#endif

      dgratio=Epsilon;


//if ((BUMP_AMP == 0.)) {
//now implement steady state solution as in Auddy et al 2022;

//	St0=STOKES1;
//	eps0=Epsilon;
//	eta0= 0.5*(smallp_eff+smallq)*pow(ASPECTRATIO,2.0);

//	vR_gas_0= 2.0*St0*eps0/(eps0*eps0+2.0*eps0+1.0+St0*St0)*eta0*omega_kep0*R0;

//	a=vR_gas_0*pow(R/R0,SIGMASLOPE-1.0);
//	b=2.0*(a-St0*eta0*omega_kep0*R0);
//	c=(1.0+St0*St0)*a;

//       if(Epsilon <= 1.){
//	dgratio=-b/(2.0*a) - 1.0/(2.0*a)*sqrt(b*b-4.0*a*c);
//       }  else {
//	dgratio=-b/(2.0*a) + 1.0/(2.0*a)*sqrt(b*b-4.0*a*c);
//       }
//
//}

        // taper off dust density at boundaries to accomodate stockholm (INNER BOUNDARY ONLY)
      DeltaR = R - R0;
       if(DeltaR <= TaperRin){
        taper=pow( (DeltaR-DeltaRin)/(TaperRin - DeltaRin), sharp);
       } else if(DeltaR >= TaperRout){
        //taper=pow( (DeltaR-DeltaRout)/(TaperRout - DeltaRout), sharp);
       } else {
        taper = 1.0;
       }

	#ifdef DTAPER
        dgratio*= taper;
	#endif


	//#ifdef DSLOPE
	//dgratio=Epsilon + 0.05*tanh(DUSTSLOPE*(R-R0));
	//#endif


	//Introduce STRATIFIED dust layer if vertical gravity for dust is activated
	z=0.;
	vert=1.0;
#ifdef DUSTVG
        z  = Zmed(k);
#ifdef EQUILIBRIUM
	vert=exp(-0.5*(STOKES1/delta)*pow(z/H,2.0));
#endif
#ifdef PRESET
	vert=exp(-0.5*pow(z/H/HDHG,2.0));
#endif
       dgratio*= vert;
#endif

#ifdef Z
      rhog = SIGMA0*pow(R/R0,-SIGMASLOPE)/sqrt(2.0*M_PI)/H; //midplane gas density
#else
      rhog = SIGMA0*pow(R/R0,-SIGMASLOPE); //RAZOR THIN DISK
#endif
      rhog*= bump; //gas surface density bump





if ((NFLUIDS >= 2)) {
       rhod = dgratio*rhog;

}


        //now do vphi

#ifdef FIXPARTICLESIZE
       tstop   = (STOKES1/omega_kep0)*cs0*rhog0/cs/rhog;
#endif
#ifdef STOKESNUMBER
          tstop  = STOKES1/omega_kep;
#endif
       St      = tstop*omega_kep;



//##################################################################################

	fg=1.0/(1.0+dgratio);  //gas fraction
	fd=dgratio/(1.0+dgratio); //dust fraction

       StPrime = St*fg;
#ifdef Z
        eta = 0.5*(smallp_eff+smallq)*pow(H/R,2.0);
#else
        eta = 0.5*(SIGMASLOPE + smallq)*pow(H/R,2.0);
#endif
        omega_gas = omega_kep*sqrt(1.0 - 2.0*eta);
        vphi_gas  = R*(omega_gas - OMEGAFRAME);
        vphi_dust = R*(omega_kep - OMEGAFRAME);

#ifdef FEEDBACK

       // omega_gas = omega_kep*sqrt(1.0 - 2.0*fg*eta);

//	vphi_gas  = R*(omega_gas - OMEGAFRAME) - fd*pow(fg,2.0)*eta*R*omega_kep*pow(St,2.0)/(1.0+pow(StPrime,2.0));
//	vphi_dust = R*(omega_gas - OMEGAFRAME)    + pow(fg,3.0)*eta*R*omega_kep*pow(St,2.0)/(1.0+pow(StPrime,2.0)) ;


#endif

        //now do vrad

        vrad_gas = 0.;
        vrad_dust = 0.;


#ifdef FEEDBACK

    //    vrad_gas=         2.0*fd*fg*eta*R*omega_kep*St/(1.0+pow(StPrime,2.0));
    //    vrad_dust= -2.0*pow(fg,2.0)*eta*R*omega_kep*St/(1.0+pow(StPrime,2.0));
#endif

//#######################################################################################


        //now do vz

          vz_gas=0.;
          vz_dust=0.;

#ifdef EQUILIBRIUM
#ifdef DUSTVG
          vz_dust = -STOKES1*z*omega_kep;
#endif
#endif




#ifdef ISOTHERMAL
        gas_energy = cs;
#else
        gas_energy = rhog*cs2/(GAMMA-1.0);
#endif




      for (i=0; i<Nx; i++) {

        if (Fluidtype == GAS) {

        v1[l] = vphi_gas;
        v2[l] = vrad_gas;
#ifdef Z
        v3[l] = vz_gas;
#endif
        rho[l] = rhog;
	e[l] = gas_energy;

        }

       if (Fluidtype == DUST) {
          rho[l]  = rhod;

          v1[l]  = vphi_dust;
          v2[l]  = vrad_dust;
          v3[l]  = vz_dust;
          e[l]   = 0.0;
        }

      }
    }
/********************************************************************************/
/* FIX 1                                     */
/* The closing brace for the k-loop is now also wrapped in an #ifdef Z block.   */
/* This prevents it from prematurely closing Init() in a 2D compile.          */
/********************************************************************************/
#ifdef Z
  }
#endif

 //       if (Fluidtype == GAS) {

  buffer = (YMAX-YMIN)/3.;

  srand48 (ArrayNb); //Same seed for all processes

  /* The following initialization loop is a bit tricky. Instead of the
     standard procedure that consists in each PE initializing its own
     patch (only), here each PE performs a global loop over the whole
     mesh (hence the indices gj, gk). The aim is that each of them
     follows the exact same sequence of random numbers, so that the
     noise initialized in the density field be independent on the
     number of processors (which allows to test that the code output
     is independent of the number of processors). The values are
     stored in _tmp variables, and used only if they fall in the
     current PE.*/

  for (gk=0; gk<NZ+2*NGHZ; gk++) {
    for (gj=0; gj<NY+2*NGHY; gj++) {
      j = gj-y0cell;
      k = gk-z0cell;
      if ((j >= 0) && (j < Ny+2*NGHY))
	r = Ymed(j);
      else
	r = Ymed(0);

      R = r;
      H = ASPECTRATIO*pow(R/R0,FLARINGINDEX)*R;

      R3 = R*R*R;
      omega_kep = sqrt(G*MSTAR/R3);
      cs    = omega_kep*H;

      cutoff = 1.0;
      if (r < YMIN+buffer)
	cutoff = pow((r-YMIN)/buffer,10.);
      if (r > YMAX-buffer)
	cutoff = pow((YMAX-r)/buffer,10.);
      for (i=NGHX; i<Nx+NGHX; i++) {
//	rho_tmp =(drand48()-.5)*2.*5.0e-5 ;

 //r = Ymed(j);

	if ((j >= 0) && (j < Ny+2*NGHY) && (k >= 0) && (k < Nz+2*NGHZ)) {
	  v1[l] += (drand48()-.5)*2.*NOISE_AMP*cs*cutoff;
	  v2[l] += (drand48()-.5)*2.*NOISE_AMP*cs*cutoff;
          #ifdef Z
          v3[l] += (drand48()-.5)*2.*NOISE_AMP*cs*cutoff;
          #endif
	}
      }
    }
  }
//} //fluid type gas
}

/********************************************************************************/
/* FIX 2                                     */
/* This is the proper end of the Init() function. The stray braces from the   */
/* original file have been removed.                                           */
/********************************************************************************/

void CondInit() {

  int id_gas = 0;
  int feedback;


#ifdef FEEDBACK
        feedback = 1;
//	printf("setting feedback=YES!\n");

#else
        feedback = 0;


//	printf("setting feedback=0!\n");

#endif


   Fluids[id_gas] = CreateFluid("gas",GAS);
   SelectFluid(id_gas);
//	printf("EXECUTING INIT!\n");
   Init();


  char dust_name[MAXNAMELENGTH];
  int id_dust;
  real InvStokes1, rhog0, H0;

   for(id_dust = 1; id_dust<NFLUIDS; id_dust++) {
    sprintf(dust_name,"dust%d",id_dust); //We assign different names to the dust fluids

    Fluids[id_dust]  = CreateFluid(dust_name, DUST);
    SelectFluid(id_dust);
    Init();

  }


#ifdef FIXPARTICLESIZE
  H0         = ASPECTRATIO*R0;
  rhog0      = SIGMA0/sqrt(2.0*M_PI)/H0;
  InvStokes1 = 1.0/(STOKES1*H0*rhog0);
#endif

#ifdef STOKESNUMBER
  InvStokes1 = 1.0/STOKES1;
#endif

masterprint("EXECUTING COLRATE!\n");
  ColRate(InvStokes1, id_gas, 1, feedback);

}
