
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
  
  real rdnum,rho_tmp,v_tmp, buffer;
  int idum;
  

  real rhog, rhod, vert, rhog_c, rhod_c, bump, dbump, bump_amp, dgmid, FR, bcoeff, Hd, Heps, Hg; 
  real omega, omega_kep, omega_gas, vphi_gas, vphi_dust, vrad_gas, vrad_dust, vtheta_gas, vtheta_dust, vR_gas, vR_dust;
  real r, R3, R, z, DeltaR, DeltaRin, DeltaRout, TaperRin, TaperRout, taper, sharp;
  real eta, St, StPrime, tstop, cs0, cs, rhog0, H0, omega_kep0, dgratio, gas_energy, bet,delta,alpha;


  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;


  sharp     = 16.0; 
  DeltaRin  = YMIN - R0;
  DeltaRout = YMAX - R0;
  TaperRin  = DUSTTAPER*DeltaRin;
  TaperRout = DUSTTAPER*DeltaRout; 

  smallp    = SIGMASLOPE + FLARINGINDEX + 1.0; 
  smallq    = 1.0 - 2.0*FLARINGINDEX; 
  
  H0         = ASPECTRATIO*R0;
  omega_kep0 = sqrt(G*MSTAR/R0/R0/R0);
  cs0        = omega_kep0*H0;
  rhog0      = SIGMA0*(1.0 + BUMP_AMP)/sqrt(2.0*M_PI)/H0;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {


//setup below is for spherical coords only (taken from Nelson et al 2013, with p -> -p and q -> -q from their notation)


      //first do densities and vphi 
      r = Ymed(j);
      R = r*sin(Zmed(k));
      H = ASPECTRATIO*pow(R/R0,FLARINGINDEX)*R;

      R3 = R*R*R;
      omega_kep = sqrt(G*MSTAR/R3);   //  cs2 = pow((sqrt(G*MSTAR/pow(Ymed(j)*sin(Zmed(k)),3.0))*ASPECTRATIO*pow(Ymed(j)*sin(Zmed(k))/R0,FLARINGINDEX)*Ymed(j)*sin(Zmed(k))),2.0)
      cs    = omega_kep*H;
      cs2   = pow(cs, 2.0);           //  cs2 = pow((sqrt(G*MSTAR/pow(Ymed(j)*sin(Zmed(k)),3.0))*ASPECTRATIO*pow(Ymed(j)*sin(Zmed(k))/R0,FLARINGINDEX)*Ymed(j)*sin(Zmed(k))),2.0)
      
      z  = r*cos(Zmed(k));
      bump_amp = BUMP_AMP;//*exp(-0.5*z*z/H0/H0);  

      bump  = 1.0 + bump_amp*exp(-0.5*pow((R-BUMP_RAD)/(BUMP_WIDTH*H0),2.0));
         dbump = -(R-BUMP_RAD)/pow(BUMP_WIDTH*H0,2.0)*(bump - 1.0);
         dbump*=  R/bump;

         smallp_eff = smallp  - dbump;

         rhog = SIGMA0*pow(R/R0,-SIGMASLOPE)/sqrt(2.0*M_PI)/H; //midplane gas density
         rhog *= bump; 
         vert = (G*MSTAR/cs2)*(1.0/r - 1.0/R);
         rhog*= exp(vert);



       //set non-uniform dust-to-gas ratio according to Chen & Lin (2018) for steady state (assume epsilon < 1 here, and "xi" parameter = 1 (locally isothermal) 
       //now also accoutning for surface density (or midplane density)  bump (CL18 assumed power law pressure profile) 
       //however, this modifiation assumes dP/dr \neq 0 everywhere, which can't be the case if we want a pressure bump. 
       //so this modis commented out 
	
	dgratio=0.;
       if(EPSILON > 0.0){  
       FR      =( EPSILON/pow(1.0+EPSILON,2.0) )*pow(R/R0, smallq/2.0);
//       FR     *=(smallq + smallp);
//       FR     /=(smallq + smallp_eff);   
       bcoeff  = 2.0 - 1.0/FR; 
       dgmid   = -bcoeff - sqrt(bcoeff*bcoeff - 4.0);
       dgmid  /= 2.0; 
       dgratio = dgmid;
				}

  
       if(FLARINGINDEX == 0.5){  //THIS CORRESPONDS TO q=0
       dgratio=EPSILON;
       }


	#ifdef DSLOPE
	dgratio=EPSILON + 0.05*tanh(DUSTSLOPE*(R-R0));
	#endif



	// Introduce initial vertical dust to gas ratio profile
	Hg=ASPECTRATIO;
	Hd = HDHG*Hg;
	Heps = pow(1./pow(Hd,2.0) - 1./pow(Hg,2.0),-0.5);
	dgratio *= exp(-pow(z,2.0)/(2. *pow(Heps,2.0)));

      DeltaR = R - R0;
       if(DeltaR <= TaperRin){
        taper=pow( (DeltaR-DeltaRin)/(TaperRin - DeltaRin), sharp);
       } else if(DeltaR >= TaperRout){
        taper=pow( (DeltaR-DeltaRout)/(TaperRout - DeltaRout), sharp);
       } else {
        taper = 1.0;
       }
	#ifdef DTAPER
        dgratio*= taper;
	#endif


      rhog_c = rhog;

 

       rhod = dgratio*rhog; 
       rhod_c = rhod; 

//printf("%2.6f", rhod_c);


        //now do vphi 

#ifdef FIXPARTICLESIZE
       tstop   = (STOKES1/omega_kep0)*cs0*rhog0/cs/rhog;
#endif
#ifdef STOKESNUMBER
          tstop  = STOKES1/omega_kep;
#endif
       St      = tstop*omega_kep;
       StPrime = St/(1.0 + dgratio);

        eta       = 0.5*( (smallp_eff+smallq)*pow(H/R,2.0)  + smallq - smallq*R/r);
        omega_gas = omega_kep*sqrt(1.0 - 2.0*eta);
        vphi_gas  = R*(omega_gas - OMEGAFRAME);
        vphi_gas += eta*R*omega_kep*( (dgratio/(1.0+dgratio))/(StPrime*StPrime + 1.0) );   

        omega     = sqrt(G*MSTAR/r/r/r);   
        vphi_dust = R*omega_kep - R*OMEGAFRAME; 
        vphi_dust-= eta*R*omega_kep*( (1.0/(1.0+dgratio))/(StPrime*StPrime + 1.0) );

        //now do vrad 

          vR_gas  = 2.0*StPrime/(StPrime*StPrime+1.0);
          vR_gas *= dgratio/(1.0+dgratio);
          vR_gas *= eta*R*omega_kep; 
          vrad_gas= vR_gas*sin(Zmed(k));


          vR_dust = -2.0*StPrime/(StPrime*StPrime+1.0); 
          vR_dust*=  1.0/(1.0+dgratio);
          vR_dust*=  eta*R*omega_kep; 
          vrad_dust= vR_dust*sin(Zmed(k));


          vtheta_gas= vR_gas*cos(Zmin(k));

          vtheta_dust= vR_dust*cos(Zmed(k));
         // vtheta_dust = -sin(Zmed(k))*bet*z*omega_kep;




#ifdef ISOTHERMAL
        gas_energy = cs;
#else
        gas_energy = rhog_c*cs2/(GAMMA-1.0);
#endif



      for (i=0; i<Nx; i++) {

        if (Fluidtype == GAS) {

        v1[l] = vphi_gas;
        v2[l] = vrad_gas;
        v3[l] = vtheta_gas;

        rho[l] = rhog_c;

          // ADD WHITE NOISE
          idum=	 rand() % 41; 
          idum -= 20;	  
          rdnum= (double) idum;
	  rdnum /= (double) 40;

//          rho[l] += 0.01*rdnum*rho[l] ;    //******************************************MODIFICATION************************************************  

	e[l] = gas_energy;

        }

       if (Fluidtype == DUST) {
          rho[l]  = rhod_c;

          // ADD WHITE NOISE
          idum=	 rand() % 41; 
          idum -= 20;	  
          rdnum= (double) idum;
	  rdnum /= (double) 40;

 //         rho[l] += 0.01*rdnum*rho[l] ;    //******************************************MODIFICATION************************************************   

          v1[l]  = vphi_dust; 
          v2[l]  = vrad_dust;
          v3[l]  = vtheta_dust;  
          e[l]   = 0.0;
        }

      }
    }
  }




        if (Fluidtype == GAS) {


 
  buffer = (YMAX-YMIN)/14.;
  
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
     


      for (i=NGHX; i<Nx+NGHX; i++) {
//	rho_tmp =(drand48()-.5)*2.*5.0e-5 ;


 r = Ymed(j);
      R = r*sin(Zmed(k));
      H = ASPECTRATIO*pow(R/R0,FLARINGINDEX)*R;

      R3 = R*R*R;
      omega_kep = sqrt(G*MSTAR/R3);
      cs    = omega_kep*H;

      

	if ((j >= 0) && (j < Ny+2*NGHY) && (k >= 0) && (k < Nz+2*NGHZ)) {
	  v1[l] += (drand48()-.5)*2.*NOISE_AMP*cs;
          v3[l] += (drand48()-.5)*2.*NOISE_AMP*cs;
	
	}
      }
    }
  }



} //fluid type gas


}


void CondInit() {

  int id_gas = 0;
  int feedback;


#ifdef FEEDBACK
        feedback = YES;
masterprint("FEEDBACK INCLUDED\n");
#else
 //       feedback = NO;
        feedback = 0;
masterprint("FEEDBACK NOT INCLUDED!\n");
#endif


   Fluids[id_gas] = CreateFluid("gas",GAS);
   SelectFluid(id_gas);
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

  ColRate(InvStokes1, id_gas, 1, feedback);

}
