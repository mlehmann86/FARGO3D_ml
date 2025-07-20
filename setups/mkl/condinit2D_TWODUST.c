#include "fargo3d.h"

void Init() {
  int i,j,k;
  real *v1;
  real *v2;
  real *v3;
  real *e;
  real *rho;
  real H, cs2, smallp, smallp_eff, smallq;
  
  real rdnum;
  int idum;


  real rhog, rhod, vert, rhog_c, rhod_c, bump, dbump, bump_amp, dgmid, FR, bcoeff; 
  real omega, omega_kep, omega_gas, vphi_gas, vphi_dust, vrad_gas, vrad_dust, vtheta_gas, vtheta_dust, vR_gas, vR_dust;
  real r, R3, R, z, DeltaR, DeltaRin, DeltaRout, TaperRin, TaperRout, taper, sharp;
  real eta, St, StPrime, tstop, cs0, cs, rhog0, H0, omega_kep0, dgratio, gas_energy;

  rho = Density->field_cpu;
  e   = Energy->field_cpu;
  v1  = Vx->field_cpu;
  v2  = Vy->field_cpu;
  v3  = Vz->field_cpu;



  sharp     = 2.0; 
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
      omega_kep = sqrt(G*MSTAR/R3);
      cs    = omega_kep*H;
      cs2   = pow(cs, 2.0); 
      
      z  = r*cos(Zmed(k));
      bump_amp = BUMP_AMP;//*exp(-0.5*z*z/H0/H0);  

      bump  = 1.0 + bump_amp*exp(-0.5*pow((R-BUMP_RAD)/(BUMP_WIDTH*H0),2.0));
      dbump = -(R-BUMP_RAD)/pow(BUMP_WIDTH*H0,2.0)*(bump - 1.0);
      dbump*=  R/bump;

      smallp_eff = smallp  - dbump; //smallp  \equiv -dln rhog_mid/dln r

      rhog = SIGMA0*pow(R/R0,-SIGMASLOPE)/sqrt(2.0*M_PI)/H; //midplane gas density 
      rhog*= bump; //gas surface density bump 
      vert = (G*MSTAR/cs2)*(1.0/r - 1.0/R); 
      rhog*= exp(vert); 
      rhog_c = rhog;

       DeltaR = R - R0;
       if(DeltaR <= TaperRin){
        taper=pow( (DeltaR-DeltaRin)/(TaperRin - DeltaRin), sharp);
       } else if(DeltaR >= TaperRout){
        taper=pow( (DeltaR-DeltaRout)/(TaperRout - DeltaRout), sharp);
       } else {
        taper = 1.0;
       }
       //set non-uniform dust-to-gas ratio according to Chen & Lin (2018) for steady state (assume epsilon < 1 here, and "xi" parameter = 1 (locally isothermal) 
       //now also accoutning for surface density (or midplane density)  bump (CL18 assumed power law pressure profile) 
       //however, this modifiation assumes dP/dr \neq 0 everywhere, which can't be the case if we want a pressure bump. 
       //so this modis commented out 
       FR      =( EPSILON/pow(1.0+EPSILON,2.0) )*pow(R/R0, smallq/2.0);
//       FR     *=(smallq + smallp);
//       FR     /=(smallq + smallp_eff);   
       bcoeff  = 2.0 - 1.0/FR; 
       dgmid   = -bcoeff - sqrt(bcoeff*bcoeff - 4.0);
       dgmid  /= 2.0; 
       dgratio = dgmid;//*taper;
       rhod = dgratio*rhog; 
       rhod_c = rhod; 




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


#ifdef ISOTHERMAL
        gas_energy = cs;
#else
        gas_energy = rhog_c*cs2/(GAMMA-1.0);
#endif


        //now do vrad 
         r = Ymin(j);
         R = r*sin(Zmed(k));
         H = ASPECTRATIO*pow(R/R0,FLARINGINDEX)*R;

         R3 = R*R*R;
         omega_kep = sqrt(G*MSTAR/R3);
         cs    = omega_kep*H;
         cs2   = pow(cs, 2.0);


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

         DeltaR = R - R0;
         if(DeltaR <= TaperRin){
         taper=pow( (DeltaR-DeltaRin)/(TaperRin - DeltaRin), sharp);
         } else if(DeltaR >= TaperRout){
         taper=pow( (DeltaR-DeltaRout)/(TaperRout - DeltaRout), sharp);
         } else {
         taper = 1.0;
         }
         //set non-uniform dust-to-gas ratio according to Chen & Lin (2018) for steady state (assume epsilon < 1 here, and "xi" parameter = 1 (locally isothermal)
         FR      =( EPSILON/pow(1.0+EPSILON,2.0) )*pow(R/R0, smallq/2.0);
//     	   FR     *=(smallq + smallp);
//         FR     /=(smallq + smallp_eff); 
         bcoeff  = 2.0 - 1.0/FR;
         dgmid   = -bcoeff - sqrt(bcoeff*bcoeff - 4.0);
         dgmid  /= 2.0;
         dgratio = dgmid;//*taper;
 
#ifdef FIXPARTICLESIZE
          tstop   = (STOKES1/omega_kep0)*cs0*rhog0/cs/rhog;
#endif
#ifdef STOKESNUMBER
          tstop  = STOKES1/omega_kep;
#endif
          St      = tstop*omega_kep;
          StPrime = St/(1.0 + dgratio);

          eta       = 0.5*( (smallp_eff+smallq)*pow(H/R,2.0)  + smallq - smallq*R/r);
        
          vR_gas  = 2.0*StPrime/(StPrime*StPrime+1.0);
          vR_gas *= dgratio/(1.0+dgratio);
          vR_gas *= eta*R*omega_kep; 
          vrad_gas= vR_gas*sin(Zmed(k));
//          vrad_gas = -3.0*NU*(1.0 + 2.0*SIGMASLOPE)/(2.0*r); //******************************************MODIFICATION************************************************

	 
 
          vR_dust = -2.0*StPrime/(StPrime*StPrime+1.0); 
          vR_dust*=  1.0/(1.0+dgratio);
          vR_dust*=  eta*R*omega_kep; 
          vrad_dust= vR_dust*sin(Zmed(k));
//          vrad_dust = 0.; //******************************************MODIFICATION************************************************

         //now do vtheta 
         r = Ymed(j);
         R = r*sin(Zmin(k));
         H = ASPECTRATIO*pow(R/R0,FLARINGINDEX)*R;

         R3 = R*R*R;
         omega_kep = sqrt(G*MSTAR/R3);
         cs    = omega_kep*H;
         cs2   = pow(cs, 2.0);

         z  = r*cos(Zmin(k));
         bump_amp = BUMP_AMP;//*exp(-0.5*z*z/H0/H0); 

         bump  = 1.0 + bump_amp*exp(-0.5*pow((R-BUMP_RAD)/(BUMP_WIDTH*H0),2.0));
         dbump = -(R-BUMP_RAD)/pow(BUMP_WIDTH*H0,2.0)*(bump - 1.0);
         dbump*=  R/bump;

         smallp_eff = smallp  - dbump;

         rhog = SIGMA0*pow(R/R0,-SIGMASLOPE)/sqrt(2.0*M_PI)/H; //midplane gas density
         rhog *= bump; 
         vert = (G*MSTAR/cs2)*(1.0/r - 1.0/R);
         rhog*= exp(vert);

         DeltaR = R - R0;
         if(DeltaR <= TaperRin){
         taper=pow( (DeltaR-DeltaRin)/(TaperRin - DeltaRin), sharp);
         } else if(DeltaR >= TaperRout){
         taper=pow( (DeltaR-DeltaRout)/(TaperRout - DeltaRout), sharp);
         } else {
         taper = 1.0;
         }
         //set non-uniform dust-to-gas ratio according to Chen & Lin (2018) for steady state (assume epsilon < 1 here, and "xi" parameter = 1 (locally isothermal)
         FR      =( EPSILON/pow(1.0+EPSILON,2.0) )*pow(R/R0, smallq/2.0);
//         FR     *=(smallq + smallp);
//         FR     /=(smallq + smallp_eff); 
         bcoeff  = 2.0 - 1.0/FR;
         dgmid   = -bcoeff - sqrt(bcoeff*bcoeff - 4.0);
         dgmid  /= 2.0;
         dgratio = dgmid;//*taper;

#ifdef FIXPARTICLESIZE
          tstop   = (STOKES1/omega_kep0)*cs0*rhog0/cs/rhog;
#endif
#ifdef STOKESNUMBER
          tstop   = STOKES1/omega_kep;
#endif
          St      = tstop*omega_kep;
          StPrime = St/(1.0 + dgratio);

          eta       = 0.5*( (smallp_eff+smallq)*pow(H/R,2.0)  + smallq - smallq*R/r);

          vR_gas  = 2.0*StPrime/(StPrime*StPrime+1.0);
          vR_gas *= dgratio/(1.0+dgratio);
          vR_gas *= eta*R*omega_kep;
          vtheta_gas= vR_gas*cos(Zmin(k));
   //       vtheta_gas *= 0.; //******************************************MODIFICATION************************************************

          vR_dust = -2.0*StPrime/(StPrime*StPrime+1.0);
          vR_dust*=  1.0/(1.0+dgratio);
          vR_dust*=  eta*R*omega_kep;
          vtheta_dust= vR_dust*cos(Zmed(k));
   //       vtheta_dust *= 0.; //******************************************MODIFICATION************************************************



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

          rho[l] += 0.01*rdnum*rho[l] ;    //******************************************MODIFICATION************************************************  

	e[l] = gas_energy;

        }

       if (Fluidtype == DUST) {
          rho[l]  = rhod_c;

          // ADD WHITE NOISE
          idum=	 rand() % 41; 
          idum -= 20;	  
          rdnum= (double) idum;
	  rdnum /= (double) 40;

          rho[l] += 0.01*rdnum*rho[l] ;    //******************************************MODIFICATION************************************************   

          v1[l]  = vphi_dust; 
          v2[l]  = vrad_dust;
          v3[l]  = vtheta_dust;  
          e[l]   = 0.0;
        }

      }
    }
  }
}

void CondInit() {

  int id_gas = 0;
  int feedback = YES;

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
