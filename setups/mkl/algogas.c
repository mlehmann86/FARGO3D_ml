#include "fargo3d.h"

TimeProcess t_Comm;
TimeProcess t_Hydro;
TimeProcess t_Mhd;
TimeProcess t_sub1;
TimeProcess t_sub1_x;
TimeProcess t_sub1_y;
TimeProcess t_sub1_z;


// Corrected version of the manual boundary fix function
void Fix_Tguess_Boundaries() {
    real *field = T_guess->field_cpu;
    int i, j, k, idx, l_src; // Renamed 'l' to 'idx'
    int pitch = Pitch_cpu;
    int stride = Stride_cpu;

    if (Gridd.bc_left) {
        for (k = 0; k < Nz + 2 * NGHZ; k++) {
            for (j = 0; j < NGHY; j++) {
                for (i = 0; i < Nx; i++) {
                    idx = i + j * pitch + k * stride;
                    l_src = i + NGHY * pitch + k * stride;
                    field[idx] = field[l_src];
                }
            }
        }
    }

    if (Gridd.bc_right) {
        for (k = 0; k < Nz + 2 * NGHZ; k++) {
            for (j = Ny + NGHY; j < Ny + 2 * NGHY; j++) {
                for (i = 0; i < Nx; i++) {
                    idx = i + j * pitch + k * stride;
                    l_src = i + (Ny + NGHY - 1) * pitch + k * stride;
                    field[idx] = field[l_src];
                }
            }
        }
    }
}



void FillGhosts (int var) {

  InitSpecificTime (&t_Comm, "MPI Communications");
  FARGO_SAFE(comm (var));
  GiveSpecificTime (t_Comm);
  FARGO_SAFE(boundaries()); // Always after a comm.

#if defined(Y)
  if (NY == 1)    /* Y dimension is mute */
    CheckMuteY();
#endif
#if defined(Z)
  if (NZ == 1)    /* Z dimension is mute */
    CheckMuteZ();
#endif

}

void Fill_Resistivity_Profiles () {

  OUTPUT2D(Eta_profile_xi);
  OUTPUT2D(Eta_profile_xizi);
  OUTPUT2D(Eta_profile_zi);

  int j,k;
  if (Resistivity_Profiles_Filled) return;
  real* eta_profile_xi = Eta_profile_xi->field_cpu;
  real* eta_profile_xizi = Eta_profile_xizi->field_cpu;
  real* eta_profile_zi = Eta_profile_zi->field_cpu;

  for (k=0; k<Nz+2*NGHZ; k++) {
    for (j=0; j<Ny+2*NGHY; j++) {
      eta_profile_xi[l2D] = Resistivity (Ymin(j),Zmed(k));
      eta_profile_xizi[l2D] = Resistivity (Ymin(j),Zmin(k));
      eta_profile_zi[l2D] = Resistivity (Ymed(j),Zmin(k));
    }
  }
  Resistivity_Profiles_Filled = YES;
}




void Fill_Viscosity_Profile () {
  masterprint("ALGOGAS_VISC: Entered Fill_Viscosity_Profile.\n");


#ifdef VISCOSITYSLOPE
masterprint("✅ VISCOSITYSLOPE is recognized inside Fill_Viscosity_Profile\n");
#else
masterprint("❌ VISCOSITYSLOPE is NOT recognized inside Fill_Viscosity_Profile\n");
#endif

  masterprint("ALGOGAS_VISC: Calling OUTPUT2D.\n");


  OUTPUT2D(Viscosity_profile);
  masterprint("ALGOGAS_VISC: Returned from OUTPUT2D.\n");


  masterprint("ALGOGAS_VISC: Checking if profile is already filled.\n");

  if (Viscosity_Profile_Filled) return;
  masterprint("ALGOGAS_VISC: Profile is not filled, proceeding.\n");

  masterprint("ALGOGAS_VISC: Getting field_cpu pointer.\n");
  real* viscosity_profile = Viscosity_profile->field_cpu;
  int j, k;

  masterprint("Fill_Viscosity_Profile DEBUG: Ny = %d, NGHY = %d\n", Ny, NGHY);
  masterprint("Fill_Viscosity_Profile DEBUG: Nz = %d, NGHZ = %d\n", Nz, NGHZ);

#ifdef VISCOSITYSLOPE
  int total_size = (Ny + 2 * NGHY) * (Nz + 2 * NGHZ);
  for (k = 0; k < Nz + 2*NGHZ; k++) {
    for (j = 0; j < Ny + 2*NGHY; j++) {
      int index2D = j + (Ny + 2*NGHY) * k;

      if (index2D < 0 || index2D >= total_size) {
        masterprint("❌ ERROR: Invalid index2D = %d at j = %d, k = %d (Ny = %d, NGHY = %d, Nz = %d, NGHZ = %d)\n",
                    index2D, j, k, Ny, NGHY, Nz, NGHZ);
        exit(1);
      }

      real bigrad_safe = fmax(ymed(j), 1e-8);
      real viscosity_raw = NU * pow(bigrad_safe / R0, SIGMASLOPE - 0.5);
      viscosity_profile[index2D] = viscosity_raw;

      if (k == Nz / 2) {
    real r_here = ymed(j);

    if (fabs(r_here - 0.4) < 1e-3 ||
        fabs(r_here - 0.6) < 1e-3 ||
        fabs(r_here - 0.8) < 1e-3 ||
        fabs(r_here - 1.0) < 1e-3 ||
        fabs(r_here - 1.2) < 1e-3 ||
        fabs(r_here - 1.4) < 1e-3 ||
        fabs(r_here - 1.6) < 1e-3) {
        
        masterprint("Fill_Viscosity_Profile DEBUG: j = %d, k = %d, r = %.3f, viscosity = %.8e\n",
                    j, k, r_here, viscosity_raw);
    }
}
    }
  }
  masterprint("Fill_Viscosity_Profile: viscosity profile filled with VISCOSITYSLOPE (power-law).\n");


#else
  int total_size = (Ny + 2 * NGHY) * (Nz + 2 * NGHZ);
  for (k = 0; k < Nz + 2*NGHZ; k++) {
    for (j = 0; j < Ny + 2*NGHY; j++) {
      int index2D = j + (Ny + 2*NGHY) * k;

      if (index2D < 0 || index2D >= total_size) {
        masterprint("❌ ERROR: Invalid index2D = %d at j = %d, k = %d\n", index2D, j, k);
        exit(1);
      }

      viscosity_profile[index2D] = NU;
    }
  }
  masterprint("Fill_Viscosity_Profile: viscosity profile filled with constant NU.\n");
#endif

  Viscosity_Profile_Filled = YES;

  // === Min/Max Check (full array) ===
  real min_visc = 1e30;
  real max_visc = -1e30;

  for (k = 0; k < Nz + 2*NGHZ; k++) {
    for (j = 0; j < Ny + 2*NGHY; j++) {
      real visc = viscosity_profile[l2D];
      if (visc < min_visc) min_visc = visc;
      if (visc > max_visc) max_visc = visc;
    }
  }

  masterprint("Fill_Viscosity_Profile: min viscosity = %.8e, max viscosity = %.8e\n", min_visc, max_visc);

// === Optional: Dump viscosity profile to a file for full inspection ===
FILE *fvisc;
char fname[256];
sprintf(fname, "viscosity_profile_rank%d.dat", CPU_Rank);
fvisc = fopen(fname, "w");

if (fvisc != NULL) {
  for (k = 0; k < Nz + 2*NGHZ; k++) {
    for (j = 0; j < Ny + 2*NGHY; j++) {
      real r_here = ymed(j);
      real visc = viscosity_profile[l2D];
      fprintf(fvisc, "%.8e %.8e\n", r_here, visc);
    }
  }
  fclose(fvisc);
  masterprint("Fill_Viscosity_Profile: full profile written to viscosity_profile_dump.dat\n");
} else {
  masterprint("Fill_Viscosity_Profile: WARNING! Could not open file for writing viscosity profile.\n");
}

}



void Fill_Chi_Profile () {
    masterprint("DEBUG: Entered Fill_Chi_Profile()\n");

    //if (Chi_profile == NULL) {
    //    Chi_profile = CreateField2D("Chi_profile", YES);
    //    masterprint("Chi_profile was NULL, created via CreateField2D()\n");
    //}

    OUTPUT2D(Chi_profile);  // Now safe

    if (Chi_Profile_Filled == YES) return;

    real* chi_profile = Chi_profile->field_cpu;
    int j, k;

    masterprint("Fill_Chi_Profile DEBUG: Ny = %d, NGHY = %d\n", Ny, NGHY);
    masterprint("Fill_Chi_Profile DEBUG: Nz = %d, NGHZ = %d\n", Nz, NGHZ);

#ifdef THERMALDIFFUSIONSLOPE
    real kslope = 1.0 - 2.0 * FLARINGINDEX;
    int total_size = (Ny + 2 * NGHY) * (Nz + 2 * NGHZ);
    for (k = 0; k < Nz + 2*NGHZ; k++) {
        for (j = 0; j < Ny + 2*NGHY; j++) {
            int index2D = j + (Ny + 2*NGHY) * k;

            if (index2D < 0 || index2D >= total_size) {
                masterprint("❌ ERROR: Invalid index2D = %d at j = %d, k = %d\n", index2D, j, k);
                exit(1);
            }

            real bigrad_safe = fmax(ymed(j), 1e-8);

if (k == Nz / 2 && j % 10 == 0) {
    masterprint("Fill DEBUG: j = %d, ymed(j) = %.6f, bigrad_safe = %.8e\n",
        j, ymed(j), bigrad_safe);
}

            real chi_raw = CHI * pow(bigrad_safe / R0, kslope);
            chi_profile[index2D] = chi_raw;

            if (k == Nz / 2) {
    real r_here = ymed(j);

    if (fabs(r_here - 0.4) < 1e-3 ||
        fabs(r_here - 0.6) < 1e-3 ||
        fabs(r_here - 0.8) < 1e-3 ||
        fabs(r_here - 1.0) < 1e-3 ||
        fabs(r_here - 1.2) < 1e-3 ||
        fabs(r_here - 1.4) < 1e-3 ||
        fabs(r_here - 1.6) < 1e-3) {
        
        masterprint("Fill_ThermalDiffusion_Profile DEBUG: j = %d, k = %d, r = %.3f, chi = %.8e\n",
            j, k, r_here, chi_raw);
    }
}
        }
    }
    masterprint("Fill_Chi_Profile: chi profile filled with THERMALDIFFUSIONSLOPE (power-law).\n");

#else
    int total_size = (Ny + 2 * NGHY) * (Nz + 2 * NGHZ);
    for (k = 0; k < Nz + 2*NGHZ; k++) {
        for (j = 0; j < Ny + 2*NGHY; j++) {
            int index2D = j + (Ny + 2*NGHY) * k;

            if (index2D < 0 || index2D >= total_size) {
                masterprint("❌ ERROR: Invalid index2D = %d at j = %d, k = %d\n", index2D, j, k);
                exit(1);
            }

            chi_profile[index2D] = CHI;
        }
    }
    masterprint("Fill_Chi_Profile: chi profile filled with constant CHI.\n");
#endif

    Chi_Profile_Filled = YES;

    // === Optional min/max diagnostic
    real min_chi = 1e30;
    real max_chi = -1e30;

    for (k = 0; k < Nz + 2*NGHZ; k++) {
        for (j = 0; j < Ny + 2*NGHY; j++) {
            int index2D = j + (Ny + 2*NGHY) * k;
            real chi_val = chi_profile[index2D];
            if (chi_val < min_chi) min_chi = chi_val;
            if (chi_val > max_chi) max_chi = chi_val;
        }
    }

    masterprint("Fill_Chi_Profile: min chi = %.8e, max chi = %.8e\n", min_chi, max_chi);
}


void Sources(real dt) {
     

  SetupHook1 (); //Setup specific hook. Defaults to empty function.
  
  //Equations of state-----------------------------------------------------------
#ifdef ADIABATIC
  //      if (Fluidtype == GAS) {
  FARGO_SAFE(ComputePressureFieldAd());
//				}
#endif
#ifdef ISOTHERMAL
  FARGO_SAFE(ComputePressureFieldIso());
#endif
#ifdef POLYTROPIC
  FARGO_SAFE(ComputePressureFieldPoly());
#endif
  //-----------------------------------------------------------------------------
    
  InitSpecificTime (&t_Hydro, "Eulerian Hydro (no transport) algorithms");
  
  // REGARDLESS OF WHETHER WE USE FARGO, Vx IS ALWAYS THE TOTAL VELOCITY IN X
 


#ifdef POTENTIAL
 FARGO_SAFE(compute_potential(dt));

  if (Corotating) {
    FARGO_SAFE(CorrectVtheta(Domega));
}
#endif
  
#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC(Vy));
#endif

#ifdef X
  FARGO_SAFE(SubStep1_x(dt));
#endif    
#ifdef Y
  FARGO_SAFE(SubStep1_y(dt));
#endif  
#ifdef Z
  FARGO_SAFE(SubStep1_z(dt));
#endif
  
#if (defined(VISCOSITY) || defined(ALPHAVISCOSITY))
  if (Fluidtype == GAS) viscosity(dt);
#endif
  
#ifndef NOSUBSTEP2

  FARGO_SAFE(SubStep2_a(dt));
//  if(Fluidtype == GAS) {
    FARGO_SAFE(SubStep2_b(dt));
//  }
#endif

  // NOW: Vx INITIAL X VELOCITY, Vx_temp UPDATED X VELOCITY FROM SOURCE TERMS + ARTIFICIAL VISCOSITY

#ifdef ADIABATIC
#if defined(BETACOOLING) || defined(THERMALDIFFUSION)

#ifdef PREDICTORCORRECTOR
	if(Fluidtype == GAS) {
	  Edamp_predict(dt);  //Calls either _cpu or _gpu, depending on its value
      // ---------------------
	  //Edamp_fillghosts(dt);  //Calls either _cpu or _gpu, depending on its value

      // --- NEW SOLUTION: Hijack the Energy communication channel ---
        Field *backup_energy_ptr = Energy; // 1. Save the real Energy pointer
        Energy = T_guess;                  // 2. Point 'Energy' to our T_guess field

        FARGO_SAFE(FillGhosts(ENERGY));    // 3. Call FillGhosts for ENERGY.
                                           //    The code now communicates T_guess data.

        Energy = backup_energy_ptr;        // 4. Restore the real Energy pointer
     
        // --- MANUALLY FIX THE PHYSICAL BOUNDARIES ---
        Fix_Tguess_Boundaries();


      // ---------------------
	  Edamp_correct(dt);  //Calls either _cpu or _gpu, depending on its value
	}
#else
	if(Fluidtype == GAS) {
	  Edamp(dt);  //Calls either _cpu or _gpu, depending on its value
	}
#endif



#endif
	if(Fluidtype == GAS) {
	  FARGO_SAFE(SubStep3(dt));
	}

#endif
    
  GiveSpecificTime (t_Hydro);
  
#ifdef MHD //-------------------------------------------------------------------
  if(Fluidtype == GAS){
    InitSpecificTime (&t_Mhd, "MHD algorithms");
    FARGO_SAFE(copy_velocities(VTEMP2V));
#ifndef STANDARD // WE USE THE FARGO ALGORITHM
    FARGO_SAFE(ComputeVmed(Vx));
    FARGO_SAFE(ChangeFrame(-1, Vx, VxMed)); //Vx becomes the residual velocity
    VxIsResidual = YES;
#endif
     
    ComputeMHD(dt);

#ifndef STANDARD
    FARGO_SAFE(ChangeFrame(+1, Vx, VxMed)); //Vx becomes the total, updated velocity
    VxIsResidual = NO;
#endif //STANDARD
    FARGO_SAFE(copy_velocities(V2VTEMP));
    // THIS COPIES Vx INTO Vx_temp
    GiveSpecificTime (t_Mhd);
  }
#endif //END MHD----------------------------------------------------------------

  InitSpecificTime (&t_Hydro, "Transport algorithms");

#if ((defined(SHEARINGSHEET2D) || defined(SHEARINGBOX3D)) && !defined(SHEARINGBC))
  FARGO_SAFE(NonReflectingBC (Vy_temp));
#endif
  
  FARGO_SAFE(copy_velocities(VTEMP2V));
  FARGO_SAFE(FillGhosts(PrimitiveVariables()));
  FARGO_SAFE(copy_velocities(V2VTEMP));

#ifdef MHD //-------------------------------------------------------------------
  if(Fluidtype == GAS){ //We do MHD only for the gaseous component
    
    FARGO_SAFE(UpdateMagneticField(dt,1,0,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,1,0));
    FARGO_SAFE(UpdateMagneticField(dt,0,0,1));

#if !defined(STANDARD)
    FARGO_SAFE(MHD_fargo (dt)); // Perform additional field update with uniform velocity
#endif

  } 
#endif //END MHD ---------------------------------------------------------------
}

void Transport(real dt) {

  //NOTE: V_temp IS USED IN TRANSPORT

#ifdef X
#ifndef STANDARD
  FARGO_SAFE(ComputeVmed(Vx_temp)); 
#endif
#endif

  transport(dt);
  
  GiveSpecificTime (t_Hydro);
  
  if (ForwardOneStep == YES) prs_exit(EXIT_SUCCESS);
  
#ifdef MHD
  if(Fluidtype == GAS) {   // We do MHD only for the gaseous component
   *(Emfx->owner) = Emfx;  // EMFs claim ownership of their storage area
   *(Emfy->owner) = Emfy;
   *(Emfz->owner) = Emfz;
 }
#endif

}
