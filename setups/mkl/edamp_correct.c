//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
//<\INCLUDES>


void Edamp_correct_cpu(real dt) {


//<USER_DEFINED>
  INPUT(Energy);
  OUTPUT(Energy);
  INPUT(Density); // Added Density as it is used
  INPUT2D(Chi_profile);
  INPUT(T_guess);
  INPUT(Div_old);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int ll;
  real bigrad;
  real Hg;
  real R3;
  real omk;
  real cs;
  real cs2;
  real tauc;
  real DeltaRin;
  real DeltaRout;
  real TaperRin;
  real TaperRout;
  real taucbcgas;
  real betas;
  real smallq;
  real dr;
  real dz;
  real r;
//<\INTERNAL>

//<EXTERNAL>
  real* e = Energy->field_cpu;
  real* rho = Density->field_cpu;
  real* chi_profile = Chi_profile->field_cpu;
  real* T_predict = T_guess->field_cpu;
  real* D_old     = Div_old->field_cpu;
  real beta = BETA;
  real gam = GAMMA;
  real r0 = R0;
  real bigg = G;
  real mstar = MSTAR;
  real asp = ASPECTRATIO;
  real flar = FLARINGINDEX;
  real a = APAR;
  real b = BPAR;
  real beta0 = BETA0PAR;
  real beta1 = BETA1PAR;
  real ymin_val = YMIN;
  real ymax_val = YMAX;
  real bnd = BETATAPER;
  int pitch = Pitch_cpu;
  int stride = Stride_cpu;
  int size_x = Nx;
  int size_y = Ny + 2 * NGHY;
  int size_z = Nz + 2 * NGHZ;
  int pitch2d = Pitch2D; 
//<\EXTERNAL>


//<CONSTANT>
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>
  for (k = 0; k < size_z; k++) {
    for (j = 0; j < size_y; j++) {
      for (i = 0; i < size_x; i++) {
//<#>
        ll = i + j * pitch + k * stride;

#ifdef THERMALDIFFUSION

        int Ny = size_y - 2 * NGHY;
        int Nz = size_z - 2 * NGHZ;
        if (j >= NGHY && j < Ny + NGHY && k >= NGHZ && k < Nz + NGHZ) {

            dr = ymed(j+1) - ymed(j);
#ifdef Z
            dz = zmed(k+1) - zmed(k);
#endif
#ifdef X
            int Nx_loc = size_x;
            real dphi = (real)(2.0 * M_PI / (double) Nx_loc);
#endif
            real r = ymed(j);

            // Step 1: Recover initial Temperature T^n and predicted Temperature T_guess
            real T_n = e[ll] * (gam - 1.0) / rho[ll];
            real T_g = T_predict[ll];

            // Step 2: Get neighbor temperatures from the predictor step
            real T_rp = T_predict[ll + pitch];
            real T_rm = T_predict[ll - pitch];
#ifdef Z
            real T_zp = T_predict[ll + stride];
            real T_zm = T_predict[ll - stride];
#endif
#ifdef X
            int ip1 = (i + 1) % Nx_loc;
            int im1 = (i - 1 + Nx_loc) % Nx_loc;
            int l_ip1 = ip1 + j * pitch + k * stride;
            int l_im1 = im1 + j * pitch + k * stride;
            real T_pp = T_predict[l_ip1];
            real T_pm = T_predict[l_im1];
#endif

#ifdef TPERTURBATION
            // Step 3: Calculate Equilibrium Temperature T0 = cs^2
#ifdef CYLINDRICAL
            bigrad = ymed(j);
#else // SPHERICAL
            bigrad = ymed(j) * sin(zmed(k));
#endif
            Hg = asp * pow(bigrad / r0, flar) * bigrad;
            R3 = bigrad * bigrad * bigrad;
            omk = sqrt(bigg * mstar / R3);
            cs2 = omk * omk * Hg * Hg;
            real T0 = cs2;

            // Step 4: Convert ALL temperatures to deviations T' = T - T0
            // --- Bug Fix: Added T_g and fixed missing semicolons ---
            T_g  = T_g  - T0;
            T_rp = T_rp - T0;
            T_rm = T_rm - T0;
  #ifdef Z
            T_zp = T_zp - T0;
            T_zm = T_zm - T0;
  #endif
  #ifdef X
            T_pp = T_pp - T0;
            T_pm = T_pm - T0;
  #endif
#endif
            // Step 5: Compute Divergence (now on T')
            real chi_c       = chi_profile[l2D];
            real chi_r_plus  = 0.5 * (chi_profile[l2D] + chi_profile[l2D+1]);
            real chi_r_minus = 0.5 * (chi_profile[l2D-1] + chi_profile[l2D]);

            real dT_dr_plus   = (T_rp - T_g) / dr;
            real dT_dr_minus  = (T_g - T_rm) / dr;
            real flux_r_plus  = chi_r_plus * dT_dr_plus;
            real flux_r_minus = chi_r_minus * dT_dr_minus;
            real div_r = ((r + 0.5 * dr) * flux_r_plus - (r - 0.5 * dr) * flux_r_minus) / (r * dr);

#ifdef X
            real dT_dphi_plus   = (T_pp - T_g) / dphi;
            real dT_dphi_minus  = (T_g - T_pm) / dphi;
            real flux_phi_plus  = chi_c * dT_dphi_plus;
            real flux_phi_minus = chi_c * dT_dphi_minus;
            real div_phi = (flux_phi_plus - flux_phi_minus) / (r * dphi);
            div_phi = div_phi / r;
#else
            real div_phi = 0.0;
#endif
#ifdef Z
            real dT_dz_plus   = (T_zp - T_g) / dz;
            real dT_dz_minus  = (T_g - T_zm) / dz;
            real flux_z_plus  = chi_c * dT_dz_plus;
            real flux_z_minus = chi_c * dT_dz_minus;
            real div_z = (flux_z_plus - flux_z_minus) / dz;
#else
            real div_z = 0.0;
#endif

            real D_guess = div_r + div_phi + div_z;
            real D_old_val = D_old[ll];

            // Step 6: Final update. This part is correct as is.
            real T_final = T_n + 0.5 * dt * (D_old_val + D_guess);
            e[ll] = rho[ll] * T_final / (gam - 1.0);
        }
#endif //End THERMALDIFFUSION

//<\#>
      }
    }
  }
//<\MAIN_LOOP>
}
