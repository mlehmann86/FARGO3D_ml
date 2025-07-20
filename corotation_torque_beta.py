#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import numpy as np
import argparse

# --- Import functions from your existing script ---
try:
    # Adjust 'planet_data' if your script file has a different name
    from planet_data import (
        determine_base_path,
        read_parameters,
        extract_planet_mass, # Assuming this exists or qp is read from params
        read_alternative_torque
        # compute_theoretical_torques # We will re-implement needed logic below
    )
    print("Successfully imported functions from planet_data.py")
except ImportError as e:
    print(f"Error: Could not import functions from 'planet_data.py': {e}")
    print("Please ensure the file exists and is in the Python path.")
    exit()
# Removed compute_theoretical_torques import, will redefine needed logic locally

# --- New Function to Calculate Specific Theoretical Torques ---

def calculate_theory_components(parameters, qp, IDEFIX=False):
    """
    Computes theoretical estimates for normalized linear corotation torque
    and non-linear horseshoe drag based on Paardekooper et al. (2010).
    Re-implements calculation logic using parameters.

    Returns tuple: (linear_CT_total_norm, nonlinear_HSD_total_norm, GAM0)
                   where torques are normalized as Gamma / Gamma0
    """
    try:
        # --- Parameter Extraction (logic adapted from original compute_theoretical_torques) ---
        if IDEFIX:
            gam = float(parameters.get("gamma", 1.666667))
            # Assuming flaringIndex determines temperature slope q = 1 - 2*fi
            fi = float(parameters.get("flaringIndex", 0.0)) # Temp slope beta_temp = q
            bet_temp_slope = 1.0 - 2.0 * fi
            h = float(parameters.get("h0", 0.05))
            alp = float(parameters.get("sigmaSlope", 1.5)) # Density slope alpha = ss
            sigma0 = float(parameters.get("Sigma0", 1.0))
            # Assume smoothing is given directly in units of H
            bh = float(parameters.get("thicknessSmoothing", 0.6)) # Smoothing b/h
        else: # FARGO
            gam = float(parameters.get('GAMMA', 1.666667))
            fi = float(parameters.get('FLARINGINDEX', 0.0))
            bet_temp_slope = 1.0 - 2.0 * fi
            h = float(parameters.get('ASPECTRATIO', 0.05))
            bh = float(parameters.get('THICKNESSSMOOTHING', 0.6))
            alp = float(parameters.get('SIGMASLOPE', 1.5))
            sigma0 = float(parameters.get('SIGMA0', 1.0))

        # --- Basic Input Validation ---
        if h == 0: raise ValueError("Aspect ratio h cannot be zero.")
        if gam == 0: raise ValueError("Gamma cannot be zero.")
        if qp == 0: raise ValueError("Planet mass qp cannot be zero.")
        if bh == 0:
            print("Warning: Zero smoothing length (bh=0). Theory results may be inf/nan.")
            fac = np.inf
        else:
            fac = 0.4 / bh

        xi = bet_temp_slope - (gam - 1) * alp # Entropy gradient index
        GAM0 = (qp / h)**2 * sigma0
        if GAM0 == 0 or not np.isfinite(GAM0):
            raise ValueError(f"GAM0 normalization is zero or invalid: {GAM0}")
        # --- End Parameter Extraction ---

        # --- Calculate Torque Components (Formulas yield gamma * Gamma / Gamma0) ---

        # Linear CT components (based on Eq. 18 logic, Paardekooper 2010)
        linear_CT_baro_gG = 0.7 * (1.5 - alp - 2.0 * xi / gam) * fac**1.26
        linear_CT_ent_gG = 2.2 * xi * fac**0.71
        linear_CT_total_gG_norm = linear_CT_baro_gG + linear_CT_ent_gG # gamma*Gamma/Gamma0

        # Nonlinear HSD components (based on Eq. 45 logic, Paardekooper 2010)
        nonlinear_HSD_baro_gG = 1.1 * fac * (1.5 - alp)
        # Handle potential issues with fac term for entropy part
        term_in_sqrt = fac
        if term_in_sqrt < 0 or np.isinf(term_in_sqrt): term_in_sqrt = 0
        nonlinear_HSD_ent_gG = (xi / gam) * fac * (10.1 * term_in_sqrt**0.5 - 2.2)
        nonlinear_HSD_total_gG_norm = nonlinear_HSD_baro_gG + nonlinear_HSD_ent_gG # gamma*Gamma/Gamma0

        # --- Final Normalization to Gamma / Gamma0 ---
        linear_CT_total_norm = linear_CT_total_gG_norm / gam
        nonlinear_HSD_total_norm = nonlinear_HSD_total_gG_norm / gam

        # Check for NaN/Inf results before returning
        if not np.isfinite(linear_CT_total_norm): linear_CT_total_norm = np.nan
        if not np.isfinite(nonlinear_HSD_total_norm): nonlinear_HSD_total_norm = np.nan

        # print(f"Debug Theory: qp={qp}, h={h}, sigma0={sigma0}, GAM0={GAM0}")
        # print(f"Debug Theory: gam={gam}, alp={alp}, bet_temp={bet_temp_slope}, bh={bh}, fac={fac}, xi={xi}")
        # print(f"Debug Theory: LinCT_Norm={linear_CT_total_norm}, NLinHSD_Norm={nonlinear_HSD_total_norm}")


        return linear_CT_total_norm, nonlinear_HSD_total_norm, GAM0

    except KeyError as e:
        print(f"Error: Missing parameter key for theoretical torque: {e}")
        raise ValueError("Missing parameters for theoretical torque.")
    except ValueError as e:
         print(f"Error: Invalid parameter value for theoretical torque: {e}")
         raise ValueError("Invalid parameters for theoretical torque.")
    except Exception as e:
         print(f"Unexpected error during theoretical torque calculation: {e}")
         raise

# --- Main Analysis Function ---

def analyze_cooling_effect_smooth_peak(simulation_list,
                                     peak_time_window_orbits=(5, 100),
                                     smoothing_time_orbits=1.0,
                                     IDEFIX=False):
    """
    Extracts smoothed peak torque for simulations with different cooling times (beta)
    and plots it vs beta, comparing with theoretical estimates.
    """
    # Imports needed within this function
    import matplotlib.pyplot as plt
    from scipy.ndimage import uniform_filter1d

    beta_values = []
    peak_smooth_torque_norm_values = []
    processed_sims = []

    theoretical_ct_lin_norm = np.nan
    theoretical_hsd_nl_norm = np.nan
    GAM0_ref = np.nan
    first_sim = True

    for sim_name in simulation_list:
        print(f"\nProcessing: {sim_name}")
        try:
            base_path = determine_base_path(sim_name, IDEFIX=IDEFIX)
            output_path = base_path

            param_file = os.path.join(output_path, "idefix.0.log") if IDEFIX else os.path.join(output_path, "summary0.dat")
            torque_file = os.path.join(output_path, "tqwk0.dat") if IDEFIX else os.path.join(output_path, "torq_planet_0.dat")

            if not os.path.exists(param_file) or not os.path.exists(torque_file):
                print(f"Skipping {sim_name}: Missing parameter or torque file.")
                continue

            parameters = read_parameters(param_file, IDEFIX=IDEFIX)
            beta_str = parameters.get("beta") if IDEFIX else parameters.get("BETA")
            if beta_str is None: raise ValueError("Cooling time not found.")
            beta = float(beta_str)

            # Get qp (adapt as needed based on where qp is stored)
            qp = extract_planet_mass(parameters, IDEFIX=IDEFIX) # Assumes this function works
            if qp == 0.0: raise ValueError("Planet mass is zero.")

            # Read TOTAL torque using imported function
            time, total_torque, orbit_numbers = read_alternative_torque(torque_file, IDEFIX=IDEFIX)
            if time.size < 2: raise ValueError("Not enough torque data.")

            time_in_orbits = time if IDEFIX else time / (2.0 * np.pi)

            # --- Calculate Theoretical Values & Gamma0 (only once) ---
            if first_sim:
                 # Use the NEW local function
                lin_ct_norm, nl_hsd_norm, GAM0_calc = calculate_theory_components(parameters, qp, IDEFIX=IDEFIX)
                if np.isnan(lin_ct_norm) or np.isnan(nl_hsd_norm) or np.isnan(GAM0_calc):
                    print("Warning: Failed to calculate theoretical values. No theory lines plotted.")
                    theoretical_ct_lin_norm = np.nan
                    theoretical_hsd_nl_norm = np.nan
                    GAM0_ref = np.nan
                else:
                    theoretical_ct_lin_norm = lin_ct_norm
                    theoretical_hsd_nl_norm = nl_hsd_norm
                    GAM0_ref = GAM0_calc # Store reference Gamma0
                first_sim = False

            # Need GAM0_ref to proceed
            if np.isnan(GAM0_ref) or GAM0_ref == 0:
                 print(f"Skipping {sim_name} due to invalid GAM0_ref.")
                 continue

            # --- Apply Smoothing ---
            dt_orbit = np.mean(np.diff(time_in_orbits))
            if dt_orbit <= 0 or not np.isfinite(dt_orbit):
                print("Warning: Invalid time step. Using default smoothing window=10.")
                smoothing_window_size = 10
            else:
                smoothing_window_size = max(1, int(smoothing_time_orbits / dt_orbit))
            # print(f"  Smoothing window size: {smoothing_window_size} points")
            smoothed_torque = uniform_filter1d(total_torque, size=smoothing_window_size, mode='nearest')

            # --- Find Peak Smoothed Torque in Window ---
            t_start, t_end = peak_time_window_orbits
            mask = (time_in_orbits >= t_start) & (time_in_orbits <= t_end)
            if not np.any(mask):
                print(f"Warning: No data in peak window [{t_start}, {t_end}]. Skipping.")
                continue

            time_in_window = time_in_orbits[mask]
            smoothed_torque_in_window = smoothed_torque[mask]
            if smoothed_torque_in_window.size == 0: continue # Skip if empty after masking

            idx_max_abs = np.argmax(np.abs(smoothed_torque_in_window))
            peak_smoothed_signed_torque = smoothed_torque_in_window[idx_max_abs]
            time_of_peak = time_in_window[idx_max_abs]

            # Normalize by Gamma0
            peak_torque_norm = peak_smoothed_signed_torque / GAM0_ref # Gamma / Gamma0

            # Store results
            beta_values.append(beta)
            peak_smooth_torque_norm_values.append(peak_torque_norm)
            processed_sims.append(sim_name)
            print(f"  Beta = {beta:.2e}, Peak Smooth Norm. Torque = {peak_torque_norm:.4f} (at {time_of_peak:.1f} orbits)")

        except FileNotFoundError as e: print(f"FNF Error for {sim_name}: {e}")
        except ValueError as e: print(f"Value Error for {sim_name}: {e}")
        except KeyError as e: print(f"Key Error processing {sim_name}: Missing key {e}")
        except Exception as e: print(f"General Error for {sim_name}: {e}")

    # --- Plotting ---
    if not beta_values:
        print("No simulations processed successfully. Cannot create plot.")
        return

    beta_values = np.array(beta_values)
    peak_smooth_torque_norm_values = np.array(peak_smooth_torque_norm_values)
    sort_indices = np.argsort(beta_values)
    beta_values_sorted = beta_values[sort_indices]
    peak_smooth_torque_norm_values_sorted = peak_smooth_torque_norm_values[sort_indices]

    plt.figure(figsize=(10, 6))
    plt.plot(beta_values_sorted, peak_smooth_torque_norm_values_sorted, marker='o', linestyle='-', color='blue',
             label=f'Peak Smoothed Sim. Torque ({peak_time_window_orbits[0]}-{peak_time_window_orbits[1]} orbits, smooth={smoothing_time_orbits:.1f} orb)')

    # Plot theoretical lines (normalized to Gamma/Gamma0)
    if np.isfinite(theoretical_ct_lin_norm):
        plt.axhline(theoretical_ct_lin_norm, color='red', linestyle='--',
                    label=f'Theory Linear CT ($\Gamma / \Gamma_0 = {theoretical_ct_lin_norm:.3f}$)')
        print(f"Plotting Theory Linear CT line at: {theoretical_ct_lin_norm:.3f}")
    if np.isfinite(theoretical_hsd_nl_norm):
        plt.axhline(theoretical_hsd_nl_norm, color='green', linestyle=':',
                    label=f'Theory Non-linear HSD ($\Gamma / \Gamma_0 = {theoretical_hsd_nl_norm:.3f}$)')
        print(f"Plotting Theory Non-linear HSD line at: {theoretical_hsd_nl_norm:.3f}")

    plt.xscale('log')
    plt.xlabel(r'Cooling Time Parameter ($\beta = t_{cool} \Omega_p$)')
    plt.ylabel(r'Normalized Peak Smoothed Torque ($\Gamma_{\mathrm{peak}} / \Gamma_0$)')
    plt.title(f'Peak Smoothed Torque vs. Cooling Time (Smoothed over {smoothing_time_orbits:.1f} orbits)')
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.tight_layout()

    plot_filename = "peak_smoothed_torque_vs_beta.pdf"
    plt.savefig(plot_filename)
    print(f"\nPlot saved to {plot_filename}")
    # plt.show() # Comment out if running non-interactively


# --- Script Execution ---
if __name__ == "__main__":
    # Assume planet_data.py is in the same directory or Python path
    # List of simulation directories provided by the user
    simulation_list = [
        "cos_bet1d2_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
        "cos_bet1d4_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
        "cos_bet1d3_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
        "cos_bet1dm2_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
        "cos_bet1dm1_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
        "cos_bet1d1_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
        "cos_bet1d0_gam53_ss15_q2_r0516_nu1dm11_COR_HR150_2D",
    ]

    parser = argparse.ArgumentParser(description="Analyze peak smoothed torque vs cooling time.")
    parser.add_argument('-s', '--simulations', nargs='+', default=simulation_list,
                        help="List of simulation directory names.")
    parser.add_argument('--idefix', action='store_true', default=False,
                        help="Flag if simulations use IDEFIX format.")
    parser.add_argument('--t_start', type=float, default=5.0,
                        help="Start time (orbits) for peak torque window.")
    parser.add_argument('--t_end', type=float, default=100.0,
                        help="End time (orbits) for peak torque window.")
    parser.add_argument('--smooth', type=float, default=1.0,
                        help="Smoothing timescale (orbits) for torque.")

    args = parser.parse_args()

    analyze_cooling_effect_smooth_peak(
        args.simulations,
        peak_time_window_orbits=(args.t_start, args.t_end),
        smoothing_time_orbits=args.smooth,
        IDEFIX=args.idefix
    )
