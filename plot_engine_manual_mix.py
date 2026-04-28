"""
plot_engine_manual_mix.py

Run the CalcShore thermal engine in adiabatic mode with manually-specified
mix kinetics + composition parameters, optionally compare against a CW
reference centerline curve.

Edit the INPUTS block below, then run:
    python plot_engine_manual_mix.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# Engine + loader. Required to be importable from the same folder.
from cw_scenario_loader import CWMixDesign
from thermal_engine_2d import build_grid_rectangular, solve_hydration_2d


# ======================================================================
# ============================  INPUTS  ================================
# ======================================================================
# EDIT THESE VALUES, then run the script. No command-line flags needed.

# ---- Kinetics (Schindler-Folliard 5-parameter model) ----
# Read these straight from CW's "Material Properties" dialog or input.dat.
# Use full precision -- a 0.5% Hu rounding can move peak T by ~0.7 deg F.
Ea_J_mol     = 25000     # activation energy
tau_hrs      = 25.999     # time constant
beta         = 0.913     # shape parameter
alpha_u      = 0.9925      # ultimate degree of hydration
Hu_J_kg      = 744.896    # ultimate heat per kg cementitious

# Optional empirical multiplier on Hu (e.g. 0.95 for SCM derating test)
Hu_factor    = 1                          # <- set to 0.95 to test the SCM derating

# ---- Mix composition (lb/yd^3) ----
cement_lb_yd3      = 1    # Type I/II cement
fly_ash_F_lb_yd3   = 0.0    # Class F fly ash
fly_ash_C_lb_yd3   = 0      # Class C fly ash
ggbfs_lb_yd3       = 0.0    # Slag
silica_fume_lb_yd3 = 575.0      # Silica fume
water_lb_yd3       = 253.0
coarse_agg_lb_yd3  = 1800.0
fine_agg_lb_yd3    = 1100.0

# ---- Mix optional (CW-typical defaults) ----
air_content_pct          = 5.0       # %
aggregate_Cp_BTU_lb_F    = 0.20467   # CW dialog displays 0.20; use full precision if available
thermal_k_BTU_hr_ft_F    = 1.55509   # CW dialog displays 1.56; only matters in non-adiabatic
CTE_microstrain_F        = 4.25
fly_ash_CaO_pct          = 10.0

# ---- Run-time settings ----
T0_F             = 73.0     # initial temperature (deg F)
duration_hrs     = 168.0    # simulation duration (hours)
output_interval_s = 300.0   # engine sampling interval (5 min matches CW grid)

# ---- I/O ----
# Set CW_CSV = None to skip the CW comparison (engine-only run).
CW_CSV   = "docs/kinetics_isolation_handoff/cw_adiabatic_reference_mix00H.csv"
# CW_CSV = None  # uncomment if you don't have a CW reference for this mix

OUT_PNG  = "engine_manual_run.png"
OUT_CSV  = "engine_manual_run.csv"
MAKE_PLOT = True

# ======================================================================
# ===========================  END INPUTS  =============================
# ======================================================================


# ----------------------------------------------------------------------
# Conversions
# ----------------------------------------------------------------------

def f_to_c(t_f: float) -> float:
    return (t_f - 32.0) * 5.0 / 9.0


def c_to_f(t_c):
    return np.asarray(t_c) * 9.0 / 5.0 + 32.0


# ----------------------------------------------------------------------
# Build mix
# ----------------------------------------------------------------------

def build_mix() -> CWMixDesign:
    Hu_effective = Hu_J_kg * Hu_factor
    return CWMixDesign(
        cement_type_I_II_lb_yd3=cement_lb_yd3,
        fly_ash_F_lb_yd3=fly_ash_F_lb_yd3,
        fly_ash_C_lb_yd3=fly_ash_C_lb_yd3,
        ggbfs_lb_yd3=ggbfs_lb_yd3,
        silica_fume_lb_yd3=silica_fume_lb_yd3,
        water_lb_yd3=water_lb_yd3,
        coarse_agg_lb_yd3=coarse_agg_lb_yd3,
        fine_agg_lb_yd3=fine_agg_lb_yd3,
        air_content_pct=air_content_pct,
        fly_ash_CaO_pct=fly_ash_CaO_pct,
        activation_energy_J_mol=Ea_J_mol,
        tau_hrs=tau_hrs,
        beta=beta,
        alpha_u=alpha_u,
        Hu_J_kg=Hu_effective,
        thermal_conductivity_BTU_hr_ft_F=thermal_k_BTU_hr_ft_F,
        aggregate_Cp_BTU_lb_F=aggregate_Cp_BTU_lb_F,
        CTE_microstrain_F=CTE_microstrain_F,
    )


# ----------------------------------------------------------------------
# Engine run
# ----------------------------------------------------------------------

def run_engine(mix: CWMixDesign):
    grid = build_grid_rectangular(Lx_m=1.0, Ly_m=1.0, nx=3, ny=3)
    T0_C = f_to_c(T0_F)
    T_init = np.full((grid.ny, grid.nx), T0_C, dtype=np.float64)
    res = solve_hydration_2d(
        grid, mix, T_init,
        duration_s=duration_hrs * 3600.0,
        output_interval_s=output_interval_s,
        boundary_mode="adiabatic",
    )
    t_h = res.t_s / 3600.0
    T_F = c_to_f(res.T_field_C[:, grid.ny // 2, grid.nx // 2])
    return t_h, T_F


# ----------------------------------------------------------------------
# Reference loader
# ----------------------------------------------------------------------

def load_cw_reference(csv_path):
    arr = np.genfromtxt(str(csv_path), delimiter=",", skip_header=1)
    return arr[:, 0], arr[:, 1]


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def make_plot_with_cw(t_eng_h, T_eng_F, t_cw_h, T_cw_F, mix, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    T_eng_on_cw = np.interp(t_cw_h, t_eng_h, T_eng_F)
    delta = T_eng_on_cw - T_cw_F
    peak_eng = float(T_eng_F[-1])
    peak_cw = float(T_cw_F[-1])
    peak_delta = peak_eng - peak_cw
    rms = float(np.sqrt(np.mean(delta ** 2)))
    max_abs = float(np.max(np.abs(delta)))

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax1.plot(t_cw_h, T_cw_F, "k-", linewidth=1.8,
             label=f"CW reference  (peak {peak_cw:.2f} F)")
    ax1.plot(t_eng_h, T_eng_F, "C0--", linewidth=1.5,
             label=f"Engine adiabatic, Hu*{Hu_factor:.3f}  "
                   f"(peak {peak_eng:.2f} F)")
    ax1.axhline(100, color="gray", linewidth=0.6, linestyle=":", alpha=0.5)
    ax1.axhline(130, color="gray", linewidth=0.6, linestyle=":", alpha=0.5)
    ax1.set_ylabel("Temperature (deg F)")
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    cm = mix.total_cementitious_lb_yd3
    scm_pct = 100.0 * (cm - mix.cement_type_I_II_lb_yd3) / cm if cm > 0 else 0.0
    ax1.set_title(
        f"Engine vs CW -- T0={T0_F:.0f}F, duration={duration_hrs:.0f}h, "
        f"adiabatic, SCM={scm_pct:.1f}%\n"
        f"Ea={mix.activation_energy_J_mol:.0f}  tau={mix.tau_hrs:.2f}  "
        f"beta={mix.beta:.3f}  alpha_u={mix.alpha_u:.4f}  "
        f"Hu_eff={mix.Hu_J_kg:.0f} (Hu*{Hu_factor:.3f})"
    )

    ax2.plot(t_cw_h, delta, "C3-", linewidth=1.0)
    ax2.axhline(0.0, color="k", linewidth=0.5)
    ax2.fill_between(t_cw_h, -1.0, 1.0, color="green", alpha=0.10,
                     label="+/- 1 deg F gate")
    ax2.set_xlabel("Time since placement (hours)")
    ax2.set_ylabel("Engine - CW (deg F)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=9)
    plt.tight_layout()
    plt.savefig(str(out_png), dpi=140)
    plt.close(fig)
    return {"peak_engine_F": peak_eng, "peak_cw_F": peak_cw,
            "peak_delta_F": peak_delta, "rms_F": rms,
            "max_abs_delta_F": max_abs}


def make_plot_engine_only(t_eng_h, T_eng_F, mix, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(t_eng_h, T_eng_F, "C0-", linewidth=1.6,
             label=f"Engine adiabatic, Hu*{Hu_factor:.3f}  "
                   f"(peak {T_eng_F[-1]:.2f} F)")
    ax1.axhline(100, color="gray", linewidth=0.6, linestyle=":", alpha=0.5)
    ax1.axhline(130, color="gray", linewidth=0.6, linestyle=":", alpha=0.5)
    ax1.set_xlabel("Time since placement (hours)")
    ax1.set_ylabel("Temperature (deg F)")
    cm = mix.total_cementitious_lb_yd3
    scm_pct = 100.0 * (cm - mix.cement_type_I_II_lb_yd3) / cm if cm > 0 else 0.0
    ax1.set_title(
        f"Engine adiabatic -- T0={T0_F:.0f}F, duration={duration_hrs:.0f}h, "
        f"SCM={scm_pct:.1f}%\n"
        f"Ea={mix.activation_energy_J_mol:.0f}  tau={mix.tau_hrs:.2f}  "
        f"beta={mix.beta:.3f}  alpha_u={mix.alpha_u:.4f}  "
        f"Hu_eff={mix.Hu_J_kg:.0f} (Hu*{Hu_factor:.3f})"
    )
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(str(out_png), dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------
# CSV writers
# ----------------------------------------------------------------------

def write_csv_with_cw(t_cw_h, T_cw_F, t_eng_h, T_eng_F, out_csv):
    T_eng_on_cw = np.interp(t_cw_h, t_eng_h, T_eng_F)
    delta = T_eng_on_cw - T_cw_F
    arr = np.column_stack([t_cw_h, T_eng_on_cw, T_cw_F, delta])
    header = "time_hrs,T_engine_F,T_cw_F,delta_F"
    np.savetxt(str(out_csv), arr, delimiter=",", header=header,
               fmt=("%.4f", "%.4f", "%.4f", "%.5f"), comments="")


def write_csv_engine_only(t_eng_h, T_eng_F, out_csv):
    arr = np.column_stack([t_eng_h, T_eng_F])
    header = "time_hrs,T_engine_F"
    np.savetxt(str(out_csv), arr, delimiter=",", header=header,
               fmt=("%.4f", "%.4f"), comments="")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main():
    mix = build_mix()
    cm = mix.total_cementitious_lb_yd3
    scm_pct = 100.0 * (cm - mix.cement_type_I_II_lb_yd3) / cm if cm > 0 else 0.0

    print("Mix composition:")
    print(f"  cementitious total = {cm:.1f} lb/yd^3  ({scm_pct:.1f}% SCM)")
    print(f"  density            = {mix.concrete_density_lb_ft3:.2f} lb/ft^3")
    print(f"  w/cm               = {mix.wcm:.3f}")
    print()
    print("Kinetics parameters:")
    print(f"  Ea       = {mix.activation_energy_J_mol:>10.2f}   J/mol")
    print(f"  tau      = {mix.tau_hrs:>10.4f}   hr")
    print(f"  beta     = {mix.beta:>10.4f}")
    print(f"  alpha_u  = {mix.alpha_u:>10.5f}")
    print(f"  Hu (raw) = {Hu_J_kg:>10.1f}   J/kg")
    print(f"  Hu_factor= {Hu_factor:>10.4f}")
    print(f"  Hu (eff) = {mix.Hu_J_kg:>10.1f}   J/kg")
    print()

    print(f"Running engine: T0={T0_F} F, "
          f"duration={duration_hrs} hr, adiabatic ...")
    t_eng, T_eng = run_engine(mix)
    print(f"  engine peak at t={duration_hrs:.0f} hr: {T_eng[-1]:.3f} deg F")

    if CW_CSV is not None and Path(CW_CSV).exists():
        print(f"\nLoading CW reference: {CW_CSV}")
        t_cw, T_cw = load_cw_reference(CW_CSV)
        print(f"  CW reference: {len(t_cw)} samples, "
              f"t in [{t_cw[0]:.2f}, {t_cw[-1]:.2f}] hr, "
              f"T in [{T_cw.min():.2f}, {T_cw.max():.2f}] F")
        if MAKE_PLOT:
            print(f"\nWriting plot -> {OUT_PNG}")
            metrics = make_plot_with_cw(t_eng, T_eng, t_cw, T_cw, mix, OUT_PNG)
            print(f"\n{'=' * 56}")
            print("Engine vs CW summary")
            print(f"{'=' * 56}")
            print(f"  Engine peak at t={duration_hrs:.0f} hr: "
                  f"{metrics['peak_engine_F']:>7.2f} deg F")
            print(f"  CW     peak at t={duration_hrs:.0f} hr: "
                  f"{metrics['peak_cw_F']:>7.2f} deg F")
            print(f"  delta peak:                "
                  f"{metrics['peak_delta_F']:>+7.2f} deg F")
            print(f"  RMS over [0, {duration_hrs:.0f}] hr:    "
                  f"{metrics['rms_F']:>7.2f} deg F")
            print(f"  max |residual|:            "
                  f"{metrics['max_abs_delta_F']:>7.2f} deg F")
            print(f"{'=' * 56}")
        print(f"Writing data -> {OUT_CSV}")
        write_csv_with_cw(t_cw, T_cw, t_eng, T_eng, OUT_CSV)
    else:
        if CW_CSV is not None:
            print(f"\nWARNING: CW reference '{CW_CSV}' not found. "
                  f"Skipping comparison.")
        if MAKE_PLOT:
            print(f"\nWriting plot -> {OUT_PNG}")
            make_plot_engine_only(t_eng, T_eng, mix, OUT_PNG)
        print(f"Writing data -> {OUT_CSV}")
        write_csv_engine_only(t_eng, T_eng, OUT_CSV)

    return 0


if __name__ == "__main__":
    sys.exit(main())
