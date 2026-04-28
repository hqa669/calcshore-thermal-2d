"""
plot_engine_vs_cw_centerline.py

Plot the CalcShore engine's adiabatic temperature trajectory against the
ConcreteWorks centerline reference curve at T0 = 73 deg F, using the
default 5 kinetics parameters from MIX-01 input.dat (no overrides):

    Ea     = 26458 J/mol      (activation energy)
    tau    = 29.40 hr         (time constant)
    beta   = 0.895            (shape parameter)
    alpha_u= 0.7585           (ultimate degree of hydration)
    Hu     = 424143 J/kg      (ultimate heat of hydration)

Files expected in the same directory (or pass paths via flags):
  - thermal_engine_2d.py
  - cw_scenario_loader.py
  - input.dat
  - cw_adiabatic_reference_mix01.csv

Outputs:
  - engine_vs_cw_centerline_73F.png   (2-panel: trajectory + residual)
  - engine_vs_cw_centerline_73F.csv   (raw data: t, T_engine, T_cw, residual)

Usage:
  python plot_engine_vs_cw_centerline.py
  python plot_engine_vs_cw_centerline.py --T0-F 73 --duration-hrs 168
  python plot_engine_vs_cw_centerline.py --input-dat path/to/input.dat \\
                                         --cw-csv path/to/reference.csv \\
                                         --out-png my_plot.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# These two imports require the engine and loader to be importable.
# Easiest way: put this script in the same folder as thermal_engine_2d.py
# and cw_scenario_loader.py, then run it from that folder.
from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_rectangular, solve_hydration_2d


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def f_to_c(t_f: float) -> float:
    return (t_f - 32.0) * 5.0 / 9.0


def c_to_f(t_c):
    return np.asarray(t_c) * 9.0 / 5.0 + 32.0


# ----------------------------------------------------------------------
# Engine run
# ----------------------------------------------------------------------

def run_engine_adiabatic(input_dat_path: str | Path,
                         T0_F: float,
                         duration_hrs: float,
                         output_interval_s: float = 300.0,
                         use_cw_calibrated_hu: bool = True):
    """
    Run solve_hydration_2d in adiabatic mode at uniform T0.

    The apr28 composition-based Hu calibration is applied by default
    via load_cw_scenario; pass use_cw_calibrated_hu=False to reproduce
    the pre-correction (raw-Hu) baseline for diagnostics.

    Returns (t_hrs, T_F, mix_kinetics_summary).
    """
    scenario = load_cw_scenario(str(input_dat_path),
                                weather_dat=None, cw_output_txt=None,
                                use_cw_calibrated_hu=use_cw_calibrated_hu)
    mix = scenario.mix

    # In adiabatic mode the field stays uniform, so geometry is irrelevant.
    # 3 x 3 is the minimum allowed by build_grid_rectangular.
    grid = build_grid_rectangular(Lx_m=1.0, Ly_m=1.0, nx=3, ny=3)

    T0_C = f_to_c(T0_F)
    T_init = np.full((grid.ny, grid.nx), T0_C, dtype=np.float64)

    res = solve_hydration_2d(
        grid, mix, T_init,
        duration_s=duration_hrs * 3600.0,
        output_interval_s=output_interval_s,
        boundary_mode="adiabatic",
    )

    t_hrs = res.t_s / 3600.0
    # Pick interior cell; in adiabatic mode every cell follows the same ODE.
    T_F = c_to_f(res.T_field_C[:, grid.ny // 2, grid.nx // 2])

    summary = {
        "Ea_J_mol":             mix.activation_energy_J_mol,
        "tau_hrs":              mix.tau_hrs,
        "beta":                 mix.beta,
        "alpha_u":              mix.alpha_u,
        "Hu_J_kg":              mix.Hu_J_kg,
        "Hu_factor_calibrated": mix.Hu_factor_calibrated,
        "Hu_J_kg_effective":    mix.Hu_J_kg_effective,
        "Hu_calibration_note":  mix.Hu_calibration_note,
        "n_inner_steps":        res.n_inner_steps,
        "dt_inner_s":           res.dt_inner_s,
    }
    return t_hrs, T_F, summary


# ----------------------------------------------------------------------
# Reference loader
# ----------------------------------------------------------------------

def load_cw_reference(csv_path: str | Path):
    """
    Load CW adiabatic centerline reference. Uses column 1
    (T_center_F_adiabatic) -- column 0 is time in hours.
    """
    arr = np.genfromtxt(str(csv_path), delimiter=",", skip_header=1)
    return arr[:, 0], arr[:, 1]


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def make_plot(t_eng_h, T_eng_F, t_cw_h, T_cw_F, summary, out_png):
    import matplotlib
    matplotlib.use("Agg")  # no display required
    import matplotlib.pyplot as plt

    # Residual on CW time grid
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
             label=f"CW centerline @ T0=73 deg F  (peak {peak_cw:.2f} F)")
    ax1.plot(t_eng_h, T_eng_F, "C0--", linewidth=1.5,
             label=f"Engine v3 adiabatic  (peak {peak_eng:.2f} F)")
    ax1.axhline(100, color="gray", linewidth=0.6, linestyle=":", alpha=0.5)
    ax1.axhline(130, color="gray", linewidth=0.6, linestyle=":", alpha=0.5)
    ax1.set_ylabel("Temperature (deg F)")
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    ax1.set_title(
        "Engine v3 vs CW centerline -- MIX-01, T0 = 73 deg F, adiabatic\n"
        f"Defaults: Ea={summary['Ea_J_mol']:.0f}  tau={summary['tau_hrs']:.2f}  "
        f"beta={summary['beta']:.3f}  alpha_u={summary['alpha_u']:.4f}  "
        f"Hu_factor={summary['Hu_factor_calibrated']:.4f}  "
        f"(Hu_eff={summary['Hu_J_kg_effective']:.0f} J/kg)"
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

    metrics = {
        "peak_engine_F": peak_eng,
        "peak_cw_F": peak_cw,
        "peak_delta_F": peak_delta,
        "rms_F": rms,
        "max_abs_delta_F": max_abs,
    }
    return metrics


def write_csv(t_cw_h, T_cw_F, t_eng_h, T_eng_F, out_csv):
    """
    Write a CSV with engine + CW + residual on the CW time grid.
    """
    T_eng_on_cw = np.interp(t_cw_h, t_eng_h, T_eng_F)
    delta = T_eng_on_cw - T_cw_F
    arr = np.column_stack([t_cw_h, T_eng_on_cw, T_cw_F, delta])
    header = "time_hrs,T_engine_F,T_cw_F,delta_F"
    np.savetxt(str(out_csv), arr, delimiter=",", header=header,
               fmt=("%.4f", "%.4f", "%.4f", "%.5f"), comments="")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dat", default="input.dat",
                   help="MIX-01 CW input.dat (default: ./input.dat)")
    p.add_argument("--cw-csv", default="cw_adiabatic_reference_mix01.csv",
                   help="CW adiabatic reference CSV.")
    p.add_argument("--T0-F", type=float, default=73.0,
                   help="Initial temperature in deg F (default 73).")
    p.add_argument("--duration-hrs", type=float, default=168.0,
                   help="Simulation duration in hours (default 168).")
    p.add_argument("--output-interval-s", type=float, default=300.0,
                   help="Engine output sample interval in seconds (default 300).")
    p.add_argument("--out-png", default="engine_vs_cw_centerline_73F.png",
                   help="Output PNG path.")
    p.add_argument("--out-csv", default="engine_vs_cw_centerline_73F.csv",
                   help="Output CSV path.")
    p.add_argument("--raw-hu", action="store_true",
                   help="Diagnostic: disable apr28 Hu calibration "
                        "(use the verbatim CW-regressed Hu_J_kg). "
                        "Reproduces the +4.32°F MIX-01 overshoot baseline.")
    args = p.parse_args(argv)

    use_calibrated = not args.raw_hu

    # ---- Run engine ----
    mode = "raw-Hu (DIAGNOSTIC)" if args.raw_hu else "apr28 calibrated"
    print(f"Running engine: T0={args.T0_F} deg F, "
          f"duration={args.duration_hrs} hr, "
          f"output_interval={args.output_interval_s} s, adiabatic, {mode} ...")
    t_eng, T_eng, summary = run_engine_adiabatic(
        args.input_dat, args.T0_F, args.duration_hrs, args.output_interval_s,
        use_cw_calibrated_hu=use_calibrated,
    )
    print("Kinetics parameters consumed (from input.dat):")
    print(f"  Ea          = {summary['Ea_J_mol']:>10.2f}   J/mol")
    print(f"  tau         = {summary['tau_hrs']:>10.4f}   hr")
    print(f"  beta        = {summary['beta']:>10.4f}")
    print(f"  alpha_u     = {summary['alpha_u']:>10.5f}")
    print(f"  Hu (raw)    = {summary['Hu_J_kg']:>10.1f}   J/kg")
    print(f"  Hu_factor   = {summary['Hu_factor_calibrated']:>10.4f}")
    print(f"  Hu (eff)    = {summary['Hu_J_kg_effective']:>10.1f}   J/kg "
          f"(consumed by solver)")
    if summary["Hu_calibration_note"]:
        print(f"  envelope    = {summary['Hu_calibration_note']}")
    print(f"  inner steps = {summary['n_inner_steps']}, "
          f"dt_inner = {summary['dt_inner_s']:.1f} s")

    # ---- Load CW reference ----
    print(f"\nLoading CW reference: {args.cw_csv}")
    t_cw, T_cw = load_cw_reference(args.cw_csv)
    print(f"  {len(t_cw)} samples, t in [{t_cw[0]:.2f}, {t_cw[-1]:.2f}] hr, "
          f"T in [{T_cw.min():.2f}, {T_cw.max():.2f}] F")

    # ---- Plot + CSV ----
    print(f"\nWriting plot -> {args.out_png}")
    metrics = make_plot(t_eng, T_eng, t_cw, T_cw, summary, args.out_png)
    print(f"Writing data -> {args.out_csv}")
    write_csv(t_cw, T_cw, t_eng, T_eng, args.out_csv)

    # ---- Summary ----
    print()
    print("=" * 56)
    print("Summary")
    print("=" * 56)
    print(f"  Engine peak at t={args.duration_hrs:.0f} hr: "
          f"{metrics['peak_engine_F']:>7.2f} deg F")
    print(f"  CW     peak at t={args.duration_hrs:.0f} hr: "
          f"{metrics['peak_cw_F']:>7.2f} deg F")
    print(f"  delta peak:                {metrics['peak_delta_F']:>+7.2f} deg F")
    print(f"  RMS over [0, {args.duration_hrs:.0f}] hr:    "
          f"{metrics['rms_F']:>7.2f} deg F")
    print(f"  max |residual|:            {metrics['max_abs_delta_F']:>7.2f} deg F")
    print("=" * 56)
    return 0


if __name__ == "__main__":
    sys.exit(main())
