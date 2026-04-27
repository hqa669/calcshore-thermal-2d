"""
thermal_engine_test.py — Kinetics-isolation validation for CalcShore engine v3.

Goal
----
Validate that the engine's Schindler-Folliard hydration ODE integration
matches ConcreteWorks' (CW) adiabatic curve for MIX-01 at T0 = 73 deg F,
with all BC physics decoupled (boundary_mode="adiabatic", uniform IC).

Procedure
---------
1. Load MIX-01 mix design from a CW input.dat via cw_scenario_loader.
2. Build a small all-concrete rectangular grid (3 x 3; geometry is
   irrelevant in adiabatic mode by construction).
3. Set initial T uniformly to 73 deg F (= 22.7778 deg C).
4. Run solve_hydration_2d for 168 hours with output_interval_s = 300 s
   (5-min sampling, matching the CW reference grid).
5. Extract any cell trajectory (in adiabatic mode, all cells follow the
   same ODE; we verify spatial uniformity as a sanity check).
6. Compare against T_center_F_adiabatic from
   cw_adiabatic_reference_mix01.csv:
   - Peak deviation at t = 168 hr  (tol  ±1.0 deg F)
   - RMS over [0, 168] hr          (tol  ≤1.0 deg F)
   - Time to T = 100 deg F         (tol  ±2 hr)
   - Time to T = 130 deg F         (tol  ±2 hr)
7. Save a 2-panel comparison plot (engine vs CW + residual).
8. Emit a pass/fail report; exit 0 on pass, 1 on fail.

If this passes, kinetics is verified clean and the v3 residuals
(0.74 deg F Centerline RMS, 2 deg F Corner RMS) are unambiguously
in the BC physics.

If this fails, the failure signature isolates the offending parameter:
  - Wrong peak amplitude  -> Hu, Cc, alpha_u, or Cp(alpha, T)
  - Wrong inflection time -> tau, beta
  - Wrong steepness       -> Ea (Arrhenius)
  - Wrong asymptote       -> alpha_u
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_rectangular, solve_hydration_2d


# ----------------------------------------------------------------------
# Tolerances (per passdown_kinetics_isolation_test.md "Proposed pass
# tolerances"). Numbers are in degrees Fahrenheit / hours.
# ----------------------------------------------------------------------

TOL_PEAK_F     = 1.0   # |delta T at t=168 hr|        (deg F)
TOL_RMS_F      = 1.0   # RMS(engine - CW) over [0,168] (deg F)
TOL_TIME_100_H = 2.0   # |t_engine(100F) - t_cw(100F)| (hr)
TOL_TIME_130_H = 2.0   # |t_engine(130F) - t_cw(130F)| (hr)

# ----------------------------------------------------------------------
# Conversions
# ----------------------------------------------------------------------

def f_to_c(t_f: float) -> float:
    return (t_f - 32.0) * 5.0 / 9.0

def c_to_f(t_c):
    return np.asarray(t_c) * 9.0 / 5.0 + 32.0


# ----------------------------------------------------------------------
# Result container
# ----------------------------------------------------------------------

@dataclass
class TestMetrics:
    peak_engine_F: float
    peak_cw_F: float
    peak_delta_F: float
    rms_F: float
    max_abs_delta_F: float
    t_engine_100F_h: float
    t_cw_100F_h: float
    delta_t_100F_h: float
    t_engine_130F_h: float
    t_cw_130F_h: float
    delta_t_130F_h: float
    spatial_uniformity_max_F: float  # max(T) - min(T) over field at t_final
    n_samples_compared: int

    @property
    def gates(self) -> dict:
        return {
            "peak_within_1F":      abs(self.peak_delta_F)    <= TOL_PEAK_F,
            "rms_within_1F":       self.rms_F                <= TOL_RMS_F,
            "t_to_100F_within_2h": abs(self.delta_t_100F_h)  <= TOL_TIME_100_H,
            "t_to_130F_within_2h": abs(self.delta_t_130F_h)  <= TOL_TIME_130_H,
        }

    @property
    def passed(self) -> bool:
        return all(self.gates.values())


# ----------------------------------------------------------------------
# Core test routine
# ----------------------------------------------------------------------

def run_engine_adiabatic(
    input_dat: str | Path,
    T_initial_F: float = 73.0,
    duration_hrs: float = 168.0,
    output_interval_s: float = 300.0,
    grid_n: int = 3,
):
    """
    Run the engine in adiabatic mode at uniform T0, return:
        t_hrs (n,), T_F (n,) for any cell, full T_field_C at final step
    """
    scenario = load_cw_scenario(str(input_dat),
                                weather_dat=None, cw_output_txt=None)
    mix = scenario.mix

    # Geometry is irrelevant in adiabatic mode (uniform field). Use the
    # smallest grid that satisfies build_grid_rectangular's nx, ny >= 3
    # constraint. Lx, Ly chosen at 1 m for numeric tameness; the dt_cfl
    # bound depends on (dx, dy_min) but we never approach it because
    # the inner step is clamped to output_interval_s anyway in adiabatic
    # mode (Q is the only RHS term, no diffusive transport across a
    # uniform field).
    grid = build_grid_rectangular(Lx_m=1.0, Ly_m=1.0, nx=grid_n, ny=grid_n)

    T0_C = f_to_c(T_initial_F)
    T_initial_C = np.full((grid.ny, grid.nx), T0_C, dtype=np.float64)

    result = solve_hydration_2d(
        grid, mix, T_initial_C,
        duration_s=duration_hrs * 3600.0,
        output_interval_s=output_interval_s,
        boundary_mode="adiabatic",
        # environment, construction not needed for adiabatic mode.
    )

    t_hrs = result.t_s / 3600.0
    # Any cell works in adiabatic mode -- pick the interior cell.
    j_pick, i_pick = grid.ny // 2, grid.nx // 2
    T_F = c_to_f(result.T_field_C[:, j_pick, i_pick])

    # Spatial uniformity check: max(T) - min(T) over the whole field at
    # t_final, in deg F. Should be ~0 to floating-point precision because
    # adiabatic mode + uniform IC = uniform field for all time.
    T_final = result.T_field_C[-1]
    spatial_spread_F = (T_final.max() - T_final.min()) * 9.0 / 5.0

    return {
        "t_hrs": t_hrs,
        "T_F": T_F,
        "spatial_spread_F": float(spatial_spread_F),
        "mix": mix,
        "n_inner_steps": result.n_inner_steps,
        "dt_inner_s": result.dt_inner_s,
    }


def load_cw_reference(csv_path: str | Path):
    """
    Load the CW adiabatic reference. Returns (t_hrs, T_F).

    Uses the T_center_F_adiabatic column (col index 1) -- the
    geometric centerline mid-depth of the 200x200x200 ft member, which
    is verified clean of BC contamination (see
    cw_adiabatic_reference_mix01_README.md).
    """
    arr = np.genfromtxt(str(csv_path), delimiter=",", skip_header=1)
    return arr[:, 0], arr[:, 1]


def compute_metrics(
    t_engine_h: np.ndarray, T_engine_F: np.ndarray,
    t_cw_h: np.ndarray, T_cw_F: np.ndarray,
    spatial_spread_F: float,
) -> TestMetrics:
    """
    Interpolate the engine onto the CW time grid and compute peak / RMS
    / time-to-threshold metrics. The CW grid is the reference, so we
    take (engine - CW) deltas at CW times.
    """
    T_engine_on_cw = np.interp(t_cw_h, t_engine_h, T_engine_F)
    delta = T_engine_on_cw - T_cw_F

    # Peak comparison at the final time (matches passdown criterion;
    # both curves are still monotonically rising at t=168 hr per the CW
    # README, so this is also the ~peak of each).
    peak_engine_F = float(T_engine_on_cw[-1])
    peak_cw_F     = float(T_cw_F[-1])

    rms = float(np.sqrt(np.mean(delta ** 2)))
    max_abs = float(np.max(np.abs(delta)))

    # Time to threshold: np.interp requires the x-data to be increasing,
    # which T_F is (both curves are monotone increasing over [0,168]).
    t_engine_100 = float(np.interp(100.0, T_engine_on_cw, t_cw_h))
    t_cw_100     = float(np.interp(100.0,        T_cw_F, t_cw_h))
    t_engine_130 = float(np.interp(130.0, T_engine_on_cw, t_cw_h))
    t_cw_130     = float(np.interp(130.0,        T_cw_F, t_cw_h))

    return TestMetrics(
        peak_engine_F=peak_engine_F,
        peak_cw_F=peak_cw_F,
        peak_delta_F=peak_engine_F - peak_cw_F,
        rms_F=rms,
        max_abs_delta_F=max_abs,
        t_engine_100F_h=t_engine_100,
        t_cw_100F_h=t_cw_100,
        delta_t_100F_h=t_engine_100 - t_cw_100,
        t_engine_130F_h=t_engine_130,
        t_cw_130F_h=t_cw_130,
        delta_t_130F_h=t_engine_130 - t_cw_130,
        spatial_uniformity_max_F=spatial_spread_F,
        n_samples_compared=int(len(t_cw_h)),
    )


def plot_comparison(
    t_cw_h: np.ndarray, T_cw_F: np.ndarray,
    t_engine_h: np.ndarray, T_engine_F: np.ndarray,
    metrics: TestMetrics,
    out_path: str | Path,
) -> None:
    """
    Two-panel plot:
      Top:  engine vs CW temperature trajectories (overlaid)
      Bot:  engine - CW residual at CW sample times
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    T_engine_on_cw = np.interp(t_cw_h, t_engine_h, T_engine_F)
    delta = T_engine_on_cw - T_cw_F

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax1.plot(t_cw_h, T_cw_F, "k-", linewidth=1.6, label="CW (T_center, adiabatic)")
    ax1.plot(t_engine_h, T_engine_F, "C0--", linewidth=1.4,
             label="Engine v3 (adiabatic)")
    ax1.set_ylabel("Temperature (deg F)")
    ax1.legend(loc="lower right")
    ax1.grid(True, alpha=0.3)
    title = (
        "Kinetics-isolation: CalcShore engine v3 vs CW (MIX-01, T0=73F)\n"
        f"peak delta = {metrics.peak_delta_F:+.3f} F   |   "
        f"RMS = {metrics.rms_F:.3f} F   |   "
        f"max |delta| = {metrics.max_abs_delta_F:.3f} F   |   "
        f"{'PASS' if metrics.passed else 'FAIL'}"
    )
    ax1.set_title(title)

    # Mark threshold-crossing points
    for thresh, color in [(100.0, "C2"), (130.0, "C3")]:
        ax1.axhline(thresh, color=color, alpha=0.25, linewidth=0.8,
                    linestyle=":")

    ax2.plot(t_cw_h, delta, "C3-", linewidth=1.0)
    ax2.axhline(0.0, color="k", linewidth=0.5)
    ax2.fill_between(t_cw_h, -TOL_RMS_F, +TOL_RMS_F,
                     color="green", alpha=0.08,
                     label=f"+/- {TOL_RMS_F:.1f} F band")
    ax2.set_xlabel("Time since placement (hours)")
    ax2.set_ylabel("Engine - CW (deg F)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=8)

    plt.tight_layout()
    plt.savefig(str(out_path), dpi=140)
    plt.close(fig)


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------

def format_report(metrics: TestMetrics, mix) -> str:
    lines = []
    lines.append("=" * 72)
    lines.append("CalcShore engine v3 -- kinetics isolation test")
    lines.append("MIX-01 adiabatic, T0 = 73 deg F, duration = 168 hr")
    lines.append("=" * 72)
    lines.append("")
    lines.append("Mix kinetics parameters consumed:")
    lines.append(f"  Ea       = {mix.activation_energy_J_mol:>10.2f}   J/mol")
    lines.append(f"  tau      = {mix.tau_hrs:>10.4f}   hr")
    lines.append(f"  beta     = {mix.beta:>10.4f}")
    lines.append(f"  alpha_u  = {mix.alpha_u:>10.5f}")
    lines.append(f"  Hu       = {mix.Hu_J_kg:>10.1f}   J/kg")
    lines.append(f"  Cc total = {mix.total_cementitious_lb_yd3:>10.1f}   lb/yd^3 "
                 f"({100.0 * (mix.total_cementitious_lb_yd3 - mix.cement_type_I_II_lb_yd3) / mix.total_cementitious_lb_yd3:.1f}% SCM)")
    lines.append("")
    lines.append("Spatial-uniformity check (adiabatic invariant):")
    lines.append(f"  max(T) - min(T) over field at t_final = "
                 f"{metrics.spatial_uniformity_max_F:.3e} deg F")
    lines.append(f"  (should be ~0 to floating-point precision)")
    lines.append("")
    lines.append(f"Comparison against CW reference ({metrics.n_samples_compared} samples):")
    lines.append("")
    lines.append(f"  Peak T at t=168 hr:")
    lines.append(f"    engine   = {metrics.peak_engine_F:>8.3f} deg F")
    lines.append(f"    CW       = {metrics.peak_cw_F:>8.3f} deg F")
    lines.append(f"    delta    = {metrics.peak_delta_F:>+8.3f} deg F")
    lines.append("")
    lines.append(f"  RMS (engine - CW) over [0, 168] hr:")
    lines.append(f"    rms      = {metrics.rms_F:>8.3f} deg F")
    lines.append(f"    max|delta|= {metrics.max_abs_delta_F:>8.3f} deg F")
    lines.append("")
    lines.append(f"  Time to T = 100 deg F:")
    lines.append(f"    engine   = {metrics.t_engine_100F_h:>8.3f} hr")
    lines.append(f"    CW       = {metrics.t_cw_100F_h:>8.3f} hr")
    lines.append(f"    delta    = {metrics.delta_t_100F_h:>+8.3f} hr")
    lines.append("")
    lines.append(f"  Time to T = 130 deg F:")
    lines.append(f"    engine   = {metrics.t_engine_130F_h:>8.3f} hr")
    lines.append(f"    CW       = {metrics.t_cw_130F_h:>8.3f} hr")
    lines.append(f"    delta    = {metrics.delta_t_130F_h:>+8.3f} hr")
    lines.append("")
    lines.append("Gate dispositions (S0-style):")
    lines.append("")
    g = metrics.gates
    fmts = [
        ("peak_within_1F",      f"|peak delta| <= {TOL_PEAK_F:.1f} F",
         f"|{metrics.peak_delta_F:+.3f}| F"),
        ("rms_within_1F",       f"RMS          <= {TOL_RMS_F:.1f} F",
         f"{metrics.rms_F:.3f} F"),
        ("t_to_100F_within_2h", f"|delta t_100|<= {TOL_TIME_100_H:.1f} h",
         f"|{metrics.delta_t_100F_h:+.3f}| h"),
        ("t_to_130F_within_2h", f"|delta t_130|<= {TOL_TIME_130_H:.1f} h",
         f"|{metrics.delta_t_130F_h:+.3f}| h"),
    ]
    for key, criterion, observed in fmts:
        status = "PASS" if g[key] else "FAIL"
        lines.append(f"  [{status}]  {criterion:<30}  observed: {observed}")
    lines.append("")
    lines.append("=" * 72)
    overall = "PASS" if metrics.passed else "FAIL"
    lines.append(f"OVERALL: {overall}")
    lines.append("=" * 72)
    if metrics.passed:
        lines.append("")
        lines.append("Kinetics ODE integration is verified clean. The validated-")
        lines.append("envelope residuals (Centerline RMS ~0.74 deg F, Corner RMS")
        lines.append("~2.0 deg F) are routed unambiguously to BC physics.")
    else:
        lines.append("")
        lines.append("Failure signature analysis (per passdown):")
        lines.append("  - Wrong peak amplitude  -> Hu, Cc, alpha_u, or Cp(alpha,T)")
        lines.append("  - Wrong inflection time -> tau, beta")
        lines.append("  - Wrong steepness       -> Ea (Arrhenius)")
        lines.append("  - Wrong asymptote       -> alpha_u")
        lines.append("Also check: input.dat parsing, T_ref=296.15 K, Imperial")
        lines.append("-> SI conversions, te integration scheme.")
    return "\n".join(lines)


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Kinetics-isolation validation for CalcShore engine v3.",
    )
    parser.add_argument(
        "--input-dat", default="input.dat",
        help="Path to MIX-01 CW input.dat (default: ./input.dat)",
    )
    parser.add_argument(
        "--cw-csv", default="cw_adiabatic_reference_mix01.csv",
        help="Path to CW adiabatic reference CSV.",
    )
    parser.add_argument(
        "--out-plot", default="kinetics_isolation_test.png",
        help="Output PNG path for the comparison plot.",
    )
    parser.add_argument(
        "--out-report", default="kinetics_isolation_report.txt",
        help="Output path for the text report.",
    )
    parser.add_argument(
        "--T0-F", type=float, default=73.0,
        help="Initial temperature in deg F (default 73, matching CW reference).",
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Skip the matplotlib plot (e.g., headless / no matplotlib).",
    )
    args = parser.parse_args(argv)

    print(f"[1/4] Loading mix from {args.input_dat} ...")
    print(f"[2/4] Running engine: 168 hr, 5-min sampling, "
          f"adiabatic, T0={args.T0_F} F ...")
    eng = run_engine_adiabatic(args.input_dat, T_initial_F=args.T0_F)
    print(f"      engine: {eng['n_inner_steps']} inner steps, "
          f"dt_inner = {eng['dt_inner_s']:.1f} s, "
          f"spatial spread at t_final = {eng['spatial_spread_F']:.2e} F")

    print(f"[3/4] Loading CW reference from {args.cw_csv} ...")
    t_cw_h, T_cw_F = load_cw_reference(args.cw_csv)
    print(f"      CW reference: {len(t_cw_h)} samples, "
          f"t in [{t_cw_h[0]:.2f}, {t_cw_h[-1]:.2f}] hr, "
          f"T in [{T_cw_F.min():.2f}, {T_cw_F.max():.2f}] F")

    print(f"[4/4] Computing metrics ...")
    metrics = compute_metrics(
        eng["t_hrs"], eng["T_F"], t_cw_h, T_cw_F, eng["spatial_spread_F"],
    )

    report = format_report(metrics, eng["mix"])
    print()
    print(report)

    Path(args.out_report).write_text(report + "\n", encoding="utf-8")
    print(f"\nReport written to {args.out_report}")

    if not args.no_plot:
        try:
            plot_comparison(
                t_cw_h, T_cw_F,
                eng["t_hrs"], eng["T_F"],
                metrics, args.out_plot,
            )
            print(f"Plot written to {args.out_plot}")
        except Exception as exc:  # pragma: no cover
            print(f"(plot skipped -- {type(exc).__name__}: {exc})")

    return 0 if metrics.passed else 1


if __name__ == "__main__":
    sys.exit(main())
