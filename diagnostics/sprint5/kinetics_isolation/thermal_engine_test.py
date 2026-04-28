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
5. Extract the engine's native centerline trajectory
   result.centerline_T_C[:, ny//2] and verify it matches any interior
   cell to floating-point precision (spatial-uniformity invariant).
6. Compare against T_center_F_adiabatic from
   cw_adiabatic_reference_mix01.csv:
   - Peak deviation at t = 168 hr  (tol  +/-1.0 deg F)
   - RMS over [0, 168] hr          (tol  <=1.0 deg F)
   - Time to T = 100 deg F         (tol  +/-2 hr)
   - Time to T = 130 deg F         (tol  +/-2 hr)
7. Save a 2-panel comparison plot (engine vs CW + residual).
8. Emit cross-section snapshot data as an npz consumable by
   visualize_xs_snapshots.py --input-npz.
9. Emit a pass/fail report; exit 0 on pass, 1 on fail.

Expected result: FAIL (engine runs ~+4.3 deg F hot vs CW). This is the
deliverable — it disambiguates v3 residuals as kinetics-driven, not
BC-driven. BC heat losses in the validated 14-mix envelope were masking
a +4.3 deg F kinetics overshoot. See findings.md for full analysis.

Failure signature:
  - Timing of early rise is correct (t_100F within 2 hr) -> tau, beta,
    Ea not the primary culprit.
  - Cumulative heat is too high (+4.3 deg F at t=168 hr) -> Hu, Cc,
    alpha_u, or Cp(alpha, T).

Usage (from repo root):
    python diagnostics/sprint5/kinetics_isolation/thermal_engine_test.py

All outputs default to the directory alongside this script.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import numpy as np

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))

from cw_scenario_loader import load_cw_scenario  # noqa: E402
from thermal_engine_2d import build_grid_rectangular, solve_hydration_2d  # noqa: E402

_HERE = Path(__file__).parent


# ---------------------------------------------------------------------------
# Tolerances (per passdown_kinetics_isolation_test.md "Proposed pass
# tolerances"). Numbers in degrees Fahrenheit / hours.
# ---------------------------------------------------------------------------

TOL_PEAK_F     = 1.0   # |delta T at t=168 hr|          (deg F)
TOL_RMS_F      = 1.0   # RMS(engine - CW) over [0,168]  (deg F)
TOL_TIME_100_H = 2.0   # |t_engine(100F) - t_cw(100F)|  (hr)
TOL_TIME_130_H = 2.0   # |t_engine(130F) - t_cw(130F)|  (hr)


# ---------------------------------------------------------------------------
# Conversions
# ---------------------------------------------------------------------------

def _f_to_c(t_f: float) -> float:
    return (t_f - 32.0) * 5.0 / 9.0


def _c_to_f(t_c):
    return np.asarray(t_c) * 9.0 / 5.0 + 32.0


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

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
    spatial_uniformity_max_F: float   # max(T) - min(T) over field at t_final
    centerline_vs_interior_max_F: float  # max |centerline - interior cell|
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


# ---------------------------------------------------------------------------
# Core engine run
# ---------------------------------------------------------------------------

def run_engine_adiabatic(
    input_dat: str | Path,
    T_initial_F: float = 73.0,
    duration_hrs: float = 168.0,
    output_interval_s: float = 300.0,
    grid_n: int = 3,
):
    """Run the engine in adiabatic mode at uniform T0.

    Returns a dict with:
        t_hrs (n_out,), T_F (n_out,) via centerline accessor,
        T_field_C (n_out, ny, nx), grid, and solver metadata.
    """
    scenario = load_cw_scenario(str(input_dat),
                                weather_dat=None, cw_output_txt=None)
    mix = scenario.mix

    # Geometry is irrelevant in adiabatic mode (uniform field). Use the
    # smallest grid that satisfies build_grid_rectangular's nx, ny >= 3
    # constraint.
    grid = build_grid_rectangular(Lx_m=1.0, Ly_m=1.0, nx=grid_n, ny=grid_n)

    T0_C = _f_to_c(T_initial_F)
    T_initial_C = np.full((grid.ny, grid.nx), T0_C, dtype=np.float64)

    result = solve_hydration_2d(
        grid, mix, T_initial_C,
        duration_s=duration_hrs * 3600.0,
        output_interval_s=output_interval_s,
        boundary_mode="adiabatic",
    )

    t_hrs = result.t_s / 3600.0

    # Primary trajectory: engine's native centerline accessor at mid-depth.
    # shape of centerline_T_C is (n_out, ny); we take ny//2 for the mid row.
    T_center_C = result.centerline_T_C[:, grid.ny // 2]
    T_F = _c_to_f(T_center_C)

    # Sanity check: centerline must match any interior cell to floating-point
    # precision (adiabatic + uniform IC = spatially uniform field always).
    j_int, i_int = grid.ny // 2, grid.nx // 2
    T_interior_C = result.T_field_C[:, j_int, i_int]
    cl_vs_int_max = float(np.max(np.abs(T_center_C - T_interior_C)))
    if not np.allclose(T_center_C, T_interior_C, atol=1e-12):
        raise RuntimeError(
            f"centerline_T_C diverges from interior cell by {cl_vs_int_max:.3e} C "
            "— adiabatic spatial-uniformity invariant violated; likely an engine bug."
        )

    # Spatial spread over the full field at t_final.
    T_final = result.T_field_C[-1]
    spatial_spread_F = float((T_final.max() - T_final.min()) * 9.0 / 5.0)

    return {
        "t_hrs": t_hrs,
        "T_F": T_F,
        "T_field_C": result.T_field_C,
        "spatial_spread_F": spatial_spread_F,
        "centerline_vs_interior_max_F": float(cl_vs_int_max * 9.0 / 5.0),
        "mix": mix,
        "grid": grid,
        "n_inner_steps": result.n_inner_steps,
        "dt_inner_s": result.dt_inner_s,
    }


# ---------------------------------------------------------------------------
# CW reference loader
# ---------------------------------------------------------------------------

def load_cw_reference(csv_path: str | Path):
    """Load the CW adiabatic reference. Returns (t_hrs, T_F).

    Uses column index 1 (T_center_F_adiabatic) — the geometric center of
    the 200x200x200 ft member, verified clean of BC contamination.
    See cw_adiabatic_reference_mix01_README.md.
    """
    arr = np.genfromtxt(str(csv_path), delimiter=",", skip_header=1)
    return arr[:, 0], arr[:, 1]


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def compute_metrics(
    t_engine_h: np.ndarray, T_engine_F: np.ndarray,
    t_cw_h: np.ndarray, T_cw_F: np.ndarray,
    spatial_spread_F: float,
    centerline_vs_interior_max_F: float,
) -> TestMetrics:
    """Interpolate the engine onto the CW time grid and compute metrics.

    Deltas are (engine - CW) at CW sample times.
    """
    T_engine_on_cw = np.interp(t_cw_h, t_engine_h, T_engine_F)
    delta = T_engine_on_cw - T_cw_F

    # Peak comparison at t=168 hr (final sample). Both curves are still
    # monotonically rising at t=168 hr per the CW README.
    peak_engine_F = float(T_engine_on_cw[-1])
    peak_cw_F     = float(T_cw_F[-1])

    rms     = float(np.sqrt(np.mean(delta ** 2)))
    max_abs = float(np.max(np.abs(delta)))

    # Time to threshold: np.interp requires increasing x-data; both curves
    # are monotone increasing over [0, 168] hr.
    t_engine_100 = float(np.interp(100.0, T_engine_on_cw, t_cw_h))
    t_cw_100     = float(np.interp(100.0, T_cw_F,         t_cw_h))
    t_engine_130 = float(np.interp(130.0, T_engine_on_cw, t_cw_h))
    t_cw_130     = float(np.interp(130.0, T_cw_F,         t_cw_h))

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
        centerline_vs_interior_max_F=centerline_vs_interior_max_F,
        n_samples_compared=int(len(t_cw_h)),
    )


# ---------------------------------------------------------------------------
# Cross-section snapshot emitter
# ---------------------------------------------------------------------------

def save_xs_snapshots(
    t_hrs: np.ndarray,
    T_field_C: np.ndarray,
    grid,
    out_path: str | Path,
) -> None:
    """Save T_field_C as an npz consumable by visualize_xs_snapshots.py.

    Schema matches CWValidationSeries: T_field_F (n_out, ny, nx) in
    degrees F; time_hrs (n_out,) in hours; depths_m (ny,) ascending;
    widths_m (nx,) descending to mirror the CW convention where index 0
    is the point furthest from the form face (centerline side).

    Note: in adiabatic mode with uniform IC the field is spatially
    uniform at every timestep, so all cross-section snapshots are
    single-colour tiles. This is physically correct.
    """
    T_field_F = _c_to_f(T_field_C)          # (n_out, ny, nx)
    depths_m  = grid.y                       # ascending 0 -> Ly
    # Flip width axis to descending (CW convention: index 0 = "far from corner")
    widths_m = grid.x[::-1]                  # descending Lx -> 0
    T_field_F = T_field_F[:, :, ::-1]        # match flipped width order

    np.savez(
        str(out_path),
        T_field_F=T_field_F.astype(np.float64),
        time_hrs=t_hrs.astype(np.float64),
        depths_m=depths_m.astype(np.float64),
        widths_m=widths_m.astype(np.float64),
    )


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_comparison(
    t_cw_h: np.ndarray, T_cw_F: np.ndarray,
    t_engine_h: np.ndarray, T_engine_F: np.ndarray,
    metrics: TestMetrics,
    out_path: str | Path,
) -> None:
    """Two-panel plot: temperature trajectories (top) and residual (bottom)."""
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

    for thresh, color in [(100.0, "C2"), (130.0, "C3")]:
        ax1.axhline(thresh, color=color, alpha=0.25, linewidth=0.8, linestyle=":")

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


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

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
    lines.append(f"  centerline vs interior cell max |diff| = "
                 f"{metrics.centerline_vs_interior_max_F:.3e} deg F")
    lines.append(f"  (both should be ~0 to floating-point precision)")
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
        lines.append("  - t_100F within 2 hr -> tau, beta, Ea not the primary issue")
        lines.append("  - Cumulative heat too high -> Hu, Cc, alpha_u, or Cp(alpha,T)")
        lines.append("  See findings.md for prime suspects and routing implications.")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Kinetics-isolation validation for CalcShore engine v3.",
    )
    parser.add_argument(
        "--input-dat",
        default=str(_HERE / "input.dat"),
        help="Path to MIX-01 CW input.dat.",
    )
    parser.add_argument(
        "--cw-csv",
        default=str(_HERE / "cw_adiabatic_reference_mix01.csv"),
        help="Path to CW adiabatic reference CSV.",
    )
    parser.add_argument(
        "--out-plot",
        default=str(_HERE / "kinetics_isolation_test.png"),
        help="Output PNG path for the comparison plot.",
    )
    parser.add_argument(
        "--out-report",
        default=str(_HERE / "kinetics_isolation_report.txt"),
        help="Output path for the text report.",
    )
    parser.add_argument(
        "--out-xs-npz",
        default=str(_HERE / "kinetics_isolation_xs_snapshots.npz"),
        help="Output npz path for cross-section snapshots (for visualize_xs_snapshots.py).",
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

    print(f"[1/5] Loading mix from {args.input_dat} ...")
    print(f"[2/5] Running engine: 168 hr, 5-min sampling, "
          f"adiabatic, T0={args.T0_F} F ...")
    eng = run_engine_adiabatic(args.input_dat, T_initial_F=args.T0_F)
    print(f"      engine: {eng['n_inner_steps']} inner steps, "
          f"dt_inner = {eng['dt_inner_s']:.1f} s, "
          f"spatial spread at t_final = {eng['spatial_spread_F']:.2e} F, "
          f"centerline vs interior = {eng['centerline_vs_interior_max_F']:.2e} F")

    print(f"[3/5] Loading CW reference from {args.cw_csv} ...")
    t_cw_h, T_cw_F = load_cw_reference(args.cw_csv)
    print(f"      CW reference: {len(t_cw_h)} samples, "
          f"t in [{t_cw_h[0]:.2f}, {t_cw_h[-1]:.2f}] hr, "
          f"T in [{T_cw_F.min():.2f}, {T_cw_F.max():.2f}] F")

    print(f"[4/5] Computing metrics ...")
    metrics = compute_metrics(
        eng["t_hrs"], eng["T_F"], t_cw_h, T_cw_F,
        eng["spatial_spread_F"], eng["centerline_vs_interior_max_F"],
    )

    report = format_report(metrics, eng["mix"])
    print()
    print(report)

    Path(args.out_report).write_text(report + "\n", encoding="utf-8")
    print(f"\nReport written to   {args.out_report}")

    print(f"[5/5] Saving cross-section snapshot data ...")
    save_xs_snapshots(
        eng["t_hrs"], eng["T_field_C"], eng["grid"], args.out_xs_npz,
    )
    print(f"XS snapshots written to  {args.out_xs_npz}")

    if not args.no_plot:
        try:
            plot_comparison(
                t_cw_h, T_cw_F,
                eng["t_hrs"], eng["T_F"],
                metrics, args.out_plot,
            )
            print(f"Plot written to          {args.out_plot}")
        except Exception as exc:  # pragma: no cover
            print(f"(plot skipped -- {type(exc).__name__}: {exc})")

    return 0 if metrics.passed else 1


if __name__ == "__main__":
    sys.exit(main())
