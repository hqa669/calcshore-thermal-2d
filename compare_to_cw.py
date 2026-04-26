#!/usr/bin/env python3
"""
CW thermal validation harness — single-scenario CLI and multi-scenario library.

Usage (CLI):
    python compare_to_cw.py <directory>

<directory> must contain input.dat, weather.dat, output.txt.

Exit 0 if ALL metrics PASS. Exit 1 if ANY FAIL.

Library use: from compare_to_cw import run_one, print_gate_table
"""

import sys
import os
import time
import argparse
import numpy as np

from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import (
    build_grid_half_mat,
    solve_hydration_2d,
)

# Sprint-0 acceptance tolerances
TOL_PEAK_MAX_F    = 1.0   # °F
TOL_PEAK_GRAD_F   = 2.0   # °F
TOL_FIELD_RMS_F   = 2.0   # °F
TOL_CENTER_RMS_F  = 1.0   # °F
TOL_CORNER_RMS_F  = 3.0   # °F

# S1 aspirational tolerances — documentary only; never drive PASS/FAIL or exit code
S1_TOL_PEAK_MAX_F    = 0.5
S1_TOL_PEAK_GRAD_F   = 1.0
S1_TOL_FIELD_RMS_F   = 1.0
S1_TOL_CENTER_RMS_F  = 0.5
S1_TOL_CORNER_RMS_F  = 2.0   # Sprint 3 tightened; Sprint 4 multi-climate target


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("scenario_dir", help="Directory containing input.dat, weather.dat, output.txt")
    return p.parse_args()


def check_files(d):
    missing = []
    for fn in ("input.dat", "weather.dat", "output.txt"):
        if not os.path.isfile(os.path.join(d, fn)):
            missing.append(fn)
    if missing:
        print(f"ERROR: missing files in {d}: {missing}", file=sys.stderr)
        sys.exit(2)


def rms(a, b):
    """RMS of (a - b)."""
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(diff**2)))


def pass_fail(ok):
    return "PASS" if ok else "FAIL"


def run_one(scenario_dir: str, *, png_path: str | None = None) -> dict:
    """Run engine vs CW comparison for a single mix scenario.

    Loads CW input.dat + weather.dat + output.txt from scenario_dir, runs the
    CalcShore thermal engine on the same construction parameters, and returns a
    structured comparison dict containing peak temperatures, gradients, RMS
    errors, and S0 / S1-aspire gate dispositions.

    Parameters
    ----------
    scenario_dir : str
        Path to a scenario directory containing CW input.dat, weather.dat, and
        output.txt. The directory name is used as the mix identifier (e.g.,
        "MIX-01"). If output.txt is absent, the scenario is treated as a
        sentinel skip and a result dict with ``skipped=True`` is returned — no
        exception is raised by the library path.
    png_path : str or None, optional
        Keyword-only. If provided, write a 4×3 comparison plot grid to this
        path. Default None means no plot is produced.

    Returns
    -------
    result : dict
        Structured comparison result. Top-level keys:

        scenario_dir : str
            Input path, echoed back (trailing slash stripped).
        scenario_name : str
            Basename of scenario_dir (e.g., "MIX-01").
        skipped : bool
            True if output.txt was missing; only skip_reason is also populated.
        skip_reason : str or None
            ``"no_cw_output"`` when skipped=True, else None.
        meta : dict
            ``scm_pct`` (float), ``placement_temp_F`` (float),
            ``total_cementitious_lb_yd3`` (float).
        engine : dict
            Engine extractions: ``peak_max_F``, ``peak_max_hr``,
            ``peak_grad_F``, ``peak_grad_hr`` (all float); ``t_wall_s``
            (solver wall time); ``n_steps``, ``n_nodes``,
            ``n_output_samples`` (int).
        cw : dict
            CW reference values: ``peak_max_F``, ``peak_max_hr``,
            ``peak_grad_F``, ``peak_grad_hr``.
        deltas : dict
            Engine minus CW: ``peak_max_F``, ``peak_grad_F``.
        rms : dict
            RMS error vs CW over the evaluation window: ``field_F``,
            ``center_F``, ``corner_F``.
        s0 : dict
            Per-gate S0 pass/fail booleans. Keys (unsuffixed): ``peak_max``,
            ``peak_grad``, ``field``, ``center``, ``corner``.
        s0_pass_count : int
        s0_pass_total : int
        s0_overall : bool
            True when all 5 S0 gates pass.
        s1_aspire : dict
            Per-gate S1-aspire booleans; same 5 keys as ``s0``.
        rms_window_hr : tuple
            ``(start_hr, end_hr)`` of the RMS evaluation window.
        rms_n_samples : int
            Number of output samples in the RMS window.

    Notes
    -----
    The engine v3 validated envelope is steel-form, half-mat geometry
    (40×60×8 ft), Austin TX summer, 60°F placement, 39.1% SCM Reference
    mixes. See docs/engine_v3_release_notes.md for the full validated
    envelope, known residuals, and out-of-envelope conditions.

    Examples
    --------
    >>> result = run_one("validation/cw_exports/MIX-01")
    >>> result["s0_overall"]
    True
    >>> result["deltas"]["peak_max_F"]
    -0.29
    """
    d = scenario_dir.rstrip("/")

    # Sentinel: missing output.txt → skip without raising, no engine run.
    # CLI path uses check_files() for hard errors; library path uses sentinel.
    if not os.path.isfile(os.path.join(d, "output.txt")):
        return {
            "scenario_dir": d,
            "skipped": True,
            "skip_reason": "no_cw_output",
            "scenario_name": os.path.basename(d),
        }

    # ------------------------------------------------------------------ #
    # 1. Load scenario
    # ------------------------------------------------------------------ #
    scn = load_cw_scenario(
        os.path.join(d, "input.dat"),
        os.path.join(d, "weather.dat"),
        os.path.join(d, "output.txt"),
    )

    # ------------------------------------------------------------------ #
    # 2. Build grid + initial conditions
    # ------------------------------------------------------------------ #
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)

    # ------------------------------------------------------------------ #
    # 3. Run solver
    # ------------------------------------------------------------------ #
    t_wall_start = time.perf_counter()
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )
    t_wall_s = time.perf_counter() - t_wall_start

    # ------------------------------------------------------------------ #
    # 4. Engine extractions (°F)
    # ------------------------------------------------------------------ #
    jslice, islice = grid.concrete_slice()
    T_conc_C = result.T_field_C[:, jslice, islice]
    T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
    t_hrs = result.t_s / 3600.0

    # Peak max T
    peak_flat = int(np.argmax(T_conc_F))
    peak_idx_3d = np.unravel_index(peak_flat, T_conc_F.shape)
    engine_peak_max_F  = float(T_conc_F.max())
    engine_peak_max_hr = float(t_hrs[peak_idx_3d[0]])

    # Peak gradient (max - min over concrete field at each time)
    grad_series = T_conc_F.max(axis=(1, 2)) - T_conc_F.min(axis=(1, 2))
    engine_peak_grad_F  = float(grad_series.max())
    engine_peak_grad_hr = float(t_hrs[int(np.argmax(grad_series))])

    # Centerline core: mid-depth concrete on centerline (i = nx-1)
    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_center_C = result.T_field_C[:, iy_mid, grid.nx - 1]
    eng_center_F = eng_center_C * 9.0 / 5.0 + 32.0

    # Top-corner surface: outer corner of concrete top surface
    iy_top = grid.iy_concrete_start
    ix_corner = grid.ix_concrete_start
    eng_corner_C = result.T_field_C[:, iy_top, ix_corner]
    eng_corner_F = eng_corner_C * 9.0 / 5.0 + 32.0

    # Ambient from solver diagnostic
    eng_amb_C = result.T_amb_C_history
    eng_amb_F = eng_amb_C * 9.0 / 5.0 + 32.0 if eng_amb_C is not None else None

    # ------------------------------------------------------------------ #
    # 5. CW extractions
    # ------------------------------------------------------------------ #
    val = scn.cw_validation
    cw_peak_max_idx = int(np.argmax(val.T_max_xs_F))
    cw_peak_max_F  = float(val.T_max_xs_F.max())
    cw_peak_max_hr = float(val.time_hrs[cw_peak_max_idx])
    cw_peak_grad_F  = float(val.T_diff_xs_F.max())
    cw_peak_grad_hr = float(val.time_hrs[int(np.argmax(val.T_diff_xs_F))])

    # Centerline and corner from T_field_F if available
    field_rms_F = None
    cw_center_F = None
    cw_corner_F = None
    center_rms_F = None
    corner_rms_F = None
    n_samples_rms = None

    # ------------------------------------------------------------------ #
    # VALIDATION WINDOW (PR 8.5)
    # ------------------------------------------------------------------ #
    # RMS metrics are computed over t ∈ [T_START_RMS_HR, 168]hr.
    # Rationale: the first ~48hr is dominated by the hydration-rise
    # transient, where engine and CW diverge in curve SHAPE (not steady-
    # state physics). Including the transient conflates two independent
    # error sources — boundary physics (Sprint 2 scope) and hydration
    # curve shape (Sprint 3 scope). Sprint 2 S0 gates validate boundary
    # physics; evaluating on the steady-state window isolates them.
    #
    # Chosen at 48hr: past the adiabatic peak for MIX-01 (~40-60hr), past
    # the first full diurnal cycle, retains 71% of CW output samples
    # (1441/2016). Diagnostic verify_pr8_floor_v2.py showed Corner RMS
    # stable at ~2.2-2.4°F across [48,168] through [96,168] windows;
    # evaluation is insensitive to the exact cutoff past 48hr.
    #
    # Peak Max T and Peak Gradient remain evaluated on the full 168hr
    # window — they are single-point metrics (max of series), not
    # integrated error metrics, so window-masking is not meaningful for
    # them. The peak in both engine and CW always occurs well after 48hr.
    T_START_RMS_HR: float = 48.0    # hours; start of RMS integration window

    if val.T_field_F is not None:
        # CW field: (n_cw_t, nD, nW); width index 0 = centerline, -1 = corner; depth 0 = top
        cw_t_s = val.time_hrs * 3600.0
        n_cw_t, n_cw_d, n_cw_w = val.T_field_F.shape

        # CW centerline at mid-depth (nD//2, width index 0)
        cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]

        # CW top corner (depth 0, width index -1)
        cw_corner_F = val.T_field_F[:, 0, -1]

        # Interpolate engine series to CW times
        eng_center_interp = np.interp(cw_t_s, result.t_s, eng_center_F)
        eng_corner_interp = np.interp(cw_t_s, result.t_s, eng_corner_F)

        # Apply steady-state window to all integrated RMS metrics.
        rms_mask = val.time_hrs >= T_START_RMS_HR
        n_samples_rms = int(rms_mask.sum())

        center_rms_F = rms(eng_center_interp[rms_mask], cw_center_F[rms_mask])
        corner_rms_F = rms(eng_corner_interp[rms_mask], cw_corner_F[rms_mask])

        # Field-wide: compare centerline column of engine to CW centerline across all depths
        n_conc_y = grid.iy_concrete_end - grid.iy_concrete_start + 1
        sq_errs = []
        for ti in range(n_cw_t):
            if val.time_hrs[ti] < T_START_RMS_HR:
                continue   # skip transient samples
            idx = int(np.searchsorted(result.t_s, cw_t_s[ti]))
            idx = min(idx, len(result.t_s) - 1)
            T_eng_col_C = result.T_field_C[idx, jslice, grid.nx - 1]  # centerline col
            T_eng_col_F = T_eng_col_C * 9.0 / 5.0 + 32.0
            n_cmp = min(len(T_eng_col_F), n_cw_d)
            sq_errs.extend((T_eng_col_F[:n_cmp] - val.T_field_F[ti, :n_cmp, 0]).tolist())
        field_rms_F = float(np.sqrt(np.mean(np.array(sq_errs) ** 2)))

    # ------------------------------------------------------------------ #
    # 6. Metrics and tolerances
    # ------------------------------------------------------------------ #
    peak_max_delta = engine_peak_max_F - cw_peak_max_F
    peak_grad_delta = engine_peak_grad_F - cw_peak_grad_F

    pass_peak_max  = abs(peak_max_delta) < TOL_PEAK_MAX_F
    pass_peak_grad = abs(peak_grad_delta) < TOL_PEAK_GRAD_F
    pass_field     = (field_rms_F  < TOL_FIELD_RMS_F)  if field_rms_F  is not None else None
    pass_center    = (center_rms_F < TOL_CENTER_RMS_F) if center_rms_F is not None else None
    pass_corner    = (corner_rms_F < TOL_CORNER_RMS_F) if corner_rms_F is not None else None

    s0_non_none = {k: v for k, v in {
        "peak_max": pass_peak_max, "peak_grad": pass_peak_grad,
        "field": pass_field, "center": pass_center, "corner": pass_corner,
    }.items() if v is not None}
    s0_overall = all(s0_non_none.values())

    # S1 aspirational (documentary only — never affects exit code)
    s1_peak_max  = abs(peak_max_delta) < S1_TOL_PEAK_MAX_F
    s1_peak_grad = abs(peak_grad_delta) < S1_TOL_PEAK_GRAD_F
    s1_field     = (field_rms_F  < S1_TOL_FIELD_RMS_F)  if field_rms_F  is not None else False
    s1_center    = (center_rms_F < S1_TOL_CENTER_RMS_F) if center_rms_F is not None else False
    s1_corner    = (corner_rms_F < S1_TOL_CORNER_RMS_F) if corner_rms_F is not None else False

    # ------------------------------------------------------------------ #
    # 7. Scenario metadata (SCM% for run_all.py markdown; computed from
    #    already-loaded scn so no extra I/O)
    # ------------------------------------------------------------------ #
    m = scn.mix
    scm_pct = 100.0 * (
        m.fly_ash_F_lb_yd3 + m.fly_ash_C_lb_yd3
        + m.ggbfs_lb_yd3 + m.silica_fume_lb_yd3
    ) / m.total_cementitious_lb_yd3

    scenario_name = (
        f"{os.path.basename(d)}"
        f" {scn.environment.location}"
        f" {scn.construction.placement_date}"
    )

    # ------------------------------------------------------------------ #
    # 8. Optional plot
    # ------------------------------------------------------------------ #
    if png_path is not None:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            _make_plot(t_hrs, engine_peak_max_F, eng_center_F, eng_corner_F, eng_amb_F,
                       val, cw_center_F, cw_corner_F, grad_series, d,
                       result=result, grid=grid, out_path=png_path,
                       soil_lag_hrs=scn.construction.soil_lag_hrs,
                       soil_damping=scn.construction.soil_damping)
        except ImportError:
            pass  # silently skip if matplotlib not installed
        except Exception:
            pass  # don't crash on plot errors

    return {
        "scenario_dir": d,
        "skipped": False,
        "skip_reason": None,
        "scenario_name": scenario_name,
        "meta": {
            "scm_pct": scm_pct,
            "placement_temp_F": scn.construction.placement_temp_F,
            "total_cementitious_lb_yd3": m.total_cementitious_lb_yd3,
        },
        "engine": {
            "peak_max_F":       engine_peak_max_F,
            "peak_max_hr":      engine_peak_max_hr,
            "peak_grad_F":      engine_peak_grad_F,
            "peak_grad_hr":     engine_peak_grad_hr,
            "t_wall_s":         t_wall_s,
            "n_steps":          result.n_inner_steps,
            "n_nodes":          grid.nx * grid.ny,
            "n_output_samples": result.n_output_samples,
        },
        "cw": {
            "peak_max_F":  cw_peak_max_F,
            "peak_max_hr": cw_peak_max_hr,
            "peak_grad_F": cw_peak_grad_F,
            "peak_grad_hr": cw_peak_grad_hr,
        },
        "deltas": {
            "peak_max_F":  peak_max_delta,
            "peak_grad_F": peak_grad_delta,
        },
        "rms": {
            "field_F":  field_rms_F,
            "center_F": center_rms_F,
            "corner_F": corner_rms_F,
        },
        "s0": {
            "peak_max": pass_peak_max,
            "peak_grad": pass_peak_grad,
            "field":    pass_field,
            "center":   pass_center,
            "corner":   pass_corner,
        },
        "s0_pass_count": sum(s0_non_none.values()),
        "s0_pass_total": len(s0_non_none),
        "s0_overall":    s0_overall,
        "s1_aspire": {
            "peak_max":  s1_peak_max,
            "peak_grad": s1_peak_grad,
            "field":     s1_field,
            "center":    s1_center,
            "corner":    s1_corner,
        },
        "rms_window_hr":  (T_START_RMS_HR, 168.0),
        "rms_n_samples":  n_samples_rms,
    }


def print_gate_table(result: dict, scenario_name: str | None = None) -> None:
    """Print the formatted gate table to stdout. Handles skipped sentinel."""
    name = scenario_name if scenario_name is not None else result.get("scenario_name", "")

    if result.get("skipped"):
        print(f"  {name}  — skipped ({result.get('skip_reason', '')})")
        return

    def s1_mark(ok): return "✓" if ok else "✗"

    eng    = result["engine"]
    cw     = result["cw"]
    deltas = result["deltas"]
    rms_d  = result["rms"]
    s0     = result["s0"]
    s1     = result["s1_aspire"]
    s1_count = sum(s1.values())

    field_rms_F  = rms_d["field_F"]
    center_rms_F = rms_d["center_F"]
    corner_rms_F = rms_d["corner_F"]

    T_START_RMS_HR = result["rms_window_hr"][0]
    n_samples_rms  = result["rms_n_samples"]

    print()
    print(f"  {name}  ({168} hr run, {eng['n_steps']} timesteps, {eng['n_nodes']} nodes)")
    if field_rms_F is not None:
        print(f"  (RMS metrics evaluated on t ∈ [{T_START_RMS_HR:.0f}, 168]hr = {n_samples_rms} samples)")
    print(f"  Peak Max T:     Engine {eng['peak_max_F']:6.1f}°F @ {eng['peak_max_hr']:5.1f} hr | CW {cw['peak_max_F']:6.1f}°F @ {cw['peak_max_hr']:.1f} hr")
    print(f"                  Δ = {deltas['peak_max_F']:+.1f}°F   [{pass_fail(s0['peak_max'])}] (S0 ±{TOL_PEAK_MAX_F:.1f}°F)  [S1-aspire ±{S1_TOL_PEAK_MAX_F:.1f}°F: {s1_mark(s1['peak_max'])}]")
    print(f"  Peak Gradient:  Engine {eng['peak_grad_F']:5.1f}°F @ {eng['peak_grad_hr']:5.1f} hr | CW {cw['peak_grad_F']:5.1f}°F @ {cw['peak_grad_hr']:.1f} hr")
    print(f"                  Δ = {deltas['peak_grad_F']:+.1f}°F   [{pass_fail(s0['peak_grad'])}] (S0 ±{TOL_PEAK_GRAD_F:.1f}°F)  [S1-aspire ±{S1_TOL_PEAK_GRAD_F:.1f}°F: {s1_mark(s1['peak_grad'])}]")

    if field_rms_F is not None:
        print(f"  Field-wide RMS: {field_rms_F:.2f}°F        [{pass_fail(s0['field'])}] (S0 tol {TOL_FIELD_RMS_F:.1f}°F)    [S1-aspire {S1_TOL_FIELD_RMS_F:.1f}°F: {s1_mark(s1['field'])}]")
    else:
        print(f"  Field-wide RMS: N/A (CW T_field_F not in fixture)")

    if center_rms_F is not None:
        print(f"  Centerline RMS: {center_rms_F:.2f}°F        [{pass_fail(s0['center'])}] (S0 tol {TOL_CENTER_RMS_F:.1f}°F)    [S1-aspire {S1_TOL_CENTER_RMS_F:.1f}°F: {s1_mark(s1['center'])}]")
    else:
        print(f"  Centerline RMS: N/A (CW T_field_F not in fixture)")

    if corner_rms_F is not None:
        print(f"  Corner RMS:     {corner_rms_F:.2f}°F        [{pass_fail(s0['corner'])}] (S0 tol {TOL_CORNER_RMS_F:.1f}°F)    [S1-aspire {S1_TOL_CORNER_RMS_F:.1f}°F: {s1_mark(s1['corner'])}]")
    else:
        print(f"  Corner RMS:     N/A (CW T_field_F not in fixture)")

    print(f"  Runtime:        {eng['t_wall_s']:.1f} s")
    overall_str = "[PASS]" if result["s0_overall"] else "[FAIL]"
    print(f"  OVERALL:        {overall_str}    S1-aspire: {s1_count}/5 metrics meet")
    print()


def main():
    args = parse_args()
    d = args.scenario_dir.rstrip("/")
    check_files(d)
    png_path = f"cw_comparison_{os.path.basename(d)}.png"
    result = run_one(d, png_path=png_path)
    print_gate_table(result)
    sys.exit(0 if result.get("s0_overall", False) else 1)


def _make_plot(t_hrs, engine_peak_max_F, eng_center_F, eng_corner_F, eng_amb_F,
               val, cw_center_F, cw_corner_F, grad_series, scenario_dir,
               result=None, grid=None, out_path: str = "cw_comparison.png",
               soil_lag_hrs=None, soil_damping=None):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 3, figsize=(15, 16))
    fig.suptitle("CalcShore Engine vs CW MIX-01 Austin", fontsize=13)

    # Row 1 — (a) Peak T center core
    ax = axes[0, 0]
    ax.plot(t_hrs, eng_center_F, label="Engine center", color="tab:blue")
    if cw_center_F is not None and val.time_hrs is not None:
        ax.plot(val.time_hrs, cw_center_F, "--", label="CW center", color="tab:orange")
    ax.axhline(engine_peak_max_F, color="tab:blue", alpha=0.3, linestyle=":")
    ax.set_xlabel("t (hr)"); ax.set_ylabel("T (°F)"); ax.set_title("(a) Peak T — centerline core")
    ax.legend(fontsize=8)

    # Row 1 — (b) Min T cross-section
    ax = axes[0, 1]
    if hasattr(val, "T_min_xs_F") and val.T_min_xs_F is not None:
        ax.plot(val.time_hrs, val.T_min_xs_F, label="CW min xs", color="tab:orange")
    ax.set_xlabel("t (hr)"); ax.set_ylabel("T (°F)"); ax.set_title("(b) Min T cross-section")
    ax.legend(fontsize=8)

    # Row 1 — (c) Peak Gradient
    ax = axes[0, 2]
    ax.plot(t_hrs, grad_series, label="Engine grad", color="tab:blue")
    if val.T_diff_xs_F is not None:
        ax.plot(val.time_hrs, val.T_diff_xs_F, "--", label="CW grad", color="tab:orange")
    ax.axhline(39.3, linestyle=":", color="gray", alpha=0.5)
    ax.set_xlabel("t (hr)"); ax.set_ylabel("ΔT (°F)"); ax.set_title("(c) Peak gradient")
    ax.legend(fontsize=8)

    # Row 2 — (d) Centerline core T(t)
    ax = axes[1, 0]
    ax.plot(t_hrs, eng_center_F, label="Engine center", color="tab:blue")
    if cw_center_F is not None:
        ax.plot(val.time_hrs, cw_center_F, "--", label="CW center", color="tab:orange")
    ax.set_xlabel("t (hr)"); ax.set_ylabel("T (°F)"); ax.set_title("(d) Centerline core T(t)")
    ax.legend(fontsize=8)

    # Row 2 — (e) Corner surface T(t)
    ax = axes[1, 1]
    ax.plot(t_hrs, eng_corner_F, label="Engine corner", color="tab:blue")
    if cw_corner_F is not None:
        ax.plot(val.time_hrs, cw_corner_F, "--", label="CW corner", color="tab:orange")
    ax.set_xlabel("t (hr)"); ax.set_ylabel("T (°F)"); ax.set_title("(e) Top-corner surface T(t)")
    ax.legend(fontsize=8)

    # Row 2 — (f) Ambient T phase check
    ax = axes[1, 2]
    if eng_amb_F is not None:
        ax.plot(t_hrs, eng_amb_F, label="Engine ambient", color="tab:blue")
    if val.T_ambient_F is not None:
        ax.plot(val.time_hrs, val.T_ambient_F, "--", label="CW ambient", color="tab:orange")
    ax.set_xlabel("t (hr)"); ax.set_ylabel("T (°F)"); ax.set_title("(f) Ambient T (phase check)")
    ax.legend(fontsize=8)

    # Row 3 — (g) Solar flux on top (PR 2)
    ax = axes[2, 0]
    if result is not None and grid is not None and result.q_solar_history is not None:
        q_cl = result.q_solar_history[:, grid.nx - 1]  # centerline column
        ax.plot(t_hrs, q_cl, color="tab:red", label="Solar flux (centerline)")
        ax.axhline(0.0, color="gray", alpha=0.3, linestyle=":")
        ax.set_xlabel("t (hr)"); ax.set_ylabel("q_sw (W/m²)")
        ax.legend(fontsize=8)
    ax.set_title("(g) Solar flux on top (W/m², negative=heat in)")

    # Row 3 — (h) PR 3: Longwave flux
    ax = axes[2, 1]
    if result is not None and grid is not None and result.q_LW_history is not None:
        q_eff = result.q_LW_history[:, grid.nx - 1]
        q_inc = result.q_LW_incident_history[:, grid.nx - 1]
        ax.plot(t_hrs, q_eff, color="tab:purple", label="q_LW effective (into concrete)")
        ax.plot(t_hrs, q_inc, "--", color="tab:purple", alpha=0.5,
                label="q_LW incident (at blanket surface)")
        ax.axhline(0.0, color="gray", alpha=0.3, linestyle=":")
        ax.set_xlabel("t (hr)"); ax.set_ylabel("q_lw (W/m²)")
        ax.legend(fontsize=8)
    ax.set_title("(h) Longwave flux (W/m², positive=heat out)")

    # Row 3 — (i) PR 4: Total top-surface flux with components
    ax = axes[2, 2]
    if (result is not None and grid is not None
            and result.q_top_total_history is not None
            and result.q_conv_history is not None):
        cl = grid.nx - 1
        ax.plot(t_hrs, result.q_conv_history[:, cl],      ls="--", color="tab:blue",   lw=0.8, label="Convection")
        ax.plot(t_hrs, result.q_evap_history[:, cl],      ls="--", color="tab:green",  lw=0.8, label="Evaporation")
        ax.plot(t_hrs, result.q_solar_history[:, cl],     ls="--", color="tab:orange", lw=0.8, label="Solar (effective)")
        ax.plot(t_hrs, result.q_LW_history[:, cl],        ls="--", color="tab:purple", lw=0.8, label="LW (effective)")
        ax.plot(t_hrs, result.q_top_total_history[:, cl], color="black", lw=1.5, label="Total")
        ax.axhline(0.0, color="gray", alpha=0.3, linestyle=":")
        ax.set_xlabel("t (hr)")
        ax.set_ylabel("W/m² (positive = heat out)")
        ax.legend(fontsize=7, loc="best")
        ax.grid(alpha=0.3)
    ax.set_title("(i) Total top-surface flux (W/m², positive=heat out)")

    # Row 4 — form-face flux panels; corner row = index 0 (topmost side row)
    # (j) Solar flux on form face — parallels (g)
    ax = axes[3, 0]
    if result is not None and result.q_side_solar_history is not None:
        q_side_sol = result.q_side_solar_history[:, 0]
        ax.plot(t_hrs, q_side_sol, color="tab:red", label="Solar flux (corner row)")
        ax.axhline(0.0, color="gray", alpha=0.3, linestyle=":")
        ax.set_xlabel("t (hr)"); ax.set_ylabel("q_sw (W/m²)")
        ax.legend(fontsize=8)
    ax.set_title("(j) Solar flux on form face (W/m², negative=heat in)")

    # (k) Longwave flux on form face — parallels (h)
    ax = axes[3, 1]
    if result is not None and result.q_side_LW_history is not None:
        q_side_lw = result.q_side_LW_history[:, 0]
        ax.plot(t_hrs, q_side_lw, color="tab:purple", label="q_LW form face (corner row)")
        ax.axhline(0.0, color="gray", alpha=0.3, linestyle=":")
        ax.set_xlabel("t (hr)"); ax.set_ylabel("q_lw (W/m²)")
        ax.legend(fontsize=8)
    _lag = f"{soil_lag_hrs:.1f}h" if soil_lag_hrs is not None else "?"
    _damp = f"{soil_damping:.2g}" if soil_damping is not None else "?"
    ax.set_title(f"(k) Form-face longwave (Barber T_ground, lag={_lag}, damp={_damp}; W/m², positive=heat out)")

    # (l) Total form-face flux decomposition — parallels (i); no evap on side face
    ax = axes[3, 2]
    if (result is not None
            and result.q_side_total_history is not None
            and result.q_side_conv_history is not None):
        ax.plot(t_hrs, result.q_side_conv_history[:, 0],  ls="--", color="tab:blue",   lw=0.8, label="Convection")
        ax.plot(t_hrs, result.q_side_solar_history[:, 0], ls="--", color="tab:orange", lw=0.8, label="Solar (effective)")
        ax.plot(t_hrs, result.q_side_LW_history[:, 0],    ls="--", color="tab:purple", lw=0.8, label="LW (effective)")
        ax.plot(t_hrs, result.q_side_total_history[:, 0], color="black", lw=1.5, label="Total")
        ax.axhline(0.0, color="gray", alpha=0.3, linestyle=":")
        ax.set_xlabel("t (hr)")
        ax.set_ylabel("W/m² (positive = heat out)")
        ax.legend(fontsize=7, loc="best")
        ax.grid(alpha=0.3)
    ax.set_title("(l) Total form-face flux (W/m², positive=heat out)")

    plt.tight_layout()
    plt.savefig(out_path, dpi=100)
    plt.close(fig)
    print(f"  Plot saved: {out_path}")


if __name__ == "__main__":
    main()
