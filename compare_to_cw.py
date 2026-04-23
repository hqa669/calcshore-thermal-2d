#!/usr/bin/env python3
"""
CW MIX-01 validation harness.

Usage:
    python compare_to_cw.py <directory>

<directory> must contain input.dat, weather.dat, output.txt.

Exit 0 if ALL metrics PASS. Exit 1 if ANY FAIL.
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

# CW reference values for MIX-01 (docs/coding_passdown_v3.md)
CW_PEAK_MAX_F    = 129.6
CW_PEAK_MAX_T_HR = 145.8
CW_PEAK_GRAD_F   = 39.3
CW_PEAK_GRAD_T_HR = 146.2


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


def main():
    args = parse_args()
    d = args.scenario_dir.rstrip("/")
    check_files(d)

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

        center_rms_F = rms(eng_center_interp, cw_center_F)
        corner_rms_F = rms(eng_corner_interp, cw_corner_F)

        # Field-wide: compare centerline column of engine to CW centerline across all depths
        n_conc_y = grid.iy_concrete_end - grid.iy_concrete_start + 1
        sq_errs = []
        for ti in range(n_cw_t):
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

    overall = pass_peak_max and pass_peak_grad
    if pass_field is not None:  overall = overall and pass_field
    if pass_center is not None: overall = overall and pass_center
    if pass_corner is not None: overall = overall and pass_corner

    # ------------------------------------------------------------------ #
    # 7. Print formatted table
    # ------------------------------------------------------------------ #
    n_out   = result.n_output_samples
    n_steps = result.n_inner_steps
    n_nodes = grid.nx * grid.ny
    scenario_name = "MIX-01 Austin 2026-07-15"

    print()
    print(f"  {scenario_name}  ({168} hr run, {n_steps} timesteps, {n_nodes} nodes)")
    print(f"  Peak Max T:     Engine {engine_peak_max_F:6.1f}°F @ {engine_peak_max_hr:5.1f} hr | CW {cw_peak_max_F:6.1f}°F @ {cw_peak_max_hr:.1f} hr")
    print(f"                  Δ = {peak_max_delta:+.1f}°F   [{pass_fail(pass_peak_max)}] (tol ±{TOL_PEAK_MAX_F:.1f}°F)")
    print(f"  Peak Gradient:  Engine {engine_peak_grad_F:5.1f}°F @ {engine_peak_grad_hr:5.1f} hr | CW {cw_peak_grad_F:5.1f}°F @ {CW_PEAK_GRAD_T_HR:.1f} hr")
    print(f"                  Δ = {peak_grad_delta:+.1f}°F   [{pass_fail(pass_peak_grad)}] (tol ±{TOL_PEAK_GRAD_F:.1f}°F)")

    if field_rms_F is not None:
        print(f"  Field-wide RMS: {field_rms_F:.2f}°F   [{pass_fail(pass_field)}] (tol {TOL_FIELD_RMS_F:.1f}°F)")
    else:
        print(f"  Field-wide RMS: N/A (CW T_field_F not in fixture)")

    if center_rms_F is not None:
        print(f"  Centerline RMS: {center_rms_F:.2f}°F   [{pass_fail(pass_center)}] (tol {TOL_CENTER_RMS_F:.1f}°F)")
    else:
        print(f"  Centerline RMS: N/A (CW T_field_F not in fixture)")

    if corner_rms_F is not None:
        print(f"  Corner RMS:     {corner_rms_F:.2f}°F   [{pass_fail(pass_corner)}] (tol {TOL_CORNER_RMS_F:.1f}°F)")
    else:
        print(f"  Corner RMS:     N/A (CW T_field_F not in fixture)")

    print(f"  Runtime:        {t_wall_s:.1f} s")
    # Determine overall status with Sprint-0 known-limitation annotation
    if overall:
        overall_str = "[PASS]"
    elif (pass_peak_max and pass_peak_grad and pass_field
          and (pass_center is None or pass_center)
          and not (pass_corner is None or pass_corner)):
        # Corner RMS is the only failure; root cause is missing solar forcing.
        overall_str = "[PARTIAL — 4/5 metrics pass; Corner RMS deferred to Sprint 1 (solar forcing)]"
    elif pass_peak_max and pass_peak_grad and pass_field:
        overall_str = "[PARTIAL — Peak, Gradient, Field pass; Ctr/Corner deferred to Sprint 1]"
    else:
        overall_str = "[FAIL]"
    print(f"  OVERALL:        {overall_str}")
    print()

    # ------------------------------------------------------------------ #
    # 8. Optional plot
    # ------------------------------------------------------------------ #
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        _make_plot(t_hrs, engine_peak_max_F, eng_center_F, eng_corner_F, eng_amb_F,
                   val, cw_center_F, cw_corner_F, grad_series, d,
                   result=result, grid=grid)
    except ImportError:
        pass  # silently skip if matplotlib not installed
    except Exception:
        pass  # don't crash on plot errors

    sys.exit(0 if overall else 1)


def _make_plot(t_hrs, engine_peak_max_F, eng_center_F, eng_corner_F, eng_amb_F,
               val, cw_center_F, cw_corner_F, grad_series, scenario_dir,
               result=None, grid=None):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
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

    # Row 3 — (i) Placeholder for PR 4: Total top flux
    ax = axes[2, 2]
    ax.set_title("(i) PR 4: Total top flux (pending)")
    ax.set_axis_off()

    plt.tight_layout()
    out = "cw_comparison_MIX-01.png"
    plt.savefig(out, dpi=100)
    plt.close(fig)
    print(f"  Plot saved: {out}")


if __name__ == "__main__":
    main()
