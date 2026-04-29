#!/usr/bin/env python3
"""
CHECK 1.2 — Instrumented MIX-01 run to confirm apr28 Hu calibration is
actually applied on the run_all.py / compare_to_cw.run_one path.

Loads MIX-01 the same way run_all.py does (load_cw_scenario with default
use_cw_calibrated_hu=True), prints the four calibration-related fields,
then runs the full-stack solver and prints the peak temperature plus the
Hu value immediately before the solve call so the value the solver
actually consumes is logged.

Read-only: makes no source changes. Writes only stdout/stderr.
"""

import os
import sys
import time

import numpy as np

# Run from repo root regardless of cwd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
os.chdir(REPO_ROOT)

from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d


SCENARIO_DIR = "validation/cw_exports/MIX-01"


def main():
    input_dat   = os.path.join(SCENARIO_DIR, "input.dat")
    weather_dat = os.path.join(SCENARIO_DIR, "weather.dat")
    output_txt  = os.path.join(SCENARIO_DIR, "output.txt")

    # Load via the same path run_all.py → run_one uses (positional, no
    # use_cw_calibrated_hu kwarg passed → default True applies).
    scn = load_cw_scenario(input_dat, weather_dat, output_txt)
    mix = scn.mix

    print("=" * 72)
    print("CHECK 1.2 — MIX-01 calibration-field dump after load_cw_scenario")
    print("=" * 72)
    print(f"Hu_J_kg              = {mix.Hu_J_kg}")
    print(f"Hu_J_kg_effective    = {mix.Hu_J_kg_effective}")
    print(f"Hu_factor_calibrated = {mix.Hu_factor_calibrated}")
    print(f"Hu_calibration_note  = {mix.Hu_calibration_note!r}")
    if mix.Hu_J_kg != 0.0:
        ratio = mix.Hu_J_kg_effective / mix.Hu_J_kg
        print(f"ratio (effective/raw) = {ratio:.6f}")
    else:
        print("ratio (effective/raw) = (Hu_J_kg is zero — cannot compute)")

    # Build grid + initial condition exactly as compare_to_cw.run_one does
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)

    print()
    print("=" * 72)
    print("Pre-solve: Hu the solver is about to consume")
    print("=" * 72)
    print(f"mix.Hu_J_kg_effective (passed via scn.mix to solver) = "
          f"{mix.Hu_J_kg_effective}")
    print(f"thermal_engine_2d.py:1438 will read: Hu = mix.Hu_J_kg_effective")

    print()
    print("=" * 72)
    print("Running solve_hydration_2d (168hr, full_2d boundary mode)")
    print("=" * 72)
    t0 = time.perf_counter()
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )
    t_wall = time.perf_counter() - t0

    # Engine peak temperature in °F over concrete domain
    jslice, islice = grid.concrete_slice()
    T_conc_C = result.T_field_C[:, jslice, islice]
    T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
    peak_F = float(T_conc_F.max())
    peak_idx = int(np.argmax(T_conc_F.max(axis=(1, 2))))
    peak_hr = float(result.t_s[peak_idx] / 3600.0)

    # CW peak from validation fixture
    cw_peak_F = None
    cw_peak_hr = None
    if scn.cw_validation is not None:
        cw_peak_F = float(scn.cw_validation.T_max_xs_F.max())
        cw_peak_idx = int(np.argmax(scn.cw_validation.T_max_xs_F))
        cw_peak_hr = float(scn.cw_validation.time_hrs[cw_peak_idx])

    print()
    print("=" * 72)
    print("Post-solve: peak temperatures")
    print("=" * 72)
    print(f"Engine peak max T = {peak_F:.2f} °F @ {peak_hr:.1f} hr")
    if cw_peak_F is not None:
        print(f"CW peak max T     = {cw_peak_F:.2f} °F @ {cw_peak_hr:.1f} hr")
        print(f"Δ (engine − CW)   = {peak_F - cw_peak_F:+.2f} °F")
    print(f"Solver wall time  = {t_wall:.2f} s")

    # Sanity comparison: what the post-apr28 run_all_output.md reported
    print()
    print("=" * 72)
    print("Sanity check vs. validation/post_apr28_full_stack/run_all_output.md")
    print("=" * 72)
    print("That file recorded MIX-01 PeakMax Δ = -3.84°F (engine 125.7°F vs CW "
          "129.6°F).")
    if cw_peak_F is not None:
        delta = peak_F - cw_peak_F
        match = abs(delta - (-3.84)) < 0.05
        print(f"This run:                         PeakMax Δ = "
              f"{delta:+.2f}°F  → {'MATCHES committed run' if match else 'DIFFERS from committed run'}")


if __name__ == "__main__":
    main()
