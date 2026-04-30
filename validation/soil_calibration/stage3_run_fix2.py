#!/usr/bin/env python3
"""Stage 3 Fix 2 — engine runs with is_submerged=True + soil_temp_F plumbed.

Both Stage 3 fixes active:
  Fix 1: T_ground_deep_C from construction.soil_temp_F; soil initialized at soil_temp_F.
  Fix 2: is_submerged=True — concrete-height rows in the soil-extension columns are
         soil (material_id=2) rather than air (material_id=3). The left face of the
         concrete now has soil contact, mirroring the bottom-soil geometry.

Uses the actual Austin TX weather (from cw_exports/MIX-01/weather.dat, shared by all 9
runs) rather than a synthetic neutral environment. All 9 runs place on July 15 at 5 AM
with the same location, so a single weather file applies to all.

Output: stage3_fix2_runs/run{A..I}_t168.csv + manifest.json + grid_info.json

Usage:
    python validation/soil_calibration/stage3_run_fix2.py
"""
import os
import sys
import json
import time
import traceback
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import parse_cw_dat, load_cw_scenario
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
OUT_DIR = os.path.join(HERE, "stage3_fix2_runs")
os.makedirs(OUT_DIR, exist_ok=True)

# All 9 runs share Austin TX, July 15 placement — load shared weather from MIX-01.
_WEATHER_SOURCE = os.path.join(
    os.path.dirname(HERE), "../validation/cw_exports/MIX-01/weather.dat"
)
# Resolve to absolute path
_WEATHER_SOURCE = os.path.normpath(
    os.path.join(HERE, "../../validation/cw_exports/MIX-01/weather.dat")
)

RUNS = [
    ("A", "runA_baseline",  73,  73),
    ("B", "runB_73_60",     73,  60),
    ("C", "runC_73_90",     73,  90),
    ("D", "runD_60_73",     60,  73),
    ("E", "runE_90_73",     90,  73),
    ("F", "runF_73_45",     73,  45),
    ("G", "runG_73_100",    73, 100),
    ("H", "runH_45_73",     45,  73),
    ("I", "runI_100_73",   100,  73),
]


def save_field_csv(path, field_2d_F, x_m, y_m):
    ny, nx = field_2d_F.shape
    with open(path, "w") as f:
        f.write("depth_m\\width_m," + ",".join(f"{x:.4f}" for x in x_m) + "\n")
        for iy in range(ny):
            row = f"{y_m[iy]:.4f}," + ",".join(f"{field_2d_F[iy, ix]:.4f}" for ix in range(nx))
            f.write(row + "\n")


def nearest_time_idx(t_s, target_hr):
    return int(np.abs(np.asarray(t_s) - target_hr * 3600.0).argmin())


def run():
    # Load shared Austin TX weather (all 9 runs placed July 15, 5 AM, Austin TX)
    if not os.path.isfile(_WEATHER_SOURCE):
        print(f"ERROR: weather.dat not found at {_WEATHER_SOURCE}")
        return
    _weather_scn = load_cw_scenario(
        os.path.join(HERE, "cw_runs/runA_baseline/input.dat"),
        _WEATHER_SOURCE,
        os.path.join(HERE, "cw_runs/runA_baseline/output.txt"),
    )
    _shared_env = _weather_scn.environment
    print(f"Loaded Austin TX weather from {_WEATHER_SOURCE}")

    manifest = {}
    grid_info = None

    for label, folder, placement, soil in RUNS:
        print(f"\nRun {label} ({folder}): placement={placement}°F, soil={soil}°F")
        dat_path = os.path.join(CW_RUNS, folder, "input.dat")
        if not os.path.isfile(dat_path):
            print(f"  SKIP: input.dat not found")
            manifest[label] = {"status": "skip", "error": "input.dat not found"}
            continue

        try:
            mix, geom, constr, _ = parse_cw_dat(dat_path)
            factor, _ = compute_hu_factor(mix)
            mix.Hu_factor_calibrated = factor
            mix.Hu_J_kg_effective = mix.Hu_J_kg * factor

            grid = build_grid_half_mat(
                geom.width_ft, geom.depth_ft,
                is_submerged=True,  # Fix 2: concrete sides contact soil
            )

            if grid_info is None:
                jslice, islice = grid.concrete_slice()
                y_conc = grid.y[jslice]
                x_conc = grid.x[islice]
                grid_info = {
                    "nx_concrete": len(x_conc),
                    "ny_concrete": len(y_conc),
                    "x_concrete_m": x_conc.tolist(),
                    "y_concrete_m": y_conc.tolist(),
                    "ix_concrete_start": grid.ix_concrete_start,
                    "iy_concrete_start": grid.iy_concrete_start,
                    "iy_concrete_end": grid.iy_concrete_end,
                    "dx_m": float(grid.dx),
                    "is_submerged": True,
                    "n_soil_cells": int(grid.is_soil.sum()),
                    "n_air_cells": int(grid.is_air.sum()),
                }
                print(f"  Grid: {grid.ny}×{grid.nx}, "
                      f"soil cells={grid.is_soil.sum()}, air cells={grid.is_air.sum()}")

            T0_C = (constr.placement_temp_F - 32.0) * 5.0 / 9.0
            T_soil_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            T_initial = np.full((grid.ny, grid.nx), T0_C)
            T_initial[grid.is_soil] = T_soil_C

            env = _shared_env  # Austin TX weather shared by all runs
            T_ground_deep_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            print(f"  soil_temp_F={constr.soil_temp_F}°F → T_ground_deep_C={T_ground_deep_C:.2f}°C")

            t0 = time.perf_counter()
            result = solve_hydration_2d(
                grid, mix, T_initial,
                duration_s=168 * 3600,
                output_interval_s=1800.0,
                boundary_mode="full_2d",
                environment=env,
                construction=constr,
                T_ground_deep_C=T_ground_deep_C,
                diagnostic_outputs=False,
            )
            t_wall = time.perf_counter() - t0
            print(f"  Solved in {t_wall:.1f}s")

            jslice, islice = grid.concrete_slice()
            T_conc_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
            t_hrs = result.t_s / 3600.0

            ti = nearest_time_idx(result.t_s, 168.0)
            csv_path = os.path.join(OUT_DIR, f"run{label}_t168.csv")
            save_field_csv(csv_path, T_conc_F[ti], grid.x[islice], grid.y[jslice])

            manifest[label] = {
                "status": "ok",
                "t_wall_s": float(t_wall),
                "T_min_F": float(T_conc_F.min()),
                "T_max_F": float(T_conc_F.max()),
                "placement_F": placement,
                "soil_F": soil,
            }
            print(f"  T range: {T_conc_F.min():.1f}–{T_conc_F.max():.1f}°F")

        except Exception as e:
            print(f"  ERROR: {e}\n{traceback.format_exc()}")
            manifest[label] = {"status": "error", "error": str(e)}

    if grid_info is not None:
        with open(os.path.join(OUT_DIR, "grid_info.json"), "w") as f:
            json.dump(grid_info, f, indent=2)

    with open(os.path.join(OUT_DIR, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

    n_ok = sum(1 for v in manifest.values() if v.get("status") == "ok")
    print(f"\nDone: {n_ok}/9 runs OK. Output: {OUT_DIR}/")


if __name__ == "__main__":
    run()
