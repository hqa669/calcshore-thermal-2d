#!/usr/bin/env python3
"""Stage 3 Fix 1 — engine runs with soil_temp_F plumbed, is_submerged=False.

Runs the engine on all 9 CW-dataset input.dat files using the production
code path after Fix 1:
  - T_ground_deep_C derived from construction.soil_temp_F (no longer CW Eq 44)
  - Soil region initialized at soil_temp_F (no IC transient)
  - is_submerged=False (concrete sides still air — Fix 2 not yet applied)
  - Neutral flat-ambient environment (isolates soil-coupling physics from top-BC)

Neutral environment is retained here to match the Stage 2 diagnostic protocol and
clearly show the side-soil residual without top-BC noise. Fix 2 uses real weather.

Outputs: stage3_fix1_runs/run{A..I}_t168.csv + manifest.json + grid_info.json

Usage:
    python validation/soil_calibration/stage3_run_fix1.py
"""
import os
import sys
import json
import time
import traceback
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import parse_cw_dat, CWEnvironment
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
OUT_DIR = os.path.join(HERE, "stage3_fix1_runs")
os.makedirs(OUT_DIR, exist_ok=True)

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


def make_neutral_env(placement_temp_F: float) -> CWEnvironment:
    n_hrs = 169
    hours = np.arange(n_hrs, dtype=float)
    T_air = np.full(n_hrs, placement_temp_F)
    T_air_C = (T_air - 32.0) * 5.0 / 9.0
    return CWEnvironment(
        hours=hours,
        T_air_F=T_air,
        RH_pct=np.full(n_hrs, 60.0),
        solar_W_m2=np.zeros(n_hrs),
        wind_m_s=np.zeros(n_hrs),
        cloud_cover=np.ones(n_hrs),
        pressure_hPa=np.full(n_hrs, 1013.25),
        T_dp_C=T_air_C - 10.0,
        T_sky_C=T_air_C - 5.0,
        daily_max_F=[placement_temp_F] * 8,
        daily_min_F=[placement_temp_F] * 8,
        cw_ave_max_daily_temp_F=placement_temp_F,
        cw_ave_min_daily_temp_F=placement_temp_F,
        cw_ave_max_solar_W_m2=0.0,
        cw_ave_max_wind_m_s=0.0,
        cw_ave_max_RH_pct=60.0,
        cw_ave_min_RH_pct=60.0,
    )


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
                is_submerged=False,  # Fix 1 only: sides still air
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
                    "is_submerged": False,
                }

            T0_C = (constr.placement_temp_F - 32.0) * 5.0 / 9.0
            T_soil_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            T_initial = np.full((grid.ny, grid.nx), T0_C)
            T_initial[grid.is_soil] = T_soil_C

            env = make_neutral_env(constr.placement_temp_F)
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
