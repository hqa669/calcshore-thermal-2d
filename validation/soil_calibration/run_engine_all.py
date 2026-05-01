#!/usr/bin/env python3
"""S2.2 — Run the engine on all 9 input.dat files (hydration suppressed).

Engine is called with:
  - Neutral synthetic environment (flat ambient = placement_temp_F, no solar/wind)
    to isolate soil-coupling physics from top-BC variability.
  - T_ground_deep_C = (soil_temp_F - 32) * 5/9 passed explicitly, to give the engine
    the correct soil Dirichlet BC for model-form comparison against CW.
    NOTE: the engine does NOT do this automatically from input.dat (known gap —
    soil_temp_F is parsed but never consumed by solve_hydration_2d); see
    STAGE2_engine_soil_audit.md for details.

Output: per-run temperature field slices at t=0, 84, 168 hr.

Usage:
    python validation/soil_calibration/run_engine_all.py

Writes:
    validation/soil_calibration/engine_runs/run{A..I}_t{0|84|168}.csv
    validation/soil_calibration/engine_runs/manifest.json
    validation/soil_calibration/engine_runs/grid_info.json
"""
import os
import sys
import json
import time
import traceback
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import (
    parse_cw_dat, CWMixDesign, CWGeometry, CWConstruction, CWEnvironment,
    load_cw_scenario,
)
from thermal_engine_2d import (
    solve_hydration_2d, build_grid_half_mat,
)

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
ENGINE_RUNS = os.path.join(HERE, "engine_runs")
os.makedirs(ENGINE_RUNS, exist_ok=True)

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

TARGET_HRS = [0.0, 84.0, 168.0]


def make_neutral_env(placement_temp_F: float) -> CWEnvironment:
    """Synthetic flat-ambient environment for top-BC calculation.

    All temperatures pinned to placement_temp_F. No solar, no wind.
    This isolates soil-coupling from top-BC variability.
    """
    n_hrs = 169  # 0..168 hr inclusive
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


def nearest_time_idx(t_s, target_hr):
    target_s = target_hr * 3600.0
    return int(np.abs(np.asarray(t_s) - target_s).argmin())


def save_field_csv(path, field_2d_F, x_m, y_m):
    """Save 2D field as CSV with header row/col coordinates.

    First row: x_m coordinates (width)
    First col: y_m coordinates (depth)
    Body: temperatures in °F
    """
    ny, nx = field_2d_F.shape
    with open(path, "w") as f:
        # Header: depth_m, then x coords
        f.write("depth_m\\width_m," + ",".join(f"{x:.4f}" for x in x_m) + "\n")
        for iy in range(ny):
            row = f"{y_m[iy]:.4f}," + ",".join(f"{field_2d_F[iy, ix]:.4f}" for ix in range(nx))
            f.write(row + "\n")


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
            # Load mix/geometry/construction from input.dat (no weather.dat, no cw output)
            mix, geom, constr, raw = parse_cw_dat(dat_path)

            # Apply Hu calibration (apr28) — preserves Hu_J_kg=1 behavior
            from kinetics_correction import compute_hu_factor
            factor, note = compute_hu_factor(mix)
            mix.Hu_factor_calibrated = factor
            mix.Hu_J_kg_effective    = mix.Hu_J_kg * factor   # 1 * factor ≈ 0

            print(f"  Hu_J_kg={mix.Hu_J_kg}, Hu_J_kg_effective={mix.Hu_J_kg_effective:.4f} "
                  f"(hydration {'suppressed' if mix.Hu_J_kg_effective < 10 else 'ACTIVE!'})")

            # Build grid
            grid = build_grid_half_mat(geom.width_ft, geom.depth_ft)
            if grid_info is None:
                jslice, islice = grid.concrete_slice()
                y_conc = grid.y[jslice]
                x_conc = grid.x[islice]
                grid_info = {
                    "nx_full": grid.nx,
                    "ny_full": grid.ny,
                    "nx_concrete": len(x_conc),
                    "ny_concrete": len(y_conc),
                    "x_concrete_m": x_conc.tolist(),
                    "y_concrete_m": y_conc.tolist(),
                    "ix_concrete_start": grid.ix_concrete_start,
                    "ix_concrete_end": grid.ix_concrete_end,
                    "iy_concrete_start": grid.iy_concrete_start,
                    "iy_concrete_end": grid.iy_concrete_end,
                    "width_ft": geom.width_ft,
                    "depth_ft": geom.depth_ft,
                    "dx_m": float(grid.dx),
                }
                print(f"  Grid: concrete {len(y_conc)}×{len(x_conc)} "
                      f"(y: {y_conc[0]:.3f}–{y_conc[-1]:.3f} m, "
                      f"x: {x_conc[0]:.3f}–{x_conc[-1]:.3f} m)")

            # Initial conditions:
            #   Concrete + blanket: placement temperature (physical: fresh pour)
            #   Soil region: soil temperature (physical: undisturbed ground)
            # NOTE: the engine has a 3m-deep soil region (n_soil_y=15, dy=0.2m).
            # If the soil is initialized at placement_temp with the Dirichlet BC
            # at the FAR soil boundary set to soil_temp, it takes t ≈ L²/α_soil
            # ≈ (3m)²/2.86e-3 m²/hr ≈ 3150 hr for the signal to reach the
            # concrete-soil interface — 20× longer than our 168hr run window.
            # Initializing soil at soil_temp_F gives the physically correct
            # starting condition and allows the comparison to isolate model-form
            # differences (spatial discretization, k/ρ/Cp) rather than initial-
            # condition transients.
            T0_C = (constr.placement_temp_F - 32.0) * 5.0 / 9.0
            T_soil_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            T_initial = np.full((grid.ny, grid.nx), T0_C)
            T_initial[grid.is_soil] = T_soil_C  # soil starts at soil temp, not placement temp

            # Synthetic neutral environment
            env = make_neutral_env(constr.placement_temp_F)

            # Soil BC: pass explicitly from input.dat soil_temp_F
            T_ground_deep_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            print(f"  T_ground_deep_C = {T_ground_deep_C:.2f}°C ({constr.soil_temp_F}°F)")

            t_wall_start = time.perf_counter()
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
            t_wall = time.perf_counter() - t_wall_start
            print(f"  Solved in {t_wall:.1f}s, {result.n_output_samples} samples")

            # Extract concrete subfield
            jslice, islice = grid.concrete_slice()
            T_conc_C = result.T_field_C[:, jslice, islice]  # (n_samples, ny_conc, nx_conc)
            T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
            t_hrs = result.t_s / 3600.0
            y_conc = grid.y[jslice]
            x_conc = grid.x[islice]

            # Save slices at target hours
            for target_hr in TARGET_HRS:
                ti = nearest_time_idx(result.t_s, target_hr)
                actual_hr = float(t_hrs[ti])
                field = T_conc_F[ti]  # (ny_conc, nx_conc)

                csv_path = os.path.join(ENGINE_RUNS, f"run{label}_t{int(target_hr)}.csv")
                save_field_csv(csv_path, field, x_conc, y_conc)
                print(f"  Saved: {csv_path} (actual t={actual_hr:.2f} hr, "
                      f"T range: {field.min():.1f}–{field.max():.1f}°F)")

            # Sanity check
            T_min = float(T_conc_F.min())
            T_max = float(T_conc_F.max())
            if T_min < 0 or T_max > 300:
                raise ValueError(f"Implausible temperature range: {T_min:.1f}–{T_max:.1f}°F")

            manifest[label] = {
                "status": "ok",
                "n_samples": int(result.n_output_samples),
                "t_wall_s": float(t_wall),
                "T_min_F": T_min,
                "T_max_F": T_max,
                "placement_F": placement,
                "soil_F": soil,
            }

        except Exception as e:
            tb = traceback.format_exc()
            print(f"  ERROR: {e}\n{tb}")
            manifest[label] = {"status": "error", "error": str(e), "traceback": tb}

    # Save grid info and manifest
    grid_path = os.path.join(ENGINE_RUNS, "grid_info.json")
    if grid_info is not None:
        with open(grid_path, "w") as f:
            json.dump(grid_info, f, indent=2)
        print(f"\nWritten: {grid_path}")

    manifest_path = os.path.join(ENGINE_RUNS, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Written: {manifest_path}")

    n_ok = sum(1 for v in manifest.values() if v.get("status") == "ok")
    n_err = sum(1 for v in manifest.values() if v.get("status") == "error")
    n_skip = sum(1 for v in manifest.values() if v.get("status") == "skip")
    print(f"\nSummary: {n_ok} OK, {n_err} errors, {n_skip} skipped")


if __name__ == "__main__":
    run()
