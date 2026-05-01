#!/usr/bin/env python3
"""Debug where the full-field max|R| lives for Run A k×1.00 vs k×0.96."""
import os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SC   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, SC)

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

CW_RUNS = os.path.join(SC, "cw_runs")

def run_engine(folder, pl_F, so_F, k_factor):
    mix, geom, constr, _ = parse_cw_dat(os.path.join(CW_RUNS, folder, "input.dat"))
    fac, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = fac
    mix.Hu_J_kg_effective = mix.Hu_J_kg * fac
    mix.thermal_conductivity_BTU_hr_ft_F *= k_factor
    constr.model_soil = False; constr.is_submerged = True
    grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                               is_submerged=True, model_soil=False,
                               blanket_thickness_m=0.0)
    T0 = (pl_F-32)*5/9; Ts = (so_F-32)*5/9
    Ti = np.full((grid.ny, grid.nx), T0); Ti[grid.is_soil] = Ts
    res = solve_hydration_2d(grid, mix, Ti, duration_s=168*3600,
                             output_interval_s=1800., boundary_mode="full_2d",
                             environment=make_neutral_env(pl_F), construction=constr,
                             T_ground_deep_C=Ts, diagnostic_outputs=False)
    jsl, isl = grid.concrete_slice()
    ti = nearest_time_idx(res.t_s, 168.0)
    TF = res.T_field_C[ti, jsl, isl]*9/5+32
    return grid.y[jsl], grid.x[isl], TF

def load_cw(folder):
    v = parse_cw_temp_output(os.path.join(CW_RUNS, folder, "output.txt"))
    ti = int(np.abs(v.time_hrs - 168).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m

cw_F, cw_widths_m, cw_depths_m = load_cw("runA_baseline")

for kfac in [1.00, 0.96]:
    ey, ex, eT = run_engine("runA_baseline", 73, 73, kfac)
    eng_on_cw = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)
    R = eng_on_cw - cw_F
    ff_max = float(np.max(np.abs(R)))
    loc = np.unravel_index(np.argmax(np.abs(R)), R.shape)
    di_max, wi_max = loc
    print(f"\nk×{kfac:.2f}: full-field max|R| = {ff_max:.4f}°F  at (di={di_max}, wi={wi_max})")
    print(f"  R[{di_max},{wi_max}] = {R[di_max, wi_max]:+.4f}°F")
    print(f"  Engine T = {eng_on_cw[di_max, wi_max]:.4f}°F,  CW T = {cw_F[di_max, wi_max]:.4f}°F")
    # Show top 5 di=0..4 max to understand top-row contribution
    top5_max = float(np.max(np.abs(R[:5, :])))
    top5_loc = np.unravel_index(np.argmax(np.abs(R[:5, :])), R[:5,:].shape)
    print(f"  Top-5-row (di=0..4) max|R| = {top5_max:.4f}°F at (di={top5_loc[0]}, wi={top5_loc[1]})")
    # Corner R3 region
    r3_max = float(np.max(np.abs(R[43:49, 8:13])))
    print(f"  R3 corner (di=43..48, wi=8..12) max|R| = {r3_max:.4f}°F")
    # Bottom row only (di=48)
    r_bot = float(np.max(np.abs(R[48, :])))
    print(f"  Bottom row (di=48) max|R| = {r_bot:.4f}°F")
    # Print the residual map for the corner window
    print(f"  Corner R map (di=44..48, wi=8..12):")
    for di in range(44, 49):
        vals = "  ".join(f"{R[di, wi]:+.3f}" for wi in range(8, 13))
        print(f"    di={di}: {vals}")
