#!/usr/bin/env python3
"""Stage 5e production gate — 9-run Structure C metric, NO k override.

k_uc calibration is now baked into the engine (K_UC_CALIBRATION_FACTOR_SPRINT7=0.96).
This script runs all 9 scenarios using the unmodified engine and checks that all pass
Structure C: R1 max|R| ≤ 0.35°F, R2 max|R| ≤ 0.35°F (R3 RMSE reported only).
"""
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

CW_RUNS    = os.path.join(SC, "cw_runs")
COMPARE_HR = 168
GATE_F     = 0.35
R1_DI      = 24
R2_WI      = 0
R3_DI      = slice(43, 49)
R3_WI      = slice(8, 13)

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


def run_engine(folder, pl_F, so_F):
    """Run engine with NO k override — uses engine's built-in calibration."""
    mix, geom, constr, _ = parse_cw_dat(os.path.join(CW_RUNS, folder, "input.dat"))
    fac, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = fac
    mix.Hu_J_kg_effective = mix.Hu_J_kg * fac
    # NO k override here — K_UC_CALIBRATION_FACTOR_SPRINT7 is applied inside engine
    constr.model_soil = False; constr.is_submerged = True
    grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                               is_submerged=True, model_soil=False,
                               blanket_thickness_m=0.0)
    T0 = (pl_F - 32) * 5/9
    Ts = (so_F - 32) * 5/9
    Ti = np.full((grid.ny, grid.nx), T0); Ti[grid.is_soil] = Ts
    res = solve_hydration_2d(grid, mix, Ti, duration_s=COMPARE_HR * 3600,
                             output_interval_s=1800., boundary_mode="full_2d",
                             environment=make_neutral_env(pl_F), construction=constr,
                             T_ground_deep_C=Ts, diagnostic_outputs=False)
    jsl, isl = grid.concrete_slice()
    ti = nearest_time_idx(res.t_s, float(COMPARE_HR))
    TF = res.T_field_C[ti, jsl, isl] * 9/5 + 32
    return grid.y[jsl], grid.x[isl], TF


def load_cw(folder):
    v = parse_cw_temp_output(os.path.join(CW_RUNS, folder, "output.txt"))
    ti = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m


print("Stage 5e — Production Gate (Structure C, k_uc calibrated in engine)")
print("=" * 70)
print(f"{'Run':>4}  {'R1 max|R|':>11}  {'R2 max|R|':>11}  {'R3 RMSE':>9}  {'R1':>5}  {'R2':>5}  {'Overall'}")
print("-" * 70)

rows = []
all_pass = True

for label, folder, pl_F, so_F in RUNS:
    ey, ex, eT = run_engine(folder, pl_F, so_F)
    cw_F, cw_widths_m, cw_depths_m = load_cw(folder)
    eng = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)
    R = eng - cw_F

    r1 = float(np.max(np.abs(R[R1_DI, :])))
    r2 = float(np.max(np.abs(R[slice(24, 49), R2_WI])))
    r3_rmse = float(np.sqrt(np.mean(R[R3_DI, R3_WI] ** 2)))

    gate_r1 = r1 <= GATE_F
    gate_r2 = r2 <= GATE_F
    overall = "PASS" if (gate_r1 and gate_r2) else "FAIL"
    if not (gate_r1 and gate_r2):
        all_pass = False

    print(f"  {label:>2}  {r1:>11.4f}  {r2:>11.4f}  {r3_rmse:>9.4f}"
          f"  {'✓' if gate_r1 else '✗':>5}  {'✓' if gate_r2 else '✗':>5}  [{overall}]")
    rows.append((label, r1, r2, r3_rmse, gate_r1, gate_r2))

print("-" * 70)
print(f"\nStructure C gate: {'ALL 9 PASS' if all_pass else 'FAILURES PRESENT'}")
print(f"  R1 ≤ {GATE_F}°F (side profile, di=24)")
print(f"  R2 ≤ {GATE_F}°F (bottom profile, wi=0, di=24..48)")
print(f"  R3 RMSE reported only (di=43..48, wi=8..12)")

# Summary stats
r1_vals = [r[1] for r in rows]
r2_vals = [r[2] for r in rows]
r3_vals = [r[3] for r in rows]
print(f"\n  R1: mean={np.mean(r1_vals):.4f}  max={max(r1_vals):.4f}")
print(f"  R2: mean={np.mean(r2_vals):.4f}  max={max(r2_vals):.4f}")
print(f"  R3 RMSE: mean={np.mean(r3_vals):.4f}  max={max(r3_vals):.4f}")
