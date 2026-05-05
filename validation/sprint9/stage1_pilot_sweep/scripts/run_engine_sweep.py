#!/usr/bin/env python3
"""Sprint 9 Stage 1-pilot — (Hu_factor, c) sensitivity sweep engine runner.

Sweeps 8 Hu_factor_override values × 5 c_multiplier values × 2 mixes = 80 runs.
No engine source modification; no .npz files written — scalars extracted in-process.

Override semantics:
  Hu_factor_override replaces compute_hu_factor(mix) entirely.
  c_eff = c_multiplier × 1.0532  (replaces C_VAL everywhere).

Writes:
  data/sweep_residual_table.csv  — 80 rows; max|R|, R(di=47/48/36), pass/fail
  data/sweep_wrapper_log.csv     — 80 rows; full provenance per run

Usage:
    cd /Users/hqa668/calcshore-thermal-2d
    python validation/sprint9/stage1_pilot_sweep/scripts/run_engine_sweep.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE  = Path(__file__).resolve().parent
SWEEP = HERE.parent                             # stage1_pilot_sweep/
PILOT = SWEEP.parent / "stage1_pilot"          # stage1_pilot/ (CW data source)
ROOT  = (SWEEP / "../../..").resolve()         # calcshore-thermal-2d/
SC    = ROOT / "validation" / "soil_calibration"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env

# Sweep grid (brief §Setup)
HU_FACTOR_GRID = [0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95]
C_MULT_GRID    = [0.95, 0.97, 1.00, 1.03, 1.05]

# Fixed constants — identical to stage1_pilot
HU_FLOOR_THRESHOLD = 10_000.0
HU_RESIDUAL        = 12_937.0
T_PL_F             = 73.0
T_SO_F             = 85.0
DURATION_HR        = 168.0
OUTPUT_INTERVAL_S  = 3600.0
C_ANCHOR           = 1.0532   # c(T_pl=73) pilot baseline; c_eff = c_multiplier × this

# Residual gate and depth indices (match analyze_pilot.py)
GATE_F  = 0.5
DI_CORE = 24
DI_BOT  = 48

DATA_DIR = SWEEP / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

MIXES = [
    ("mix01", PILOT / "cw_data/mix01"),
    ("mix07", PILOT / "cw_data/mix07"),
]


def load_cw_once(cw_folder: Path) -> dict:
    v = parse_cw_temp_output(str(cw_folder / "output.txt"))
    return {
        "t_hrs":    v.time_hrs,
        "T_cw_F":   v.T_field_F,    # (n_t_cw, nD, nW)
        "depths_m": v.depths_m,
        "widths_m": v.widths_m,
    }


def align_cw_to_engine(eng_t, cw_t, cw_T):
    n_eng = len(eng_t)
    nD, nW = cw_T.shape[1], cw_T.shape[2]
    T_aligned = np.zeros((n_eng, nD, nW), dtype=np.float32)
    for i, th in enumerate(eng_t):
        ti_cw = int(np.abs(cw_t - th).argmin())
        T_aligned[i] = cw_T[ti_cw]
    return T_aligned


def ti_near(t_arr, target_hr):
    return int(np.abs(np.asarray(t_arr) - target_hr).argmin())


def compute_scalars(T_engine_F, t_hrs, cw):
    """Extract residual scalars in-process; no .npz needed."""
    T_cw_aln = align_cw_to_engine(t_hrs, cw["t_hrs"], cw["T_cw_F"])
    ti168 = ti_near(t_hrs, 168.0)

    R_t168_wi0 = T_engine_F[ti168, :, 0] - T_cw_aln[ti168, :, 0]
    R_bottom   = R_t168_wi0[DI_CORE:DI_BOT + 1]   # [24..48], 25 pts
    max_R      = float(np.max(np.abs(R_bottom)))
    di_rel     = int(np.argmax(np.abs(R_bottom)))
    di_at_max  = DI_CORE + di_rel

    return {
        "max_R_F":   max_R,
        "di_at_max": di_at_max,
        "R_di47_F":  float(R_t168_wi0[47]),
        "R_di48_F":  float(R_t168_wi0[48]),
        "R_di36_F":  float(R_t168_wi0[36]),
        "pass":      "PASS" if max_R < GATE_F else "FAIL",
    }


def run_one(mix_name, cw_folder, Hu_fac, c_mult, cw):
    """Run engine for one (Hu_factor_override, c_multiplier, mix) cell.

    Returns (scalars_dict, log_dict). No .npz written.
    """
    mix, geom, constr, _ = parse_cw_dat(str(cw_folder / "input.dat"))

    alpha_u_raw = mix.alpha_u
    Hu_raw      = mix.Hu_J_kg

    # Override 1: replace compute_hu_factor with the direct Hu_factor value
    fac = Hu_fac

    # Override 2: c_eff = c_multiplier × c(T_pl=73)
    c_eff = c_mult * C_ANCHOR

    # Hu_residual conditional (passthrough for both pilot mixes; log it)
    if Hu_raw < HU_FLOOR_THRESHOLD:
        Hu_base  = HU_RESIDUAL
        cond_txt = f"APPLIED (Hu_raw={Hu_raw:.0f} < {HU_FLOOR_THRESHOLD:.0f})"
    else:
        Hu_base  = Hu_raw
        cond_txt = "passthrough"

    Hu_factored = Hu_base * fac
    alpha_u_eff = c_eff * alpha_u_raw
    Hu_eff      = Hu_factored / c_eff

    # Mutate mix in-place (engine reads mix.alpha_u:1513, mix.Hu_J_kg_effective:1510)
    mix.alpha_u           = alpha_u_eff
    mix.Hu_J_kg_effective = Hu_eff
    constr.model_soil     = False
    constr.is_submerged   = True

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )

    T0_C = (T_PL_F - 32.0) * 5.0 / 9.0
    Ts_C = (T_SO_F - 32.0) * 5.0 / 9.0
    Ti   = np.full((grid.ny, grid.nx), T0_C)
    if hasattr(grid, "is_soil") and grid.is_soil is not None:
        Ti[grid.is_soil] = Ts_C

    res = solve_hydration_2d(
        grid, mix, Ti,
        duration_s=DURATION_HR * 3600.0,
        output_interval_s=OUTPUT_INTERVAL_S,
        boundary_mode="full_2d",
        environment=make_neutral_env(T_PL_F),
        construction=constr,
        T_ground_deep_C=Ts_C,
        diagnostic_outputs=False,
    )
    t_hrs = np.asarray(res.t_s) / 3600.0

    # Resample engine to CW grid (same as run_engine_pilot.py lines 196–209)
    jsl, isl = grid.concrete_slice()
    eng_y = grid.y[jsl]
    eng_x = grid.x[isl]
    nD = len(cw["depths_m"])
    nW = len(cw["widths_m"])
    T_engine_F = np.zeros((len(t_hrs), nD, nW), dtype=np.float32)
    for ti in range(len(t_hrs)):
        slab_T_F = res.T_field_C[ti, jsl, isl] * 9.0 / 5.0 + 32.0
        T_engine_F[ti] = resample_engine_to_cw(
            eng_y, eng_x, slab_T_F, cw["depths_m"], cw["widths_m"]
        )

    scalars = compute_scalars(T_engine_F, t_hrs, cw)

    log = {
        "mix":          mix_name,
        "Hu_factor":    f"{fac:.4f}",
        "c_multiplier": f"{c_mult:.4f}",
        "c_eff":        f"{c_eff:.6f}",
        "alpha_u_raw":  f"{alpha_u_raw:.6f}",
        "Hu_raw":       f"{Hu_raw:.0f}",
        "alpha_u_eff":  f"{alpha_u_eff:.6f}",
        "Hu_eff":       f"{Hu_eff:.1f}",
        "conditional":  cond_txt,
        "t_hrs_count":  len(t_hrs),
    }

    return scalars, log


def main():
    total = len(MIXES) * len(HU_FACTOR_GRID) * len(C_MULT_GRID)
    print("Sprint 9 Stage 1-pilot — (Hu_factor, c) sensitivity sweep")
    print(f"Grid: {len(HU_FACTOR_GRID)} Hu_factor × {len(C_MULT_GRID)} c_multiplier × {len(MIXES)} mixes = {total} runs")
    print(f"C_ANCHOR = {C_ANCHOR:.4f}  (c_eff = c_multiplier × C_ANCHOR)")
    print(f"CW data: {PILOT}/cw_data/{{mix01,mix07}}/")
    print(f"Output:  {DATA_DIR}/")
    print()

    # Pre-load CW T-field once per mix — invariant across sweep
    cw_cache = {name: load_cw_once(folder) for name, folder in MIXES}

    residual_rows = []
    log_rows      = []
    run_n         = 0

    for mix_name, cw_folder in MIXES:
        cw = cw_cache[mix_name]
        for Hu_fac in HU_FACTOR_GRID:
            for c_mult in C_MULT_GRID:
                run_n += 1
                c_eff = c_mult * C_ANCHOR
                print(f"[{run_n:02d}/{total}] {mix_name.upper()}  "
                      f"Hu_factor={Hu_fac:.2f}  c_mult={c_mult:.2f}  c_eff={c_eff:.4f}",
                      flush=True)

                scalars, log = run_one(mix_name, cw_folder, Hu_fac, c_mult, cw)

                print(f"         max|R|={scalars['max_R_F']:.4f}°F @di={scalars['di_at_max']}  "
                      f"R47={scalars['R_di47_F']:+.4f}  R48={scalars['R_di48_F']:+.4f}  "
                      f"R36={scalars['R_di36_F']:+.4f}  → {scalars['pass']}", flush=True)

                residual_rows.append({
                    "mix":          mix_name,
                    "Hu_factor":    Hu_fac,
                    "c_multiplier": c_mult,
                    **scalars,
                })
                log_rows.append(log)

    # Write sweep_residual_table.csv
    res_path   = DATA_DIR / "sweep_residual_table.csv"
    res_fields = ["mix", "Hu_factor", "c_multiplier", "max_R_F", "di_at_max",
                  "R_di47_F", "R_di48_F", "R_di36_F", "pass"]
    with open(res_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=res_fields)
        w.writeheader()
        w.writerows(residual_rows)
    print(f"\nWrote {res_path}  ({len(residual_rows)} rows)")

    # Write sweep_wrapper_log.csv
    log_path   = DATA_DIR / "sweep_wrapper_log.csv"
    log_fields = ["mix", "Hu_factor", "c_multiplier", "c_eff", "alpha_u_raw",
                  "Hu_raw", "alpha_u_eff", "Hu_eff", "conditional", "t_hrs_count"]
    with open(log_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=log_fields)
        w.writeheader()
        w.writerows(log_rows)
    print(f"Wrote {log_path}  ({len(log_rows)} rows)")

    # Headline diagnostic
    print("\n" + "=" * 70)
    print("Headline diagnostic — R(di=47) statistics per mix across full 8×5 grid:")
    for mix_name, _ in MIXES:
        vals = [r["R_di47_F"] for r in residual_rows if r["mix"] == mix_name]
        arr  = np.array(vals)
        print(f"  {mix_name.upper()}: mean={arr.mean():+.4f}°F  std={arr.std():.4f}°F  "
              f"range=[{arr.min():+.4f}, {arr.max():+.4f}]°F")
    print("=" * 70)
    print(f"\nAll {run_n} runs complete. Next: run analyze_sweep.py")


if __name__ == "__main__":
    main()
