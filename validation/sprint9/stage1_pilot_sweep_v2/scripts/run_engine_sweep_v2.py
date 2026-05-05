#!/usr/bin/env python3
"""Sprint 9 Stage 1-pilot — refined sensitivity sweep v2 (composition-centered Hu_factor).

Mix-specific Hu_factor grids (±1% around each mix's composition-correct value):
  MIX-01: Hu_factor ∈ {0.9414, 0.9464, 0.9514, 0.9564, 0.9614}
  MIX-07: Hu_factor ∈ {0.8846, 0.8896, 0.8946, 0.8996, 0.9046}

Shared c_multiplier grid: {0.95, 0.97, 1.00, 1.03, 1.05}

Total: 5 × 5 × 2 mixes = 50 runs. No .npz files written.

Expanded scalar set vs. v1: adds R_di42, stencil_drop (= R_di47 − R_di42),
T_core_eng, T_core_cw, and c_eff column.

Writes:
  data/sweep_v2_residual_table.csv  — 50 rows
  data/sweep_v2_wrapper_log.csv     — 50 rows; provenance

Usage:
    cd /Users/hqa668/calcshore-thermal-2d
    python validation/sprint9/stage1_pilot_sweep_v2/scripts/run_engine_sweep_v2.py
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE  = Path(__file__).resolve().parent
SWEEP = HERE.parent                               # stage1_pilot_sweep_v2/
PILOT = SWEEP.parent / "stage1_pilot"            # stage1_pilot/ (CW data source)
ROOT  = (SWEEP / "../../..").resolve()           # calcshore-thermal-2d/
SC    = ROOT / "validation" / "soil_calibration"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env

# Mix-specific Hu_factor grids — ±1% around each mix's composition-correct value
HU_FACTOR_BY_MIX = {
    "mix01": [0.9414, 0.9464, 0.9514, 0.9564, 0.9614],
    "mix07": [0.8846, 0.8896, 0.8946, 0.8996, 0.9046],
}
# Composition-correct reference points (centre of each grid)
HU_CORRECT = {"mix01": 0.9514, "mix07": 0.8946}

C_MULT_GRID = [0.95, 0.97, 1.00, 1.03, 1.05]

# Fixed constants — identical to stage1_pilot and sweep v1
HU_FLOOR_THRESHOLD = 10_000.0
HU_RESIDUAL        = 12_937.0
T_PL_F             = 73.0
T_SO_F             = 85.0
DURATION_HR        = 168.0
OUTPUT_INTERVAL_S  = 3600.0
C_ANCHOR           = 1.0532

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
        "T_cw_F":   v.T_field_F,
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
    T_cw_aln = align_cw_to_engine(t_hrs, cw["t_hrs"], cw["T_cw_F"])
    ti168 = ti_near(t_hrs, 168.0)
    ti0   = ti_near(t_hrs, 0.0)

    R_t168 = T_engine_F[ti168, :, 0] - T_cw_aln[ti168, :, 0]
    R_bottom  = R_t168[DI_CORE:DI_BOT + 1]
    max_R     = float(np.max(np.abs(R_bottom)))
    di_at_max = DI_CORE + int(np.argmax(np.abs(R_bottom)))

    R_di36 = float(R_t168[36])
    R_di42 = float(R_t168[42])
    R_di47 = float(R_t168[47])
    R_di48 = float(R_t168[48])

    T_core_eng = float(T_engine_F[ti168, DI_CORE, 0])
    T_core_cw  = float(T_cw_aln[ti168, DI_CORE, 0])
    T_ic_eng   = float(T_engine_F[ti0, DI_CORE, 0])
    T_bot_eng  = float(T_engine_F[ti168, DI_BOT, 0])

    return {
        "max_R":        max_R,
        "di_at_max":    di_at_max,
        "R_di36":       R_di36,
        "R_di42":       R_di42,
        "R_di47":       R_di47,
        "R_di48":       R_di48,
        "stencil_drop": R_di47 - R_di42,
        "T_core_eng":   T_core_eng,
        "T_core_cw":    T_core_cw,
        # Internals for sanity (not written to CSV)
        "_T_ic_eng":    T_ic_eng,
        "_T_bot_eng":   T_bot_eng,
        "pass":         "PASS" if max_R < GATE_F else "FAIL",
    }


def run_sanity_check(scalars, mix_name):
    """Halt-on-failure sanity check. Called on first run per mix only."""
    ic_ok   = abs(scalars["_T_ic_eng"] - T_PL_F) <= 0.05
    bot_ok  = abs(scalars["_T_bot_eng"] - T_SO_F) <= 0.5
    core_ok = scalars["T_core_eng"] > 110.0

    print(f"  Sanity check ({mix_name.upper()}):")
    print(f"    T_IC  (di=24,t=0)   = {scalars['_T_ic_eng']:.4f}°F  "
          f"(expect {T_PL_F}±0.05°F) → {'OK' if ic_ok else 'FAIL'}")
    print(f"    T_bot (di=48,t=168) = {scalars['_T_bot_eng']:.4f}°F  "
          f"(expect {T_SO_F}±0.5°F)  → {'OK' if bot_ok else 'FAIL'}")
    print(f"    T_core(di=24,t=168) = {scalars['T_core_eng']:.4f}°F  "
          f"(expect > 110°F)         → {'OK' if core_ok else 'FAIL'}")

    if not ic_ok:
        raise RuntimeError(
            f"HALT: IC sanity failed for {mix_name}: "
            f"T_IC={scalars['_T_ic_eng']:.4f}°F vs {T_PL_F}±0.05°F"
        )
    if not bot_ok:
        raise RuntimeError(
            f"HALT: bottom Dirichlet sanity failed for {mix_name}: "
            f"T_bot={scalars['_T_bot_eng']:.4f}°F vs {T_SO_F}±0.5°F"
        )
    if not core_ok:
        raise RuntimeError(
            f"HALT: core-temperature sanity failed for {mix_name}: "
            f"T_core={scalars['T_core_eng']:.4f}°F (heat path may not have fired)"
        )
    print(f"    All sanity checks PASSED.")


def run_one(mix_name, cw_folder, Hu_fac, c_mult, cw, do_sanity=False):
    mix, geom, constr, _ = parse_cw_dat(str(cw_folder / "input.dat"))

    alpha_u_raw = mix.alpha_u
    Hu_raw      = mix.Hu_J_kg

    # Hu_residual conditional (passthrough for both pilot mixes)
    if Hu_raw < HU_FLOOR_THRESHOLD:
        Hu_base  = HU_RESIDUAL
        cond_txt = f"APPLIED (Hu_raw={Hu_raw:.0f} < {HU_FLOOR_THRESHOLD:.0f})"
    else:
        Hu_base  = Hu_raw
        cond_txt = "passthrough"

    c_eff       = c_mult * C_ANCHOR
    Hu_factored = Hu_base * Hu_fac
    alpha_u_eff = c_eff * alpha_u_raw
    Hu_eff      = Hu_factored / c_eff

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

    if do_sanity:
        run_sanity_check(scalars, mix_name)

    log = {
        "mix":          mix_name,
        "Hu_factor":    f"{Hu_fac:.4f}",
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
    total = sum(len(v) for v in HU_FACTOR_BY_MIX.values()) * len(C_MULT_GRID)
    print("Sprint 9 Stage 1-pilot — refined (Hu_factor, c) sensitivity sweep v2")
    print(f"Grid: mix-specific 5-point Hu_factor × {len(C_MULT_GRID)} c_multiplier × 2 mixes = {total} runs")
    print(f"C_ANCHOR = {C_ANCHOR:.4f}")
    print(f"MIX-01 Hu grid: {HU_FACTOR_BY_MIX['mix01']}  (correct = {HU_CORRECT['mix01']})")
    print(f"MIX-07 Hu grid: {HU_FACTOR_BY_MIX['mix07']}  (correct = {HU_CORRECT['mix07']})")
    print()

    cw_cache = {name: load_cw_once(folder) for name, folder in MIXES}

    residual_rows = []
    log_rows      = []
    run_n         = 0
    t_wall_start  = time.perf_counter()

    RES_FIELDS = ["mix", "Hu_factor", "c_multiplier", "c_eff", "max_R", "di_at_max",
                  "R_di36", "R_di42", "R_di47", "R_di48", "stencil_drop",
                  "T_core_eng", "T_core_cw"]

    for mix_name, cw_folder in MIXES:
        cw = cw_cache[mix_name]
        for i, Hu_fac in enumerate(HU_FACTOR_BY_MIX[mix_name]):
            for j, c_mult in enumerate(C_MULT_GRID):
                run_n += 1
                do_sanity = (i == 0 and j == 0)
                c_eff = c_mult * C_ANCHOR
                print(f"[{run_n:02d}/{total}] {mix_name.upper()}  "
                      f"Hu_factor={Hu_fac:.4f}  c_mult={c_mult:.2f}  c_eff={c_eff:.4f}",
                      flush=True)

                scalars, log = run_one(mix_name, cw_folder, Hu_fac, c_mult, cw,
                                       do_sanity=do_sanity)

                print(f"         max|R|={scalars['max_R']:.4f}°F @di={scalars['di_at_max']}  "
                      f"R36={scalars['R_di36']:+.4f}  R42={scalars['R_di42']:+.4f}  "
                      f"R47={scalars['R_di47']:+.4f}  R48={scalars['R_di48']:+.4f}  "
                      f"drop={scalars['stencil_drop']:+.4f}  "
                      f"T_core={scalars['T_core_eng']:.2f}°F  → {scalars['pass']}", flush=True)

                row = {
                    "mix":          mix_name,
                    "Hu_factor":    Hu_fac,
                    "c_multiplier": c_mult,
                    "c_eff":        round(c_eff, 6),
                    "max_R":        scalars["max_R"],
                    "di_at_max":    scalars["di_at_max"],
                    "R_di36":       scalars["R_di36"],
                    "R_di42":       scalars["R_di42"],
                    "R_di47":       scalars["R_di47"],
                    "R_di48":       scalars["R_di48"],
                    "stencil_drop": scalars["stencil_drop"],
                    "T_core_eng":   scalars["T_core_eng"],
                    "T_core_cw":    scalars["T_core_cw"],
                }
                residual_rows.append(row)
                log_rows.append(log)

    elapsed = time.perf_counter() - t_wall_start

    # Write sweep_v2_residual_table.csv
    res_path = DATA_DIR / "sweep_v2_residual_table.csv"
    with open(res_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=RES_FIELDS)
        w.writeheader()
        w.writerows(residual_rows)
    print(f"\nWrote {res_path}  ({len(residual_rows)} rows)")

    # Write sweep_v2_wrapper_log.csv
    log_path   = DATA_DIR / "sweep_v2_wrapper_log.csv"
    log_fields = ["mix", "Hu_factor", "c_multiplier", "c_eff", "alpha_u_raw",
                  "Hu_raw", "alpha_u_eff", "Hu_eff", "conditional", "t_hrs_count"]
    with open(log_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=log_fields)
        w.writeheader()
        w.writerows(log_rows)
    print(f"Wrote {log_path}  ({len(log_rows)} rows)")

    # Quick headline diagnostic
    print("\n" + "=" * 70)
    print(f"Runtime: {elapsed:.1f} s total  ({elapsed/run_n:.2f} s/run)")
    print("Headline stencil_drop statistics per mix (all 25 grid points each):")
    for mix_name, _ in MIXES:
        drops = np.array([r["stencil_drop"] for r in residual_rows if r["mix"] == mix_name])
        print(f"  {mix_name.upper()}: mean={drops.mean():+.4f}°F  "
              f"std={drops.std():.4f}°F  range=[{drops.min():+.4f}, {drops.max():+.4f}]°F")
    print("=" * 70)
    print(f"\nAll {run_n} runs complete. Next: run analyze_sweep_v2.py")


if __name__ == "__main__":
    main()
