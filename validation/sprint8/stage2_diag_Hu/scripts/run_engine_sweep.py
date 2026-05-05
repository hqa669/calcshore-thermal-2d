#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-Hu §4.9.1–§4.9.2 — engine A-scenario sweep.

Runs the engine on the same 15 (Hu, α_u) points as §4.1 in the A scenario:
T_pl = T_soil = 73°F, ambient = 73°F, model_soil=False, is_submerged=True,
blanket_thickness=0, k_uc factor = 0.96.

The engine reads `mix.Hu_J_kg_effective` (NOT mix.Hu_J_kg) and `mix.alpha_u`.
We bypass the apr28 `compute_hu_factor` correction and assign these directly
per sweep point.

Writes:
  data/T_core_trajectories_engine.npz  — raw (time_hrs, T_core_F) per run
  data/Hu_engine_table.csv             — engine ΔT metrics, one row per run
  data/Hu_residual_table_full.csv      — CW + engine + disagreement, joined

Run from validation/sprint8/stage2_diag_Hu/.
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "soil_calibration"))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat
from stage4b_run import make_neutral_env, nearest_time_idx

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"

# Use any one of the 15 input files as the base mix; α_u and Hu are overridden per point.
BASE_INPUT = CW_DATA / "thermal_Hu_res_au02_Hu1000" / "input.dat"
if not BASE_INPUT.exists():
    BASE_INPUT = CW_DATA / "thermal_Hu_res_au02_Hu1000" / "input_test.dat"

T_PL_F   = 73.0
T_SOIL_F = 73.0
DURATION_HR     = 168.0
OUTPUT_INTERVAL = 1800.0   # 30-min cadence → 337 samples
K_UC_FACTOR     = 0.96     # Sprint 7 calibrated value (matches engine default)

# 15 sweep points — match the 15 cw_data folders
SWEEP = [
    ("au02_Hu001",   0.20, 0.01),
    ("au02_Hu1",     0.20, 1.0),
    ("au02_Hu10",    0.20, 10.0),
    ("au02_Hu100",   0.20, 100.0),
    ("au02_Hu1000",  0.20, 1000.0),
    ("au02_Hu10000", 0.20, 10000.0),
    ("au02_Hu50000", 0.20, 50000.0),
    ("au02_Hu100000",0.20, 100000.0),
    ("au08_Hu1",     0.80, 1.0),
    ("au08_Hu10",    0.80, 10.0),
    ("au08_Hu100",   0.80, 100.0),
    ("au08_Hu1000",  0.80, 1000.0),
    ("au08_Hu10000", 0.80, 10000.0),
    ("au08_Hu50000", 0.80, 50000.0),
    ("au08_Hu100000",0.80, 100000.0),
]


def run_one(label: str, alpha_u: float, Hu_J_kg: float):
    """Run engine for one (α_u, Hu) point. Returns (time_hrs, T_core_F)."""
    mix, geom, constr, _ = parse_cw_dat(str(BASE_INPUT))

    # Override sweep parameters. The solver consumes mix.Hu_J_kg_effective
    # and mix.alpha_u; mix.Hu_J_kg is not consumed by the solver but we
    # update it for logging consistency.
    mix.alpha_u = alpha_u
    mix.Hu_J_kg = Hu_J_kg
    mix.Hu_J_kg_effective = Hu_J_kg
    mix.Hu_factor_calibrated = 1.0

    # A-scenario BC discipline
    constr.model_soil   = False
    constr.is_submerged = True

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )

    T0_C = (T_PL_F   - 32.0) * 5.0 / 9.0
    Ts_C = (T_SOIL_F - 32.0) * 5.0 / 9.0
    Ti_C = np.full((grid.ny, grid.nx), T0_C)
    Ti_C[grid.is_soil] = Ts_C

    # Apply k_uc factor (default is 0.96 — set explicitly for clarity)
    original = te2d.K_UC_CALIBRATION_FACTOR_SPRINT7
    try:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = K_UC_FACTOR
        res = solve_hydration_2d(
            grid, mix, Ti_C,
            duration_s=DURATION_HR * 3600.0,
            output_interval_s=OUTPUT_INTERVAL,
            boundary_mode="full_2d",
            environment=make_neutral_env(T_PL_F),
            construction=constr,
            T_ground_deep_C=Ts_C,
            diagnostic_outputs=False,
        )
    finally:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = original

    # Extract centerline mid-depth trajectory
    jsl, isl = grid.concrete_slice()
    T_conc_C = res.T_field_C[:, jsl, isl]      # (n_t, nD, nW)
    nD = T_conc_C.shape[1]
    # In the concrete sub-grid, x runs soil-contact (i=0) → centerline (i=-1).
    # Centerline mid-depth = [:, nD//2, -1].
    T_core_C = T_conc_C[:, nD // 2, -1]
    T_core_F = T_core_C * 9.0 / 5.0 + 32.0
    time_hrs = np.asarray(res.t_s) / 3600.0
    return time_hrs, T_core_F


def metrics(time_hrs: np.ndarray, T_F: np.ndarray) -> dict:
    mask = time_hrs <= DURATION_HR + 0.01
    t = time_hrs[mask]
    T = T_F[mask]
    T0   = float(T[0])
    T168 = float(np.interp(DURATION_HR, t, T))
    dT_168 = T168 - T0
    dT_max = float(np.max(T)) - T0
    t_peak = float(t[np.argmax(T)])
    return {
        "T0_F":      T0,
        "T168_F":    T168,
        "dT_168_F":  dT_168,
        "dT_max_F":  dT_max,
        "t_peak_hr": t_peak,
        "T0_C":      (T0 - 32) * 5/9,
        "T168_C":    (T168 - 32) * 5/9,
        "dT_168_C":  dT_168 * 5/9,
        "dT_max_C":  dT_max * 5/9,
    }


def main():
    if not BASE_INPUT.exists():
        print(f"ERROR: base input not found at {BASE_INPUT}", file=sys.stderr)
        sys.exit(1)
    print(f"Base mix loaded from: {BASE_INPUT}")
    print(f"k_uc factor: {K_UC_FACTOR}")
    print(f"Default engine k_uc factor: {te2d.K_UC_CALIBRATION_FACTOR_SPRINT7}")
    print(f"Output cadence: {OUTPUT_INTERVAL}s ({OUTPUT_INTERVAL/60:.0f} min)\n")

    rows = []
    trajs = {}
    for label, au, hu in SWEEP:
        print(f"  Engine run {label:<16} α_u={au:.2f} Hu={hu:>9.2g} ...", end="", flush=True)
        time_hrs, T_F = run_one(label, au, hu)
        m = metrics(time_hrs, T_F)
        row = {
            "label":   label,
            "alpha_u": au,
            "Hu_J_kg": hu,
            **m,
        }
        rows.append(row)
        trajs[label] = {"time_hrs": time_hrs, "T_core_F": T_F}
        print(f"  T0={m['T0_F']:.3f}°F  dT_168={m['dT_168_F']:>7.4f}°F "
              f"({m['dT_168_C']:>6.4f}°C)  t_peak={m['t_peak_hr']:.1f}hr")

    # Sanity checks (§5)
    print("\nRunning §5 sanity checks ...")
    sane = True
    for r in rows:
        if abs(r["T0_F"] - T_PL_F) > 0.05:
            print(f"  WARN: {r['label']} T_core(t=0)={r['T0_F']:.4f}°F (expected {T_PL_F})", file=sys.stderr)
            sane = False
    # Hu near 0 should give near-zero ΔT
    for r in rows:
        if r["Hu_J_kg"] < 1.0 and r["dT_168_C"] > 0.05:
            print(f"  WARN: {r['label']} Hu={r['Hu_J_kg']} but engine dT_168={r['dT_168_C']:.4f}°C "
                  "(expected ≤ 0.05°C — check engine kinetics path)", file=sys.stderr)
    # Monotonicity per α_u
    for au in (0.20, 0.80):
        sub = sorted([r for r in rows if abs(r["alpha_u"] - au) < 0.05], key=lambda r: r["Hu_J_kg"])
        for i in range(1, len(sub)):
            if sub[i]["dT_168_F"] < sub[i-1]["dT_168_F"] - 0.01:
                print(f"  WARN: monotonicity α_u≈{au}: Hu={sub[i-1]['Hu_J_kg']} dT={sub[i-1]['dT_168_F']:.3f} "
                      f"> Hu={sub[i]['Hu_J_kg']} dT={sub[i]['dT_168_F']:.3f}", file=sys.stderr)
                sane = False
    if sane:
        print("  All sanity checks passed.")

    # Save engine trajectories
    npz_data = {}
    for label, d in trajs.items():
        npz_data[label + "_time_hrs"] = d["time_hrs"]
        npz_data[label + "_T_core_F"] = d["T_core_F"]
    npz_path = DATA / "T_core_trajectories_engine.npz"
    np.savez(npz_path, **npz_data)
    print(f"\nWrote {npz_path}")

    # Save engine table
    eng_path = DATA / "Hu_engine_table.csv"
    fields = ["label", "alpha_u", "Hu_J_kg",
              "T0_F", "T168_F", "dT_168_F", "dT_max_F", "t_peak_hr",
              "T0_C", "T168_C", "dT_168_C", "dT_max_C"]
    with open(eng_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["alpha_u"], r["Hu_J_kg"])):
            w.writerow({k: r[k] for k in fields})
    print(f"Wrote {eng_path}")

    # Build joined CW + engine table
    cw_rows = {}
    with open(DATA / "Hu_residual_table.csv", newline="") as f:
        for r in csv.DictReader(f):
            au = float(r["alpha_u"])
            hu = float(r["Hu_J_kg"])
            cw_rows[(round(au, 3), hu)] = r

    full_path = DATA / "Hu_residual_table_full.csv"
    full_fields = ["alpha_u", "Hu_J_kg",
                   "dT_CW_168_C", "dT_eng_168_C", "disagree_C",
                   "dT_CW_168_F", "dT_eng_168_F", "disagree_F",
                   "dT_CW_max_C", "dT_eng_max_C",
                   "t_peak_CW_hr", "t_peak_eng_hr"]
    with open(full_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=full_fields)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["alpha_u"], r["Hu_J_kg"])):
            cw = cw_rows.get((round(r["alpha_u"], 3), r["Hu_J_kg"]))
            if cw is None:
                continue
            cw_dT_C = float(cw["dT_168_C"])
            cw_dT_F = float(cw["dT_168_F"])
            row = {
                "alpha_u":       r["alpha_u"],
                "Hu_J_kg":       r["Hu_J_kg"],
                "dT_CW_168_C":   cw_dT_C,
                "dT_eng_168_C":  r["dT_168_C"],
                "disagree_C":    cw_dT_C - r["dT_168_C"],
                "dT_CW_168_F":   cw_dT_F,
                "dT_eng_168_F":  r["dT_168_F"],
                "disagree_F":    cw_dT_F - r["dT_168_F"],
                "dT_CW_max_C":   float(cw["dT_max_C"]),
                "dT_eng_max_C":  r["dT_max_C"],
                "t_peak_CW_hr":  float(cw["t_peak_hr"]),
                "t_peak_eng_hr": r["t_peak_hr"],
            }
            w.writerow(row)
    print(f"Wrote {full_path}")

    # Print side-by-side summary
    print("\n§4.9.3 Side-by-side comparison:")
    print(f"{'α_u':>5} {'Hu':>9}  {'CW dT°C':>8}  {'eng dT°C':>9}  {'disagree°C':>10}  {'disagree°F':>10}")
    print("-" * 65)
    for r in sorted(rows, key=lambda r: (r["alpha_u"], r["Hu_J_kg"])):
        cw = cw_rows.get((round(r["alpha_u"], 3), r["Hu_J_kg"]))
        if cw is None:
            continue
        cw_dT_C = float(cw["dT_168_C"])
        cw_dT_F = float(cw["dT_168_F"])
        d_C     = cw_dT_C - r["dT_168_C"]
        d_F     = cw_dT_F - r["dT_168_F"]
        print(f"{r['alpha_u']:>5.2f} {r['Hu_J_kg']:>9.2g}  {cw_dT_C:>8.4f}  {r['dT_168_C']:>9.4f}  "
              f"{d_C:>10.4f}  {d_F:>10.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
