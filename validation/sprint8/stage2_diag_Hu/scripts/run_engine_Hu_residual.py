#!/usr/bin/env python3
"""Sprint 8 Stage 2-floor-test §2.2–§2.3 — engine reruns at Hu_residual.

Runs the engine on the 12 Sprint 8 Stage 2 datasets with
mix.Hu_J_kg_effective overridden to Hu_residual (the floor-equivalent
Hu computed in §2.1), bypassing the apr28 compute_hu_factor path.
k_uc factor stays at 0.96 (Sprint 7 calibration).

Computes Structure C residuals (R1, R2, R3) at t=168 hr against the
matching CW output.txt and writes:

  data/sprint8_corrected_residuals.csv          (R1/R2/R3 per dataset)
  data/T_fields_engine_Hu_residual/<label>.npz  (T(z, x, t) at 24/84/168 hr)

Run from validation/sprint8/stage2_diag_Hu/.
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE  = Path(__file__).resolve().parents[1]
REPO  = Path(__file__).resolve().parents[4]
STAGE = REPO / "validation" / "sprint8" / "stage2_calibration"
SC    = REPO / "validation" / "soil_calibration"

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STAGE))
sys.path.insert(0, str(SC))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx
from engine_runner import (
    R1_DI, R2_WI, R2_DI, R3_DI, R3_WI, GATE_F, COMPARE_HR,
)

CW_DATA = STAGE / "cw_data"
DATA    = HERE / "data"
TFIELDS = DATA / "T_fields_engine_Hu_residual"

HU_RESIDUAL_J_KG = 13641.5      # geometric mean of α_u=0.20 (13691) and α_u=0.80 (13592)
K_UC_FACTOR      = 0.96
DURATION_HR      = 168.0
OUTPUT_INTERVAL  = 1800.0
SNAPSHOT_HRS     = (24.0, 84.0, 168.0)

RUNS = [
    ("alpha02_A", "thermal_alpha02_A_73_73",  0.20,  73,  73),
    ("alpha02_F", "thermal_alpha02_F_73_45",  0.20,  73,  45),
    ("alpha02_I", "thermal_alpha02_I_100_73", 0.20, 100,  73),
    ("alpha04_A", "thermal_alpha04_A_73_73",  0.40,  73,  73),
    ("alpha04_F", "thermal_alpha04_F_73_45",  0.40,  73,  45),
    ("alpha04_I", "thermal_alpha04_I_100_73", 0.40, 100,  73),
    ("alpha06_A", "thermal_alpha06_A_73_73",  0.60,  73,  73),
    ("alpha06_F", "thermal_alpha06_F_73_45",  0.60,  73,  45),
    ("alpha06_I", "thermal_alpha06_I_100_73", 0.60, 100,  73),
    ("alpha08_A", "thermal_alpha08_A_73_73",  0.80,  73,  73),
    ("alpha08_F", "thermal_alpha08_F_73_45",  0.80,  73,  45),
    ("alpha08_I", "thermal_alpha08_I_100_73", 0.80, 100,  73),
]


def run_one(folder: Path, T_pl_F: float, T_so_F: float, factor: float = K_UC_FACTOR):
    """Engine solve for one dataset with mix.Hu_J_kg_effective = HU_RESIDUAL_J_KG.

    Returns: t_s, grid, T_field_C (n_t, ny, nx), eng_y, eng_x,
             eng_T_F_at_168, cw_T_F_at_168, cw_widths_m, cw_depths_m.
    """
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))

    # Override the apr28-corrected effective Hu with the floor-equivalent value.
    # mix.alpha_u, tau, beta, etc. stay at the dataset's native values.
    mix.Hu_factor_calibrated = 1.0
    mix.Hu_J_kg_effective    = HU_RESIDUAL_J_KG

    constr.model_soil   = False
    constr.is_submerged = True

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0 = (T_pl_F - 32.0) * 5.0 / 9.0
    Ts = (T_so_F - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0)
    Ti[grid.is_soil] = Ts

    original = te2d.K_UC_CALIBRATION_FACTOR_SPRINT7
    try:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = factor
        res = solve_hydration_2d(
            grid, mix, Ti,
            duration_s=DURATION_HR * 3600.0,
            output_interval_s=OUTPUT_INTERVAL,
            boundary_mode="full_2d",
            environment=make_neutral_env(T_pl_F),
            construction=constr,
            T_ground_deep_C=Ts,
            diagnostic_outputs=False,
        )
    finally:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = original

    jsl, isl = grid.concrete_slice()
    ti168 = nearest_time_idx(res.t_s, COMPARE_HR)
    eng_T_F_168 = res.T_field_C[ti168, jsl, isl] * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jsl]
    eng_x = grid.x[isl]

    # Load CW reference
    v  = parse_cw_temp_output(str(folder / "output.txt"))
    ti_cw = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    cw_T_F_168 = v.T_field_F[ti_cw]

    return {
        "t_s":        np.asarray(res.t_s),
        "T_field_C":  res.T_field_C,
        "grid":       grid,
        "eng_y":      eng_y,
        "eng_x":      eng_x,
        "eng_T_F_168":eng_T_F_168,
        "cw_T_F_168": cw_T_F_168,
        "cw_widths_m":v.widths_m,
        "cw_depths_m":v.depths_m,
        "cw_series":  v,
    }


def residuals_F(d, eng_y, eng_x, eng_T_F):
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_T_F,
                                      d["cw_depths_m"], d["cw_widths_m"])
    R = eng_on_cw - d["cw_T_F_168"]
    r1 = float(np.max(np.abs(R[R1_DI, :])))
    r2 = float(np.max(np.abs(R[R2_DI, R2_WI])))
    r3 = float(np.sqrt(np.mean(R[R3_DI, R3_WI] ** 2)))
    return r1, r2, r3, R, eng_on_cw


def save_t_field_snapshots(label: str, t_s: np.ndarray, T_field_C: np.ndarray,
                           grid):
    """Save concrete-region T (°C) at t = 24, 84, 168 hr."""
    jsl, isl = grid.concrete_slice()
    snaps = {}
    for hr in SNAPSHOT_HRS:
        ti = nearest_time_idx(t_s, hr)
        snaps[f"T_C_t{int(hr)}hr"] = T_field_C[ti, jsl, isl]
    snaps["depths_m"] = grid.y[jsl]
    snaps["widths_m"] = grid.x[isl]
    snaps["snapshot_hrs"] = np.array(SNAPSHOT_HRS)
    out = TFIELDS / f"{label}.npz"
    np.savez(out, **snaps)
    return out


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    TFIELDS.mkdir(parents=True, exist_ok=True)

    print(f"\n§2.2 Engine reruns at Hu_residual = {HU_RESIDUAL_J_KG} J/kg")
    print(f"     k_uc factor = {K_UC_FACTOR}")
    print("=" * 100)
    print(f"{'Dataset':<14} {'α_u':>5} {'T_pl':>5} {'T_so':>5}"
          f"  {'R1':>9}  {'R2':>9}  {'R3 RMSE':>9}"
          f"  {'gateR1':>6}  {'gateR2':>6}  {'pass'}")
    print("-" * 100)

    rows = []
    a_scenario_check_data = {}

    for label, folder, alpha_u, T_pl, T_so in RUNS:
        dpath = CW_DATA / folder
        print(f"  {label:<14} ...", end="", flush=True)
        d = run_one(dpath, T_pl, T_so)
        r1, r2, r3, _, _ = residuals_F(d, d["eng_y"], d["eng_x"], d["eng_T_F_168"])
        gate_r1 = r1 <= GATE_F
        gate_r2 = r2 <= GATE_F
        passes  = gate_r1 and gate_r2

        # Save snapshots
        save_t_field_snapshots(label, d["t_s"], d["T_field_C"], d["grid"])

        # For §2.6 A-scenario sanity (track per α_u)
        if "_A" in label:
            jsl, isl = d["grid"].concrete_slice()
            T_conc_C = d["T_field_C"][:, jsl, isl]
            nD = T_conc_C.shape[1]
            T_core_C = T_conc_C[:, nD // 2, -1]   # centerline mid-depth (engine convention)
            T_core_F = T_core_C * 9.0 / 5.0 + 32.0
            t_hr = d["t_s"] / 3600.0
            ti168 = nearest_time_idx(d["t_s"], COMPARE_HR)
            dT_eng_168_F = float(T_core_F[ti168] - T_core_F[0])
            # CW core warming for matching dataset
            cw = d["cw_series"]
            cw_t = cw.time_hrs
            cw_core_F = cw.T_core_center_F  # mid-depth × wi=0 (CW centerline)
            cw_ti168 = int(np.abs(cw_t - COMPARE_HR).argmin())
            dT_cw_168_F = float(cw_core_F[cw_ti168] - cw_core_F[0])
            a_scenario_check_data[label] = {
                "alpha_u":     alpha_u,
                "T_eng_0_F":   float(T_core_F[0]),
                "T_eng_168_F": float(T_core_F[ti168]),
                "dT_eng_168_F":dT_eng_168_F,
                "dT_cw_168_F": dT_cw_168_F,
                "dT_eng_C":    dT_eng_168_F * 5/9,
                "dT_cw_C":     dT_cw_168_F * 5/9,
            }

        verdict = "PASS" if passes else "FAIL"
        print(f"\r  {label:<14} {alpha_u:>5.2f} {T_pl:>5.0f} {T_so:>5.0f}"
              f"  {r1:>9.4f}  {r2:>9.4f}  {r3:>9.4f}"
              f"  {'✓' if gate_r1 else '✗':>6}  {'✓' if gate_r2 else '✗':>6}"
              f"  [{verdict}]")
        rows.append({
            "label": label, "folder": folder, "alpha_u": alpha_u,
            "T_pl_F": T_pl, "T_so_F": T_so,
            "Hu_J_kg_effective": HU_RESIDUAL_J_KG,
            "factor": K_UC_FACTOR,
            "r1_max_F": r1, "r2_max_F": r2, "r3_rmse_F": r3,
            "gate_r1": gate_r1, "gate_r2": gate_r2,
            "passes": passes, "minimax": max(r1, r2),
        })

    # Sanity §5: T(t=0) ≈ T_pl
    print("\n§5 sanity: engine T_core(t=0) per dataset")
    for label, info in a_scenario_check_data.items():
        T_pl = next(r["T_pl_F"] for r in rows if r["label"] == label)
        ok = abs(info["T_eng_0_F"] - T_pl) < 0.05
        print(f"  {label}: T_core(0)={info['T_eng_0_F']:.4f}  "
              f"(expected {T_pl}) {'✓' if ok else '✗ FAIL'}")

    # Side-by-side summary by α target
    print("\nGrouped summary by α_u:")
    for au in [0.20, 0.40, 0.60, 0.80]:
        group = [r for r in rows if abs(r["alpha_u"] - au) < 0.01]
        r1s = [r["r1_max_F"] for r in group]
        r2s = [r["r2_max_F"] for r in group]
        n_pass = sum(1 for r in group if r["passes"])
        print(f"  α_u={au:.2f}:  R1 max={max(r1s):.4f}  R2 max={max(r2s):.4f}"
              f"  ({n_pass}/{len(group)} pass gate)")

    # Write CSV
    out_csv = DATA / "sprint8_corrected_residuals.csv"
    fields = ["label", "folder", "alpha_u", "T_pl_F", "T_so_F",
              "Hu_J_kg_effective", "factor",
              "r1_max_F", "r2_max_F", "r3_rmse_F",
              "gate_r1", "gate_r2", "passes", "minimax"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_csv}")

    # Write A-scenario sanity table
    sanity_csv = DATA / "sprint8_floor_test_A_sanity.csv"
    sf = ["label", "alpha_u", "T_eng_0_F", "T_eng_168_F",
          "dT_eng_168_F", "dT_cw_168_F", "dT_eng_C", "dT_cw_C",
          "ratio_eng_over_cw", "abs_pct_disagreement"]
    with open(sanity_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sf)
        w.writeheader()
        for label, info in sorted(a_scenario_check_data.items()):
            ratio = info["dT_eng_168_F"] / info["dT_cw_168_F"] if info["dT_cw_168_F"] else float("nan")
            pct   = abs(ratio - 1.0) * 100
            row = {"label": label, **info,
                   "ratio_eng_over_cw": ratio,
                   "abs_pct_disagreement": pct}
            w.writerow(row)
    print(f"Wrote {sanity_csv}")
    print("\nDone.")


if __name__ == "__main__":
    main()
