#!/usr/bin/env python3
"""§4.11.8 Task 2 — Apply c(T_pl) correction to 12 Sprint 8 Stage 2 datasets.

For each of the 12 (α_u_nominal × {A, F, I}) datasets in
`validation/sprint8/stage2_calibration/cw_data/`, runs the engine at:
  - α_u_engine = c(T_pl) × α_u_nominal      (c from §4.11.8 chosen quadratic fit)
  - Hu_J_kg_effective = HU_EFF (12,936.96 J/kg, §4.11 calibration anchor)
  - K_UC_CALIBRATION_FACTOR_SPRINT7 = 0.96  (Sprint 7 calibration, no engine source change)
  - model_soil = False, is_submerged = True, blanket_thickness_m = 0.0

Computes Sprint 7 Structure C residuals (R1, R2, R3) at t = 168 hr.

Compares to:
  - Sprint 7 / original Stage 2 (Hu=1, c=1)         → from existing literature in report
  - Stage 2-floor-test (Hu_residual, c=1)            → from
        validation/sprint8/stage2_diag_Hu/data/sprint8_corrected_residuals.csv
  - This stage (Hu_eff = 12,937 J/kg, c=c(T_pl))     → produced here

Outputs:
  data/sprint8_corrected_v3_residuals.csv           — 12-dataset before/after table
  figures/sprint8_residuals_three_way.png           — bar chart R1/R2 across 12 datasets
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE  = Path(__file__).resolve().parents[1]
DATA  = HERE / "data"
FIGS  = HERE / "figures"
REPO  = Path(__file__).resolve().parents[4]
STAGE = REPO / "validation" / "sprint8" / "stage2_calibration"
SC    = REPO / "validation" / "soil_calibration"
DIAG  = REPO / "validation" / "sprint8" / "stage2_diag_Hu"

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STAGE))
sys.path.insert(0, str(SC))
sys.path.insert(0, str(DIAG / "scripts"))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx
from engine_runner import R1_DI, R2_WI, R2_DI, R3_DI, R3_WI, GATE_F, COMPARE_HR

CW_DATA = STAGE / "cw_data"

HU_EFF_J_KG     = 12936.96     # §4.11 calibration anchor (70°F α_u=0.20 probe)
K_UC_FACTOR     = 0.96
DURATION_HR     = 168.0
OUTPUT_INTERVAL = 1800.0

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


def load_chosen_fit():
    """Read the chosen c(T_pc) form and parameters from §4.11.8 Task 1."""
    with open(DATA / "c_T_pc_fit_chosen.csv") as f:
        rows = list(csv.reader(f))
    info = {r[0]: r[1] for r in rows if len(r) >= 2 and r[0]}
    form   = info["form"]
    params = json.loads(info["params_json"])
    if form != "B_quadratic":
        raise RuntimeError(f"Expected B_quadratic chosen form; got {form}")
    a0, a1, a2 = params
    def c_of_T(T_F):
        T = np.asarray(T_F, dtype=float)
        return a0 + a1 * T + a2 * T ** 2
    return form, params, c_of_T


def run_one_corrected(folder: Path, alpha_u_nom: float, c_factor: float,
                      T_pl_F: float, T_so_F: float):
    """Engine solve with mix.alpha_u = alpha_u_nom × c, mix.Hu = HU_EFF."""
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))

    mix.Hu_factor_calibrated = 1.0
    mix.Hu_J_kg_effective    = HU_EFF_J_KG
    mix.alpha_u              = alpha_u_nom * c_factor

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
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = K_UC_FACTOR
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
    eng_T_F = res.T_field_C[ti168, jsl, isl] * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jsl]
    eng_x = grid.x[isl]

    v = parse_cw_temp_output(str(folder / "output.txt"))
    ti_cw = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    cw_T_F = v.T_field_F[ti_cw]

    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_T_F, v.depths_m, v.widths_m)
    R = eng_on_cw - cw_T_F
    r1 = float(np.max(np.abs(R[R1_DI, :])))
    r2 = float(np.max(np.abs(R[R2_DI, R2_WI])))
    r3 = float(np.sqrt(np.mean(R[R3_DI, R3_WI] ** 2)))
    return r1, r2, r3


def load_floor_test():
    """Pull R1/R2 from §2.2 Stage 2-floor-test (Hu_residual, c=1) for comparison."""
    csv_path = DIAG / "data" / "sprint8_corrected_residuals.csv"
    out = {}
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            out[r["label"]] = (float(r["r1_max_F"]), float(r["r2_max_F"]),
                               float(r["r3_rmse_F"]))
    return out


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    form, params, c_of_T = load_chosen_fit()
    floor = load_floor_test()

    print(f"§4.11.8 Task 2 — Apply c(T_pl) to 12 Sprint 8 datasets")
    print(f"  Chosen fit: {form}, params = {params}")
    print(f"  Hu_J_kg_effective = {HU_EFF_J_KG} J/kg (§4.11 calibration anchor)")
    print(f"  k_uc factor = {K_UC_FACTOR}")
    print(f"  c(73°F) = {float(c_of_T(73)):.4f}")
    print(f"  c(100°F) = {float(c_of_T(100)):.4f}")
    print("=" * 100)

    rows = []
    print(f"\n{'Dataset':<14}  {'α_u_nom':>7}  {'T_pl':>4}  {'c(T_pl)':>8}"
          f"  {'α_u_eng':>8}  {'R1':>9}  {'R2':>9}  {'R3':>9}  pass")
    print("-" * 100)
    for label, folder, alpha_u_nom, T_pl, T_so in RUNS:
        c_val = float(c_of_T(T_pl))
        alpha_u_eng = alpha_u_nom * c_val
        r1, r2, r3 = run_one_corrected(CW_DATA / folder, alpha_u_nom, c_val,
                                       T_pl, T_so)
        passes = (r1 <= GATE_F) and (r2 <= GATE_F)
        verdict = "PASS" if passes else "FAIL"
        print(f"  {label:<14}  {alpha_u_nom:>7.2f}  {T_pl:>4.0f}  {c_val:>8.4f}"
              f"  {alpha_u_eng:>8.4f}  {r1:>9.4f}  {r2:>9.4f}  {r3:>9.4f}  [{verdict}]")
        floor_r1, floor_r2, floor_r3 = floor.get(label, (np.nan, np.nan, np.nan))
        rows.append({
            "label":          label,
            "folder":         folder,
            "alpha_u_nom":    alpha_u_nom,
            "T_pl_F":         T_pl,
            "T_so_F":         T_so,
            "c_T_pl":         c_val,
            "alpha_u_eng":    alpha_u_eng,
            "Hu_J_kg":        HU_EFF_J_KG,
            "factor":         K_UC_FACTOR,
            "r1_v3_F":        r1,
            "r2_v3_F":        r2,
            "r3_v3_F":        r3,
            "r1_floor_F":     floor_r1,
            "r2_floor_F":     floor_r2,
            "r3_floor_F":     floor_r3,
            "delta_r1":       r1 - floor_r1,
            "delta_r2":       r2 - floor_r2,
            "gate_r1":        r1 <= GATE_F,
            "gate_r2":        r2 <= GATE_F,
            "passes":         passes,
        })

    # Per-scenario summary
    print("\n" + "=" * 100)
    print("Per-scenario summary (gate = R1, R2 ≤ 0.35°F):")
    print("=" * 100)
    print(f"{'Scenario':>10}  {'#datasets':>9}  {'#pass v3':>8}  "
          f"{'#pass floor':>11}  {'worst_R1 v3':>12}  {'worst_R2 v3':>12}")
    for sc in ["A", "F", "I"]:
        group = [r for r in rows if f"_{sc}" in r["label"]]
        n_pass_v3 = sum(1 for r in group if r["passes"])
        n_pass_fl = sum(1 for r in group if r["r1_floor_F"] <= GATE_F
                        and r["r2_floor_F"] <= GATE_F)
        worst_r1 = max(r["r1_v3_F"] for r in group)
        worst_r2 = max(r["r2_v3_F"] for r in group)
        print(f"  {sc:>8}-scenario  {len(group):>9}  {n_pass_v3:>8}  "
              f"{n_pass_fl:>11}  {worst_r1:>12.4f}  {worst_r2:>12.4f}")
    n_pass_total = sum(1 for r in rows if r["passes"])
    print(f"\nOverall: {n_pass_total}/12 pass gate (R1 & R2 ≤ {GATE_F}°F)")

    # Improvement vs floor-test for I-scenarios
    print("\nI-scenario improvement vs Stage 2-floor-test (Hu_res, c=1):")
    for r in rows:
        if "_I" in r["label"]:
            d1 = r["delta_r1"]
            d2 = r["delta_r2"]
            pct1 = 100.0 * (-d1) / r["r1_floor_F"] if r["r1_floor_F"] > 0 else 0
            pct2 = 100.0 * (-d2) / r["r2_floor_F"] if r["r2_floor_F"] > 0 else 0
            print(f"  {r['label']}: R1 {r['r1_floor_F']:.3f} → {r['r1_v3_F']:.3f} "
                  f"({pct1:+.1f}%)  R2 {r['r2_floor_F']:.3f} → {r['r2_v3_F']:.3f} "
                  f"({pct2:+.1f}%)")

    # CSV
    out_csv = DATA / "sprint8_corrected_v3_residuals.csv"
    fields = ["label", "folder", "alpha_u_nom", "T_pl_F", "T_so_F",
              "c_T_pl", "alpha_u_eng", "Hu_J_kg", "factor",
              "r1_v3_F", "r2_v3_F", "r3_v3_F",
              "r1_floor_F", "r2_floor_F", "r3_floor_F",
              "delta_r1", "delta_r2",
              "gate_r1", "gate_r2", "passes"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_csv}")

    # Figure: three-way bar chart
    labels = [r["label"] for r in rows]
    r1_v3  = [r["r1_v3_F"] for r in rows]
    r2_v3  = [r["r2_v3_F"] for r in rows]
    r1_fl  = [r["r1_floor_F"] for r in rows]
    r2_fl  = [r["r2_floor_F"] for r in rows]

    fig, axes = plt.subplots(2, 1, figsize=(13, 7), sharex=True)
    x = np.arange(len(labels))
    w = 0.40

    ax = axes[0]
    ax.bar(x - w / 2, r1_fl, w, color="#888", label="Stage 2-floor-test (Hu_res, c=1)")
    ax.bar(x + w / 2, r1_v3, w, color="#d62728",
           label="§4.11.8 (Hu_eff, c(T_pl))")
    ax.axhline(GATE_F, color="black", linestyle="--", linewidth=1.0,
               label=f"gate = {GATE_F}°F")
    ax.set_ylabel("R1 max|R| (°F)", fontsize=11)
    ax.set_title("§4.11.8 Task 2 — R1/R2 residuals across 12 Sprint 8 Stage 2 datasets",
                 fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)

    ax = axes[1]
    ax.bar(x - w / 2, r2_fl, w, color="#888", label="Stage 2-floor-test (Hu_res, c=1)")
    ax.bar(x + w / 2, r2_v3, w, color="#d62728",
           label="§4.11.8 (Hu_eff, c(T_pl))")
    ax.axhline(GATE_F, color="black", linestyle="--", linewidth=1.0)
    ax.set_ylabel("R2 max|R| (°F)", fontsize=11)
    ax.set_xlabel("Dataset", fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.grid(True, axis="y", linestyle="--", alpha=0.4)

    fig.tight_layout()
    fig.savefig(FIGS / "sprint8_residuals_three_way.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGS / 'sprint8_residuals_three_way.png'}")
    print("\nDone.")


if __name__ == "__main__":
    main()
