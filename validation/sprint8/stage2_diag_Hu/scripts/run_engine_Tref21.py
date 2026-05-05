#!/usr/bin/env python3
"""Sprint 8 Stage 2-floor-test follow-up — T_ref hypothesis test.

Re-runs the 12 Sprint 8 Stage 2 datasets with the engine's Schindler-Folliard
Arrhenius reference temperature changed from 23°C → 21°C (T_REF_K: 296.15 → 294.15).
Wrapper override only (monkey-patch te2d.T_REF_K). Hu_J_kg_effective stays at
13,641.5 J/kg, k_uc multiplier stays at 0.96.

Strategy per the brief: run the 4 I scenarios first (these are the failing
datasets from §4.10.3). If residuals drop substantially, re-run A and F too
to check for regression.

Reads:
  ../stage2_calibration/cw_data/<folder>/{input.dat,output.txt}

Writes:
  data/sprint8_Tref21_residuals_I.csv         (always)
  data/sprint8_Tref21_residuals_AF.csv        (only if I improves)
  data/sprint8_Tref21_residuals_full.csv      (combined, if AF run)

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
sys.path.insert(0, str(HERE / "scripts"))

import thermal_engine_2d as te2d
from run_engine_Hu_residual import (
    run_one, residuals_F, HU_RESIDUAL_J_KG, K_UC_FACTOR, GATE_F, CW_DATA,
)

DATA = HERE / "data"
T_REF_K_NEW = 294.15      # 21°C — hypothesis test
T_REF_K_OLD = 296.15      # 23°C — committed default

# Prior residuals from §4.10.3 (Hu_residual, k_uc=0.96, T_ref=23°C)
PRIOR = {
    "alpha02_A": (0.0177, 0.0163),
    "alpha02_F": (0.1374, 0.1393),
    "alpha02_I": (0.2019, 0.2296),
    "alpha04_A": (0.0208, 0.0193),
    "alpha04_F": (0.1612, 0.1218),
    "alpha04_I": (0.4094, 0.4099),
    "alpha06_A": (0.0119, 0.0089),
    "alpha06_F": (0.1721, 0.1343),
    "alpha06_I": (0.6147, 0.6150),
    "alpha08_A": (0.0129, 0.0168),
    "alpha08_F": (0.1944, 0.1499),
    "alpha08_I": (0.8202, 0.8240),
}

I_RUNS = [
    ("alpha02_I", "thermal_alpha02_I_100_73", 0.20, 100, 73),
    ("alpha04_I", "thermal_alpha04_I_100_73", 0.40, 100, 73),
    ("alpha06_I", "thermal_alpha06_I_100_73", 0.60, 100, 73),
    ("alpha08_I", "thermal_alpha08_I_100_73", 0.80, 100, 73),
]
AF_RUNS = [
    ("alpha02_A", "thermal_alpha02_A_73_73",  0.20,  73, 73),
    ("alpha02_F", "thermal_alpha02_F_73_45",  0.20,  73, 45),
    ("alpha04_A", "thermal_alpha04_A_73_73",  0.40,  73, 73),
    ("alpha04_F", "thermal_alpha04_F_73_45",  0.40,  73, 45),
    ("alpha06_A", "thermal_alpha06_A_73_73",  0.60,  73, 73),
    ("alpha06_F", "thermal_alpha06_F_73_45",  0.60,  73, 45),
    ("alpha08_A", "thermal_alpha08_A_73_73",  0.80,  73, 73),
    ("alpha08_F", "thermal_alpha08_F_73_45",  0.80,  73, 45),
]


def run_set(runs, label_set: str):
    rows = []
    print(f"\n  T_ref = {T_REF_K_NEW} K = {T_REF_K_NEW-273.15:.1f}°C "
          f"(was {T_REF_K_OLD} K = {T_REF_K_OLD-273.15:.1f}°C)")
    print(f"  Hu_J_kg_effective = {HU_RESIDUAL_J_KG} J/kg, k_uc factor = {K_UC_FACTOR}")
    print("  " + "-" * 90)
    print(f"  {'Dataset':<14} {'α_u':>5}  "
          f"{'R1 prior':>9}  {'R1 new':>9}  {'ΔR1':>8}  "
          f"{'R2 prior':>9}  {'R2 new':>9}  {'ΔR2':>8}  verdict")
    print("  " + "-" * 90)

    original_T_ref = te2d.T_REF_K
    try:
        te2d.T_REF_K = T_REF_K_NEW
        for label, folder, alpha_u, T_pl, T_so in runs:
            d = run_one(CW_DATA / folder, T_pl, T_so, factor=K_UC_FACTOR)
            r1, r2, r3, _, _ = residuals_F(d, d["eng_y"], d["eng_x"], d["eng_T_F_168"])
            r1p, r2p = PRIOR[label]
            d_r1 = r1 - r1p
            d_r2 = r2 - r2p
            gate_r1 = r1 <= GATE_F
            gate_r2 = r2 <= GATE_F
            passes  = gate_r1 and gate_r2
            v = "PASS" if passes else "FAIL"
            print(f"  {label:<14} {alpha_u:>5.2f}  "
                  f"{r1p:>9.4f}  {r1:>9.4f}  {d_r1:>+8.4f}  "
                  f"{r2p:>9.4f}  {r2:>9.4f}  {d_r2:>+8.4f}  [{v}]")
            rows.append({
                "label": label, "folder": folder, "alpha_u": alpha_u,
                "T_pl_F": T_pl, "T_so_F": T_so,
                "T_REF_K": T_REF_K_NEW,
                "Hu_J_kg_effective": HU_RESIDUAL_J_KG,
                "factor": K_UC_FACTOR,
                "r1_prior_F": r1p, "r2_prior_F": r2p,
                "r1_max_F":   r1,  "r2_max_F":   r2, "r3_rmse_F": r3,
                "delta_r1_F": d_r1, "delta_r2_F": d_r2,
                "gate_r1": gate_r1, "gate_r2": gate_r2,
                "passes": passes, "minimax": max(r1, r2),
            })
    finally:
        te2d.T_REF_K = original_T_ref
    return rows


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    print("\n§4.10 follow-up — T_ref hypothesis test")
    print("=" * 100)
    print("Phase A: 4 I-scenario datasets")
    print("=" * 100)

    rows_I = run_set(I_RUNS, "I")

    # Decision: did I residuals drop "substantially"?
    # Define substantial = mean ΔR1 ≤ -0.05°F (i.e. avg drop > 0.05°F across the four I)
    mean_dR1_I = sum(r["delta_r1_F"] for r in rows_I) / len(rows_I)
    mean_dR2_I = sum(r["delta_r2_F"] for r in rows_I) / len(rows_I)
    n_pass_new = sum(1 for r in rows_I if r["passes"])
    n_pass_pri = 0  # all 4 I failed in §4.10.3
    print(f"\n  I-scenario mean ΔR1 = {mean_dR1_I:+.4f}°F, "
          f"mean ΔR2 = {mean_dR2_I:+.4f}°F")
    print(f"  I-scenario gate pass count: prior {n_pass_pri}/4 → new {n_pass_new}/4")

    substantial = (mean_dR1_I <= -0.05) or (n_pass_new > 0)
    print(f"  Substantial improvement? {'YES' if substantial else 'NO'} "
          f"(threshold: mean ΔR1 ≤ -0.05°F or new pass count > 0)")

    # Write I CSV
    fields = ["label", "folder", "alpha_u", "T_pl_F", "T_so_F",
              "T_REF_K", "Hu_J_kg_effective", "factor",
              "r1_prior_F", "r2_prior_F",
              "r1_max_F", "r2_max_F", "r3_rmse_F",
              "delta_r1_F", "delta_r2_F",
              "gate_r1", "gate_r2", "passes", "minimax"]
    out_I = DATA / "sprint8_Tref21_residuals_I.csv"
    with open(out_I, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows_I)
    print(f"\nWrote {out_I}")

    rows_AF = []
    if substantial:
        print("\n" + "=" * 100)
        print("Phase B: 8 A and F datasets (regression check)")
        print("=" * 100)
        rows_AF = run_set(AF_RUNS, "AF")

        # Regression check
        regressions = [r for r in rows_AF
                       if (r["delta_r1_F"] >  0.05) or (r["delta_r2_F"] >  0.05)
                       or (not r["passes"] and PRIOR[r["label"]][0] <= GATE_F
                           and PRIOR[r["label"]][1] <= GATE_F)]
        if regressions:
            print(f"\n  ⚠ {len(regressions)} A/F dataset(s) regressed (ΔR > +0.05°F or new gate FAIL):")
            for r in regressions:
                print(f"    {r['label']}: ΔR1={r['delta_r1_F']:+.4f}, ΔR2={r['delta_r2_F']:+.4f}, "
                      f"gate {'PASS' if r['passes'] else 'FAIL'}")
        else:
            print("\n  No A/F regressions (all stayed within ±0.05°F of prior, all still pass gate).")

        out_AF = DATA / "sprint8_Tref21_residuals_AF.csv"
        with open(out_AF, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows_AF)
        print(f"Wrote {out_AF}")

        out_full = DATA / "sprint8_Tref21_residuals_full.csv"
        with open(out_full, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows_I + rows_AF)
        print(f"Wrote {out_full}")
    else:
        print("\nSkipping A/F regression check (I-scenario change not substantial).")

    print("\nDone.")


if __name__ == "__main__":
    main()
