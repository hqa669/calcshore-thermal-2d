#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.4 — engine baseline residuals at Sprint 7 calibration.

Runs all 12 new datasets with the engine's current k_uc factor (0.96) and
reports R1/R2/R3 residuals.  This answers: does Sprint 7's calibration hold
at the new α targets without modification?

Output: findings/baseline_residuals.csv + printed table
"""

import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
from engine_runner import run_and_compare, GATE_F

CW_DATA  = HERE / "cw_data"
FINDINGS = HERE / "findings"

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

SPRINT7_FACTOR = 0.96


def main():
    FINDINGS.mkdir(parents=True, exist_ok=True)

    print("\n§3.4 Baseline residuals at Sprint 7 k_uc factor = 0.96")
    print("=" * 75)
    print(f"{'Dataset':<14} {'α_u':>5} {'T_pl':>5} {'T_so':>5}"
          f"  {'R1 max|R|':>10}  {'R2 max|R|':>10}  {'R3 RMSE':>8}"
          f"  {'R1':>4}  {'R2':>4}  {'Overall'}")
    print("-" * 75)

    rows = []
    all_pass = True
    for label, folder, alpha_u, T_pl, T_so in RUNS:
        dpath = CW_DATA / folder
        print(f"  {label:<14} ...", end="", flush=True)
        rpt = run_and_compare(dpath, T_pl, T_so, factor=SPRINT7_FACTOR)
        verdict = "PASS" if rpt.passes else "FAIL"
        if not rpt.passes:
            all_pass = False
        print(f"\r  {label:<14} {alpha_u:>5.2f} {T_pl:>5.0f} {T_so:>5.0f}"
              f"  {rpt.r1_max_F:>10.4f}  {rpt.r2_max_F:>10.4f}  {rpt.r3_rmse_F:>8.4f}"
              f"  {'✓' if rpt.gate_r1 else '✗':>4}  {'✓' if rpt.gate_r2 else '✗':>4}"
              f"  [{verdict}]")
        rows.append({
            "label": label, "folder": folder, "alpha_u": alpha_u,
            "T_pl_F": T_pl, "T_so_F": T_so,
            "factor": SPRINT7_FACTOR,
            "r1_max_F": rpt.r1_max_F, "r2_max_F": rpt.r2_max_F,
            "r3_rmse_F": rpt.r3_rmse_F,
            "gate_r1": rpt.gate_r1, "gate_r2": rpt.gate_r2,
            "passes": rpt.passes, "minimax": rpt.minimax,
        })

    print("-" * 75)
    print(f"\nSprint 7 factor (0.96) at new α targets: "
          f"{'ALL 12 PASS' if all_pass else 'FAILURES PRESENT'}")

    # Grouped summary by α_u
    for au in [0.20, 0.40, 0.60, 0.80]:
        group = [r for r in rows if abs(r["alpha_u"] - au) < 0.01]
        r1s = [r["r1_max_F"] for r in group]
        r2s = [r["r2_max_F"] for r in group]
        print(f"  α_u={au:.2f}:  R1 max={max(r1s):.4f}  R2 max={max(r2s):.4f}"
              f"  {'PASS' if all(r['passes'] for r in group) else 'FAIL'}")

    csv_path = FINDINGS / "baseline_residuals.csv"
    fieldnames = ["label","folder","alpha_u","T_pl_F","T_so_F","factor",
                  "r1_max_F","r2_max_F","r3_rmse_F","gate_r1","gate_r2",
                  "passes","minimax"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote: {csv_path}")


if __name__ == "__main__":
    main()
