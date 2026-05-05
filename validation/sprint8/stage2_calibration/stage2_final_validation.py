#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.8 + §3.9 — final 21-point validation and Sprint 7 regression.

§3.8  Re-runs all 21 data points (9 Sprint 7 + 12 Sprint 8) using the optimal
      calibration found by §3.5/§3.6.

      If Outcome A (constant factor): all 21 runs use that single factor.
      If Outcome B (α-dependent): reads a function defined below; user must
      update CALIBRATION_MODE and f_alpha() before running.

§3.9  Reports whether Sprint 7's 9-run gate still passes with the new
      calibration (regression check).

Usage:
  python stage2_final_validation.py [--factor FLOAT]
    --factor  override: apply this constant factor to all 21 runs
              (use after determining Outcome A optimal factor)

Output:
  findings/final_validation.csv
  findings/sprint7_regression.csv
"""

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = (HERE / "../../..").resolve()
SC   = (HERE / "../../soil_calibration").resolve()

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from engine_runner import run_and_compare, GATE_F

CW_DATA  = HERE / "cw_data"
FINDINGS = HERE / "findings"

S7_RUNS  = SC / "cw_runs"

SPRINT7_RUNS = [
    ("S7_A", S7_RUNS / "runA_baseline",  0.036,  73,  73),
    ("S7_B", S7_RUNS / "runB_73_60",     0.036,  73,  60),
    ("S7_C", S7_RUNS / "runC_73_90",     0.036,  73,  90),
    ("S7_D", S7_RUNS / "runD_60_73",     0.036,  60,  73),
    ("S7_E", S7_RUNS / "runE_90_73",     0.036,  90,  73),
    ("S7_F", S7_RUNS / "runF_73_45",     0.036,  73,  45),
    ("S7_G", S7_RUNS / "runG_73_100",    0.036,  73, 100),
    ("S7_H", S7_RUNS / "runH_45_73",     0.036,  45,  73),
    ("S7_I", S7_RUNS / "runI_100_73",    0.036, 100,  73),
]

SPRINT8_RUNS = [
    ("alpha02_A", CW_DATA / "thermal_alpha02_A_73_73",  0.20,  73,  73),
    ("alpha02_F", CW_DATA / "thermal_alpha02_F_73_45",  0.20,  73,  45),
    ("alpha02_I", CW_DATA / "thermal_alpha02_I_100_73", 0.20, 100,  73),
    ("alpha04_A", CW_DATA / "thermal_alpha04_A_73_73",  0.40,  73,  73),
    ("alpha04_F", CW_DATA / "thermal_alpha04_F_73_45",  0.40,  73,  45),
    ("alpha04_I", CW_DATA / "thermal_alpha04_I_100_73", 0.40, 100,  73),
    ("alpha06_A", CW_DATA / "thermal_alpha06_A_73_73",  0.60,  73,  73),
    ("alpha06_F", CW_DATA / "thermal_alpha06_F_73_45",  0.60,  73,  45),
    ("alpha06_I", CW_DATA / "thermal_alpha06_I_100_73", 0.60, 100,  73),
    ("alpha08_A", CW_DATA / "thermal_alpha08_A_73_73",  0.80,  73,  73),
    ("alpha08_F", CW_DATA / "thermal_alpha08_F_73_45",  0.80,  73,  45),
    ("alpha08_I", CW_DATA / "thermal_alpha08_I_100_73", 0.80, 100,  73),
]


def load_optimal_factor():
    """Read single-factor optimal from calibration_decision.md or kuc_sweep_optimal.txt."""
    for fname in ("kuc_sweep_optimal_refined.txt", "kuc_sweep_optimal.txt"):
        p = FINDINGS / fname
        if p.exists():
            factors = []
            for line in p.read_text().splitlines():
                if line.startswith("#") or "anchor" in line or "alpha_u" in line or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        factors.append(float(parts[1]))
                    except ValueError:
                        pass
            if factors:
                mean_f = float(np.mean(factors))
                # Also include Sprint 7 anchor
                mean_f = float(np.mean(factors + [0.96]))
                return mean_f
    return 0.96   # fallback


def f_alpha(alpha_u: float, constant_factor: float) -> float:
    """Return k_uc factor for a given α_u.

    For Outcome A: returns constant_factor for all α.
    For Outcome B: replace this with the fitted function after stage2_calibration_decision.py.
    """
    return constant_factor


def run_batch(runs, constant_factor, label_prefix=""):
    rows = []
    all_pass = True
    for label, folder, alpha_u, T_pl, T_so in runs:
        factor = f_alpha(alpha_u, constant_factor)
        print(f"  {label:<14} α_u={alpha_u:.3f}  factor={factor:.4f} ...", end="", flush=True)
        rpt = run_and_compare(folder, T_pl, T_so, factor=factor)
        verdict = "PASS" if rpt.passes else "FAIL"
        if not rpt.passes:
            all_pass = False
        print(f"\r  {label:<14} α_u={alpha_u:.3f}  factor={factor:.4f}"
              f"  R1={rpt.r1_max_F:.4f}  R2={rpt.r2_max_F:.4f}"
              f"  R3={rpt.r3_rmse_F:.4f}  [{verdict}]")
        rows.append({
            "label": label, "folder": str(folder), "alpha_u": alpha_u,
            "T_pl_F": T_pl, "T_so_F": T_so, "factor": factor,
            "r1_max_F": rpt.r1_max_F, "r2_max_F": rpt.r2_max_F,
            "r3_rmse_F": rpt.r3_rmse_F,
            "gate_r1": rpt.gate_r1, "gate_r2": rpt.gate_r2,
            "passes": rpt.passes, "minimax": rpt.minimax,
        })
    return rows, all_pass


def print_section_header(title):
    print(f"\n{'='*70}")
    print(f" {title}")
    print(f"{'='*70}")
    print(f"  {'Label':<14} {'α_u':>5}  {'factor':>7}"
          f"  {'R1 max|R|':>10}  {'R2 max|R|':>10}  {'R3 RMSE':>8}  {'Verdict':>8}")
    print("-" * 70)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--factor", type=float, default=None,
                    help="override: apply this constant factor to all 21 runs")
    args = ap.parse_args()

    FINDINGS.mkdir(parents=True, exist_ok=True)

    if args.factor is not None:
        constant_factor = args.factor
        print(f"\nUsing user-supplied factor: {constant_factor:.4f}")
    else:
        constant_factor = load_optimal_factor()
        print(f"\nLoaded mean optimal factor from findings: {constant_factor:.4f}")

    fieldnames = ["label","folder","alpha_u","T_pl_F","T_so_F","factor",
                  "r1_max_F","r2_max_F","r3_rmse_F","gate_r1","gate_r2",
                  "passes","minimax"]

    # ── §3.9 Sprint 7 regression ──────────────────────────────────────────
    print_section_header("§3.9 Sprint 7 regression check (9 original datasets)")
    s7_rows, s7_all_pass = run_batch(SPRINT7_RUNS, constant_factor)

    csv_path_s7 = FINDINGS / "sprint7_regression.csv"
    with open(csv_path_s7, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(s7_rows)
    print(f"\nSprint 7 regression: {'ALL 9 PASS ✓' if s7_all_pass else 'FAILURES — see above ✗'}")

    # ── §3.8 Sprint 8 validation ──────────────────────────────────────────
    print_section_header("§3.8 Sprint 8 Stage 2 validation (12 new datasets)")
    s8_rows, s8_all_pass = run_batch(SPRINT8_RUNS, constant_factor)

    csv_path_s8 = FINDINGS / "final_validation.csv"
    with open(csv_path_s8, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(s8_rows)
    print(f"\nSprint 8 Stage 2 validation: {'ALL 12 PASS ✓' if s8_all_pass else 'FAILURES — see above ✗'}")

    # ── Combined 21-point summary ─────────────────────────────────────────
    all_rows  = s7_rows + s8_rows
    all_pass  = s7_all_pass and s8_all_pass
    n_pass    = sum(r["passes"] for r in all_rows)
    print(f"\n{'='*70}")
    print(f" Combined 21-point Structure C gate (R1 ≤ {GATE_F}°F, R2 ≤ {GATE_F}°F)")
    print(f"  Passing: {n_pass}/21  ({'ALL PASS ✓' if all_pass else 'FAILURES PRESENT ✗'})")

    if not all_pass:
        print("\n  Failing runs:")
        for r in all_rows:
            if not r["passes"]:
                print(f"    {r['label']:<14}  R1={r['r1_max_F']:.4f}  R2={r['r2_max_F']:.4f}"
                      f"  R1_gate={'✓' if r['gate_r1'] else '✗'}  R2_gate={'✓' if r['gate_r2'] else '✗'}")

    print(f"\nWrote: {csv_path_s7}")
    print(f"Wrote: {csv_path_s8}")


if __name__ == "__main__":
    main()
