#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.5 — bidirectional k_uc sensitivity sweep.

For each of the 4 new α targets (0.20, 0.40, 0.60, 0.80), sweeps k_uc factor
across {0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.04} and reports R1/R2/R3 at
each factor.  Identifies the optimal factor (minimising minimax of R1 max,
R2 max) for each α target.

Output: findings/kuc_sweep.csv + printed tables
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

ALPHA_GROUPS = {
    0.20: [
        ("alpha02_A", "thermal_alpha02_A_73_73",  0.20,  73,  73),
        ("alpha02_F", "thermal_alpha02_F_73_45",  0.20,  73,  45),
        ("alpha02_I", "thermal_alpha02_I_100_73", 0.20, 100,  73),
    ],
    0.40: [
        ("alpha04_A", "thermal_alpha04_A_73_73",  0.40,  73,  73),
        ("alpha04_F", "thermal_alpha04_F_73_45",  0.40,  73,  45),
        ("alpha04_I", "thermal_alpha04_I_100_73", 0.40, 100,  73),
    ],
    0.60: [
        ("alpha06_A", "thermal_alpha06_A_73_73",  0.60,  73,  73),
        ("alpha06_F", "thermal_alpha06_F_73_45",  0.60,  73,  45),
        ("alpha06_I", "thermal_alpha06_I_100_73", 0.60, 100,  73),
    ],
    0.80: [
        ("alpha08_A", "thermal_alpha08_A_73_73",  0.80,  73,  73),
        ("alpha08_F", "thermal_alpha08_F_73_45",  0.80,  73,  45),
        ("alpha08_I", "thermal_alpha08_I_100_73", 0.80, 100,  73),
    ],
}

K_FACTORS = [0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.04]


def sweep_alpha_group(alpha_u, runs):
    """Run all factor × scenario combinations for one α target.

    Returns list of dicts and the optimal factor (minimising worst-case minimax
    over the 3 scenarios at each factor).
    """
    all_rows = []
    factor_minimax = {}   # factor → worst minimax across the 3 scenarios

    for kf in K_FACTORS:
        print(f"  k_uc={kf:.2f} ...", end="", flush=True)
        worst = 0.0
        for label, folder, au, T_pl, T_so in runs:
            rpt = run_and_compare(CW_DATA / folder, T_pl, T_so, factor=kf)
            all_rows.append({
                "alpha_u": alpha_u, "factor": kf, "label": label,
                "T_pl_F": T_pl, "T_so_F": T_so,
                "r1_max_F": rpt.r1_max_F, "r2_max_F": rpt.r2_max_F,
                "r3_rmse_F": rpt.r3_rmse_F, "minimax": rpt.minimax,
                "gate_r1": rpt.gate_r1, "gate_r2": rpt.gate_r2,
                "passes": rpt.passes,
            })
            worst = max(worst, rpt.minimax)
        factor_minimax[kf] = worst
        print(f"  worst minimax={worst:.4f}")

    # Optimal factor = minimises worst-case minimax
    best_factor = min(factor_minimax, key=factor_minimax.get)

    # Check convexity (warn if monotonic)
    minimaxes = [factor_minimax[kf] for kf in K_FACTORS]
    min_idx   = minimaxes.index(min(minimaxes))
    if min_idx == 0:
        print(f"  *** WARNING: optimum at left edge (factor={K_FACTORS[0]}) — consider extending sweep ***")
    elif min_idx == len(K_FACTORS) - 1:
        print(f"  *** WARNING: optimum at right edge (factor={K_FACTORS[-1]}) — consider extending sweep ***")

    return all_rows, best_factor, factor_minimax


def print_alpha_table(alpha_u, rows, factor_minimax):
    best = min(factor_minimax, key=factor_minimax.get)
    print(f"\n  α_u={alpha_u:.2f}  factor sweep summary (worst-case minimax across 3 scenarios):")
    print(f"  {'Factor':>8}  {'Worst minimax':>14}  {'Sprint7 factor':>14}  {'Best?':>6}")
    for kf in K_FACTORS:
        mm = factor_minimax[kf]
        s7_flag = "← S7" if abs(kf - 0.96) < 0.001 else ""
        best_flag = "← BEST" if abs(kf - best) < 0.001 else ""
        gate_sym  = "✓" if mm <= GATE_F else "✗"
        print(f"  {kf:>8.3f}  {mm:>14.4f} {gate_sym}  {s7_flag:>14}  {best_flag}")
    print(f"  → Optimal factor at α_u={alpha_u:.2f}: {best:.3f}  (minimax={factor_minimax[best]:.4f}°F)")


def main():
    FINDINGS.mkdir(parents=True, exist_ok=True)

    print("\n§3.5 k_uc bidirectional sweep (4 α targets × 7 factors × 3 scenarios)")
    print("=" * 70)

    all_rows = []
    optimal_factors = {}   # alpha_u → best factor

    for alpha_u in sorted(ALPHA_GROUPS):
        runs = ALPHA_GROUPS[alpha_u]
        print(f"\n── α_u = {alpha_u:.2f} ──────────────────────────────────────────")
        rows, best, fm = sweep_alpha_group(alpha_u, runs)
        all_rows.extend(rows)
        optimal_factors[alpha_u] = best
        print_alpha_table(alpha_u, rows, fm)

    # ── Summary table ──────────────────────────────────────────────────────
    print("\n\n§3.5 Optimal factor summary")
    print("-" * 60)
    print(f"  {'α_u':>6}  {'Optimal factor':>14}  {'S7 factor (0.96)':>16}  "
          f"{'passes at best':>14}  {'passes at 0.96':>14}")
    for alpha_u in sorted(optimal_factors):
        bf = optimal_factors[alpha_u]
        rows_au = [r for r in all_rows if abs(r["alpha_u"] - alpha_u) < 0.01]
        rows_best = [r for r in rows_au if abs(r["factor"] - bf) < 0.001]
        rows_s7   = [r for r in rows_au if abs(r["factor"] - 0.96) < 0.001]
        pass_best = sum(r["passes"] for r in rows_best)
        pass_s7   = sum(r["passes"] for r in rows_s7)
        print(f"  {alpha_u:>6.2f}  {bf:>14.3f}  {'0.960':>16}"
              f"  {pass_best}/3{' PASS' if pass_best==3 else ' FAIL':>7}"
              f"  {pass_s7}/3{' PASS' if pass_s7==3 else ' FAIL':>7}")

    # Sprint 7 anchor: α=0.036 → factor=0.96 (known)
    print(f"  {'0.036':>6}  {'0.960':>14}  {'0.960':>16}  {'3/3 PASS':>14}  {'3/3 PASS':>14}  (Sprint 7)")

    # ── Write CSV ──────────────────────────────────────────────────────────
    fieldnames = ["alpha_u","factor","label","T_pl_F","T_so_F",
                  "r1_max_F","r2_max_F","r3_rmse_F","minimax",
                  "gate_r1","gate_r2","passes"]
    csv_path = FINDINGS / "kuc_sweep.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nWrote: {csv_path}")

    # ── Write optimal factor summary ───────────────────────────────────────
    opt_path = FINDINGS / "kuc_sweep_optimal.txt"
    lines = ["# Sprint 8 Stage 2 §3.5 — Optimal k_uc factor per α target\n"]
    lines.append("alpha_u  optimal_factor")
    lines.append("0.036    0.960  (Sprint 7 anchor)")
    for au in sorted(optimal_factors):
        lines.append(f"{au:.3f}    {optimal_factors[au]:.3f}")
    opt_path.write_text("\n".join(lines))
    print(f"Wrote: {opt_path}")

    # ── Hint for refinement ────────────────────────────────────────────────
    print("\nRefinement hint:")
    for au in sorted(optimal_factors):
        bf = optimal_factors[au]
        rows_au = [r for r in all_rows if abs(r["alpha_u"] - au) < 0.01]
        rows_bf = [r for r in rows_au if abs(r["factor"] - bf) < 0.001]
        mm_best = max(r["minimax"] for r in rows_bf)
        # Find adjacent factor values
        idx = K_FACTORS.index(bf)
        left  = K_FACTORS[idx - 1] if idx > 0 else None
        right = K_FACTORS[idx + 1] if idx < len(K_FACTORS) - 1 else None
        if left is not None and right is not None:
            rows_l = [r for r in rows_au if abs(r["factor"] - left) < 0.001]
            rows_r = [r for r in rows_au if abs(r["factor"] - right) < 0.001]
            mm_l   = max(r["minimax"] for r in rows_l)
            mm_r   = max(r["minimax"] for r in rows_r)
            gap    = min(mm_l, mm_r) - mm_best
            if gap < 0.005:
                print(f"  α_u={au:.2f}: gap to neighbour={gap:.4f}°F — run stage2_kuc_sweep_refined.py to pin to ±0.5%")
            else:
                print(f"  α_u={au:.2f}: gap to neighbour={gap:.4f}°F — factor={bf:.3f} is clear minimum, no refinement needed")
        else:
            print(f"  α_u={au:.2f}: optimal at edge — extend sweep range")


if __name__ == "__main__":
    main()
