#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.5 refined sweep — 0.005-step around coarse minimum.

Run only if stage2_kuc_sweep.py reports that the gap to the neighbouring
factor is < 0.005°F (i.e. the optimum is not clearly resolved).

Reads findings/kuc_sweep_optimal.txt to determine starting factors, then
sweeps ±0.03 in steps of 0.005 around each coarse optimum.

Output: findings/kuc_sweep_refined.csv  (appends new rows)
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


def load_coarse_optima():
    opt_path = FINDINGS / "kuc_sweep_optimal.txt"
    if not opt_path.exists():
        print("ERROR: findings/kuc_sweep_optimal.txt not found. Run stage2_kuc_sweep.py first.")
        sys.exit(1)
    result = {}
    for line in opt_path.read_text().splitlines():
        if line.startswith("#") or not line.strip() or "anchor" in line or "alpha_u" in line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                au = float(parts[0])
                bf = float(parts[1])
                result[au] = bf
            except ValueError:
                pass
    return result


def refined_sweep(alpha_u, center, runs):
    lo = max(0.85, center - 0.030)
    hi = min(1.10, center + 0.030)
    factors = sorted(set([round(lo + i * 0.005, 3) for i in range(int((hi - lo) / 0.005) + 2)]
                         + [center]))
    # Remove already-run factors from the coarse sweep (0.92..1.04 step 0.02)
    coarse = {round(0.92 + i * 0.02, 3) for i in range(7)}
    factors = [f for f in factors if round(f, 3) not in coarse]

    all_rows = []
    factor_minimax = {}

    for kf in factors:
        print(f"  k_uc={kf:.3f} ...", end="", flush=True)
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

    best = min(factor_minimax, key=factor_minimax.get)
    return all_rows, best, factor_minimax


def main():
    optima = load_coarse_optima()
    print(f"\n§3.5 Refined sweep (0.005-step) around coarse optima")
    print("=" * 60)

    fieldnames = ["alpha_u","factor","label","T_pl_F","T_so_F",
                  "r1_max_F","r2_max_F","r3_rmse_F","minimax",
                  "gate_r1","gate_r2","passes"]

    csv_path = FINDINGS / "kuc_sweep_refined.csv"
    updated_optima = {}

    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()

        for alpha_u in sorted(ALPHA_GROUPS):
            if alpha_u not in optima:
                print(f"  α_u={alpha_u:.2f}: no coarse optimum found, skipping")
                updated_optima[alpha_u] = 0.96
                continue
            center = optima[alpha_u]
            runs   = ALPHA_GROUPS[alpha_u]
            print(f"\n── α_u = {alpha_u:.2f}  (coarse optimum = {center:.3f}) ─────────")
            rows, best, fm = refined_sweep(alpha_u, center, runs)
            w.writerows(rows)

            print(f"\n  Refined sweep summary for α_u={alpha_u:.2f}:")
            for kf in sorted(fm):
                star = " ← BEST" if abs(kf - best) < 0.001 else ""
                gate = "✓" if fm[kf] <= GATE_F else "✗"
                print(f"    factor={kf:.3f}  minimax={fm[kf]:.4f} {gate}{star}")
            print(f"  → Refined optimal: {best:.3f}")
            updated_optima[alpha_u] = best

    print(f"\nWrote: {csv_path}")

    # Update optimal file with refined values
    opt_path = FINDINGS / "kuc_sweep_optimal_refined.txt"
    lines = ["# Sprint 8 Stage 2 §3.5 — Refined optimal k_uc factor per α target\n"]
    lines.append("alpha_u  optimal_factor  source")
    lines.append("0.036    0.960           Sprint 7 anchor")
    for au in sorted(updated_optima):
        lines.append(f"{au:.3f}    {updated_optima[au]:.3f}           refined sweep")
    opt_path.write_text("\n".join(lines))
    print(f"Wrote: {opt_path}")


if __name__ == "__main__":
    main()
