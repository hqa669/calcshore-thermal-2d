#!/usr/bin/env python3
"""Sprint 8 Stage 2-floor-test §2.5 — k_uc sweep at Hu_residual.

For each α target (0.20, 0.40, 0.60, 0.80), sweeps k_uc multiplier across
{0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.04} with mix.Hu_J_kg_effective fixed
at HU_RESIDUAL_J_KG.  Reports R1/R2/R3 per scenario and the optimal factor
that minimises minimax(R1, R2) across the 3 scenarios per α target.

Computes spread = max(best_factor) − min(best_factor) across the 4 α targets
plus the Sprint 7 anchor (0.96 at α=0.036).

Run from validation/sprint8/stage2_diag_Hu/.
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE   = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from run_engine_Hu_residual import (
    run_one, residuals_F, HU_RESIDUAL_J_KG, GATE_F, RUNS, CW_DATA,
)

DATA = HERE.parent / "data"

K_FACTORS = [0.92, 0.94, 0.96, 0.98, 1.00, 1.02, 1.04]
SPRINT7_ANCHOR = (0.036, 0.96)   # (α_u, factor) — Sprint 7 baseline


def sweep():
    DATA.mkdir(parents=True, exist_ok=True)

    # Group RUNS by α target
    groups = {}
    for label, folder, au, T_pl, T_so in RUNS:
        groups.setdefault(au, []).append((label, folder, au, T_pl, T_so))

    print(f"\n§2.5 k_uc sweep at Hu_residual = {HU_RESIDUAL_J_KG} J/kg")
    print(f"     factors: {K_FACTORS}")
    print("=" * 100)

    all_rows = []
    best_per_alpha = {}

    for au in sorted(groups):
        print(f"\nα_u = {au:.2f}")
        print("-" * 100)
        # factor → worst minimax across the 3 scenarios
        per_factor = {}
        for kf in K_FACTORS:
            print(f"  k_uc={kf:.2f}", end="", flush=True)
            sc_results = []
            worst_minimax = 0.0
            r1_max_alpha  = 0.0
            r2_max_alpha  = 0.0
            for label, folder, au_, T_pl, T_so in groups[au]:
                d = run_one(CW_DATA / folder, T_pl, T_so, factor=kf)
                r1, r2, r3, _, _ = residuals_F(d, d["eng_y"], d["eng_x"], d["eng_T_F_168"])
                mini = max(r1, r2)
                sc_results.append({
                    "alpha_u": au, "factor": kf, "label": label,
                    "T_pl_F": T_pl, "T_so_F": T_so,
                    "r1_max_F": r1, "r2_max_F": r2, "r3_rmse_F": r3,
                    "minimax": mini,
                    "gate_r1": r1 <= GATE_F, "gate_r2": r2 <= GATE_F,
                    "passes":  (r1 <= GATE_F and r2 <= GATE_F),
                })
                worst_minimax = max(worst_minimax, mini)
                r1_max_alpha  = max(r1_max_alpha, r1)
                r2_max_alpha  = max(r2_max_alpha, r2)
            all_rows.extend(sc_results)
            per_factor[kf] = {
                "worst_minimax": worst_minimax,
                "r1_max":        r1_max_alpha,
                "r2_max":        r2_max_alpha,
                "rows":          sc_results,
            }
            n_pass = sum(1 for r in sc_results if r["passes"])
            print(f"  worst_minimax={worst_minimax:.4f}  R1_max={r1_max_alpha:.4f}  "
                  f"R2_max={r2_max_alpha:.4f}  ({n_pass}/{len(sc_results)} pass)")

        # Optimal factor = min worst_minimax
        best_kf = min(per_factor, key=lambda k: per_factor[k]["worst_minimax"])
        b = per_factor[best_kf]
        n_pass_at_best = sum(1 for r in b["rows"] if r["passes"])
        gate_pass = (b["r1_max"] <= GATE_F and b["r2_max"] <= GATE_F)
        edge_pegged = best_kf in (K_FACTORS[0], K_FACTORS[-1])
        best_per_alpha[au] = {
            "factor":        best_kf,
            "worst_minimax": b["worst_minimax"],
            "r1_max":        b["r1_max"],
            "r2_max":        b["r2_max"],
            "n_pass":        n_pass_at_best,
            "gate_pass":     gate_pass,
            "edge_pegged":   edge_pegged,
        }
        print(f"  → best k_uc = {best_kf}  (R1_max={b['r1_max']:.4f}, "
              f"R2_max={b['r2_max']:.4f}, gate {'PASS' if gate_pass else 'FAIL'}"
              f"{'  EDGE-PEGGED' if edge_pegged else ''})")

    # Spread analysis
    print("\n" + "=" * 100)
    print("§2.5 Best-factor summary")
    print("=" * 100)
    print(f"{'α target':>8}  {'best factor':>12}  {'R1 max':>9}  {'R2 max':>9}  "
          f"{'gate':>5}  {'edge?':>5}")
    for au in sorted(best_per_alpha):
        bp = best_per_alpha[au]
        print(f"{au:>8.2f}  {bp['factor']:>12.2f}  {bp['r1_max']:>9.4f}  "
              f"{bp['r2_max']:>9.4f}  {'PASS' if bp['gate_pass'] else 'FAIL':>5}  "
              f"{'YES' if bp['edge_pegged'] else 'no':>5}")

    factors_incl_anchor = [SPRINT7_ANCHOR[1]] + [best_per_alpha[au]["factor"]
                                                  for au in sorted(best_per_alpha)]
    spread = max(factors_incl_anchor) - min(factors_incl_anchor)
    print(f"\nFactors including Sprint 7 anchor (0.96 at α≈0.036): {factors_incl_anchor}")
    print(f"Spread = max − min = {max(factors_incl_anchor):.2f} − {min(factors_incl_anchor):.2f} "
          f"= {spread:.2f}")
    if spread <= 0.01:
        verdict = "single-factor calibration (≤0.01)"
    elif spread <= 0.04:
        verdict = "marginal α-dependence (0.01 < spread ≤ 0.04)"
    else:
        verdict = "clear α-dependent factor (> 0.04) — would need f(α) fit"
    print(f"Spread verdict: {verdict}")

    # Write CSV
    out_csv = DATA / "sprint8_kuc_sweep_Hu_residual.csv"
    fields = ["alpha_u", "factor", "label", "T_pl_F", "T_so_F",
              "r1_max_F", "r2_max_F", "r3_rmse_F", "minimax",
              "gate_r1", "gate_r2", "passes"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(all_rows)
    print(f"\nWrote {out_csv}")

    # Optimal-factor summary CSV
    opt_csv = DATA / "sprint8_kuc_sweep_optimal.csv"
    with open(opt_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["alpha_u", "best_factor",
                                          "r1_max_F", "r2_max_F",
                                          "n_pass_of_3", "gate_pass",
                                          "edge_pegged"])
        w.writeheader()
        for au in sorted(best_per_alpha):
            bp = best_per_alpha[au]
            w.writerow({"alpha_u": au, "best_factor": bp["factor"],
                        "r1_max_F": bp["r1_max"], "r2_max_F": bp["r2_max"],
                        "n_pass_of_3": bp["n_pass"], "gate_pass": bp["gate_pass"],
                        "edge_pegged": bp["edge_pegged"]})
    print(f"Wrote {opt_csv}")
    print("Done.")


if __name__ == "__main__":
    sweep()
