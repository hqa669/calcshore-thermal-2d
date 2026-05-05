#!/usr/bin/env python3
"""§4.11.7 §4 — v3 c_optimal table + per-α_u linear fit + classification.

Combines:
  - T_pc=40°F: from Task 1 (`sweep_results_t1_t40_extended.npz`)
  - T_pc ∈ {50,60,70,80,90,100,110}°F: from Task 2 (`sweep_results_t2_refined.npz`)

Outputs:
  data/table_c_optimal_v3.csv             — final 8×4 c_optimal table
  data/linear_fits_per_alpha.csv          — slope, intercept, R², residuals per α_u
  figures/c_optimal_vs_Tpc_v3.png         — c_opt curve with linear fit overlays
  figures/linearity_residuals_v3.png      — c_opt - c_linear_fit per α_u
  figures/L2_landscape_examples.png       — L2 vs c at two representative points
"""
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE    = Path(__file__).resolve().parents[1]
DATA    = HERE / "data"
FIGURES = HERE / "figures"

T_PCS_ALL = [40, 50, 60, 70, 80, 90, 100, 110]
ALPHAS    = [0.20, 0.40, 0.60, 0.80]
ALPHA_COLORS = {0.20: "#1f77b4", 0.40: "#ff7f0e", 0.60: "#2ca02c", 0.80: "#d62728"}


def lookup(au_arr, T_arr, vals, au, T_pc):
    n = len(au_arr)
    idx = next(i for i in range(n)
               if abs(au_arr[i] - au) < 0.01 and abs(T_arr[i] - T_pc) < 0.5)
    return float(vals[idx])


def merged_c(au, T_pc, t1, t2):
    if T_pc == 40:
        return lookup(t1["alpha_u"], t1["T_pc_F"], t1["c_optimal_t1"], au, T_pc)
    return lookup(t2["alpha_u"], t2["T_pc_F"], t2["c_optimal_t2"], au, T_pc)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    t1 = np.load(DATA / "sweep_results_t1_t40_extended.npz", allow_pickle=True)
    t2 = np.load(DATA / "sweep_results_t2_refined.npz",      allow_pickle=True)

    # Build 8×4 v3 table
    table = []
    for T_pc in T_PCS_ALL:
        row = {"T_pc_F": T_pc}
        for au in ALPHAS:
            row[f"α_u={au:.2f}"] = merged_c(au, T_pc, t1, t2)
        table.append(row)

    fields = ["T_pc_F"] + [f"α_u={au:.2f}" for au in ALPHAS]
    out = DATA / "table_c_optimal_v3.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(table)
    print(f"Wrote {out}")

    # Print v3 table
    print("\nv3 c_optimal table (8×4, all interior optima):")
    print(f"  {'T_pc':>5}" + "".join(f"  {'α='+str(au):>9}" for au in ALPHAS))
    for row in table:
        cells = [f"  {row[f'α_u={au:.2f}']:>9.4f}" for au in ALPHAS]
        print(f"  {int(row['T_pc_F']):>5}" + "".join(cells))

    # Per-α_u linear fits
    print("\nPer-α_u linear fits c(T_pc) = a + b·T_pc :")
    print(f"  {'α_u':>5}  {'a (intercept)':>13}  {'b (slope/°F)':>13}  {'R²':>7}  {'max|res|':>9}")
    fit_rows = []
    fits_by_au = {}
    for au in ALPHAS:
        T = np.array(T_PCS_ALL, dtype=float)
        c = np.array([merged_c(au, t, t1, t2) for t in T_PCS_ALL])
        b, a = np.polyfit(T, c, 1)            # slope, intercept
        c_pred = a + b * T
        res = c - c_pred
        ss_res = float(np.sum(res ** 2))
        ss_tot = float(np.sum((c - c.mean()) ** 2))
        R2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
        max_abs_res = float(np.max(np.abs(res)))
        print(f"  {au:>5.2f}  {a:>13.5f}  {b:>13.5f}  {R2:>7.5f}  {max_abs_res:>9.5f}")

        fits_by_au[au] = (a, b, R2, c, res)
        for i, T_pc in enumerate(T_PCS_ALL):
            fit_rows.append({
                "alpha_u":    au,
                "T_pc_F":     T_pc,
                "c_actual":   c[i],
                "c_linear":   c_pred[i],
                "residual":   res[i],
                "residual_pct_of_c": 100.0 * res[i] / c[i] if c[i] != 0 else 0.0,
                "intercept_a": a,
                "slope_b":    b,
                "R2":         R2,
            })

    out_fits = DATA / "linear_fits_per_alpha.csv"
    with open(out_fits, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fit_rows[0].keys()))
        w.writeheader()
        w.writerows(fit_rows)
    print(f"\nWrote {out_fits}")

    # Classification
    print("\nLinearity classification per α_u row:")
    for au, (a, b, R2, c, res) in fits_by_au.items():
        max_abs = float(np.max(np.abs(res)))
        if max_abs < 0.01:
            cls = "LINEAR (within ±0.01)"
        elif max_abs < 0.03:
            cls = "MILDLY NONLINEAR (max|res| in [0.01, 0.03])"
        else:
            cls = "STRONGLY NONLINEAR (max|res| > 0.03)"
        print(f"  α_u={au:.2f}: {cls}  (max|res|={max_abs:.4f}, R²={R2:.4f})")

    # Cross-α_u slope consistency
    slopes = np.array([fits_by_au[au][1] for au in ALPHAS])
    intercepts = np.array([fits_by_au[au][0] for au in ALPHAS])
    print("\nCross-α_u slope consistency:")
    print(f"  slopes:     {slopes}")
    print(f"  mean:       {slopes.mean():.5f}")
    print(f"  std/mean:   {slopes.std() / abs(slopes.mean()) * 100.0:.2f}%")
    print(f"  intercepts: {intercepts}")
    print(f"  intercept std/mean: {intercepts.std() / abs(intercepts.mean()) * 100.0:.2f}%")

    # --- Figure 1: c_optimal vs T_pc with linear fit overlays ---
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for au in ALPHAS:
        a, b, R2, c, _ = fits_by_au[au]
        T = np.array(T_PCS_ALL, dtype=float)
        ax.plot(T, c, "o-", color=ALPHA_COLORS[au], linewidth=1.6,
                label=f"α_u={au:.2f}", markersize=7)
        T_dense = np.linspace(T.min(), T.max(), 100)
        ax.plot(T_dense, a + b * T_dense, "--", color=ALPHA_COLORS[au],
                linewidth=1.0, alpha=0.5)
    ax.axhline(1.0, color="black", linestyle=":", linewidth=0.8, alpha=0.5)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("c_optimal", fontsize=11)
    ax.set_title("c_optimal vs T_pc (v3 — all interior optima, refined to ±0.005)\n"
                 "Solid = data; dashed = per-α_u linear fit",
                 fontsize=10)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "c_optimal_vs_Tpc_v3.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nWrote {out}")

    # --- Figure 2: linearity residuals ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for au in ALPHAS:
        _, _, _, _, res = fits_by_au[au]
        ax.plot(T_PCS_ALL, res, "o-", color=ALPHA_COLORS[au],
                linewidth=1.6, label=f"α_u={au:.2f}", markersize=6)
    ax.axhline(0.0,   color="black", linewidth=0.8)
    ax.axhline(0.01,  color="black", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.axhline(-0.01, color="black", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.axhline(0.03,  color="red",   linestyle=":", linewidth=0.6, alpha=0.5)
    ax.axhline(-0.03, color="red",   linestyle=":", linewidth=0.6, alpha=0.5)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("c_optimal − c_linear_fit", fontsize=11)
    ax.set_title("Linearity residuals per α_u (v3)\n"
                 "Black dotted = ±0.01 (linear threshold); red dotted = ±0.03 (mildly nonlinear)",
                 fontsize=10)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "linearity_residuals_v3.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Figure 3: L2 landscape examples ---
    # Re-derive the L2 vs c scan from Task 2 npz isn't possible (we only stored c_opt_t2),
    # so plot L2_min, L2_at_v2, L2_max from the diagnostic CSV at two example points.
    # Show simple bar comparison: L2 at v2 anchor vs at refined optimum.
    diag_csv = DATA / "resolution_floor_diagnostic.csv"
    if diag_csv.exists():
        with open(diag_csv) as f:
            rows = list(csv.DictReader(f))
        # Pick two example points for a "before vs after" L2 view
        examples = [(0.20, 70), (0.80, 110)]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
        for ax, (au_x, T_x) in zip(axes, examples):
            row = next(r for r in rows
                       if abs(float(r["alpha_u"]) - au_x) < 0.01
                       and abs(float(r["T_pc_F"]) - T_x) < 0.5)
            xs    = ["L2 at c_v2", "L2 at c_opt_t2", "L2 max in window"]
            vals  = [float(row["L2_at_v2"]), float(row["L2_min"]), float(row["L2_max"])]
            ax.bar(xs, vals, color=["#888", "#2ca02c", "#d62728"])
            for x, v in zip(xs, vals):
                ax.text(x, v, f"{v:.4f}", ha="center", va="bottom", fontsize=9)
            ax.set_ylabel("L2 residual (°F RMSE)", fontsize=10)
            ax.set_title(f"T_pc={T_x}°F, α_u={au_x:.2f}\n"
                         f"c_v2={float(row['c_v2']):.4f} → c_opt_t2={float(row['c_optimal_t2']):.4f}",
                         fontsize=10)
            ax.grid(True, axis="y", linestyle="--", alpha=0.4)
        fig.suptitle("L2 landscape examples — refinement effect on residual",
                     fontsize=11)
        fig.tight_layout()
        out = FIGURES / "L2_landscape_examples.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"Wrote {out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
