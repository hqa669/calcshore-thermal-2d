#!/usr/bin/env python3
"""§4.11.8 Task 1 — Fit closed-form c(T_pc) to v3 α_u-averaged c_optimal data.

Loads `data/table_c_optimal_v3.csv`, computes α_u-averaged c at each T_pc, and
fits five candidate functional forms:

  A.  Arrhenius (1 param):  c = exp[B × (1/T_ref_K − 1/T_K)],  T_ref = 70°F → 294.261 K
  B.  Quadratic in T_pc (3 params)
  C.  Cubic in T_pc (4 params)
  D.  Saturating exponential (2 params): c = c_inf − (c_inf − 1) × exp[−k × (T − T_ref)]
  E.  Logistic (3 params): c = L / (1 + exp(−k × (T − T0)))

Reports best-fit parameters, R², max|residual|, and selects the simplest form
that meets the ~0.01 absolute residual threshold. Anchored at T_pc=70°F → c=1.

Outputs:
  data/c_T_pc_fit_candidates.csv   — all candidates with parameters and quality
  data/c_T_pc_fit_chosen.csv       — chosen form with parameters
  figures/c_T_pc_fit_comparison.png — data + all candidate curves
  figures/c_T_pc_fit_chosen.png     — chosen fit and residuals
"""
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

HERE    = Path(__file__).resolve().parents[1]
DATA    = HERE / "data"
FIGURES = HERE / "figures"

T_PCS_F = np.array([40, 50, 60, 70, 80, 90, 100, 110], dtype=float)
ALPHAS  = ["α_u=0.20", "α_u=0.40", "α_u=0.60", "α_u=0.80"]
T_REF_F = 70.0
T_REF_K = (T_REF_F - 32.0) * 5.0 / 9.0 + 273.15  # 294.261 K
THRESHOLD_ABS = 0.01  # selection threshold: max|residual| in c

def F_to_K(T_F):
    return (T_F - 32.0) * 5.0 / 9.0 + 273.15


def load_c_avg():
    rows = []
    with open(DATA / "table_c_optimal_v3.csv") as f:
        for r in csv.DictReader(f):
            rows.append(r)
    c_grid = np.array([[float(r[a]) for a in ALPHAS] for r in rows])
    c_avg  = c_grid.mean(axis=1)
    c_std  = c_grid.std(axis=1)
    return c_grid, c_avg, c_std


# Candidate forms ----------------------------------------------------------------

def fA_arrhenius(T_F, B):
    return np.exp(B * (1.0 / T_REF_K - 1.0 / F_to_K(T_F)))

def fB_quadratic(T_F, a0, a1, a2):
    return a0 + a1 * T_F + a2 * T_F ** 2

def fC_cubic(T_F, a0, a1, a2, a3):
    return a0 + a1 * T_F + a2 * T_F ** 2 + a3 * T_F ** 3

def fD_saturating(T_F, c_inf, k):
    return c_inf - (c_inf - 1.0) * np.exp(-k * (T_F - T_REF_F))

def fE_logistic(T_F, L, k, T0):
    # No exact c=1 anchor; constrained via additional "anchor" residual in fit
    return L / (1.0 + np.exp(-k * (T_F - T0)))


def fit_and_score(name, fn, T, c, p0, n_params, anchor_T=None, anchor_val=1.0,
                  anchor_weight=10.0):
    """Fit fn to (T, c) data; optionally enforce a soft anchor at (anchor_T, anchor_val)."""
    if anchor_T is not None:
        # Concatenate anchor as an extra data point with high weight (sigma small)
        T_fit = np.concatenate([T, [anchor_T]])
        c_fit = np.concatenate([c, [anchor_val]])
        sigma = np.concatenate([np.ones_like(c), [1.0 / anchor_weight]])
        popt, _ = curve_fit(fn, T_fit, c_fit, p0=p0, sigma=sigma, absolute_sigma=False)
    else:
        popt, _ = curve_fit(fn, T, c, p0=p0)

    c_pred = fn(T, *popt)
    res    = c - c_pred
    ss_res = float(np.sum(res ** 2))
    ss_tot = float(np.sum((c - c.mean()) ** 2))
    R2     = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    max_abs = float(np.max(np.abs(res)))
    max_pct = float(np.max(np.abs(res / c)) * 100.0)
    anchor_dev = float(abs(fn(np.array([T_REF_F]), *popt)[0] - 1.0))
    return {
        "name": name,
        "n_params": n_params,
        "params": popt.tolist(),
        "R2": R2,
        "max_abs_res": max_abs,
        "max_pct_res": max_pct,
        "anchor_dev_at_70F": anchor_dev,
        "residuals": res.tolist(),
        "predict": lambda Tarr, p=popt, fn=fn: fn(Tarr, *p),
    }


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    c_grid, c_avg, c_std = load_c_avg()

    print("§4.11.8 Task 1 — Fit c(T_pc) closed-form")
    print(f"  Anchor: T_ref = {T_REF_F}°F → {T_REF_K:.3f} K, c(T_ref) ≡ 1 by construction\n")
    print(f"  α_u-averaged c_optimal:")
    print(f"  {'T_pc °F':>7}  {'c_avg':>8}  {'c_std':>8}  {'std/mean %':>10}")
    for T, ca, cs in zip(T_PCS_F, c_avg, c_std):
        pct = 100.0 * cs / ca if ca > 0 else 0
        print(f"  {T:>7.0f}  {ca:>8.4f}  {cs:>8.4f}  {pct:>10.3f}%")
    max_au_spread = float(c_grid.max(axis=1).max() - c_grid.min(axis=1).min())
    print(f"  Max α_u spread (any T_pc): {(c_grid.max(axis=1) - c_grid.min(axis=1)).max():.4f}")
    print()

    # Fit candidates
    cands = []
    cands.append(fit_and_score("A_arrhenius", fA_arrhenius, T_PCS_F, c_avg,
                                p0=[1500.0], n_params=1))
    cands.append(fit_and_score("B_quadratic", fB_quadratic, T_PCS_F, c_avg,
                                p0=[-1.0, 0.02, 0.0], n_params=3,
                                anchor_T=T_REF_F))
    cands.append(fit_and_score("C_cubic",     fC_cubic,    T_PCS_F, c_avg,
                                p0=[-1.0, 0.02, 0.0, 0.0], n_params=4,
                                anchor_T=T_REF_F))
    cands.append(fit_and_score("D_saturating", fD_saturating, T_PCS_F, c_avg,
                                p0=[1.5, 0.05], n_params=2))
    cands.append(fit_and_score("E_logistic",  fE_logistic,  T_PCS_F, c_avg,
                                p0=[1.5, 0.05, 65.0], n_params=3,
                                anchor_T=T_REF_F))

    # Print + write CSV
    print(f"  {'Candidate':>14}  {'#p':>2}  {'R²':>8}  {'max|res|':>9}  "
          f"{'max %':>7}  {'@70°F':>9}")
    for cd in cands:
        print(f"  {cd['name']:>14}  {cd['n_params']:>2}  {cd['R2']:>8.5f}  "
              f"{cd['max_abs_res']:>9.5f}  {cd['max_pct_res']:>7.2f}  "
              f"{cd['anchor_dev_at_70F']:>9.5f}")

    out = DATA / "c_T_pc_fit_candidates.csv"
    fields = ["name", "n_params", "params_json", "R2",
              "max_abs_res", "max_pct_res", "anchor_dev_at_70F"]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for cd in cands:
            w.writerow({
                "name": cd["name"],
                "n_params": cd["n_params"],
                "params_json": json.dumps(cd["params"]),
                "R2": cd["R2"],
                "max_abs_res": cd["max_abs_res"],
                "max_pct_res": cd["max_pct_res"],
                "anchor_dev_at_70F": cd["anchor_dev_at_70F"],
            })
    print(f"\nWrote {out}")

    # Selection: simplest meeting the threshold
    preferred_order = ["A_arrhenius", "D_saturating", "B_quadratic",
                       "E_logistic", "C_cubic"]
    chosen = None
    for name in preferred_order:
        cd = next(c for c in cands if c["name"] == name)
        if cd["max_abs_res"] <= THRESHOLD_ABS:
            chosen = cd
            break
    if chosen is None:
        chosen = min(cands, key=lambda c: c["max_abs_res"])
        print(f"\nNO candidate meets max|res| ≤ {THRESHOLD_ABS}; "
              f"selecting best-fit: {chosen['name']} (max|res|={chosen['max_abs_res']:.5f})")
    else:
        print(f"\nSelected: {chosen['name']}  (max|res|={chosen['max_abs_res']:.5f}, "
              f"R²={chosen['R2']:.5f}, {chosen['n_params']} params)")

    out_chosen = DATA / "c_T_pc_fit_chosen.csv"
    with open(out_chosen, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["form", chosen["name"]])
        w.writerow(["n_params", chosen["n_params"]])
        w.writerow(["params_json", json.dumps(chosen["params"])])
        w.writerow(["R2", chosen["R2"]])
        w.writerow(["max_abs_res", chosen["max_abs_res"]])
        w.writerow(["max_pct_res", chosen["max_pct_res"]])
        w.writerow(["anchor_dev_at_70F", chosen["anchor_dev_at_70F"]])
        w.writerow(["T_ref_F", T_REF_F])
        w.writerow([])
        w.writerow(["T_pc_F", "c_avg_data", "c_predicted", "residual"])
        for T, ca, res in zip(T_PCS_F, c_avg, chosen["residuals"]):
            w.writerow([T, ca, ca - res, res])
    print(f"Wrote {out_chosen}")

    # Figure 1: all candidates overlaid
    fig, ax = plt.subplots(figsize=(9, 5.5))
    T_dense = np.linspace(40, 110, 200)
    ax.errorbar(T_PCS_F, c_avg, yerr=c_std, fmt="o", color="black",
                markersize=8, capsize=4, label="data (mean ± α_u std)", zorder=5)
    colors = {"A_arrhenius": "#1f77b4", "B_quadratic": "#ff7f0e",
              "C_cubic": "#2ca02c", "D_saturating": "#d62728",
              "E_logistic": "#9467bd"}
    for cd in cands:
        ax.plot(T_dense, cd["predict"](T_dense), color=colors[cd["name"]],
                linewidth=1.6,
                label=f"{cd['name']} (max|res|={cd['max_abs_res']:.3f}, "
                      f"R²={cd['R2']:.3f})")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.axvline(T_REF_F, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("c_optimal (α_u-averaged)", fontsize=11)
    ax.set_title("§4.11.8 Task 1 — c(T_pc) candidate fits", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(FIGURES / "c_T_pc_fit_comparison.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGURES / 'c_T_pc_fit_comparison.png'}")

    # Figure 2: chosen + residuals
    fig, axes = plt.subplots(2, 1, figsize=(8, 7),
                             gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    ax = axes[0]
    ax.errorbar(T_PCS_F, c_avg, yerr=c_std, fmt="o", color="black",
                markersize=8, capsize=4, label="data (mean ± α_u std)")
    ax.plot(T_dense, chosen["predict"](T_dense), color="#d62728",
            linewidth=2.0, label=f"{chosen['name']} fit")
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.axvline(T_REF_F, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
    ax.set_ylabel("c_optimal", fontsize=11)
    title = f"§4.11.8 Task 1 — chosen fit: {chosen['name']}"
    if chosen["n_params"] == 1:
        title += f" (B={chosen['params'][0]:.1f} K)"
    elif chosen["name"] == "D_saturating":
        title += f" (c_inf={chosen['params'][0]:.4f}, k={chosen['params'][1]:.5f})"
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)

    ax2 = axes[1]
    ax2.bar(T_PCS_F, chosen["residuals"], width=4.0, color="#d62728",
            edgecolor="black", linewidth=0.5)
    ax2.axhline(0.0, color="black", linewidth=0.6)
    ax2.axhline(THRESHOLD_ABS, color="black", linestyle=":",
                linewidth=0.6, alpha=0.5)
    ax2.axhline(-THRESHOLD_ABS, color="black", linestyle=":",
                linewidth=0.6, alpha=0.5)
    ax2.set_xlabel("T_pc (°F)", fontsize=11)
    ax2.set_ylabel("residual (c)", fontsize=11)
    ax2.grid(True, linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(FIGURES / "c_T_pc_fit_chosen.png", dpi=150)
    plt.close(fig)
    print(f"Wrote {FIGURES / 'c_T_pc_fit_chosen.png'}")

    print("\nDone.")


if __name__ == "__main__":
    main()
