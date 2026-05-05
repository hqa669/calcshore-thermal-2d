#!/usr/bin/env python3
"""§3.7 — Trend plots from sweep_results.npz.

Figures:
  figures/c_optimal_vs_Tpc.png          — c_optimal vs T_pc (4 lines, one per α_u)
  figures/max_dT_at_c1_vs_Tpc.png       — baseline gap (c=1) vs T_pc
  figures/max_dT_at_optimal_vs_Tpc.png  — residual after best-fit c vs T_pc
  figures/example_trajectories.png      — T_pc=110°F, α_u=0.80: T_CW, T_eng(c=1), T_eng(c_opt) vs t
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE    = Path(__file__).resolve().parents[1]
DATA    = HERE / "data"
FIGURES = HERE / "figures"

T_PCS       = [40, 50, 60, 70, 80, 90, 100, 110]
ALPHAS      = [0.20, 0.40, 0.60, 0.80]
ALPHA_COLORS = {0.20: "#1f77b4", 0.40: "#ff7f0e", 0.60: "#2ca02c", 0.80: "#d62728"}
ALPHA_LABELS = {0.20: "α_u=0.20", 0.40: "α_u=0.40", 0.60: "α_u=0.60", 0.80: "α_u=0.80"}


def get_val(alpha_u_arr, T_pc_F_arr, values_arr, au, T_pc):
    n = len(alpha_u_arr)
    idx = next((i for i in range(n)
                if abs(alpha_u_arr[i] - au) < 0.01 and abs(T_pc_F_arr[i] - T_pc) < 0.5), None)
    return float(values_arr[idx]) if idx is not None else float("nan")


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)

    npz = np.load(DATA / "sweep_results.npz", allow_pickle=True)
    alpha_u_arr   = npz["alpha_u"]
    T_pc_F_arr    = npz["T_pc_F"]
    c_optimal     = npz["c_optimal"]
    max_dT_opt    = npz["max_dT_at_optimal_F"]
    max_dT_c1     = npz["max_dT_at_c1_F"]
    T_eng_opt_F   = npz["T_core_eng_at_optimal_F"]
    T_eng_c1_F    = npz["T_core_eng_at_c1_F"]
    t_cw_hr       = npz["t_cw_hr"]
    Hu_eff        = float(npz["Hu_eff_J_kg"])

    cw_npz = np.load(DATA / "cw_trajectories.npz", allow_pickle=True)
    T_cw_F = cw_npz["T_core_CW_F"]   # (32, N_t), sorted by (alpha_u, T_pc)

    # --- Figure 1: c_optimal vs T_pc ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for au in ALPHAS:
        vals = [get_val(alpha_u_arr, T_pc_F_arr, c_optimal, au, T) for T in T_PCS]
        ax.plot(T_PCS, vals, "o-", color=ALPHA_COLORS[au], label=ALPHA_LABELS[au], linewidth=1.8)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("c_optimal", fontsize=11)
    ax.set_title(f"c_optimal vs T_pc — engine α_u multiplier to match CW\n"
                 f"(Hu_eff = {Hu_eff:.0f} J/kg, anchored at c=1 for T_pc=70°F, α_u=0.20)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "c_optimal_vs_Tpc.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Figure 2: max_dT_at_c1 vs T_pc ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for au in ALPHAS:
        vals = [get_val(alpha_u_arr, T_pc_F_arr, max_dT_c1, au, T) for T in T_PCS]
        ax.plot(T_PCS, vals, "o-", color=ALPHA_COLORS[au], label=ALPHA_LABELS[au], linewidth=1.8)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("max|T_eng − T_CW| at c=1 (°F)", fontsize=11)
    ax.set_title("Baseline engine–CW gap (c=1) vs T_pc\n(max absolute deviation over t ∈ [0, 168] hr)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "max_dT_at_c1_vs_Tpc.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Figure 3: max_dT_at_optimal vs T_pc ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for au in ALPHAS:
        vals = [get_val(alpha_u_arr, T_pc_F_arr, max_dT_opt, au, T) for T in T_PCS]
        ax.plot(T_PCS, vals, "o-", color=ALPHA_COLORS[au], label=ALPHA_LABELS[au], linewidth=1.8)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("max|T_eng − T_CW| at c_optimal (°F)", fontsize=11)
    ax.set_title("Residual after best-fit c vs T_pc\n(max absolute deviation at c_optimal)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "max_dT_at_optimal_vs_Tpc.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Figure 4: example trajectories at T_pc=110°F, α_u=0.80 ---
    ex_T_pc, ex_au = 110.0, 0.80
    ex_idx = next(i for i in range(len(alpha_u_arr))
                  if abs(alpha_u_arr[i] - ex_au) < 0.01 and abs(T_pc_F_arr[i] - ex_T_pc) < 0.5)
    ex_c_opt = float(c_optimal[ex_idx])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_cw_hr, T_cw_F[ex_idx],      color="#2166ac", linestyle="-",  linewidth=2.0,
            label=f"CW")
    ax.plot(t_cw_hr, T_eng_c1_F[ex_idx],  color="#d6604d", linestyle="--", linewidth=1.8,
            label=f"Engine (c=1.00)")
    ax.plot(t_cw_hr, T_eng_opt_F[ex_idx], color="#4dac26", linestyle="-",  linewidth=1.8,
            label=f"Engine (c_opt={ex_c_opt:.3f})")
    ax.set_xlabel("Time (hr)", fontsize=11)
    ax.set_ylabel("T_core (°F)", fontsize=11)
    ax.set_title(f"T_core(t) trajectories — T_pc=110°F, α_u=0.80\n"
                 f"(Worst-case scenario from §4.10; Hu_eff={Hu_eff:.0f} J/kg)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "example_trajectories.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    print("Done.")


if __name__ == "__main__":
    main()
