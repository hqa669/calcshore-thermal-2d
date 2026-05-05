#!/usr/bin/env python3
"""§4.11.6 §5 — v2 trend plots using merged §4.11 + extended-range data.

Figures:
  figures/c_optimal_vs_Tpc_v2.png            — merged c_optimal vs T_pc, edge points marked
  figures/max_dT_at_optimal_vs_Tpc_v2.png    — residual after best-fit c (merged)
  figures/example_trajectories_T110_au08_v2.png — CW vs engine(c=1) vs engine(c_opt_ext)
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE    = Path(__file__).resolve().parents[1]
DATA    = HERE / "data"
FIGURES = HERE / "figures"

T_PCS_INTERIOR = [60, 70, 80]
T_PCS_EDGE     = [40, 50, 90, 100, 110]
T_PCS_ALL      = sorted(T_PCS_INTERIOR + T_PCS_EDGE)
ALPHAS         = [0.20, 0.40, 0.60, 0.80]
ALPHA_COLORS   = {0.20: "#1f77b4", 0.40: "#ff7f0e", 0.60: "#2ca02c", 0.80: "#d62728"}
ALPHA_LABELS   = {au: f"α_u={au:.2f}" for au in ALPHAS}


def lookup(au_arr, T_arr, vals, au, T_pc):
    n = len(au_arr)
    idx = next((i for i in range(n)
                if abs(au_arr[i] - au) < 0.01 and abs(T_arr[i] - T_pc) < 0.5), None)
    return float(vals[idx]) if idx is not None else float("nan")


def merged_curve(key_v1, key_ext, au, sweep0, sweep1):
    vals = []
    for T_pc in T_PCS_ALL:
        if T_pc in T_PCS_INTERIOR:
            vals.append(lookup(sweep0["alpha_u"], sweep0["T_pc_F"], sweep0[key_v1], au, T_pc))
        else:
            vals.append(lookup(sweep1["alpha_u"], sweep1["T_pc_F"], sweep1[key_ext], au, T_pc))
    return vals


def edge_flags(au, sweep1):
    return [bool(lookup(sweep1["alpha_u"], sweep1["T_pc_F"], sweep1["edge_flag"], au, T_pc))
            if T_pc in T_PCS_EDGE else False
            for T_pc in T_PCS_ALL]


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    sweep0 = np.load(DATA / "sweep_results.npz",          allow_pickle=True)
    sweep1 = np.load(DATA / "sweep_results_extended.npz", allow_pickle=True)
    Hu_eff = float(sweep0["Hu_eff_J_kg"])

    # --- Figure 1: merged c_optimal vs T_pc ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for au in ALPHAS:
        cvals = merged_curve("c_optimal", "c_optimal_extended", au, sweep0, sweep1)
        eflag = edge_flags(au, sweep1)
        # Solid line
        ax.plot(T_PCS_ALL, cvals, "-", color=ALPHA_COLORS[au],
                label=ALPHA_LABELS[au], linewidth=1.8, zorder=2)
        # Interior points: filled circles
        ix_int = [i for i, T in enumerate(T_PCS_ALL) if T in T_PCS_INTERIOR]
        ax.plot([T_PCS_ALL[i] for i in ix_int], [cvals[i] for i in ix_int],
                "o", color=ALPHA_COLORS[au], markersize=7, zorder=3)
        # Extended points (not edge-clamped): squares
        ix_ext = [i for i, T in enumerate(T_PCS_ALL) if T in T_PCS_EDGE and not eflag[i]]
        ax.plot([T_PCS_ALL[i] for i in ix_ext], [cvals[i] for i in ix_ext],
                "s", color=ALPHA_COLORS[au], markersize=7, zorder=3)
        # Edge-clamped (extended): X markers
        ix_clp = [i for i, T in enumerate(T_PCS_ALL) if T in T_PCS_EDGE and eflag[i]]
        if ix_clp:
            ax.plot([T_PCS_ALL[i] for i in ix_clp], [cvals[i] for i in ix_clp],
                    "x", color=ALPHA_COLORS[au], markersize=10,
                    markeredgewidth=2.5, zorder=4)
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("c_optimal", fontsize=11)
    ax.set_title(f"c_optimal vs T_pc (v2 — extended search range [0.40, 1.80])\n"
                 f"○ interior (§4.11) · □ extended (§4.11.6) · ✕ edge-clamped at v2 bounds\n"
                 f"(Hu_eff = {Hu_eff:.0f} J/kg, anchored at c=1 for T_pc=70°F, α_u=0.20)",
                 fontsize=10)
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "c_optimal_vs_Tpc_v2.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Figure 2: merged max_dT_at_optimal vs T_pc ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for au in ALPHAS:
        vals = merged_curve("max_dT_at_optimal_F", "max_dT_at_optimal_F", au, sweep0, sweep1)
        ax.plot(T_PCS_ALL, vals, "o-", color=ALPHA_COLORS[au],
                label=ALPHA_LABELS[au], linewidth=1.8)
    ax.set_xlabel("T_pc (°F)", fontsize=11)
    ax.set_ylabel("max|T_eng − T_CW| at c_optimal (°F)", fontsize=11)
    ax.set_title("Residual after best-fit c vs T_pc (v2)\n"
                 "(extended range resolves edge clamping at T_pc ∈ {40,50,90,100,110})",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "max_dT_at_optimal_vs_Tpc_v2.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # --- Figure 3: trajectory at T_pc=110°F, α_u=0.80 (v2 c_opt_extended) ---
    ex_T_pc, ex_au = 110.0, 0.80
    # CW
    traj = np.load(DATA / "cw_trajectories.npz", allow_pickle=True)
    cw_au = traj["alpha_u"]; cw_T = traj["T_pc_F"]
    n_cw = len(cw_au)
    cw_idx = next(i for i in range(n_cw)
                  if abs(cw_au[i] - ex_au) < 0.01 and abs(cw_T[i] - ex_T_pc) < 0.5)
    t_cw = traj["t_hrs"]
    T_cw = traj["T_core_CW_F"][cw_idx]

    # extended c_opt (and trajectories)
    ext_au = sweep1["alpha_u"]; ext_T = sweep1["T_pc_F"]
    n_ext = len(ext_au)
    ex_idx = next(i for i in range(n_ext)
                  if abs(ext_au[i] - ex_au) < 0.01 and abs(ext_T[i] - ex_T_pc) < 0.5)
    c_opt_ext = float(sweep1["c_optimal_extended"][ex_idx])
    edge      = bool(sweep1["edge_flag"][ex_idx])
    T_opt_ext = sweep1["T_core_eng_at_optimal_F"][ex_idx]
    T_c1_ext  = sweep1["T_core_eng_at_c1_F"][ex_idx]
    t_grid    = sweep1["t_cw_hr"]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t_cw, T_cw,       color="#2166ac", linestyle="-",  linewidth=2.0, label="CW")
    ax.plot(t_grid, T_c1_ext, color="#d6604d", linestyle="--", linewidth=1.8,
            label="Engine (c=1.00)")
    label_opt = (f"Engine (c_opt_ext={c_opt_ext:.3f}"
                 + (", EDGE-CLAMPED" if edge else "") + ")")
    ax.plot(t_grid, T_opt_ext, color="#4dac26", linestyle="-",  linewidth=1.8,
            label=label_opt)
    ax.set_xlabel("Time (hr)", fontsize=11)
    ax.set_ylabel("T_core (°F)", fontsize=11)
    ax.set_title(f"T_core(t) trajectories — T_pc=110°F, α_u=0.80 (v2)\n"
                 f"Extended-range c_opt (vs §4.11 c_opt=1.198 from clamped sweep)",
                 fontsize=10)
    ax.legend(fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.4)
    out = FIGURES / "example_trajectories_T110_au08_v2.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    print("Done.")


if __name__ == "__main__":
    main()
