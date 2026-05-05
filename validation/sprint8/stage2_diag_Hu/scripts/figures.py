#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-Hu §4.4–§4.5 — plots.

Reads data/Hu_residual_table.csv and data/T_core_trajectories.npz.

Writes:
  figures/dT_core_vs_Hu.png        — §4.4 ΔT vs Hu (log-x)
  figures/T_core_trajectories.png  — §4.5 T_core(t) by Hu, two panels

Run from validation/sprint8/stage2_diag_Hu/.
"""
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE    = Path(__file__).resolve().parents[1]
DATA    = HERE / "data"
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)


def load_residual_table() -> list[dict]:
    rows = []
    with open(DATA / "Hu_residual_table.csv", newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "folder":    r["folder"],
                "alpha_u":   float(r["alpha_u"]),
                "Hu_J_kg":   float(r["Hu_J_kg"]),
                "dT_168_C":  float(r["dT_168_C"]),
                "dT_168_F":  float(r["dT_168_F"]),
                "dT_max_C":  float(r["dT_max_C"]),
            })
    return rows


def load_trajectories() -> dict:
    npz = np.load(DATA / "T_core_trajectories.npz")
    folders = sorted({k.replace("_time_hrs", "").replace("_T_core_F", "")
                      for k in npz.files})
    result = {}
    for name in folders:
        result[name] = {
            "time_hrs":  npz[name + "_time_hrs"],
            "T_core_F":  npz[name + "_T_core_F"],
        }
    return result


def plot_dT_vs_Hu(rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))

    for au, color, marker in [(0.20, "#2166ac", "o"), (0.80, "#d6604d", "s")]:
        subset = sorted([r for r in rows if abs(r["alpha_u"] - au) < 0.05],
                        key=lambda r: r["Hu_J_kg"])
        hu_arr  = np.array([r["Hu_J_kg"]  for r in subset])
        dT_arr  = np.array([r["dT_168_C"] for r in subset])
        ax.plot(hu_arr, dT_arr, color=color, marker=marker, linewidth=1.5,
                markersize=7, label=f"α_u = {au:.2f}", zorder=3)
        # annotate Hu=0.01 / Hu=1 if present
        for r in subset:
            if r["Hu_J_kg"] < 2:
                ax.annotate(f"Hu={r['Hu_J_kg']:.2g}",
                            (r["Hu_J_kg"], r["dT_168_C"]),
                            textcoords="offset points", xytext=(6, 4),
                            fontsize=7, color=color)

    ax.set_xscale("log")
    ax.set_xlabel("Hu input (J/kg)", fontsize=11)
    ax.set_ylabel("ΔT_core(t=168 hr)  [°C]", fontsize=11)
    ax.set_title("CW core warming vs input Hu — A scenario (T_pl = T_soil = 73°F)", fontsize=10)
    ax.legend(fontsize=10)
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.set_xlim(0.008, 2e5)

    out = FIGURES / "dT_core_vs_Hu.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


def plot_trajectories(rows: list[dict], trajs: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    cmap = plt.get_cmap("viridis")

    for idx, (au, ax) in enumerate([(0.20, axes[0]), (0.80, axes[1])]):
        subset = sorted([r for r in rows if abs(r["alpha_u"] - au) < 0.05],
                        key=lambda r: r["Hu_J_kg"])
        n = len(subset)
        for i, r in enumerate(subset):
            color = cmap(i / max(n - 1, 1))
            traj  = trajs[r["folder"]]
            t = traj["time_hrs"]
            T = traj["T_core_F"]
            mask = t <= 168.01
            label = f"Hu = {r['Hu_J_kg']:.5g}"
            ax.plot(t[mask], T[mask], color=color, linewidth=1.4, label=label)

        ax.set_xlabel("Time (hr)", fontsize=10)
        ax.set_title(f"α_u = {au:.2f}", fontsize=11)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=7, loc="upper left")
        if idx == 0:
            ax.set_ylabel("T_core (°F)", fontsize=10)

    fig.suptitle("CW core temperature trajectories — A scenario (T_pl = T_soil = 73°F)", fontsize=11)
    out = FIGURES / "T_core_trajectories.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


def main():
    rows  = load_residual_table()
    trajs = load_trajectories()
    plot_dT_vs_Hu(rows)
    plot_trajectories(rows, trajs)
    print("Done.")


if __name__ == "__main__":
    main()
