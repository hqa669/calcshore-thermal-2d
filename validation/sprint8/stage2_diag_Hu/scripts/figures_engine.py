#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-Hu §4.9.4–§4.9.5 — engine comparison plots.

Reads:
  data/Hu_residual_table.csv
  data/Hu_engine_table.csv
  data/T_core_trajectories.npz
  data/T_core_trajectories_engine.npz

Writes:
  figures/dT_core_vs_Hu_with_engine.png       — §4.9.4 CW vs engine, both α_u
  figures/T_core_trajectories_engine.png      — §4.9.5 engine-only trajectories
  figures/T_core_trajectories_4panel.png      — §4.9.5 4-panel CW vs engine

Run from validation/sprint8/stage2_diag_Hu/.
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


def load_csv(path: Path, hu_key: str = "Hu_J_kg") -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "alpha_u":   float(r["alpha_u"]),
                "Hu_J_kg":   float(r[hu_key]),
                "dT_168_C":  float(r["dT_168_C"]),
                "dT_168_F":  float(r["dT_168_F"]),
            })
    return rows


def load_traj(npz_path: Path) -> dict:
    npz = np.load(npz_path)
    folders = sorted({k.replace("_time_hrs", "").replace("_T_core_F", "")
                      for k in npz.files})
    out = {}
    for name in folders:
        out[name] = {
            "time_hrs": npz[name + "_time_hrs"],
            "T_core_F": npz[name + "_T_core_F"],
        }
    return out


def plot_dT_with_engine(cw_rows, eng_rows):
    fig, ax = plt.subplots(figsize=(8, 5.5))

    series = [
        (0.20, cw_rows,  "CW α_u=0.20",     "#2166ac", "o", "-"),
        (0.80, cw_rows,  "CW α_u=0.80",     "#d6604d", "s", "-"),
        (0.20, eng_rows, "engine α_u=0.20", "#2166ac", "o", "--"),
        (0.80, eng_rows, "engine α_u=0.80", "#d6604d", "s", "--"),
    ]

    for au, rows, label, color, marker, linestyle in series:
        sub = sorted([r for r in rows if abs(r["alpha_u"] - au) < 0.05],
                     key=lambda r: r["Hu_J_kg"])
        hu = np.array([r["Hu_J_kg"]  for r in sub])
        dT = np.array([r["dT_168_C"] for r in sub])
        # Replace nonpositive dT with a small floor for log-x display only
        ax.plot(hu, dT, color=color, marker=marker, linewidth=1.5,
                markersize=6, linestyle=linestyle, label=label,
                markerfacecolor=color if linestyle == "-" else "white",
                markeredgecolor=color)

    ax.set_xscale("log")
    ax.set_xlabel("Hu input (J/kg)", fontsize=11)
    ax.set_ylabel("ΔT_core(t=168 hr)  [°C]", fontsize=11)
    ax.set_title("CW vs engine core warming — A scenario (T_pl = T_soil = 73°F, k_uc=0.96)",
                 fontsize=10)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.set_xlim(0.008, 2e5)
    ax.axhline(0, color="black", linewidth=0.5, alpha=0.5)

    out = FIGURES / "dT_core_vs_Hu_with_engine.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


def plot_engine_trajectories(eng_rows, eng_traj):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    cmap = plt.get_cmap("viridis")

    # eng_rows holds (alpha_u, Hu) — labels use the engine key directly
    for idx, (au, ax) in enumerate([(0.20, axes[0]), (0.80, axes[1])]):
        sub = sorted([r for r in eng_rows if abs(r["alpha_u"] - au) < 0.05],
                     key=lambda r: r["Hu_J_kg"])
        n = len(sub)
        for i, r in enumerate(sub):
            color = cmap(i / max(n - 1, 1))
            # Engine trajectory keys are like "au02_Hu1000" not the full folder name
            key = f"au{int(au*10):02d}_Hu{int(r['Hu_J_kg']) if r['Hu_J_kg'] >= 1 else '001'}"
            traj = eng_traj.get(key)
            if traj is None:
                continue
            t = traj["time_hrs"]
            T = traj["T_core_F"]
            mask = t <= 168.01
            ax.plot(t[mask], T[mask], color=color, linewidth=1.4,
                    label=f"Hu = {r['Hu_J_kg']:.5g}")
        ax.set_xlabel("Time (hr)", fontsize=10)
        ax.set_title(f"engine α_u = {au:.2f}", fontsize=11)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=7, loc="upper left")
        if idx == 0:
            ax.set_ylabel("T_core (°F)", fontsize=10)

    fig.suptitle("Engine core temperature trajectories — A scenario", fontsize=11)
    out = FIGURES / "T_core_trajectories_engine.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


def plot_4panel(cw_rows, eng_rows, cw_traj, eng_traj):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    cmap = plt.get_cmap("viridis")

    # Top row: CW (map folder names → use full names like 'thermal_Hu_res_au02_Hu1000')
    cw_folder_for = {
        (0.20, 0.01):    "thermal_Hu_res_au02_Hu001",
        (0.20, 1.0):     "thermal_Hu_res_au02_Hu1",
        (0.20, 10.0):    "thermal_Hu_res_au02_Hu10",
        (0.20, 100.0):   "thermal_Hu_res_au02_Hu100",
        (0.20, 1000.0):  "thermal_Hu_res_au02_Hu1000",
        (0.20, 10000.0): "thermal_Hu_res_au02_Hu10000",
        (0.20, 50000.0): "thermal_Hu_res_au02_Hu50000",
        (0.20, 100000.0):"thermal_Hu_res_au02_Hu100000",
        (0.80, 1.0):     "thermal_Hu_res_au08_Hu1",
        (0.80, 10.0):    "thermal_Hu_res_au08_Hu10",
        (0.80, 100.0):   "thermal_Hu_res_au08_Hu100",
        (0.80, 1000.0):  "thermal_Hu_res_au08_Hu1000",
        (0.80, 10000.0): "thermal_Hu_res_au08_Hu10000",
        (0.80, 50000.0): "thermal_Hu_res_au08_Hu50000",
        (0.80, 100000.0):"thermal_Hu_res_au08_Hu100000",
    }

    eng_label_for = {
        (0.20, 0.01):    "au02_Hu001",
        (0.20, 1.0):     "au02_Hu1",
        (0.20, 10.0):    "au02_Hu10",
        (0.20, 100.0):   "au02_Hu100",
        (0.20, 1000.0):  "au02_Hu1000",
        (0.20, 10000.0): "au02_Hu10000",
        (0.20, 50000.0): "au02_Hu50000",
        (0.20, 100000.0):"au02_Hu100000",
        (0.80, 1.0):     "au08_Hu1",
        (0.80, 10.0):    "au08_Hu10",
        (0.80, 100.0):   "au08_Hu100",
        (0.80, 1000.0):  "au08_Hu1000",
        (0.80, 10000.0): "au08_Hu10000",
        (0.80, 50000.0): "au08_Hu50000",
        (0.80, 100000.0):"au08_Hu100000",
    }

    panels = [
        (axes[0, 0], 0.20, "CW α_u=0.20",     cw_rows,  cw_traj,  cw_folder_for),
        (axes[0, 1], 0.80, "CW α_u=0.80",     cw_rows,  cw_traj,  cw_folder_for),
        (axes[1, 0], 0.20, "engine α_u=0.20", eng_rows, eng_traj, eng_label_for),
        (axes[1, 1], 0.80, "engine α_u=0.80", eng_rows, eng_traj, eng_label_for),
    ]

    for ax, au, title, rows, traj, key_for in panels:
        sub = sorted([r for r in rows if abs(r["alpha_u"] - au) < 0.05],
                     key=lambda r: r["Hu_J_kg"])
        n = len(sub)
        for i, r in enumerate(sub):
            color = cmap(i / max(n - 1, 1))
            key = key_for[(au, r["Hu_J_kg"])]
            d = traj.get(key)
            if d is None:
                continue
            t = d["time_hrs"]
            T = d["T_core_F"]
            mask = t <= 168.01
            ax.plot(t[mask], T[mask], color=color, linewidth=1.3,
                    label=f"Hu={r['Hu_J_kg']:.5g}")
        ax.set_title(title, fontsize=11)
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.legend(fontsize=6, loc="upper left", ncol=1)

    for ax in axes[1, :]:
        ax.set_xlabel("Time (hr)", fontsize=10)
    for ax in axes[:, 0]:
        ax.set_ylabel("T_core (°F)", fontsize=10)

    # Share y across the top row, and share y across the bottom row, but allow rows to differ
    # Actually share y-axis between matched columns for easier comparison
    for col in range(2):
        ymin = min(axes[0, col].get_ylim()[0], axes[1, col].get_ylim()[0])
        ymax = max(axes[0, col].get_ylim()[1], axes[1, col].get_ylim()[1])
        for row in range(2):
            axes[row, col].set_ylim(ymin, ymax)

    fig.suptitle("CW vs engine core trajectories — A scenario, T_pl=T_soil=73°F", fontsize=11)
    out = FIGURES / "T_core_trajectories_4panel.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")


def main():
    cw_rows  = load_csv(DATA / "Hu_residual_table.csv")
    eng_rows = load_csv(DATA / "Hu_engine_table.csv")
    cw_traj  = load_traj(DATA / "T_core_trajectories.npz")
    eng_traj = load_traj(DATA / "T_core_trajectories_engine.npz")

    plot_dT_with_engine(cw_rows, eng_rows)
    plot_engine_trajectories(eng_rows, eng_traj)
    plot_4panel(cw_rows, eng_rows, cw_traj, eng_traj)
    print("Done.")


if __name__ == "__main__":
    main()
