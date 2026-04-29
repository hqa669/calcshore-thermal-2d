#!/usr/bin/env python3
"""S1.2 — Generate CW temperature distribution plots for all 9 runs.

Per-run 3-panel figure (t≈0, 84 hr, 168 hr) and 3×3 composite at t=168 hr.

Usage:
    python validation/soil_calibration/plot_cw_distributions.py

Writes:
    validation/soil_calibration/plots/run{A..I}_distribution.png
    validation/soil_calibration/plots/all_runs_t168.png
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
from cw_scenario_loader import parse_cw_temp_output

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
PLOTS = os.path.join(HERE, "plots")
os.makedirs(PLOTS, exist_ok=True)

RUNS = [
    ("A", "runA_baseline",  73,  73),
    ("B", "runB_73_60",     73,  60),
    ("C", "runC_73_90",     73,  90),
    ("D", "runD_60_73",     60,  73),
    ("E", "runE_90_73",     90,  73),
    ("F", "runF_73_45",     73,  45),
    ("G", "runG_73_100",    73, 100),
    ("H", "runH_45_73",     45,  73),
    ("I", "runI_100_73",   100,  73),
]

VMIN, VMAX = 45.0, 100.0
TARGET_HRS = [0.0, 84.0, 168.0]


def nearest_time_idx(time_hrs, target):
    return int(np.abs(time_hrs - target).argmin())


def plot_cw_field(field_F, widths_m, depths_m, title, ax, vmin=VMIN, vmax=VMAX):
    """Draw one 2D CW temperature slice on ax.

    field_F : (nD, nW) — depth rows, width columns; widths_m is descending
    Depth increases downward (depth index 0 = top surface).
    Width index 0 = centerline (max width), index -1 = edge (0 m).
    """
    # widths_m is descending: e.g. [6.1, 5.59, ..., 0]. Flip for imshow extent.
    w_min = float(widths_m[-1])   # 0 m (edge)
    w_max = float(widths_m[0])    # 6.1 m (centerline)
    d_min = float(depths_m[0])    # 0 m (top)
    d_max = float(depths_m[-1])   # 24.38 m (bottom)

    # imshow expects (row=y, col=x); field_F is (nD, nW) with nD rows = depth axis
    # widths_m is descending, so field_F[:, 0] = centerline.
    # Flip width axis so image x-axis goes from edge (left) to centerline (right).
    field_plot = field_F[:, ::-1]   # (nD, nW) with col 0 = edge, col -1 = centerline

    extent = [w_min, w_max, d_max, d_min]   # [left, right, bottom, top]
    im = ax.imshow(
        field_plot,
        aspect="equal",
        origin="upper",
        extent=extent,
        cmap="inferno",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )

    # 1°F isotherm contour overlay
    # Build the meshgrid in the flipped-width space
    w_asc = widths_m[::-1]           # ascending
    W, D = np.meshgrid(w_asc, depths_m)
    levels = np.arange(np.floor(vmin / 1) * 1, vmax + 1, 1.0)
    cs = ax.contour(W, D, field_plot, levels=levels, colors="white", linewidths=0.3, alpha=0.5)
    # Label only every 5°F to avoid clutter
    label_levels = [l for l in levels if l % 5 == 0]
    ax.clabel(cs, levels=label_levels, fmt="%g°F", fontsize=5, inline=True, inline_spacing=1)

    ax.set_xlabel("Width (m, 0=edge, max=CL)", fontsize=7)
    ax.set_ylabel("Depth (m)", fontsize=7)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=6)
    return im


def make_run_figure(label, folder, placement, soil):
    dat_path = os.path.join(CW_RUNS, folder, "output.txt")
    if not os.path.isfile(dat_path):
        print(f"  SKIP Run {label}: output.txt not found")
        return None

    v = parse_cw_temp_output(dat_path)
    # v.T_field_F: (n_time, nD, nW), widths_m descending, depths_m ascending

    t_indices = [nearest_time_idx(v.time_hrs, t) for t in TARGET_HRS]
    t_labels = ["t≈0 min", "t=84 hr", "t=168 hr"]

    delta_T = soil - placement
    run_title = f"Run {label}: placement={placement}°F, soil={soil}°F, ΔT={delta_T:+d}°F"

    fig, axes = plt.subplots(1, 3, figsize=(12, 5))
    fig.suptitle(run_title, fontsize=10, fontweight="bold")

    im = None
    for ax, ti, tlab in zip(axes, t_indices, t_labels):
        actual_t = v.time_hrs[ti]
        title = f"{tlab} (actual {actual_t:.2f} hr)"
        im = plot_cw_field(v.T_field_F[ti], v.widths_m, v.depths_m, title, ax)

    # Shared colorbar
    fig.subplots_adjust(right=0.88, wspace=0.35)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Temperature (°F)", fontsize=8)

    out_path = os.path.join(PLOTS, f"run{label}_distribution.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written: {out_path}")
    return v


def make_composite(run_data):
    """3×3 composite of all 9 runs at t=168 hr."""
    fig, axes = plt.subplots(3, 3, figsize=(15, 13))
    fig.suptitle("CW Temperature Distribution at t=168 hr — All 9 Runs", fontsize=12)

    run_order = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    im = None
    for idx, label in enumerate(run_order):
        ax = axes[idx // 3][idx % 3]
        if label not in run_data:
            ax.set_visible(False)
            continue

        v, placement, soil = run_data[label]
        ti = nearest_time_idx(v.time_hrs, 168.0)
        delta_T = soil - placement
        title = f"Run {label}: {placement}°F / {soil}°F (ΔT={delta_T:+d})"
        im = plot_cw_field(v.T_field_F[ti], v.widths_m, v.depths_m, title, ax)

    if im is not None:
        fig.subplots_adjust(right=0.88, hspace=0.4, wspace=0.4)
        cbar_ax = fig.add_axes([0.90, 0.1, 0.015, 0.8])
        cbar = fig.colorbar(im, cax=cbar_ax)
        cbar.set_label("Temperature (°F)", fontsize=9)

    out_path = os.path.join(PLOTS, "all_runs_t168.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Written composite: {out_path}")


def run():
    run_data = {}
    for label, folder, placement, soil in RUNS:
        print(f"Processing Run {label} ({folder})…")
        v = make_run_figure(label, folder, placement, soil)
        if v is not None:
            run_data[label] = (v, placement, soil)

    if run_data:
        print("Building composite…")
        make_composite(run_data)

    print("Done.")


if __name__ == "__main__":
    run()
