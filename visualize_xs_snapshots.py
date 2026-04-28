"""
Visualize 2D temperature distribution at fixed time intervals.

Generates a grid of subplots showing the half-mat cross-section temperature
field at every N hours. Each subplot is a heatmap with depth on the y-axis
and width on the x-axis (width=0 is the form face / corner edge).

Two input modes (mutually exclusive):
  --input-cw-output PATH   parse a CW temp.txt export via parse_cw_temp_output
  --input-npz PATH         load an npz saved by thermal_engine_test.py
                           (arrays: T_field_F, time_hrs, depths_m, widths_m)

Usage examples:

  # CW output (existing workflow):
  python visualize_xs_snapshots.py \
      --input-cw-output ~/Downloads/HydrationCenter_mix01/output.txt \
      --out-xs-png ~/Downloads/HydrationCenter_mix01/xs_temperature_snapshots.png \
      --out-centerline-png ~/Downloads/HydrationCenter_mix01/centerline_column_profiles.png

  # Engine npz (kinetics-isolation test):
  python visualize_xs_snapshots.py \
      --input-npz diagnostics/sprint5/kinetics_isolation/kinetics_isolation_xs_snapshots.npz \
      --out-xs-png diagnostics/sprint5/kinetics_isolation/kinetics_isolation_xs_snapshots.png \
      --out-centerline-png diagnostics/sprint5/kinetics_isolation/kinetics_isolation_centerline_profiles.png

Run from the repo root so that cw_scenario_loader is importable for --input-cw-output mode.
"""
import argparse
import os
import sys
from types import SimpleNamespace

import numpy as np
import matplotlib.pyplot as plt


def _load_cw_output(path: str):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from cw_scenario_loader import parse_cw_temp_output
    return parse_cw_temp_output(path)


def _load_npz(path: str):
    data = np.load(path)
    missing = [k for k in ("T_field_F", "time_hrs", "depths_m", "widths_m")
               if k not in data]
    if missing:
        raise ValueError(f"npz missing required arrays: {missing}")
    return SimpleNamespace(
        T_field_F=data["T_field_F"],
        time_hrs=data["time_hrs"],
        depths_m=data["depths_m"],
        widths_m=data["widths_m"],
    )


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Visualize 2D temperature cross-section snapshots.",
    )
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--input-cw-output", metavar="PATH",
                     help="CW temp.txt export file.")
    src.add_argument("--input-npz", metavar="PATH",
                     help="Engine npz (T_field_F, time_hrs, depths_m, widths_m).")
    parser.add_argument("--out-xs-png", metavar="PATH", default=None,
                        help="Output path for cross-section snapshot grid.")
    parser.add_argument("--out-centerline-png", metavar="PATH", default=None,
                        help="Output path for centerline column profiles.")
    parser.add_argument("--snapshot-interval-hr", type=float, default=5.0,
                        metavar="HR",
                        help="Interval between snapshots in hours (default 5).")
    parser.add_argument("--zoom-corner-fraction", type=float, default=None,
                        metavar="F",
                        help="Show only the corner-most fraction of the width "
                             "(e.g. 0.05). Default: show full width.")
    args = parser.parse_args(argv)

    # ---- LOAD DATA ----
    if args.input_cw_output:
        print(f"Loading CW output from {args.input_cw_output}...")
        v = _load_cw_output(args.input_cw_output)
        out_dir = os.path.dirname(os.path.abspath(args.input_cw_output))
        label = "CW 2D Temperature Distribution"
    else:
        print(f"Loading engine npz from {args.input_npz}...")
        v = _load_npz(args.input_npz)
        out_dir = os.path.dirname(os.path.abspath(args.input_npz))
        label = "Engine 2D Temperature Distribution"

    nt, nd, nw = v.T_field_F.shape
    print(f"  Shape: {v.T_field_F.shape}  (n_time, n_depth, n_width)")
    print(f"  Time:  {v.time_hrs[0]:.2f} to {v.time_hrs[-1]:.2f} hrs")
    print(f"  Depth: {v.depths_m[0]:.2f} to {v.depths_m[-1]:.2f} m  ({nd} nodes)")
    print(f"  Width: {v.widths_m[0]:.2f} to {v.widths_m[-1]:.2f} m  ({nw} nodes)")
    print(f"  T:     {v.T_field_F.min():.1f} to {v.T_field_F.max():.1f} F")

    # Default output paths
    out_xs = args.out_xs_png or os.path.join(out_dir, "xs_temperature_snapshots.png")
    out_cl = args.out_centerline_png or os.path.join(out_dir, "centerline_column_profiles.png")
    snap_hr = args.snapshot_interval_hr

    # ---- PICK SNAPSHOT TIMES ----
    t_targets = np.arange(0.0, v.time_hrs[-1] + 0.01, snap_hr)
    snapshot_indices = [int(np.abs(v.time_hrs - tt).argmin()) for tt in t_targets]
    n_snap = len(snapshot_indices)
    print(f"\nGenerating {n_snap} snapshots at {snap_hr} hr intervals")

    # ---- ZOOM (OPTIONAL) ----
    if args.zoom_corner_fraction is not None:
        n_zoom_cols = max(2, int(nw * args.zoom_corner_fraction))
        w_slice = slice(nw - n_zoom_cols, nw)
        print(f"  Zoom: showing last {n_zoom_cols} of {nw} width nodes (corner)")
    else:
        w_slice = slice(0, nw)
    field_slice = v.T_field_F[:, :, w_slice]
    widths_plot = v.widths_m[w_slice]

    # ---- PLOT 1: XS grid ----
    T_min_global = field_slice.min()
    T_max_global = field_slice.max()
    print(f"  Color scale: {T_min_global:.1f} to {T_max_global:.1f} F\n")

    n_cols = int(np.ceil(np.sqrt(n_snap)))
    n_rows = int(np.ceil(n_snap / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.5 * n_cols, 2.8 * n_rows),
                             squeeze=False)

    widths_ft = widths_plot * 3.28084
    depths_ft = v.depths_m * 3.28084
    extent = [widths_ft[0], widths_ft[-1], depths_ft[-1], depths_ft[0]]

    for i_ax, (t_idx, ax) in enumerate(zip(snapshot_indices, axes.flat)):
        field = field_slice[t_idx]
        im = ax.imshow(field, aspect="auto", extent=extent,
                       vmin=T_min_global, vmax=T_max_global,
                       cmap="inferno", origin="upper")
        ax.set_title(f"t = {v.time_hrs[t_idx]:.1f} hr   "
                     f"(min={field.min():.1f}, max={field.max():.1f} F)",
                     fontsize=9)
        if i_ax % n_cols == 0:
            ax.set_ylabel("Depth (ft)\n(top -> bottom)", fontsize=8)
        if i_ax >= (n_rows - 1) * n_cols:
            ax.set_xlabel("Width (ft)\n(0 = centerline, -> form face)", fontsize=8)
        ax.tick_params(labelsize=7)

    for j in range(n_snap, n_rows * n_cols):
        axes.flat[j].set_visible(False)

    fig.subplots_adjust(right=0.90, hspace=0.45, wspace=0.3, top=0.94, bottom=0.08)
    cbar_ax = fig.add_axes([0.92, 0.10, 0.015, 0.80])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Temperature (°F)", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    zoom_label = (f" (zoomed to corner, last {len(widths_plot)}/{nw} cols)"
                  if args.zoom_corner_fraction else "")
    fig.suptitle(f"{label} — every {snap_hr} hr" + zoom_label,
                 fontsize=11, fontweight="bold")

    print(f"Saving XS snapshot grid to {out_xs}...")
    plt.savefig(out_xs, dpi=120, bbox_inches="tight")
    plt.close(fig)

    # ---- PLOT 2: centerline column profiles ----
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
    cmap = plt.cm.viridis
    for k, t_idx in enumerate(snapshot_indices):
        color = cmap(k / max(1, n_snap - 1))
        centerline_col = v.T_field_F[t_idx, :, 0]  # all depths, width index 0
        ax2.plot(centerline_col, depths_ft, color=color,
                 label=f"t={v.time_hrs[t_idx]:.0f} hr", linewidth=1.5)
    ax2.invert_yaxis()
    ax2.set_xlabel("Temperature (°F)")
    ax2.set_ylabel("Depth (ft)  (top -> bottom)")
    ax2.set_title("Centerline column temperature profile vs depth\n"
                  f"(width index 0, every {snap_hr} hr)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="center right", fontsize=8, ncol=2)
    fig2.tight_layout()

    print(f"Saving centerline profiles to {out_cl}...")
    plt.savefig(out_cl, dpi=120, bbox_inches="tight")
    plt.close(fig2)

    print("\nDone. Open the PNG files to inspect.")


if __name__ == "__main__":
    main()
