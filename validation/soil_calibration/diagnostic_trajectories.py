#!/usr/bin/env python3
"""S1.4 — Cross-run cooling/heating profiles at three diagnostic points.

Usage:
    python validation/soil_calibration/diagnostic_trajectories.py

Writes:
    validation/soil_calibration/plots/trajectory_side_near.png
    validation/soil_calibration/plots/trajectory_bottom_near.png
    validation/soil_calibration/plots/trajectory_deep_interior.png
    validation/soil_calibration/STAGE1_diagnostic_trajectories.md
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

TARGET_WIDTH_M_SIDE  = 0.51   # just inside concrete from side (~1.67 ft)
TARGET_DEPTH_M_MID   = 12.19  # mid depth (~40 ft in 80 ft deep mat)
TARGET_DEPTH_M_BOT   = 23.88  # near bottom (~78.35 ft)
TARGET_WIDTH_M_CL    = 6.10   # centerline

COLORS = plt.cm.tab10(np.linspace(0, 1, 9))
LINESTYLES = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-"]


def find_nearest_idx(arr, val):
    return int(np.abs(np.asarray(arr) - val).argmin())


def run():
    print("Loading all 9 CW runs…")
    data = {}
    for label, folder, placement, soil in RUNS:
        path = os.path.join(CW_RUNS, folder, "output.txt")
        if not os.path.isfile(path):
            print(f"  SKIP {label}: output.txt not found")
            continue
        v = parse_cw_temp_output(path)
        data[label] = (v, placement, soil)

    if not data:
        print("No data loaded — aborting.")
        return 1

    # Determine grid indices once from first loaded run
    sample_v = next(iter(data.values()))[0]
    widths_m = sample_v.widths_m   # descending: [6.10, 5.59, …, 0.51, 0]
    depths_m = sample_v.depths_m   # ascending:  [0, 0.51, …, 24.38]

    # width index 0 = centerline (max width); index -1 = edge (0 m)
    wi_side = find_nearest_idx(widths_m, TARGET_WIDTH_M_SIDE)   # near edge
    wi_cl   = find_nearest_idx(widths_m, TARGET_WIDTH_M_CL)      # centerline
    di_mid  = find_nearest_idx(depths_m, TARGET_DEPTH_M_MID)
    di_bot  = find_nearest_idx(depths_m, TARGET_DEPTH_M_BOT)

    actual_w_side = widths_m[wi_side]
    actual_w_cl   = widths_m[wi_cl]
    actual_d_mid  = depths_m[di_mid]
    actual_d_bot  = depths_m[di_bot]

    print(f"  Side-near index: width={actual_w_side:.2f}m (target {TARGET_WIDTH_M_SIDE}m), "
          f"depth={actual_d_mid:.2f}m")
    print(f"  Bottom-near index: width={actual_w_cl:.2f}m (CL), depth={actual_d_bot:.2f}m")
    print(f"  Deep interior index: width={actual_w_cl:.2f}m (CL), depth={actual_d_mid:.2f}m")

    POINTS = [
        ("side_near",      wi_side, di_mid, "Side-near (w≈0.51m, d≈12m)"),
        ("bottom_near",    wi_cl,   di_bot, "Bottom-near (centerline, d≈24m)"),
        ("deep_interior",  wi_cl,   di_mid, "Deep interior (centerline, d≈12m)"),
    ]

    figures = {}
    for point_name, wi, di, point_label in POINTS:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.set_title(f"Temperature Trajectory — {point_label}", fontsize=11)
        ax.set_xlabel("Time (hr)", fontsize=10)
        ax.set_ylabel("Temperature (°F)", fontsize=10)
        ax.grid(True, alpha=0.3)

        for i, (label, folder, placement, soil) in enumerate(RUNS):
            if label not in data:
                continue
            v, p, s = data[label]
            T_traj = v.T_field_F[:, di, wi]  # (n_time,)
            delta_T = s - p
            line_label = f"Run {label}: {p}°F/{s}°F (ΔT={delta_T:+d})"
            ax.plot(
                v.time_hrs, T_traj,
                color=COLORS[i], linestyle=LINESTYLES[i],
                linewidth=1.2, label=line_label, alpha=0.85,
            )

        ax.legend(fontsize=7, loc="best")
        fig.tight_layout()
        figures[point_name] = fig

        out_path = os.path.join(PLOTS, f"trajectory_{point_name}.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {out_path}")

    # --- Analysis ---
    # Check symmetry: B mirrors D, C mirrors E
    sym_analysis = []
    for lp, ln, pair_name in [("B", "D", "B↔D"), ("C", "E", "C↔E"),
                               ("F", "H", "F↔H"), ("G", "I", "G↔I")]:
        if lp not in data or ln not in data:
            sym_analysis.append(f"- {pair_name}: data missing")
            continue
        vp, pp, sp = data[lp]
        vn, pn, sn = data[ln]
        va, pa, sa = data["A"]

        for point_name, wi, di, point_label in POINTS:
            Tp = vp.T_field_F[:, di, wi] - pa   # deviation from placement A
            Tn = vn.T_field_F[:, di, wi] - pa
            # For symmetric pair: Tp(t) + Tn(t) should equal 2*baseline(t)
            # => Tp deviation from baseline = -Tn deviation from baseline
            dTp = vp.T_field_F[:, di, wi] - va.T_field_F[:, di, wi]
            dTn = vn.T_field_F[:, di, wi] - va.T_field_F[:, di, wi]
            max_asym = float(np.abs(dTp + dTn).max())
            sym_analysis.append(
                f"- {pair_name} at {point_label}: max asymmetry = {max_asym:.3f}°F"
            )

    # Saturation temperature analysis (final value at t=168)
    sat_analysis = []
    for label, folder, placement, soil in RUNS:
        if label not in data:
            continue
        v, p, s = data[label]
        for point_name, wi, di, point_label in POINTS:
            T_final = float(v.T_field_F[-1, di, wi])
            sat_analysis.append(f"Run {label} ({p}°F/{s}°F) at {point_label}: T_final={T_final:.2f}°F")

    # Linearity check: does delta-T at t=168 scale linearly with |ΔT|?
    lin_analysis = []
    baseline_v = data["A"][0]
    for point_name, wi, di, point_label in POINTS:
        T_base_168 = float(baseline_v.T_field_F[-1, di, wi])
        lin_analysis.append(f"\n### {point_label}")
        lin_analysis.append("| Run | ΔT (soil-placement) | ΔT_actual at t=168 | Ratio |")
        lin_analysis.append("| --- | ------------------- | ------------------ | ----- |")
        dt_vals = []
        delta_T_inputs = []
        for label, folder, placement, soil in RUNS:
            if label == "A" or label not in data:
                continue
            v, p, s = data[label]
            T_168 = float(v.T_field_F[-1, di, wi])
            dT_actual = T_168 - T_base_168
            dT_input = s - p
            ratio = dT_actual / dT_input if dT_input != 0 else float("nan")
            lin_analysis.append(
                f"| {label} | {dT_input:+d}°F | {dT_actual:+.3f}°F | {ratio:.4f} |"
            )
            dt_vals.append(dT_actual)
            delta_T_inputs.append(dT_input)
        # Fit linear
        if len(dt_vals) >= 2:
            x = np.array(delta_T_inputs, dtype=float)
            y = np.array(dt_vals, dtype=float)
            slope = float(np.polyfit(x, y, 1)[0])
            residuals = y - slope * x
            lin_analysis.append(f"\nLinear fit slope: {slope:.4f} (1.0 = perfect linearity)")
            lin_analysis.append(f"Max residual from linear fit: {np.abs(residuals).max():.4f}°F")

    # Write MD report
    md = [
        "# STAGE1 — Diagnostic Temperature Trajectories",
        "",
        "## Diagnostic Points",
        "",
        f"| Point | Width (m) | Depth (m) | Description |",
        f"| ----- | --------- | --------- | ----------- |",
        f"| side_near | {widths_m[wi_side]:.2f} (idx {wi_side}) | {depths_m[di_mid]:.2f} (idx {di_mid}) | Just inside concrete from side soil interface, mid-depth |",
        f"| bottom_near | {widths_m[wi_cl]:.2f} (idx {wi_cl}) | {depths_m[di_bot]:.2f} (idx {di_bot}) | Centerline, near bottom soil interface |",
        f"| deep_interior | {widths_m[wi_cl]:.2f} (idx {wi_cl}) | {depths_m[di_mid]:.2f} (idx {di_mid}) | Centerline, mid-depth |",
        "",
        "## Figures",
        "",
        "- `plots/trajectory_side_near.png`",
        "- `plots/trajectory_bottom_near.png`",
        "- `plots/trajectory_deep_interior.png`",
        "",
        "## Symmetry Analysis",
        "",
        "Max asymmetry `dT_pos + dT_neg` (should be 0 if linear):",
        "",
    ] + sym_analysis + [
        "",
        "## Saturation Temperature at t=168 hr",
        "",
        "Expected: if system fully equilibrates, T_final → soil_temp_F.",
        "Deviation indicates the thermal mass is still relaxing.",
        "",
    ] + sat_analysis + [
        "",
        "## Linearity of ΔT Response",
        "",
        "If cooling/heating rate scales linearly with |ΔT|, the ratio ΔT_actual/ΔT_input",
        "should be constant across runs at the same diagnostic point.",
        "",
    ] + lin_analysis

    out_path = os.path.join(HERE, "STAGE1_diagnostic_trajectories.md")
    with open(out_path, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"\nWritten: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
