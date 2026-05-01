#!/usr/bin/env python3
"""S1.5 — Penetration depth analysis: 1°F front vs √t.

Usage:
    python validation/soil_calibration/penetration_analysis.py

Writes:
    validation/soil_calibration/plots/penetration_side.png
    validation/soil_calibration/plots/penetration_bottom.png
    validation/soil_calibration/STAGE1_penetration_analysis.md
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import linregress

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

FRONT_THRESHOLD_F = 1.0
LINEAR_FIT_END_HR = 168.0   # use full run; CW's 0.51m grid makes early-window fit unreliable
COLORS = plt.cm.tab10(np.linspace(0, 1, 9))
LINESTYLES = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-"]


def find_nearest_idx(arr, val):
    return int(np.abs(np.asarray(arr) - val).argmin())


def compute_penetration_side(T_field_F, widths_m, depths_m, placement_F, di_mid):
    """Penetration depth from side interface into concrete.

    Walk from width index 12 (edge, 0 m) inward (decreasing idx / increasing width).
    Find the first width where |T - placement_F| < FRONT_THRESHOLD_F.
    Returns array of penetration distances (m) per time step.
    """
    n_time = T_field_F.shape[0]
    penetration = np.zeros(n_time)

    edge_idx = len(widths_m) - 1   # edge = index 12 (0 m), CL = index 0 (6.1 m)
    for ti in range(n_time):
        T_row = T_field_F[ti, di_mid, :]   # (nW,) — index 0 = centerline, index 12 = edge
        # widths_m is descending: [6.1, 5.59, ..., 0.51, 0]
        # Walk from edge (ii=12) toward CL (ii=0).
        # The 1°F front is the last point (from edge) where |T-T_place| >= threshold.
        # widths_m[ii] = distance from edge (0m) for that node.
        pen_m = 0.0
        for ii in range(edge_idx, -1, -1):  # 12, 11, 10, ..., 0
            if abs(T_row[ii] - placement_F) < FRONT_THRESHOLD_F:
                # This node is unaffected; last affected was ii+1
                if ii < edge_idx:
                    pen_m = float(widths_m[ii + 1])  # widths_m[13..edge+1] ascending from 0
                # else pen_m = 0.0 (even edge is unaffected = no penetration yet)
                break
        else:
            # All nodes exceed threshold: front has passed CL
            pen_m = float(widths_m[0])
        penetration[ti] = pen_m

    return penetration


def compute_penetration_bottom(T_field_F, widths_m, depths_m, placement_F, wi_cl):
    """Penetration depth from bottom interface upward into concrete.

    Walk from depth index 48 (max depth) upward (decreasing index).
    Find the first depth where |T - placement_F| < FRONT_THRESHOLD_F.
    Returns penetration from bottom as m.
    """
    n_time = T_field_F.shape[0]
    penetration = np.zeros(n_time)
    bottom_idx = len(depths_m) - 1  # 48 = deepest point

    # depths_m is ascending: [0, 0.51, ..., 24.38]
    # Bottom = index 48 (24.38m), top = index 0 (0m)
    # Walk from bottom (ii=48) upward (decreasing ii).
    # The 1°F front = last affected node from bottom.
    # penetration = depths_m[bottom_idx] - depths_m[ii_front]  (distance from bottom)
    for ti in range(n_time):
        T_col = T_field_F[ti, :, wi_cl]   # (nD,) — depth column at centerline
        pen_m = 0.0
        for ii in range(bottom_idx, -1, -1):  # 48, 47, ..., 0
            if abs(T_col[ii] - placement_F) < FRONT_THRESHOLD_F:
                # Unaffected here; last affected was ii+1
                if ii < bottom_idx:
                    pen_m = float(depths_m[bottom_idx] - depths_m[ii + 1])
                break
        else:
            # All nodes affected: front has passed top
            pen_m = float(depths_m[bottom_idx] - depths_m[0])
        penetration[ti] = pen_m

    return penetration


def fit_alpha_eff(time_hrs, penetration_m, end_hr=LINEAR_FIT_END_HR):
    """Fit penetration ∝ √t to extract effective diffusivity.

    slope = 2 * √α_eff  (half-space approximation: x_front = 2 * erf_inv(1 - δ) * √(α t))
    For a 1°F front on a semi-infinite slab with step BC, the leading factor is
    approximately 2·erfinv(1 - 1/|ΔT_applied|), but for simplicity we just
    fit pen vs √t and report the slope (units: m/√hr).

    Returns (slope_m_per_sqrth, r_value)
    """
    mask = (time_hrs > 0) & (time_hrs <= end_hr)
    if mask.sum() < 5:
        return float("nan"), float("nan")
    sqrt_t = np.sqrt(time_hrs[mask])
    pen = penetration_m[mask]
    # Only fit where penetration > 0
    valid = pen > 0
    if valid.sum() < 5:
        return float("nan"), float("nan")
    slope, intercept, r, p, se = linregress(sqrt_t[valid], pen[valid])
    return float(slope), float(r)


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

    sample_v = next(iter(data.values()))[0]
    widths_m = sample_v.widths_m
    depths_m = sample_v.depths_m

    wi_cl = find_nearest_idx(widths_m, 6.10)   # centerline
    di_mid = find_nearest_idx(depths_m, 12.19)  # mid depth

    fig_side, ax_side = plt.subplots(figsize=(10, 6))
    ax_side.set_title(f"Side Penetration of {FRONT_THRESHOLD_F}°F Front vs √t", fontsize=11)
    ax_side.set_xlabel("√(time) (√hr)", fontsize=10)
    ax_side.set_ylabel("Penetration from side edge (m)", fontsize=10)
    ax_side.grid(True, alpha=0.3)

    fig_bot, ax_bot = plt.subplots(figsize=(10, 6))
    ax_bot.set_title(f"Bottom Penetration of {FRONT_THRESHOLD_F}°F Front vs √t", fontsize=11)
    ax_bot.set_xlabel("√(time) (√hr)", fontsize=10)
    ax_bot.set_ylabel("Penetration from bottom (m)", fontsize=10)
    ax_bot.grid(True, alpha=0.3)

    md_rows = []
    md_rows.append(
        "| Run | Placement°F | Soil°F | ΔT | "
        "side slope (m/√hr) | side R² | side pen@168hr (m) | "
        "bot slope (m/√hr) | bot R² | bot pen@168hr (m) |"
    )
    md_rows.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")

    for i, (label, folder, placement, soil) in enumerate(RUNS):
        if label not in data:
            continue
        v, p, s = data[label]
        delta_T = s - p
        if delta_T == 0:
            # Baseline: no driving force, expect no penetration
            md_rows.append(f"| {label} | {p} | {s} | 0 | — (baseline) | — | — | — | — | — |")
            continue

        pen_side = compute_penetration_side(v.T_field_F, widths_m, depths_m, float(p), di_mid)
        pen_bot  = compute_penetration_bottom(v.T_field_F, widths_m, depths_m, float(p), wi_cl)

        sqrt_t = np.sqrt(v.time_hrs)
        lbl = f"Run {label}: {p}°F/{s}°F (ΔT={delta_T:+d})"

        ax_side.plot(sqrt_t, pen_side, color=COLORS[i], linestyle=LINESTYLES[i],
                     linewidth=1.2, label=lbl, alpha=0.85)
        ax_bot.plot(sqrt_t, pen_bot, color=COLORS[i], linestyle=LINESTYLES[i],
                    linewidth=1.2, label=lbl, alpha=0.85)

        slope_side, r_side = fit_alpha_eff(v.time_hrs, pen_side)
        slope_bot, r_bot   = fit_alpha_eff(v.time_hrs, pen_bot)

        r2_side = r_side**2 if not np.isnan(r_side) else float("nan")
        r2_bot  = r_bot**2  if not np.isnan(r_bot)  else float("nan")

        # Max penetration at t=168hr
        pen_side_168 = float(pen_side[-1])
        pen_bot_168  = float(pen_bot[-1])

        md_rows.append(
            f"| {label} | {p} | {s} | {delta_T:+d} "
            f"| {slope_side:.4f} | {r2_side:.4f} | {pen_side_168:.2f} "
            f"| {slope_bot:.4f} | {r2_bot:.4f} | {pen_bot_168:.2f} |"
        )

        print(f"  Run {label}: side slope={slope_side:.4f} m/√hr (R²={r2_side:.3f}), "
              f"bot slope={slope_bot:.4f} m/√hr (R²={r2_bot:.3f})")

    for ax in (ax_side, ax_bot):
        ax.legend(fontsize=7, loc="upper left")

    out_side = os.path.join(PLOTS, "penetration_side.png")
    out_bot  = os.path.join(PLOTS, "penetration_bottom.png")
    fig_side.tight_layout()
    fig_bot.tight_layout()
    fig_side.savefig(out_side, dpi=150)
    fig_bot.savefig(out_bot, dpi=150)
    plt.close(fig_side)
    plt.close(fig_bot)
    print(f"  Saved: {out_side}")
    print(f"  Saved: {out_bot}")

    # Write MD
    md = [
        "# STAGE1 — Penetration Depth Analysis",
        "",
        f"## Method",
        "",
        f"1°F thermal front penetration from each soil-concrete interface vs √t.",
        "Linear fit across full 168hr run → effective slope (m/√hr).",
        "In a half-space with step BC: penetration ≈ C·√(α_eff·t).",
        "",
        f"**Grid resolution caveat**: CW's output grid has 0.51m node spacing.",
        "The 1°F front position is quantized to 0.51m steps. This makes the step",
        "curves jagged and degrades R² even when the underlying physics is √t-like.",
        "The slope and R² values should be treated as indicative, not precise.",
        "",
        f"- Threshold: {FRONT_THRESHOLD_F}°F from placement temperature",
        f"- Linear fit window: 0 – {LINEAR_FIT_END_HR} hr (full run)",
        f"- Side: walk from edge (0 m) inward at mid-depth ({depths_m[di_mid]:.2f} m)",
        f"- Bottom: walk from bottom upward on centerline (w={widths_m[wi_cl]:.2f} m)",
        "",
        "## Results",
        "",
    ] + md_rows + [
        "",
        "## Interpretation",
        "",
        "- If slope values cluster across runs: single effective diffusivity, model is",
        "  consistent (property is not ΔT-dependent).",
        "- If slopes scatter: model-form non-linearity or material property mismatch.",
        "- R² < 0.95 on early-time linear fit: penetration is not √t-like even at early time",
        "  (e.g., BC-dominated or grid-resolution artifact).",
        "- R² drops at long time (deviation from linearity) is expected when the finite",
        "  soil domain saturates to Dirichlet BC temperature.",
        "",
        "## Plots",
        "",
        "- `plots/penetration_side.png`",
        "- `plots/penetration_bottom.png`",
    ]

    out_md = os.path.join(HERE, "STAGE1_penetration_analysis.md")
    with open(out_md, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"\nWritten: {out_md}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
