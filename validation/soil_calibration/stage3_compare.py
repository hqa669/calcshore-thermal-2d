#!/usr/bin/env python3
"""Stage 3 — Engine-vs-CW comparison: residual tables + plots.

Reads engine CSVs from a stage3_fix{1|2}_runs/ directory and CW output.txt files.
Resamples engine field onto CW grid (bilinear), computes residuals, generates plots.

Usage:
    python validation/soil_calibration/stage3_compare.py --fix 1
    python validation/soil_calibration/stage3_compare.py --fix 2

Outputs for --fix 1:
    validation/soil_calibration/STAGE3_fix1_residuals.md

Outputs for --fix 2:
    validation/soil_calibration/STAGE3_fix2_residuals.md
    validation/soil_calibration/plots/STAGE3_run{A..I}_engine_vs_cw.png
"""
import argparse
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator

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

COMPARE_HR = 168.0
VMIN, VMAX = 45.0, 100.0


def load_engine_csv(path):
    with open(path) as f:
        lines = f.readlines()
    header = lines[0].strip().split(",")
    x_m = np.array([float(h) for h in header[1:]])
    y_m_list, rows = [], []
    for line in lines[1:]:
        parts = line.strip().split(",")
        y_m_list.append(float(parts[0]))
        rows.append([float(v) for v in parts[1:]])
    return np.array(y_m_list), x_m, np.array(rows)


def load_cw_slice(folder, target_hr):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    if not os.path.isfile(path):
        return None, None, None, None
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - target_hr).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m, float(v.time_hrs[ti])


def resample_engine_to_cw(eng_y_m, eng_x_m, eng_field_F, cw_depths_m, cw_widths_m):
    interp = RegularGridInterpolator(
        (eng_y_m, eng_x_m), eng_field_F,
        method="linear", bounds_error=False, fill_value=None,
    )
    D, W = np.meshgrid(cw_depths_m, cw_widths_m, indexing="ij")
    pts = np.column_stack([D.ravel(), W.ravel()])
    return interp(pts).reshape(len(cw_depths_m), len(cw_widths_m))


def plot_comparison(label, placement, soil, fix_num,
                    cw_field, eng_field, residual,
                    cw_depths_m, cw_widths_m):
    delta_T = soil - placement
    w_min = float(cw_widths_m[-1])
    w_max = float(cw_widths_m[0])
    d_max = float(cw_depths_m[-1])
    extent = [w_min, w_max, d_max, 0.0]

    cw_plot  = cw_field[:, ::-1]
    eng_plot = eng_field[:, ::-1]
    res_plot = residual[:, ::-1]
    sym = max(0.01, float(np.abs(residual).max()))

    fig, axes = plt.subplots(3, 1, figsize=(6, 14))
    fig.suptitle(
        f"Stage 3 Fix {fix_num} — Run {label}: placement={placement}°F, "
        f"soil={soil}°F, ΔT={delta_T:+d}°F\nt={COMPARE_HR:.0f} hr",
        fontsize=10, fontweight="bold",
    )
    panel_data = [
        (cw_plot,  "CW (reference)",             "inferno", (VMIN, VMAX)),
        (eng_plot, "Engine (resampled)",          "inferno", (VMIN, VMAX)),
        (res_plot, f"Engine−CW (max|R|={sym:.2f}°F)", "RdBu_r", (-sym, sym)),
    ]
    for ax, (fdata, title, cmap, (vmin, vmax)) in zip(axes, panel_data):
        im = ax.imshow(fdata, aspect="equal", origin="upper",
                       extent=extent, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Width (m, 0=edge, max=CL)", fontsize=7)
        ax.set_ylabel("Depth (m)", fontsize=7)
        ax.tick_params(labelsize=6)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=6)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    out_path = os.path.join(PLOTS, f"STAGE3_run{label}_engine_vs_cw.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run(fix_num: int):
    eng_dir = os.path.join(HERE, f"stage3_fix{fix_num}_runs")
    if not os.path.isdir(eng_dir):
        print(f"ERROR: {eng_dir} not found. Run stage3_run_fix{fix_num}.py first.")
        return 1

    manifest_path = os.path.join(eng_dir, "manifest.json")
    with open(manifest_path) as f:
        manifest = json.load(f)

    do_plots = (fix_num == 2)

    rows = []
    rows.append("| Run | Placement/Soil | max\\|R\\| (°F) | mean\\|R\\| (°F) | RMS R (°F) | Location of max |")
    rows.append("| --- | -------------- | -------------- | --------------- | ------------ | --------------- |")

    for label, folder, placement, soil in RUNS:
        if manifest.get(label, {}).get("status") != "ok":
            rows.append(f"| {label} | {placement}/{soil} | — | — | — | skip |")
            continue

        eng_csv = os.path.join(eng_dir, f"run{label}_t168.csv")
        if not os.path.isfile(eng_csv):
            rows.append(f"| {label} | {placement}/{soil} | — | — | — | csv missing |")
            continue

        eng_y, eng_x, eng_F = load_engine_csv(eng_csv)
        cw_field, cw_widths_m, cw_depths_m, cw_actual_hr = load_cw_slice(folder, COMPARE_HR)
        if cw_field is None:
            rows.append(f"| {label} | {placement}/{soil} | — | — | — | no CW output |")
            continue

        eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
        residual = eng_on_cw - cw_field

        max_abs_R = float(np.abs(residual).max())
        mean_abs_R = float(np.abs(residual).mean())
        rms_R = float(np.sqrt(np.mean(residual**2)))

        flat_idx = int(np.argmax(np.abs(residual)))
        di_max, wi_max = np.unravel_index(flat_idx, residual.shape)
        d_at_max = float(cw_depths_m[di_max])
        w_at_max = float(cw_widths_m[wi_max])

        print(f"  Run {label}: max|R|={max_abs_R:.3f}°F  mean|R|={mean_abs_R:.3f}°F  "
              f"RMS={rms_R:.3f}°F  @ depth={d_at_max:.2f}m, w={w_at_max:.2f}m")

        if do_plots:
            out_path = plot_comparison(
                label, placement, soil, fix_num,
                cw_field, eng_on_cw, residual,
                cw_depths_m, cw_widths_m,
            )
            print(f"    Plot: {out_path}")

        rows.append(
            f"| {label} | {placement}/{soil} | {max_abs_R:.3f} | {mean_abs_R:.3f} | "
            f"{rms_R:.3f} | depth={d_at_max:.2f}m, w={w_at_max:.2f}m |"
        )

    # Write markdown report
    fix_tag = "fix1" if fix_num == 1 else "fix2"
    fix_desc = (
        "soil_temp_F plumbed; is_submerged=False (concrete sides still air)"
        if fix_num == 1 else
        "soil_temp_F plumbed + is_submerged=True (concrete sides contact soil)"
    )
    md = [
        f"# STAGE3 Fix {fix_num} — Residual Summary",
        "",
        f"**Fix applied:** {fix_desc}  ",
        f"**Comparison time:** t=168 hr  ",
        "**Method:** bilinear interpolation of engine (21×13) onto CW (49×13) grid.",
        "",
        "## Residual Table",
        "",
    ] + rows

    if do_plots:
        md += [
            "",
            "## Figures",
            "",
        ]
        for label, *_ in RUNS:
            md.append(f"- `plots/STAGE3_run{label}_engine_vs_cw.png`")

    out_md = os.path.join(HERE, f"STAGE3_{fix_tag}_residuals.md")
    with open(out_md, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"\nWritten: {out_md}")
    return 0


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--fix", type=int, choices=[1, 2], required=True,
                   help="Which fix stage to evaluate (1 or 2)")
    args = p.parse_args()
    sys.exit(run(args.fix))
