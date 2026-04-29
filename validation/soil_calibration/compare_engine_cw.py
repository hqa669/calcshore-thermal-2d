#!/usr/bin/env python3
"""S2.3 — Engine-vs-CW comparison plots and residual summary.

Reads engine slice CSVs from engine_runs/ and CW output.txt files.
Resamples engine (21×13 concrete grid) onto CW's (13×49) grid using bilinear
interpolation. Generates 3-row comparison figures and populates summary table.

Usage:
    python validation/soil_calibration/compare_engine_cw.py

Writes:
    validation/soil_calibration/plots/run{A..I}_engine_vs_cw.png
    validation/soil_calibration/STAGE2_engine_cw_comparison.md
"""
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
ENGINE_RUNS = os.path.join(HERE, "engine_runs")
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
    """Load engine CSV written by run_engine_all.py.

    Returns (y_m, x_m, field_F) where field_F.shape = (ny, nx).
    """
    with open(path) as f:
        lines = f.readlines()

    header = lines[0].strip().split(",")
    # header[0] = "depth_m\\width_m", rest = x coords
    x_m = np.array([float(h) for h in header[1:]])

    y_m_list = []
    rows = []
    for line in lines[1:]:
        parts = line.strip().split(",")
        y_m_list.append(float(parts[0]))
        rows.append([float(v) for v in parts[1:]])

    y_m = np.array(y_m_list)
    field_F = np.array(rows)   # (ny, nx)
    return y_m, x_m, field_F


def load_cw_slice(folder, target_hr):
    """Load CW output and return slice nearest target_hr."""
    path = os.path.join(CW_RUNS, folder, "output.txt")
    if not os.path.isfile(path):
        return None, None, None, None
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - target_hr).argmin())
    actual_hr = float(v.time_hrs[ti])
    return v.T_field_F[ti], v.widths_m, v.depths_m, actual_hr


def resample_engine_to_cw(eng_y_m, eng_x_m, eng_field_F, cw_depths_m, cw_widths_m):
    """Bilinearly interpolate engine field onto CW grid.

    Engine field: (ny_conc, nx_conc) with y ascending, x ascending
    CW grid: depths_m ascending (0→24.38), widths_m descending (6.1→0)

    Returns resampled_F with shape (nD, nW) matching CW convention.
    """
    # Build interpolator: axes must be strictly increasing
    # eng_y_m: ascending, eng_x_m: ascending
    interp = RegularGridInterpolator(
        (eng_y_m, eng_x_m),
        eng_field_F,
        method="linear",
        bounds_error=False,
        fill_value=None,   # extrapolate at boundaries
    )

    # CW query points: depths_m ascending, widths_m descending
    # For each CW (depth, width) point, query the engine
    D, W = np.meshgrid(cw_depths_m, cw_widths_m, indexing="ij")  # (nD, nW)
    pts = np.column_stack([D.ravel(), W.ravel()])
    resampled_F = interp(pts).reshape(len(cw_depths_m), len(cw_widths_m))
    return resampled_F


def plot_comparison(label, placement, soil,
                    cw_field, eng_field, residual,
                    cw_depths_m, cw_widths_m, target_hr):
    """3-row comparison figure: CW | engine | residual."""
    delta_T = soil - placement
    w_min = float(cw_widths_m[-1])   # 0 m
    w_max = float(cw_widths_m[0])    # ~6.1 m
    d_max = float(cw_depths_m[-1])   # ~24.38 m
    extent = [w_min, w_max, d_max, 0.0]

    # For imshow: need to flip width axis (widths_m is descending)
    cw_plot = cw_field[:, ::-1]
    eng_plot = eng_field[:, ::-1]
    res_plot = residual[:, ::-1]

    sym = max(0.01, float(np.abs(residual).max()))

    fig, axes = plt.subplots(3, 1, figsize=(6, 14))
    fig.suptitle(
        f"Run {label}: placement={placement}°F, soil={soil}°F, ΔT={delta_T:+d}°F\n"
        f"t={target_hr:.0f} hr",
        fontsize=10, fontweight="bold",
    )

    panel_data = [
        (cw_plot,  "CW (reference)",       "inferno",  (VMIN, VMAX)),
        (eng_plot, "Engine (resampled)",    "inferno",  (VMIN, VMAX)),
        (res_plot, f"Engine−CW (max|R|={sym:.2f}°F)", "RdBu_r", (-sym, sym)),
    ]

    for ax, (field_data, title, cmap_name, (vmin, vmax)) in zip(axes, panel_data):
        im = ax.imshow(
            field_data, aspect="equal", origin="upper",
            extent=extent, cmap=cmap_name, vmin=vmin, vmax=vmax,
            interpolation="nearest",
        )
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("Width (m, 0=edge, max=CL)", fontsize=7)
        ax.set_ylabel("Depth (m)", fontsize=7)
        ax.tick_params(labelsize=6)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=6)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out_path = os.path.join(PLOTS, f"run{label}_engine_vs_cw.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run():
    # Load manifest to know which runs succeeded
    manifest_path = os.path.join(ENGINE_RUNS, "manifest.json")
    if not os.path.isfile(manifest_path):
        print("ERROR: engine_runs/manifest.json not found. Run run_engine_all.py first.")
        return 1

    with open(manifest_path) as f:
        manifest = json.load(f)

    grid_path = os.path.join(ENGINE_RUNS, "grid_info.json")
    with open(grid_path) as f:
        grid_info = json.load(f)
    print(f"Engine concrete grid: {grid_info['ny_concrete']}×{grid_info['nx_concrete']} nodes")

    summary_rows = []
    summary_rows.append(
        "| Run | max\\|R\\| (°F) | mean\\|R\\| (°F) | RMS R (°F) | Location of max |"
    )
    summary_rows.append("| --- | ------------ | ------------ | ---------- | --------------- |")

    for label, folder, placement, soil in RUNS:
        status = manifest.get(label, {}).get("status", "missing")
        if status != "ok":
            print(f"  SKIP Run {label}: status={status}")
            summary_rows.append(f"| {label} | — | — | — | {status} |")
            continue

        print(f"Processing Run {label} ({folder})…")

        # Load engine slice at t=168
        eng_csv = os.path.join(ENGINE_RUNS, f"run{label}_t168.csv")
        if not os.path.isfile(eng_csv):
            print(f"  SKIP: {eng_csv} not found")
            continue

        eng_y, eng_x, eng_F = load_engine_csv(eng_csv)
        # eng_F: (ny_conc, nx_conc), y ascending, x ascending (0→6.1 m)

        # Load CW slice at t=168
        cw_field, cw_widths_m, cw_depths_m, cw_actual_hr = load_cw_slice(folder, COMPARE_HR)
        if cw_field is None:
            print(f"  SKIP: CW output.txt not found for {folder}")
            continue

        print(f"  CW t={cw_actual_hr:.2f} hr, engine t=168 hr")

        # Resample engine → CW grid
        eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
        # eng_on_cw: (nD, nW) matching CW

        residual = eng_on_cw - cw_field  # engine − CW

        max_abs_R = float(np.abs(residual).max())
        mean_abs_R = float(np.abs(residual).mean())
        rms_R = float(np.sqrt(np.mean(residual**2)))

        flat_idx = int(np.argmax(np.abs(residual)))
        di_max, wi_max = np.unravel_index(flat_idx, residual.shape)
        d_at_max = float(cw_depths_m[di_max])
        w_at_max = float(cw_widths_m[wi_max])

        print(f"  max|R|={max_abs_R:.3f}°F  mean|R|={mean_abs_R:.3f}°F  "
              f"RMS={rms_R:.3f}°F  @ depth={d_at_max:.2f}m, width={w_at_max:.2f}m")

        # Plot
        out_path = plot_comparison(
            label, placement, soil,
            cw_field, eng_on_cw, residual,
            cw_depths_m, cw_widths_m, COMPARE_HR,
        )
        print(f"  Saved: {out_path}")

        summary_rows.append(
            f"| {label} | {max_abs_R:.3f} | {mean_abs_R:.3f} | {rms_R:.3f} "
            f"| depth={d_at_max:.2f}m, w={w_at_max:.2f}m |"
        )

    # Write MD report
    md = [
        "# STAGE2 — Engine vs CW Comparison",
        "",
        "## Method",
        "",
        "Engine temperature fields (concrete subgrid: "
        f"{grid_info['ny_concrete']}y × {grid_info['nx_concrete']}x nodes, "
        f"dx={grid_info['dx_m']:.4f} m) are bilinearly resampled onto CW's 49×13 grid "
        "using `scipy.interpolate.RegularGridInterpolator`. Engine output is in °C; "
        "converted to °F before residual computation.",
        "",
        "Engine runs use explicit `T_ground_deep_C = soil_temp_F` (from input.dat), bypassing",
        "the default `compute_T_gw_C(env)` path. This allows the comparison to isolate",
        "model-form differences in how soil is spatially represented, independent of the",
        "known soil-BC plumbing gap. See STAGE2_engine_soil_audit.md §Engine-vs-CW protocol.",
        "",
        "Comparison time: t=168 hr.",
        "",
        "**Interpolation note**: ~0.1°F noise is expected from grid-resolution mismatch.",
        "",
        "## Summary Table",
        "",
    ] + summary_rows + [
        "",
        "## Figures",
        "",
    ]

    for label, *_ in RUNS:
        md.append(f"- `plots/run{label}_engine_vs_cw.png` — CW | engine | residual at t=168 hr")

    out_path = os.path.join(HERE, "STAGE2_engine_cw_comparison.md")
    with open(out_path, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"\nWritten: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(run())
