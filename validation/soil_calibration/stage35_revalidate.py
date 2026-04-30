#!/usr/bin/env python3
"""Stage 3.5 — Masked-region re-validation against the correct gate.

Builds a Run-A-derived validity mask (cells where |residual_RunA| < 0.3°F)
that identifies cells dominated by soil-coupling physics (vs. top-surface
and corner artifacts). Applies mask to all 9 Stage 3 Fix-2 runs at t=168 hr
and recomputes residual stats over the masked region only.

No engine code changed. Reads stage3_fix2_runs/ CSVs + CW output.txt files.

Outputs:
    validation/soil_calibration/STAGE35_validity_mask.npy
    validation/soil_calibration/plots/STAGE35_mask.png
    validation/soil_calibration/STAGE35_masked_residuals.md
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

# Reuse IO + resampling helpers from Stage 3 — do not reimplement
from stage3_compare import (
    load_engine_csv,
    load_cw_slice,
    resample_engine_to_cw,
    COMPARE_HR,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ENG_DIR = os.path.join(HERE, "stage3_fix2_runs")
PLOTS = os.path.join(HERE, "plots")
os.makedirs(PLOTS, exist_ok=True)

MASK_THRESHOLD_F = 0.3   # do NOT tune — specified in Stage 3.5 brief

# RUNS: (label, cw_folder, placement_F, soil_F)
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

# ΔT = soil - placement (used in bottom-corner ΔT-scaling check)
DELTA_T = {lbl: (soil - pl) for lbl, _, pl, soil in RUNS}


def load_run(label, cw_folder):
    """Return (residual_2d, cw_depths_m, cw_widths_m) for one run at t=168 hr."""
    eng_csv = os.path.join(ENG_DIR, f"run{label}_t168.csv")
    eng_y, eng_x, eng_F = load_engine_csv(eng_csv)
    cw_field, cw_widths_m, cw_depths_m, _ = load_cw_slice(cw_folder, COMPARE_HR)
    if cw_field is None:
        raise RuntimeError(f"CW data missing for run {label} / folder {cw_folder}")
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    return eng_on_cw - cw_field, cw_depths_m, cw_widths_m


def build_mask(residual_A):
    """Build boolean validity mask from Run A residual field."""
    mask = np.abs(residual_A) < MASK_THRESHOLD_F
    pct = 100.0 * mask.sum() / mask.size
    n_true = int(mask.sum())
    n_false = int(mask.size - n_true)
    print(f"  Mask: {n_true} True / {n_false} False  ({pct:.1f}% kept)")
    if pct < 30.0:
        raise RuntimeError(
            f"Mask keeps only {pct:.1f}% of cells — threshold may be too tight. "
            "Do not tune; stop and report."
        )
    if pct > 99.0:
        raise RuntimeError(
            f"Mask keeps {pct:.1f}% of cells — threshold may be too loose. "
            "Do not tune; stop and report."
        )
    return mask, pct, n_true, n_false


def plot_mask(mask, cw_depths_m, cw_widths_m, pct, n_true, n_false):
    """Save binary mask heatmap to plots/STAGE35_mask.png."""
    # Mirror column convention of plot_comparison in stage3_compare.py:
    # x-axis = width (0=edge, max=CL); y-axis = depth increasing downward.
    # stage3_compare uses imshow with extent=[w_min, w_max, d_max, 0.0] and
    # flips columns (::-1) so that column 0 (w=0, outer edge) is on the left.
    w_min = float(cw_widths_m[-1])   # smallest width (CL side)
    w_max = float(cw_widths_m[0])    # largest width (outer edge)
    d_max = float(cw_depths_m[-1])
    extent = [w_min, w_max, d_max, 0.0]

    mask_plot = mask.astype(float)[:, ::-1]   # flip columns: edge on left

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(
        mask_plot, aspect="equal", origin="upper",
        extent=extent, cmap="RdYlGn", vmin=0.0, vmax=1.0,
        interpolation="nearest",
    )
    ax.set_title(
        f"Mask covers {pct:.1f}% of grid; "
        "excludes top-surface-influenced and corner regions",
        fontsize=9, fontweight="bold",
    )
    ax.set_xlabel("Width (m, 0=edge, max=CL)", fontsize=8)
    ax.set_ylabel("Depth (m)", fontsize=8)
    ax.tick_params(labelsize=7)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04,
                 label="1=kept (True), 0=excluded (False)").ax.tick_params(labelsize=7)
    ax.text(
        0.02, 0.98,
        f"True (kept): {n_true}   False (excluded): {n_false}\n"
        f"Threshold: |residual_RunA| < {MASK_THRESHOLD_F}°F",
        transform=ax.transAxes, fontsize=7, va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8),
    )
    fig.tight_layout()
    out = os.path.join(PLOTS, "STAGE35_mask.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Mask plot: {out}")
    return out


def masked_stats(label, residual, mask, cw_depths_m, cw_widths_m):
    """Compute masked and unmasked residual stats for one run."""
    abs_R = np.abs(residual)
    unmasked_max = float(abs_R.max())

    masked_abs = abs_R[mask]
    masked_max = float(masked_abs.max())
    masked_mean = float(masked_abs.mean())

    # Location of masked max in grid coordinates
    masked_full = np.where(mask, abs_R, -np.inf)
    flat_idx = int(np.argmax(masked_full))
    di, wi = np.unravel_index(flat_idx, residual.shape)
    width_ft = float(cw_widths_m[wi]) * 3.28084
    depth_ft = float(cw_depths_m[di]) * 3.28084

    return {
        "unmasked_max": unmasked_max,
        "masked_max": masked_max,
        "masked_mean": masked_mean,
        "max_width_ft": width_ft,
        "max_depth_ft": depth_ft,
        "max_depth_m": float(cw_depths_m[di]),
        "max_width_m": float(cw_widths_m[wi]),
        "di": di,
        "wi": wi,
    }


def bottom_corner_disposition(mask, run_data, cw_depths_m, cw_widths_m):
    """Analyse whether the bottom-corner artifact survives the mask."""
    bottom_row = len(cw_depths_m) - 1  # depth ≈ 24.38 m
    bottom_mask = mask[bottom_row, :]
    n_bottom_kept = int(bottom_mask.sum())

    lines = ["## Bottom-corner artifact disposition", ""]

    if n_bottom_kept == 0:
        lines += [
            f"All {len(cw_widths_m)} bottom-row cells (depth ≈ "
            f"{cw_depths_m[bottom_row]:.2f} m) are **excluded by the Run-A mask**.",
            "",
            "**Disposition: scope-out.** The bottom-row artifact does not appear in "
            "the soil-coupling validation region. The Dirichlet-pinning / corner "
            "extrapolation noise is confined to cells where Run A already disagreed "
            "with CW (|residual_A| ≥ 0.3°F). No Stage 4 action required for this "
            "artifact.",
        ]
        return "\n".join(lines)

    # Some bottom cells survived — dig deeper
    lines += [
        f"{n_bottom_kept} of {len(cw_widths_m)} bottom-row cells "
        f"(depth ≈ {cw_depths_m[bottom_row]:.2f} m) survive the mask.",
        "",
        "### Bottom-row max|R| per run (masked cells only)",
        "",
        "| Run | ΔT (soil−pl, °F) | Bottom-row max\\|R\\| (°F) | Location (width_ft) |",
        "| --- | ----------------- | ------------------------- | ------------------- |",
    ]

    bottom_maxes = {}
    for label, cw_folder, placement, soil in RUNS:
        if label == "A":
            continue
        residual = run_data[label]["residual"]
        abs_R = np.abs(residual)
        bottom_abs = np.where(bottom_mask, abs_R[bottom_row, :], -np.inf)
        bmax = float(bottom_abs.max()) if n_bottom_kept > 0 else 0.0
        bmax_wi = int(np.argmax(bottom_abs)) if n_bottom_kept > 0 else 0
        bmax_w_ft = float(cw_widths_m[bmax_wi]) * 3.28084
        dt = DELTA_T[label]
        lines.append(
            f"| {label} | {dt:+d} | {bmax:.3f} | {bmax_w_ft:.2f} |"
        )
        bottom_maxes[label] = bmax

    # ΔT-scaling check: compare paired runs
    pairs = [
        ("B", "D", 13,  "B (soil=60,pl=73,ΔT=−13) vs D (soil=73,pl=60,ΔT=+13)"),
        ("C", "E", 17,  "C (soil=90,pl=73,ΔT=+17) vs E (soil=73,pl=90,ΔT=−17)"),
        ("F", "H", 28,  "F (soil=45,pl=73,ΔT=−28) vs H (soil=73,pl=45,ΔT=+28)"),
        ("G", "I", 27,  "G (soil=100,pl=73,ΔT=+27) vs I (soil=73,pl=100,ΔT=−27)"),
    ]
    lines += [
        "",
        "### ΔT-scaling check",
        "",
        "Paired runs share the same |ΔT|; if the bottom-row residual tracks |ΔT| "
        "monotonically the artifact is soil-driven.",
        "",
    ]
    scale_ok = True
    prev_abs_dt, prev_max = None, None
    ordered = [(13, "B/D"), (17, "C/E"), (27, "G/I"), (28, "F/H")]
    for abs_dt, pair_lbl in ordered:
        lbl1, lbl2 = pair_lbl.split("/")
        avg = (bottom_maxes.get(lbl1, 0.0) + bottom_maxes.get(lbl2, 0.0)) / 2.0
        lines.append(f"- |ΔT|={abs_dt}°F ({pair_lbl}): avg bottom-row max|R|={avg:.3f}°F")
        if prev_max is not None and avg < prev_max - 0.01:
            scale_ok = False
        prev_max = avg

    lines += ["", f"**Monotonically scales with |ΔT|: {'yes' if scale_ok else 'no'}**", ""]

    # Symmetry check: outer column (j=0) vs centerline (j=-1)
    # cw_widths_m[0] = max width = CL (6.096 m); cw_widths_m[-1] = min width = outer edge (0 m)
    lines += [
        "### Symmetry check (CL w≈6.10m vs outer edge w≈0m at bottom row)",
        "",
        "Note: in the CW/engine arrays, column 0 = CL (max width = 6.10 m); column -1 = outer edge (w ≈ 0 m).",
        "",
    ]
    for label, cw_folder, placement, soil in RUNS:
        if label == "A":
            continue
        residual = run_data[label]["residual"]
        cl_val = float(np.abs(residual[bottom_row, 0]))       # CL: largest width
        outer_val = float(np.abs(residual[bottom_row, -1]))   # outer edge: w≈0
        lines.append(
            f"- Run {label}: CL (w={cw_widths_m[0]:.2f}m) |R|={cl_val:.3f}°F, "
            f"outer edge (w={cw_widths_m[-1]:.2f}m) |R|={outer_val:.3f}°F"
        )

    lines += [
        "",
        "**Note:** Do not attempt a fix in this session. "
        "If the artifact is confirmed soil-coupled (large magnitude, scales with |ΔT|), "
        "flag as a Stage 4 task.",
    ]

    return "\n".join(lines)


def main():
    print("Stage 3.5 — Building Run-A validity mask ...")
    residual_A, cw_depths_m, cw_widths_m = load_run("A", "runA_baseline")
    mask, pct, n_true, n_false = build_mask(residual_A)

    np.save(os.path.join(HERE, "STAGE35_validity_mask.npy"), mask)
    print(f"  Saved: STAGE35_validity_mask.npy  shape={mask.shape}  dtype={mask.dtype}")

    plot_mask(mask, cw_depths_m, cw_widths_m, pct, n_true, n_false)

    # Sanity: Run A masked max must be < threshold by construction
    run_A_masked_max = float(np.abs(residual_A[mask]).max())
    assert run_A_masked_max < MASK_THRESHOLD_F, (
        f"BUG: Run A masked max = {run_A_masked_max:.4f}°F ≥ {MASK_THRESHOLD_F}°F"
    )

    print("\nComputing masked residuals for all 9 runs ...")
    stats = {}
    run_data = {}
    for label, cw_folder, placement, soil in RUNS:
        print(f"  Run {label} (placement={placement}°F, soil={soil}°F) ...")
        if label == "A":
            residual = residual_A
        else:
            residual, _, _ = load_run(label, cw_folder)
        s = masked_stats(label, residual, mask, cw_depths_m, cw_widths_m)
        stats[label] = s
        run_data[label] = {"residual": residual, "placement": placement, "soil": soil}
        print(
            f"    unmasked max|R|={s['unmasked_max']:.3f}°F  "
            f"masked max|R|={s['masked_max']:.3f}°F  "
            f"masked mean|R|={s['masked_mean']:.3f}°F  "
            f"@ ({s['max_width_ft']:.1f}ft, {s['max_depth_ft']:.1f}ft)"
        )

    # Build markdown report
    lines = [
        "# STAGE 3.5 — Masked-Region Residual Report",
        "",
        "**Mask definition:** cells where |residual_RunA| < 0.3°F at t=168 hr.  ",
        f"**Mask coverage:** {pct:.1f}% of grid ({n_true} of {mask.size} cells kept).  ",
        "**Engine source:** `stage3_fix2_runs/run{A..I}_t168.csv`  ",
        "**Comparison time:** t=168 hr  ",
        "**Method:** bilinear interpolation of engine (21×13) onto CW (49×13) grid.",
        "",
        "## Residual Table",
        "",
        "| Run | Placement/Soil (°F) | Unmasked max\\|R\\| (°F) | "
        "Masked max\\|R\\| (°F) | Masked mean\\|R\\| (°F) | "
        "Masked max location (width_ft, depth_ft) |",
        "| --- | ------------------- | ----------------------- | "
        "--------------------- | ---------------------- | "
        "--------------------------------------- |",
    ]
    for label, _, placement, soil in RUNS:
        s = stats[label]
        lines.append(
            f"| {label} | {placement}/{soil} | "
            f"{s['unmasked_max']:.3f} | "
            f"{s['masked_max']:.3f} | "
            f"{s['masked_mean']:.3f} | "
            f"({s['max_width_ft']:.2f}, {s['max_depth_ft']:.2f}) |"
        )

    # Gate verdict
    gate_a_pass = stats["A"]["masked_max"] <= 0.5
    non_a_max = max(stats[lbl]["masked_max"] for lbl in "BCDEFGHI")
    worst_run = max("BCDEFGHI", key=lambda lbl: stats[lbl]["masked_max"])
    gate_bi_pass = non_a_max <= 2.0
    overall_pass = gate_a_pass and gate_bi_pass

    lines += [
        "",
        "## Gate Verdict (masked region)",
        "",
        f"- **Gate A:** Run A masked max\\|R\\| = {stats['A']['masked_max']:.3f}°F ≤ 0.5°F → "
        f"{'**PASS**' if gate_a_pass else '**FAIL**'}",
        f"- **Gate B–I:** worst-case Run {worst_run} masked max\\|R\\| = {non_a_max:.3f}°F ≤ 2.0°F → "
        f"{'**PASS**' if gate_bi_pass else '**FAIL**'}",
        f"- **Overall:** {'**PASS** — Stage 3 soil-coupling validated in the correct gate region.' if overall_pass else '**FAIL**'}",
        "",
        f"Stage 4 calibration needed: **{'no' if overall_pass else 'yes'}**",
    ]

    # Bottom-corner disposition
    lines += [
        "",
        bottom_corner_disposition(mask, run_data, cw_depths_m, cw_widths_m),
    ]

    out_md = os.path.join(HERE, "STAGE35_masked_residuals.md")
    with open(out_md, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\nWritten: {out_md}")

    # Summary printout for the deliverable
    print("\n" + "=" * 60)
    print("DELIVERABLE SUMMARY")
    print("=" * 60)
    print(f"Run A masked max|R|:         {stats['A']['masked_max']:.3f}°F  (was {stats['A']['unmasked_max']:.3f}°F unmasked)")
    print(f"Worst B-I masked max|R|:     Run {worst_run}: {non_a_max:.3f}°F")
    print(f"Gate A (≤0.5°F):             {'PASS' if gate_a_pass else 'FAIL'}")
    print(f"Gate B-I (≤2.0°F):           {'PASS' if gate_bi_pass else 'FAIL'}")
    print(f"Overall gate:                {'PASS' if overall_pass else 'FAIL'}")
    print(f"Stage 4 needed:              {'no' if overall_pass else 'yes'}")
    print(f"Mask coverage:               {pct:.1f}%  ({n_true} true / {n_false} false)")


if __name__ == "__main__":
    main()
