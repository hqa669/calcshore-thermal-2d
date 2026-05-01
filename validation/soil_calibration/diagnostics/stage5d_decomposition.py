#!/usr/bin/env python3
"""Stage 5d-prep: Residual decomposition diagnostic.

Decomposes Run B–I residuals by subtracting Run A bulk baseline.
Tests linear magnitude scaling and sign symmetry across all 4 matched pairs.
Generates heatmap and wi=9 profile figures.  No engine re-runs for t=168;
Run A re-run only for t=24/84 integrated-heat trajectory (§2.6c).

No engine source changes, no commits.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SC   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, SC)

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import load_engine_csv, resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx, RUNS

ENG_DIR  = os.path.join(SC, "stage5c_runs")
CW_RUNS  = os.path.join(SC, "cw_runs")
MASK_PATH = os.path.join(SC, "STAGE35_validity_mask.npy")
OUT_DIR  = os.path.join(SC, "diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

COMPARE_HR = 168

# Matched pairs: (label_cool, label_warm, |ΔT|)
PAIRS = [
    ("B", "D", 13),
    ("E", "C", 17),
    ("I", "G", 27),
    ("F", "H", 28),
]

# ΔT (placement - soil) per run: positive = cooling wave
DT_MAP = {lbl: pl - so for lbl, _, pl, so in RUNS}


# ---------------------------------------------------------------------------
# Step 1: Load all residual fields at t=168
# ---------------------------------------------------------------------------
def load_residual_field(label, folder, cw_depths_m, cw_widths_m):
    csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
    eng_y, eng_x, eng_F = load_engine_csv(csv_path)
    cw_field, _, _, _ = _load_cw(folder)
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    return eng_on_cw - cw_field   # (49, 13) Engine - CW


def _load_cw(folder):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m, float(v.time_hrs[ti])


# ---------------------------------------------------------------------------
# Step 2: §2.2 magnitude scaling table
# ---------------------------------------------------------------------------
def magnitude_table(residuals, run_a_residual, mask):
    print("\n=== §2.2 Magnitude Scaling Table ===")
    header = (f"{'Run':>4}  {'|ΔT|':>6}  {'sign':>5}  "
              f"{'max|R_X|':>10}  {'R_A@same':>10}  "
              f"{'max|R_dec|':>11}  {'loc(di,wi)':>12}  {'scale':>8}")
    print(header)
    print("-" * len(header))
    scales = []
    for lbl, folder, pl, so in RUNS:
        if lbl == "A":
            continue
        R_X   = residuals[lbl]
        R_A   = run_a_residual
        R_dec = np.where(mask, R_X - R_A, np.nan)
        dT    = pl - so
        idx   = np.unravel_index(np.nanargmax(np.abs(R_dec)), R_dec.shape)
        maxR  = float(np.nanmax(np.abs(R_dec)))
        rA_at = float(R_A[idx]) if not np.isnan(R_A[idx]) else float("nan")
        rX_at = float(np.nanmax(np.abs(np.where(mask, R_X, np.nan))))
        sign  = "cool" if dT > 0 else "warm"
        scale = maxR / abs(dT) if abs(dT) > 0 else float("nan")
        scales.append(scale)
        print(f"{lbl:>4}  {abs(dT):>6}  {sign:>5}  "
              f"{rX_at:>10.4f}  {rA_at:>10.4f}  "
              f"{maxR:>11.4f}  ({idx[0]:>2},{idx[1]:>2})         {scale:>8.5f}")
    scales = np.array(scales)
    print(f"\n  scale_factor: mean={scales.mean():.5f}  std={scales.std():.5f}  "
          f"cv={scales.std()/scales.mean():.3f}  "
          f"({'CLEAN (<0.10)' if scales.std()/scales.mean() < 0.10 else 'NON-LINEAR (>0.10)'})")
    return scales


# ---------------------------------------------------------------------------
# Step 3: §2.3 sign symmetry check
# ---------------------------------------------------------------------------
def symmetry_table(residuals, run_a_residual, mask, cw_widths_m):
    print("\n=== §2.3 Sign Symmetry Table ===")
    print(f"{'Pair':>8}  {'|ΔT|':>5}  {'max|R_cool_dec|':>16}  "
          f"{'max|R_warm_dec|':>16}  {'asymmetry':>10}  {'verdict':>15}")
    print("-" * 80)
    asym_vals = []
    for lbl_c, lbl_w, dT in PAIRS:
        _, folder_c, pl_c, so_c = next(r for r in RUNS if r[0] == lbl_c)
        _, folder_w, pl_w, so_w = next(r for r in RUNS if r[0] == lbl_w)
        R_cool = np.where(mask, residuals[lbl_c] - run_a_residual, np.nan)
        R_warm = np.where(mask, residuals[lbl_w] - run_a_residual, np.nan)
        asym   = float(np.nanmax(np.abs(R_cool + R_warm)))
        asym_vals.append(asym)
        verdict = "CLEAN (<0.05)" if asym < 0.05 else ("BORDERLINE" if asym < 0.10 else "ASYMMETRIC (>0.10)")
        print(f"{lbl_c+'/'+lbl_w:>8}  {dT:>5}  {float(np.nanmax(np.abs(R_cool))):>16.4f}  "
              f"{float(np.nanmax(np.abs(R_warm))):>16.4f}  {asym:>10.4f}  {verdict:>15}")
    print(f"\n  Max asymmetry across pairs: {max(asym_vals):.4f}°F")
    return asym_vals


# ---------------------------------------------------------------------------
# Step 4: §2.4 4×3 decomposed residual heatmaps
# ---------------------------------------------------------------------------
def plot_heatmaps(residuals, run_a_residual, mask, cw_widths_m, cw_depths_m, out_path):
    extent = [cw_widths_m[-1] - 0.25, cw_widths_m[0] + 0.25,
              cw_depths_m[-1] + 0.25, cw_depths_m[0] - 0.25]

    fig, axes = plt.subplots(4, 3, figsize=(15, 18))
    fig.suptitle(
        "Stage 5d-prep: Decomposed Residual Fields (R_X − R_A) at t=168 hr\n"
        "Columns: −R_cool_dec (sign-flipped)  |  R_warm_dec  |  Pair sum (asymmetry)\n"
        "blanket=0.0m, model_soil=False, is_submerged=True",
        fontsize=10, fontweight="bold",
    )

    col_titles = ["−R_cool_decomposed", "R_warm_decomposed", "Sum (asymmetry)"]
    for ci, title in enumerate(col_titles):
        axes[0, ci].set_title(title, fontsize=9, fontweight="bold")

    for row_idx, (lbl_c, lbl_w, dT) in enumerate(PAIRS):
        R_cool_dec = np.where(mask, residuals[lbl_c] - run_a_residual, np.nan)
        R_warm_dec = np.where(mask, residuals[lbl_w] - run_a_residual, np.nan)
        pair_sum   = R_cool_dec + R_warm_dec  # would be 0 if perfect mirrors

        # Shared scale for cool/warm columns
        sym_scale = max(
            float(np.nanmax(np.abs(R_cool_dec))),
            float(np.nanmax(np.abs(R_warm_dec))),
        ) * 1.1

        data_cols = [-R_cool_dec, R_warm_dec, pair_sum]
        vmins     = [-sym_scale, -sym_scale, -0.08]
        vmaxs     = [ sym_scale,  sym_scale,  0.08]
        subtitles = [
            f"Run {lbl_c} (|ΔT|={dT}°F, cool, negated)",
            f"Run {lbl_w} (|ΔT|={dT}°F, warm)",
            f"Sum |ΔT|={dT}°F pair",
        ]

        for ci, (data, vmin, vmax, subtitle) in enumerate(
                zip(data_cols, vmins, vmaxs, subtitles)):
            ax = axes[row_idx, ci]
            im = ax.imshow(
                data, aspect="auto", cmap="RdBu_r",
                vmin=vmin, vmax=vmax, extent=extent,
            )
            plt.colorbar(im, ax=ax, label="°F", fraction=0.046, pad=0.04)
            ax.set_title(subtitle, fontsize=8)
            ax.set_xlabel("Width from edge (m)", fontsize=7)
            ax.set_ylabel("Depth (m)", fontsize=7)
            ax.tick_params(labelsize=6)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nHeatmap figure saved: {out_path}")


# ---------------------------------------------------------------------------
# Step 5: §2.5 wi=9 decomposed profile (all 9 runs, di=4..48)
# ---------------------------------------------------------------------------
def plot_wi9_profiles(residuals, run_a_residual, mask, cw_depths_m, out_path):
    wi9 = 9
    di_range = range(4, 49)
    depths = [cw_depths_m[di] for di in di_range]

    # Build label/color/style map
    cool_runs  = {"B": "#1f77b4", "E": "#17becf", "F": "#0000cc", "I": "#7f7f7f"}
    warm_runs  = {"C": "#d62728", "D": "#e377c2", "G": "#8c564b", "H": "#ff7f0e"}
    run_order  = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.axhline(0, color="black", lw=0.8, zorder=1)
    ax.axhspan(-0.35, 0.35, color="green", alpha=0.06, zorder=0, label="±0.35°F gate")
    ax.axhline( 0.35, color="green", lw=1.0, ls="--", zorder=1)
    ax.axhline(-0.35, color="green", lw=1.0, ls="--", zorder=1)

    for lbl in run_order:
        R_A_dec = np.zeros(len(di_range))
        if lbl == "A":
            vals = R_A_dec  # zero by construction
            color, ls, lw, zord = "black", "-", 1.5, 5
            label = "Run A (baseline, = 0 by construction)"
        else:
            R_X = residuals[lbl]
            R_dec = R_X - run_a_residual  # full (49,13)
            vals = np.array([
                R_dec[di, wi9] if mask[di, wi9] else np.nan
                for di in di_range
            ])
            dT = DT_MAP[lbl]
            if lbl in cool_runs:
                color, ls, lw, zord = cool_runs[lbl], "-", 1.5, 3
            else:
                color, ls, lw, zord = warm_runs[lbl], "--", 1.5, 3
            label = f"Run {lbl} (|ΔT|={abs(dT)}°F, {'cool' if dT>0 else 'warm'})"

        ax.plot(vals, depths, color=color, ls=ls, lw=lw, zorder=zord, label=label)

    ax.set_xlabel("Decomposed residual (Engine − CW) − R_A (°F)", fontsize=10)
    ax.set_ylabel("Depth (m)", fontsize=10)
    ax.set_ylim(cw_depths_m[-1] + 0.5, cw_depths_m[4] - 0.5)
    ax.set_title(
        "Stage 5d-prep: Decomposed Residual at wi=9 (1.52 m from edge), all runs, t=168 hr\n"
        "R_X_decomposed = R_X − R_A  (bulk baseline removed)",
        fontsize=10, fontweight="bold",
    )
    ax.legend(fontsize=8, loc="lower right")
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wi=9 profile figure saved: {out_path}")


# ---------------------------------------------------------------------------
# Step 6: §2.6 Run A bulk characterization
# ---------------------------------------------------------------------------
def run_a_bulk(run_a_residual, cw_widths_m, cw_depths_m, mask):
    print("\n=== §2.6 Run A Bulk Characterization ===")

    # (a) Spatial profile at di=24
    di24 = 24
    print(f"\n(a) R_A at di=24 (depth={cw_depths_m[di24]:.2f} m):")
    print(f"    {'wi':>3}  {'x(m)':>6}  {'R_A(°F)':>10}")
    for wi in range(13):
        print(f"    {wi:>3}  {cw_widths_m[wi]:>6.3f}  {run_a_residual[di24, wi]:>+10.4f}")

    # (b) We need the CW and engine peak T at t=168 — load from CSV and CW output
    csv_a = os.path.join(ENG_DIR, "runA_t168.csv")
    eng_y, eng_x, eng_F = load_engine_csv(csv_a)
    cw_f168, _, _, _ = _load_cw("runA_baseline")
    print(f"\n(b) Peak temperature at t=168 hr:")
    print(f"    Engine max T = {eng_F.max():.4f}°F  at concrete grid node")
    print(f"    CW     max T = {cw_f168.max():.4f}°F")
    print(f"    Engine−CW peak ΔT = {eng_F.max() - cw_f168.max():+.4f}°F")

    # (c) Integrated heat trajectory — need t=24 and t=84; re-run Run A
    print("\n(c) Integrated heat trajectory (re-running Run A to capture t=24/84/168)...")
    _run_a_heat_trajectory(cw_depths_m, cw_widths_m)


def _run_a_heat_trajectory(cw_depths_m, cw_widths_m):
    mix, geom, constr, _ = parse_cw_dat(os.path.join(CW_RUNS, "runA_baseline", "input.dat"))
    factor, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = factor
    mix.Hu_J_kg_effective = mix.Hu_J_kg * factor
    constr.model_soil = False
    constr.is_submerged = True
    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0_C = (73.0 - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    env = make_neutral_env(73.0)
    print("    Running Run A (for heat trajectory)...")
    result = solve_hydration_2d(
        grid, mix, T_initial,
        duration_s=168 * 3600, output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=env, construction=constr,
        T_ground_deep_C=T0_C, diagnostic_outputs=False,
    )
    jslice, islice = grid.concrete_slice()
    v = parse_cw_temp_output(os.path.join(CW_RUNS, "runA_baseline", "output.txt"))
    print(f"\n    {'t(hr)':>6}  {'Eng sum(T)':>12}  {'CW sum(T)':>11}  {'Δ sum':>10}  "
          f"{'Eng maxT':>10}  {'CW maxT':>10}  {'Δ maxT':>9}")
    for hr in [24, 84, 168]:
        ti_e = nearest_time_idx(result.t_s, float(hr))
        ti_c = int(np.abs(v.time_hrs - hr).argmin())
        eng_F_hr = result.T_field_C[ti_e, jslice, islice] * 9.0 / 5.0 + 32.0
        cw_F_hr  = v.T_field_F[ti_c]
        eng_on_cw = resample_engine_to_cw(
            grid.y[jslice], grid.x[islice], eng_F_hr, cw_depths_m, cw_widths_m
        )
        eng_sum = float(eng_on_cw.sum())
        cw_sum  = float(cw_F_hr.sum())
        print(f"    {hr:>6}  {eng_sum:>12.2f}  {cw_sum:>11.2f}  "
              f"{eng_sum-cw_sum:>+10.2f}  "
              f"{float(eng_on_cw.max()):>10.4f}  {float(cw_F_hr.max()):>10.4f}  "
              f"{float(eng_on_cw.max()-cw_F_hr.max()):>+9.4f}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== Stage 5d-prep: Residual Decomposition Diagnostic ===\n")

    # Reference CW geometry (identical across all runs)
    _, cw_widths_m, cw_depths_m, _ = _load_cw("runA_baseline")
    mask = np.load(MASK_PATH)

    # Load all residual fields
    print("Loading residual fields from cached CSVs...")
    residuals = {}
    for lbl, folder, pl, so in RUNS:
        residuals[lbl] = load_residual_field(lbl, folder, cw_depths_m, cw_widths_m)
        print(f"  Run {lbl}: max|R|={float(np.nanmax(np.abs(np.where(mask, residuals[lbl], np.nan)))):.4f}°F")
    run_a_residual = np.where(mask, residuals["A"], np.nan)

    # §2.2 magnitude table
    scale_factors = magnitude_table(residuals, run_a_residual, mask)

    # §2.3 symmetry table
    asym_vals = symmetry_table(residuals, run_a_residual, mask, cw_widths_m)

    # §2.4 heatmap
    plot_heatmaps(
        residuals, run_a_residual, mask, cw_widths_m, cw_depths_m,
        os.path.join(OUT_DIR, "stage5d_prep_decomposed_residuals.png"),
    )

    # §2.5 wi=9 profile
    plot_wi9_profiles(
        residuals, run_a_residual, mask, cw_depths_m,
        os.path.join(OUT_DIR, "stage5d_prep_wi9_decomposed_profiles.png"),
    )

    # §2.6 bulk characterization
    run_a_bulk(run_a_residual, cw_widths_m, cw_depths_m, mask)


if __name__ == "__main__":
    main()
