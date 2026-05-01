#!/usr/bin/env python3
"""Stage 5c diagnostic — IC alignment test (Run I only).

Hypothesis: the wi=9 lateral residual hump is caused (or contaminated) by an
IC initialization mismatch: CW pre-applies the side Dirichlet at t=0, but the
engine's BC writes only fire at the end of step 1.  This script runs Run I
twice — baseline (unmodified IC) and experimental (T_soil pre-applied to the
concrete edge column and bottom face row before solve_hydration_2d) — and
compares the resulting residuals to classify into Bucket 1/2/3.

No engine source changes.  No commits.  Observational diagnostic only.
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
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
CW_RUNS  = os.path.join(SC, "cw_runs")
MASK_PATH = os.path.join(SC, "STAGE35_validity_mask.npy")
OUT_DIR  = os.path.join(SC, "diagnostics")
OUT_FILE = os.path.join(OUT_DIR, "stage5c_diagnostic_ic_alignment_runI.png")
os.makedirs(OUT_DIR, exist_ok=True)

FOLDER       = "runI_100_73"
PLACEMENT_F  = 100.0
SOIL_F       = 73.0
DI_MID       = 24
TARGET_HRS   = [0, 84, 168]
GATE_TOL_F   = 0.35
BASELINE_MAX_R_REF = 0.546   # Stage 5c gate run value; reproduce within 0.01°F

BLANKET_M    = 0.0
MODEL_SOIL   = False
IS_SUBMERGED = True


# ---------------------------------------------------------------------------
def run_engine(placement_F, soil_F, experimental_ic):
    dat_path = os.path.join(CW_RUNS, FOLDER, "input.dat")
    mix, geom, constr, _ = parse_cw_dat(dat_path)
    factor, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = factor
    mix.Hu_J_kg_effective = mix.Hu_J_kg * factor

    constr.model_soil = MODEL_SOIL
    constr.is_submerged = IS_SUBMERGED

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=IS_SUBMERGED,
        model_soil=MODEL_SOIL,
        blanket_thickness_m=BLANKET_M,
    )

    T0_C     = (placement_F - 32.0) * 5.0 / 9.0
    T_soil_C = (soil_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_soil] = T_soil_C   # no-op for model_soil=False

    if experimental_ic:
        # Pre-apply the same Dirichlet writes that thermal_engine_2d.py:2256-2259
        # normally applies at the end of step 1 (model_soil=False, is_submerged=True).
        T_initial[grid.iy_concrete_end, grid.ix_concrete_start:] = T_soil_C
        T_initial[grid.iy_concrete_start:grid.iy_concrete_end + 1, grid.ix_concrete_start] = T_soil_C

    label = "experimental" if experimental_ic else "baseline"
    env = make_neutral_env(placement_F)
    print(f"  Running Run I ({label}) ...")
    result = solve_hydration_2d(
        grid, mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=env,
        construction=constr,
        T_ground_deep_C=T_soil_C,
        diagnostic_outputs=False,
    )
    print(f"    Done: {result.n_inner_steps} steps, dt={result.dt_inner_s:.0f}s")

    jslice, islice = grid.concrete_slice()
    T_conc_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jslice]
    eng_x = grid.x[islice]

    snapshots = {}
    for hr in TARGET_HRS:
        ti = nearest_time_idx(result.t_s, float(hr))
        snapshots[hr] = (eng_y.copy(), eng_x.copy(), T_conc_F[ti].copy())
        t_actual = result.t_s[ti] / 3600.0
        print(f"    t={hr} hr → ti={ti}, t_actual={t_actual:.2f} hr, "
              f"Tmin={T_conc_F[ti].min():.2f}°F Tmax={T_conc_F[ti].max():.2f}°F")

    return grid, snapshots


# ---------------------------------------------------------------------------
def sanity_check(grid, snaps_base, snaps_exp, cw_by_hr, cw_widths_m):
    print("\n=== Sanity Check (t=0 edge column, §3.3) ===")
    # Concrete edge column in engine: ix = ix_concrete_start, all concrete rows
    # In the CW comparison grid: wi=12 is the edge column (widths_m[12]=0.0)
    ics = grid.ix_concrete_start
    iy_start = grid.iy_concrete_start
    iy_end   = grid.iy_concrete_end

    # (a) Baseline edge column at t=0 should be at placement temp (100°F)
    _, _, eng_base_F = snaps_base[0]
    edge_col_base = eng_base_F[iy_start - iy_start : iy_end - iy_start + 1, 0]
    base_min, base_max = edge_col_base.min(), edge_col_base.max()
    print(f"  Baseline edge col (ix_concrete_start) at t=0: "
          f"min={base_min:.2f}°F  max={base_max:.2f}°F  (expect ≈{PLACEMENT_F}°F)")
    base_ok = abs(base_min - PLACEMENT_F) < 0.1 and abs(base_max - PLACEMENT_F) < 0.1
    print(f"  → {'PASS' if base_ok else 'FAIL'}")

    # (b) Experimental edge column at t=0 should be at T_soil (73°F)
    _, _, eng_exp_F = snaps_exp[0]
    edge_col_exp = eng_exp_F[iy_start - iy_start : iy_end - iy_start + 1, 0]
    exp_min, exp_max = edge_col_exp.min(), edge_col_exp.max()
    print(f"  Experimental edge col at t=0: "
          f"min={exp_min:.2f}°F  max={exp_max:.2f}°F  (expect ≈{SOIL_F}°F)")
    exp_ok = abs(exp_min - SOIL_F) < 0.1 and abs(exp_max - SOIL_F) < 0.1
    print(f"  → {'PASS' if exp_ok else 'FAIL'}")

    # (c) CW edge column (wi=12) at t=0 should be at T_soil
    cw_t0 = cw_by_hr[0]
    cw_edge = cw_t0[:, 12]   # wi=12 = edge, all di
    cw_min, cw_max = float(cw_edge.min()), float(cw_edge.max())
    print(f"  CW wi=12 (edge) at t=0: min={cw_min:.2f}°F  max={cw_max:.2f}°F  (expect ≈{SOIL_F}°F)")
    cw_ok = abs(cw_min - SOIL_F) < 0.5 and abs(cw_max - SOIL_F) < 0.5
    print(f"  → {'PASS' if cw_ok else 'FAIL'}")

    all_ok = base_ok and exp_ok and cw_ok
    if not all_ok:
        raise RuntimeError("Sanity check FAILED — do not proceed to analysis.")
    print("  All sanity checks PASS.\n")


# ---------------------------------------------------------------------------
def compute_masked_max(eng_y, eng_x, eng_F, cw_field, cw_depths_m, cw_widths_m, mask):
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    residual = eng_on_cw - cw_field
    masked_r = np.where(mask, residual, np.nan)
    return float(np.nanmax(np.abs(masked_r)))


# ---------------------------------------------------------------------------
def print_comparison_table(snaps_base, snaps_exp, cw_by_hr, cw_depths_m, cw_widths_m):
    print("\n=== Comparison Table: di=24, wi=8/9/10, t=0/84/168 (§3.1) ===")
    print(f"{'t':>5}  {'wi':>3}  {'x(m)':>6}  {'CW':>7}  {'Base':>7}  {'Exp':>7}  "
          f"{'R_base':>8}  {'R_exp':>8}  {'Δ(exp-base)':>12}")
    print("-" * 80)
    for hr in TARGET_HRS:
        eng_y_b, eng_x_b, eng_F_b = snaps_base[hr]
        eng_y_e, eng_x_e, eng_F_e = snaps_exp[hr]
        cw_field = cw_by_hr[hr]
        b_on_cw = resample_engine_to_cw(eng_y_b, eng_x_b, eng_F_b, cw_depths_m, cw_widths_m)
        e_on_cw = resample_engine_to_cw(eng_y_e, eng_x_e, eng_F_e, cw_depths_m, cw_widths_m)
        for wi in [8, 9, 10]:
            cw_T   = cw_field[DI_MID, wi]
            base_T = b_on_cw[DI_MID, wi]
            exp_T  = e_on_cw[DI_MID, wi]
            R_b    = base_T - cw_T
            R_e    = exp_T  - cw_T
            delta  = R_e - R_b
            x_m    = cw_widths_m[wi]
            print(f"{hr:>5}  {wi:>3}  {x_m:>6.3f}  {cw_T:>7.3f}  {base_T:>7.3f}  "
                  f"{exp_T:>7.3f}  {R_b:>+8.3f}  {R_e:>+8.3f}  {delta:>+12.3f}")
        print()


# ---------------------------------------------------------------------------
def make_figure(snaps_base, snaps_exp, cw_by_hr, cw_widths_m, cw_depths_m):
    x_cw  = cw_widths_m[::-1]    # edge→CL order
    wi9_x = cw_widths_m[9]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(
        "Stage 5c Diagnostic: IC Alignment Test — Run I (|ΔT|=27°F)\n"
        "Residual (Engine − CW) at di=24 (mid-depth ≈12.19 m) — "
        "blanket_thickness_m=0.0, model_soil=False, is_submerged=True",
        fontsize=10, fontweight="bold",
    )

    # Pre-compute for shared y-limits
    all_res = []
    res_store = {}
    for hr in TARGET_HRS:
        eng_y_b, eng_x_b, eng_F_b = snaps_base[hr]
        eng_y_e, eng_x_e, eng_F_e = snaps_exp[hr]
        cw_field = cw_by_hr[hr]
        b_on_cw = resample_engine_to_cw(eng_y_b, eng_x_b, eng_F_b, cw_depths_m, cw_widths_m)
        e_on_cw = resample_engine_to_cw(eng_y_e, eng_x_e, eng_F_e, cw_depths_m, cw_widths_m)
        R_base = (b_on_cw[DI_MID, :] - cw_field[DI_MID, :])[::-1]
        R_exp  = (e_on_cw[DI_MID, :] - cw_field[DI_MID, :])[::-1]
        diff   = R_exp - R_base
        res_store[hr] = (R_base, R_exp, diff)
        all_res.extend([R_base, R_exp])

    all_vals = np.concatenate(all_res)
    rng = max(float(np.abs(all_vals).max()) * 1.2, GATE_TOL_F * 1.5)
    ylim = (-rng, rng)

    for col_idx, hr in enumerate(TARGET_HRS):
        ax = axes[col_idx]
        R_base, R_exp, diff = res_store[hr]

        ax.axhspan(-GATE_TOL_F, GATE_TOL_F, color="green", alpha=0.07, zorder=0)
        ax.axhline(0.0,          color="black", lw=0.8, zorder=1)
        ax.axhline( GATE_TOL_F, color="green", lw=1.0, ls="--", zorder=1)
        ax.axhline(-GATE_TOL_F, color="green", lw=1.0, ls="--", zorder=1)

        l1, = ax.plot(x_cw, R_base, "r-o",  ms=5, lw=1.5, zorder=3, label="Baseline")
        l2, = ax.plot(x_cw, R_exp,  "g-o",  ms=5, lw=1.5, zorder=3, label="Experimental (IC-aligned)")
        l3, = ax.plot(x_cw, diff,   "k--s", ms=3, lw=1.0, zorder=2, alpha=0.6, label="Δ (exp − base)")

        ax.axvline(wi9_x, color="gray", ls="--", lw=0.8, zorder=1)
        ymin, ymax = ylim
        ax.text(wi9_x + 0.05, ymin + (ymax - ymin) * 0.03,
                "wi=9", fontsize=6, color="gray", va="bottom")

        # Annotate wi=9 values for both curves
        wi9_orig = 9
        R_b9 = res_store[hr][0][::-1][wi9_orig]
        R_e9 = res_store[hr][1][::-1][wi9_orig]
        ax.annotate(f"base={R_b9:+.3f}°F", xy=(wi9_x, R_b9),
                    xytext=(wi9_x + 0.6, R_b9 - (ymax - ymin) * 0.1),
                    fontsize=7, color="darkred",
                    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.7))
        ax.annotate(f"exp={R_e9:+.3f}°F", xy=(wi9_x, R_e9),
                    xytext=(wi9_x + 0.6, R_e9 + (ymax - ymin) * 0.1),
                    fontsize=7, color="darkgreen",
                    arrowprops=dict(arrowstyle="->", color="darkgreen", lw=0.7))

        ax.set_xlim(-0.2, x_cw[-1] + 0.3)
        ax.set_ylim(ylim)
        ax.set_title(f"t={hr} hr", fontsize=9)
        ax.set_xlabel("Distance from edge (m)", fontsize=8)
        ax.set_ylabel("Residual: Engine − CW (°F)", fontsize=8)
        ax.tick_params(labelsize=7)
        if col_idx == 2:
            ax.legend(fontsize=7, loc="upper right")

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    return fig


# ---------------------------------------------------------------------------
def main():
    print("=== Stage 5c IC Alignment Test — Run I ===\n")
    print(f"Config: blanket={BLANKET_M}m, model_soil={MODEL_SOIL}, "
          f"is_submerged={IS_SUBMERGED}\n")

    # Run both configurations
    grid, snaps_base = run_engine(PLACEMENT_F, SOIL_F, experimental_ic=False)
    _,    snaps_exp  = run_engine(PLACEMENT_F, SOIL_F, experimental_ic=True)

    # Load CW
    cw_path = os.path.join(CW_RUNS, FOLDER, "output.txt")
    v = parse_cw_temp_output(cw_path)
    cw_by_hr = {}
    for hr in TARGET_HRS:
        ti = int(np.abs(v.time_hrs - hr).argmin())
        cw_by_hr[hr] = v.T_field_F[ti]
    cw_widths_m = v.widths_m
    cw_depths_m = v.depths_m
    print(f"\nCW loaded: di=24 at {cw_depths_m[DI_MID]:.4f} m, "
          f"wi=9 at {cw_widths_m[9]:.4f} m from edge")

    # Step 2 — sanity check
    sanity_check(grid, snaps_base, snaps_exp, cw_by_hr, cw_widths_m)

    # Step 3 — quantitative comparison table
    print_comparison_table(snaps_base, snaps_exp, cw_by_hr, cw_depths_m, cw_widths_m)

    # Masked max|R|
    mask = np.load(MASK_PATH)
    eng_y_b168, eng_x_b168, eng_F_b168 = snaps_base[168]
    eng_y_e168, eng_x_e168, eng_F_e168 = snaps_exp[168]
    cw_168 = cw_by_hr[168]
    maxR_base = compute_masked_max(eng_y_b168, eng_x_b168, eng_F_b168, cw_168, cw_depths_m, cw_widths_m, mask)
    maxR_exp  = compute_masked_max(eng_y_e168, eng_x_e168, eng_F_e168, cw_168, cw_depths_m, cw_widths_m, mask)
    print(f"Masked max|R| at t=168 hr:")
    print(f"  Baseline:     {maxR_base:.3f}°F  (Stage 5c reference: {BASELINE_MAX_R_REF}°F)")
    if abs(maxR_base - BASELINE_MAX_R_REF) > 0.01:
        print(f"  WARNING: baseline max|R| drifts {abs(maxR_base - BASELINE_MAX_R_REF):.3f}°F "
              f"from reference — investigate before trusting experimental result.")
    print(f"  Experimental: {maxR_exp:.3f}°F")
    print(f"  Δ (exp − base): {maxR_exp - maxR_base:+.3f}°F")

    # Step 4 — figure
    print("\n=== Generating Figure ===")
    fig = make_figure(snaps_base, snaps_exp, cw_by_hr, cw_widths_m, cw_depths_m)
    fig.savefig(OUT_FILE, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_FILE}")

    # Step 5 — bucket classification
    wi9_orig = 9
    b_on_cw = resample_engine_to_cw(eng_y_b168, eng_x_b168, eng_F_b168, cw_depths_m, cw_widths_m)
    e_on_cw = resample_engine_to_cw(eng_y_e168, eng_x_e168, eng_F_e168, cw_depths_m, cw_widths_m)
    R_base_wi9 = float(b_on_cw[DI_MID, wi9_orig] - cw_by_hr[168][DI_MID, wi9_orig])
    R_exp_wi9  = float(e_on_cw[DI_MID, wi9_orig] - cw_by_hr[168][DI_MID, wi9_orig])
    reduction_pct = (1.0 - abs(R_exp_wi9) / abs(R_base_wi9)) * 100.0 if abs(R_base_wi9) > 0 else 0.0

    print(f"\n=== Bucket Classification (§4) ===")
    print(f"  wi=9 residual at t=168 hr:  baseline={R_base_wi9:+.3f}°F  "
          f"experimental={R_exp_wi9:+.3f}°F  reduction={reduction_pct:.1f}%")
    print(f"  Masked max|R|: baseline={maxR_base:.3f}°F  experimental={maxR_exp:.3f}°F")

    if abs(R_exp_wi9) < 0.10 and maxR_exp < GATE_TOL_F:
        bucket = 1
        msg = ("Bucket 1 — IC was the dominant cause. The wi=9 hump collapses to <0.10°F "
               "when the engine's IC matches CW's t=0 frame. The Stage 5c gate failure is "
               "primarily a BC initialization protocol mismatch: pre-applying T_soil at "
               "the concrete edge column and bottom face before step 1 is sufficient to "
               "close the gate.")
    elif reduction_pct >= 30.0 and reduction_pct < 70.0:
        bucket = 2
        msg = ("Bucket 2 — IC contributed but is not the dominant cause. The wi=9 hump "
               "shrinks by {:.0f}% under IC alignment, but a residual hump remains. Two "
               "mechanisms are at play: the IC mismatch and a genuine lateral-diffusion-rate "
               "difference (likely cell-centered vs vertex-centered BC discretization). IC "
               "alignment is part of any eventual fix, but does not by itself close the gate.".format(reduction_pct))
    else:
        bucket = 3
        msg = ("Bucket 3 — IC is independent of the hump. The wi=9 residual is essentially "
               "unchanged ({:+.3f}°F vs baseline {:+.3f}°F, {:.1f}% change). The lateral "
               "diffusion hump is a genuine model-form mismatch — the IC alignment hypothesis "
               "is ruled out. Stage 5d must investigate the lateral BC discretization "
               "directly.".format(R_exp_wi9, R_base_wi9, reduction_pct))

    print(f"\n  ** BUCKET {bucket} **")
    print(f"  {msg}\n")


if __name__ == "__main__":
    main()
