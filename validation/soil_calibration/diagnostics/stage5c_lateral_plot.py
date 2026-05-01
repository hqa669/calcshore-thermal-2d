#!/usr/bin/env python3
"""Stage 5c diagnostic: lateral temperature profile at di=24 (mid-depth ~12 m).

Generates a 2×3 figure (Runs A and I × t=0, 84, 168 hr) showing CW reference
vs engine lateral profile, to characterize the wi=9 comparison artifact.

No code changes to engine, harness, or mask. Visualization only.
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
OUT_DIR  = os.path.join(SC, "diagnostics")
OUT_FILE = os.path.join(OUT_DIR, "stage5c_diagnostic_lateral_profiles_A_I.png")
os.makedirs(OUT_DIR, exist_ok=True)

DI_MID    = 24          # mid-depth comparison row (~12.19 m)
TARGET_HRS = [0, 84, 168]

RUN_SPECS = [
    ("A", "runA_baseline",  73, 73),
    ("I", "runI_100_73",   100, 73),
]

# Stage 5c engine config (same as gate run)
BLANKET_M = 0.0
MODEL_SOIL = False
IS_SUBMERGED = True

# Sanity threshold at t=0
SANITY_TOL_F = 0.1


# ---------------------------------------------------------------------------
# Step 1: Run engine for A and I, capturing snapshots at needed times
# ---------------------------------------------------------------------------
def run_engine(label, folder, placement_F, soil_F):
    dat_path = os.path.join(CW_RUNS, folder, "input.dat")
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

    T0_C = (placement_F - 32.0) * 5.0 / 9.0
    T_soil_C = (soil_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_soil] = T_soil_C

    env = make_neutral_env(placement_F)

    print(f"  Running engine Run {label} (placement={placement_F}°F, soil={soil_F}°F) ...")
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
    print(f"    Solved: {result.n_inner_steps} steps, dt={result.dt_inner_s:.0f}s")

    jslice, islice = grid.concrete_slice()
    T_conc_C = result.T_field_C[:, jslice, islice]
    T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jslice]   # (n_concrete_y,) in meters
    eng_x = grid.x[islice]   # (n_concrete_x,) in meters

    snapshots = {}  # target_hr -> (eng_y, eng_x, T_field_F at that time)
    for hr in TARGET_HRS:
        ti = nearest_time_idx(result.t_s, float(hr))  # nearest_time_idx takes hours
        snapshots[hr] = (eng_y.copy(), eng_x.copy(), T_conc_F[ti].copy())
        print(f"    t={hr} hr: engine ti={ti}, t_actual={result.t_s[ti]/3600:.2f} hr, "
              f"Tmin={T_conc_F[ti].min():.2f}°F Tmax={T_conc_F[ti].max():.2f}°F")

    return snapshots


# ---------------------------------------------------------------------------
# Step 2: Load CW data for A and I
# ---------------------------------------------------------------------------
def load_cw_all(folder):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    cw_by_hr = {}
    for hr in TARGET_HRS:
        ti = int(np.abs(v.time_hrs - hr).argmin())
        cw_by_hr[hr] = {
            "field": v.T_field_F[ti],       # (49, 13)
            "actual_hr": float(v.time_hrs[ti]),
        }
        print(f"  CW Run {folder}: t={hr} hr → ti={ti}, actual={v.time_hrs[ti]:.2f} hr, "
              f"Tmin={v.T_field_F[ti].min():.2f}°F Tmax={v.T_field_F[ti].max():.2f}°F")
    return cw_by_hr, v.widths_m, v.depths_m


# ---------------------------------------------------------------------------
# Step 3: Sanity check at t=0 for a given run
# ---------------------------------------------------------------------------
def sanity_check(label, eng_snap, cw_by_hr, cw_widths_m, cw_depths_m):
    eng_y, eng_x, eng_F = eng_snap[0]
    cw_field_t0 = cw_by_hr[0]["field"]
    cw_t0_hr    = cw_by_hr[0]["actual_hr"]
    eng_on_cw   = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    row_cw  = cw_field_t0[DI_MID, :]    # CW di=24, all wi
    row_eng = eng_on_cw[DI_MID, :]      # Engine resampled to CW nodes, di=24
    diff    = row_eng - row_cw
    print(f"\nSanity check Run {label} at t~0 hr (CW actual={cw_t0_hr:.2f} hr):")
    for wi in range(len(cw_widths_m)):
        print(f"  wi={wi:2d} x={cw_widths_m[wi]:.2f}m: CW={row_cw[wi]:.3f}°F  "
              f"Eng={row_eng[wi]:.3f}°F  diff={diff[wi]:+.3f}°F")
    # Sanity check excludes wi=11,12 (edge columns): the side BC sets T_soil
    # there immediately, but the engine t=0 IC is still at placement temp.
    # wi=0..10 are interior and should be near the placement temperature at t~0.
    max_diff = float(np.abs(diff[:11]).max())
    ok = max_diff <= SANITY_TOL_F
    print(f"  → max|diff| = {max_diff:.4f}°F  "
          f"({'PASS' if ok else 'FAIL — stopping'})")
    if not ok:
        raise RuntimeError(
            f"Sanity FAIL Run {label} t~0: max|diff|={max_diff:.4f}°F > {SANITY_TOL_F}°F"
        )
    return True


# ---------------------------------------------------------------------------
# Step 4: Build plot
# ---------------------------------------------------------------------------
def make_figure(eng_all, cw_all, cw_widths_m, cw_depths_m):
    # x-axis: "distance from edge" = cw_widths_m (already in edge=0, CL=max convention).
    # CW array order is [CL-side ... edge-side] i.e. widths_m[0]=6.1, widths_m[-1]=0.
    # Reverse so x increases left-to-right from edge (0 m) to CL (~6.1 m).
    x_cw = cw_widths_m[::-1]    # [0.0, 0.51, ..., 6.1] — distance from edge

    wi9_x = cw_widths_m[9]      # 1.52 m from edge
    cl_x  = cw_widths_m[0]      # 6.1 m from edge (CL)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(
        "Stage 5c Diagnostic: Lateral Temperature Profile, di=24 (mid-depth ≈12.19 m)\n"
        "Runs A (|ΔT|=0°F) and I (|ΔT|=27°F) — blanket_thickness_m=0.0, model_soil=False, is_submerged=True",
        fontsize=10, fontweight="bold",
    )

    run_labels   = ["A", "I"]
    placements   = {"A": 73, "I": 100}
    soils        = {"A": 73, "I": 73}

    # Collect y-limits per column (same timestep) for consistent y-axis
    ylims_per_col = {}   # col_idx -> (ymin, ymax)

    for row_idx, label in enumerate(run_labels):
        eng_snap = eng_all[label]
        cw_by_hr, _, _ = cw_all[label]

        for col_idx, hr in enumerate(TARGET_HRS):
            ax = axes[row_idx, col_idx]
            eng_y, eng_x, eng_F = eng_snap[hr]
            cw_field = cw_by_hr[hr]["field"]    # (49, 13)

            # CW di=24 lateral profile (reversed to edge→CL order)
            cw_row = cw_field[DI_MID, :][::-1]

            # Engine full-resolution profile at y = cw_depths_m[DI_MID]
            y_target = float(cw_depths_m[DI_MID])
            eng_at_y = np.array([
                float(np.interp(y_target, eng_y, eng_F[:, xi]))
                for xi in range(eng_F.shape[1])
            ])
            # Engine x goes from 0 (edge) to ~6.096 m (CL)

            # Engine downsampled to CW comparison nodes (matches harness)
            eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
            eng_cw_row = eng_on_cw[DI_MID, :][::-1]  # reversed to edge→CL order

            # Residual at wi=9 (engine minus CW, at the artifact column)
            # wi=9 in original CW order; in reversed x_cw order it's index 12-9=3
            wi9_orig_idx = 9
            res_wi9 = eng_on_cw[DI_MID, wi9_orig_idx] - cw_field[DI_MID, wi9_orig_idx]

            # Plot
            l1, = ax.plot(x_cw, cw_row, "b-o", ms=5, lw=1.5, label="CW")
            l2, = ax.plot(eng_x, eng_at_y, "r-", lw=1.0, alpha=0.7, label="Engine (full)")
            l3, = ax.plot(x_cw, eng_cw_row, "r--o", ms=5, lw=1.5, label="Engine (CW nodes)")

            # Mark wi=9 and CL
            ax.axvline(wi9_x, color="gray", linestyle="--", lw=0.8)
            ax.axvline(cl_x,  color="gray", linestyle=":",  lw=0.8)
            ax.text(wi9_x, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else 0,
                    " wi=9", fontsize=6, color="gray", va="bottom")
            ax.text(cl_x, ax.get_ylim()[0] if ax.get_ylim()[0] != 0 else 0,
                    " CL", fontsize=6, color="gray", va="bottom")

            # Annotate wi=9 residual for Run I
            if label == "I" or abs(res_wi9) > 0.02:
                ax.annotate(
                    f"R(wi=9)={res_wi9:+.3f}°F",
                    xy=(wi9_x, eng_on_cw[DI_MID, wi9_orig_idx]),
                    xytext=(wi9_x + 0.4, eng_on_cw[DI_MID, wi9_orig_idx]),
                    fontsize=7, color="darkred",
                    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.7),
                )

            ax.set_title(
                f"Run {label}, t={hr} hr (|ΔT|={abs(soils[label]-placements[label])}°F)",
                fontsize=9,
            )
            ax.set_xlabel("Distance from edge (m)", fontsize=8)
            ax.set_ylabel("Temperature (°F)", fontsize=8)
            ax.tick_params(labelsize=7)
            if row_idx == 0 and col_idx == 0:
                ax.legend(fontsize=7, loc="upper right")

            # Collect y-range for matching within column
            all_vals = np.concatenate([cw_row, eng_at_y, eng_cw_row])
            mn, mx = float(all_vals.min()), float(all_vals.max())
            pad = max(0.5, (mx - mn) * 0.1)
            if col_idx not in ylims_per_col:
                ylims_per_col[col_idx] = [mn - pad, mx + pad]
            else:
                ylims_per_col[col_idx][0] = min(ylims_per_col[col_idx][0], mn - pad)
                ylims_per_col[col_idx][1] = max(ylims_per_col[col_idx][1], mx + pad)

    # Apply matched y-limits per column and re-draw vertical markers at correct positions
    for col_idx, hr in enumerate(TARGET_HRS):
        ymin, ymax = ylims_per_col[col_idx]
        for row_idx in range(2):
            ax = axes[row_idx, col_idx]
            ax.set_ylim(ymin, ymax)
            # Re-draw vertical markers after ylim is set
            for line in ax.lines:
                pass  # already plotted; ylim change is enough
            # Re-annotate wi=9 and CL text at correct position
            for txt in ax.texts:
                if "wi=9" in txt.get_text():
                    txt.set_y(ymin + (ymax - ymin) * 0.02)
                elif txt.get_text().strip() == "CL":
                    txt.set_y(ymin + (ymax - ymin) * 0.02)

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=== Stage 5c Lateral Profile Diagnostic ===\n")
    print(f"Config: blanket={BLANKET_M}m, model_soil={MODEL_SOIL}, "
          f"is_submerged={IS_SUBMERGED}, grid_refinement=6")
    print(f"Output: {OUT_FILE}\n")

    # Run engine
    eng_all = {}
    for label, folder, pl, so in RUN_SPECS:
        print(f"--- Run {label} ---")
        eng_all[label] = run_engine(label, folder, pl, so)

    # Load CW
    cw_all = {}
    cw_widths_m = None
    cw_depths_m = None
    for label, folder, pl, so in RUN_SPECS:
        print(f"\n--- CW data Run {label} ---")
        cw_by_hr, cw_widths_m, cw_depths_m = load_cw_all(folder)
        cw_all[label] = (cw_by_hr, cw_widths_m, cw_depths_m)

    print(f"\ndi=24: {cw_depths_m[DI_MID]:.4f} m, "
          f"wi=9: {cw_widths_m[9]:.4f} m from edge")

    # Sanity checks
    print("\n=== t~0 Sanity Checks ===")
    for label, folder, pl, so in RUN_SPECS:
        cw_by_hr, cw_w, cw_d = cw_all[label]
        sanity_check(label, eng_all[label], cw_by_hr, cw_w, cw_d)

    # Generate figure
    print("\n=== Generating Figure ===")
    cw_by_hr_A, cw_w, cw_d = cw_all["A"]
    cw_by_hr_I, _, _       = cw_all["I"]
    cw_all_for_plot = {
        "A": (cw_by_hr_A, cw_w, cw_d),
        "I": (cw_by_hr_I, cw_w, cw_d),
    }

    fig = make_figure(eng_all, cw_all_for_plot, cw_w, cw_d)
    fig.savefig(OUT_FILE, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_FILE}")

    # Print residuals at wi=9 for both runs at all timestamps
    print("\n=== wi=9 Residuals (Engine − CW) at di=24 ===")
    for label, folder, pl, so in RUN_SPECS:
        cw_by_hr, cw_w, cw_d = cw_all[label]
        print(f"Run {label} (|ΔT|={abs(so-pl)}°F):")
        for hr in TARGET_HRS:
            eng_y, eng_x, eng_F = eng_all[label][hr]
            cw_field = cw_by_hr[hr]["field"]
            eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_d, cw_w)
            res = eng_on_cw[DI_MID, :] - cw_field[DI_MID, :]
            print(f"  t={hr:3d} hr:  wi=8={res[8]:+.3f}  wi=9={res[9]:+.3f}  "
                  f"wi=10={res[10]:+.3f}  wi=12(edge)={res[12]:+.3f}  wi=0(CL)={res[0]:+.3f}")


if __name__ == "__main__":
    main()
