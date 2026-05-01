#!/usr/bin/env python3
"""Stage 5c diagnostic: lateral residual (Engine − CW) at di=24 (mid-depth ~12 m).

2×3 figure (Runs A and I × t=0, 84, 168 hr). Gate threshold ±0.35°F shown.
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
OUT_FILE = os.path.join(OUT_DIR, "stage5c_diagnostic_residuals_all.png")
os.makedirs(OUT_DIR, exist_ok=True)

DI_MID     = 24
TARGET_HRS = [0, 84, 168]
GATE_TOL_F = 0.35

RUN_SPECS = [
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

BLANKET_M    = 0.0
MODEL_SOIL   = False
IS_SUBMERGED = True


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

    T0_C     = (placement_F - 32.0) * 5.0 / 9.0
    T_soil_C = (soil_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_soil] = T_soil_C

    env = make_neutral_env(placement_F)

    print(f"  Running engine Run {label} ...")
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
    return snapshots


def load_cw(folder):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    cw_by_hr = {}
    for hr in TARGET_HRS:
        ti = int(np.abs(v.time_hrs - hr).argmin())
        cw_by_hr[hr] = v.T_field_F[ti]
    return cw_by_hr, v.widths_m, v.depths_m


# ---------------------------------------------------------------------------
def make_figure(eng_all, cw_all, cw_widths_m, cw_depths_m):
    # Reversed so x goes edge→CL (0 → ~6.1 m)
    x_cw   = cw_widths_m[::-1]
    wi9_x  = cw_widths_m[9]    # 1.52 m from edge
    cl_x   = cw_widths_m[0]    # 6.1 m (CL)

    run_labels = [s[0] for s in RUN_SPECS]
    dT_map     = {lbl: abs(so - pl) for lbl, _, pl, so in RUN_SPECS}

    n_runs = len(run_labels)
    fig, axes = plt.subplots(n_runs, 3, figsize=(15, 3.5 * n_runs), sharey=False)
    fig.suptitle(
        "Stage 5c Diagnostic: Lateral Residual (Engine − CW), di=24 (mid-depth ≈12.19 m)\n"
        "Runs A–I — blanket_thickness_m=0.0, model_soil=False, is_submerged=True",
        fontsize=10, fontweight="bold",
    )

    # Compute all residuals first so we can set shared y-limits per column
    res_data = {}  # (label, hr) -> residual array (13,) in edge→CL order
    for label in run_labels:
        cw_by_hr = cw_all[label]
        for hr in TARGET_HRS:
            eng_y, eng_x, eng_F = eng_all[label][hr]
            cw_field = cw_by_hr[hr]
            eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
            res = eng_on_cw[DI_MID, :] - cw_field[DI_MID, :]
            res_data[(label, hr)] = res[::-1]  # reversed to edge→CL

    # Per-column y-limits (shared across all runs at same timestep)
    col_ylims = {}
    for col_idx, hr in enumerate(TARGET_HRS):
        vals = np.concatenate([res_data[(lbl, hr)] for lbl in run_labels])
        rng = max(float(np.abs(vals).max()) * 1.2, GATE_TOL_F * 1.4)
        col_ylims[col_idx] = (-rng, rng)

    for row_idx, label in enumerate(run_labels):
        cw_by_hr = cw_all[label]
        for col_idx, hr in enumerate(TARGET_HRS):
            ax = axes[row_idx, col_idx]
            res = res_data[(label, hr)]

            # Gate threshold bands
            ax.axhspan(-GATE_TOL_F, GATE_TOL_F, color="green", alpha=0.07, zorder=0)
            ax.axhline(0.0,          color="black",  lw=0.8, zorder=1)
            ax.axhline( GATE_TOL_F, color="green",  lw=1.0, ls="--", zorder=1,
                        label=f"±{GATE_TOL_F}°F gate")
            ax.axhline(-GATE_TOL_F, color="green",  lw=1.0, ls="--", zorder=1)

            # Residual
            ax.plot(x_cw, res, "r-o", ms=5, lw=1.5, zorder=3)
            ax.fill_between(x_cw, res, 0, where=(np.abs(res) > GATE_TOL_F),
                            color="red", alpha=0.20, zorder=2, label="Outside gate")

            # Vertical markers
            ymin, ymax = col_ylims[col_idx]
            ax.axvline(wi9_x, color="gray", ls="--", lw=0.8, zorder=1)
            ax.axvline(cl_x,  color="gray", ls=":",  lw=0.8, zorder=1)
            ax.text(wi9_x + 0.05, ymin + (ymax - ymin) * 0.03,
                    "wi=9", fontsize=6, color="gray", va="bottom")
            ax.text(cl_x  + 0.05, ymin + (ymax - ymin) * 0.03,
                    "CL",   fontsize=6, color="gray", va="bottom")

            # Annotate wi=9 value
            wi9_orig = 9
            res_wi9_val = res_data[(label, hr)][::-1][wi9_orig]  # back to original order
            ax.annotate(
                f"{res_wi9_val:+.3f}°F",
                xy=(wi9_x, res_wi9_val),
                xytext=(wi9_x + 0.5, res_wi9_val + (ymax - ymin) * 0.08),
                fontsize=7, color="darkred",
                arrowprops=dict(arrowstyle="->", color="darkred", lw=0.7),
            )

            ax.set_xlim(-0.2, x_cw[-1] + 0.3)
            ax.set_ylim(ymin, ymax)
            ax.set_title(
                f"Run {label}, t={hr} hr (|ΔT|={dT_map[label]}°F)",
                fontsize=9,
            )
            ax.set_xlabel("Distance from edge (m)", fontsize=8)
            ax.set_ylabel("Residual: Engine − CW (°F)", fontsize=8)
            ax.tick_params(labelsize=7)
            if row_idx == 0 and col_idx == 2:
                ax.legend(fontsize=7, loc="upper right")

    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


# ---------------------------------------------------------------------------
def main():
    print("=== Stage 5c Lateral Residual Diagnostic ===\n")

    eng_all = {}
    for label, folder, pl, so in RUN_SPECS:
        print(f"--- Run {label} ---")
        eng_all[label] = run_engine(label, folder, pl, so)

    cw_all = {}
    cw_widths_m = cw_depths_m = None
    for label, folder, pl, so in RUN_SPECS:
        cw_by_hr, cw_widths_m, cw_depths_m = load_cw(folder)
        cw_all[label] = cw_by_hr

    print(f"\ndi=24: {cw_depths_m[DI_MID]:.4f} m,  wi=9: {cw_widths_m[9]:.4f} m from edge")

    print("\n=== Generating Figure ===")
    fig = make_figure(eng_all, cw_all, cw_widths_m, cw_depths_m)
    fig.savefig(OUT_FILE, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT_FILE}")

    print("\n=== wi=9 Residuals (Engine − CW) at di=24 ===")
    print(f"{'Run':>4}  {'t(hr)':>6}  {'wi=8':>7}  {'wi=9':>7}  {'wi=10':>7}  {'CL':>7}")
    for label, folder, pl, so in RUN_SPECS:
        for hr in TARGET_HRS:
            eng_y, eng_x, eng_F = eng_all[label][hr]
            cw_field = cw_all[label][hr]
            eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
            res = eng_on_cw[DI_MID, :] - cw_field[DI_MID, :]
            print(f"  {label:>2}  {hr:>6d}  {res[8]:>+7.3f}  {res[9]:>+7.3f}  "
                  f"{res[10]:>+7.3f}  {res[0]:>+7.3f}")


if __name__ == "__main__":
    main()
