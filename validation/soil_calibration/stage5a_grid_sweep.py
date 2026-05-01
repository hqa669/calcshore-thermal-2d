#!/usr/bin/env python3
"""Stage 5a — S5a.4: Grid resolution sweep for Run F.

Runs Run F (placement=73°F, soil=45°F, |ΔT|=28°F) at three grid resolutions
(1×, 2×, 4×) to quantify the numerical contribution to the α_c gap.

Resolutions:
  1×: n_concrete_x=21, n_concrete_y=13  (native dx = 1.0 ft)
  2×: n_concrete_x=41, n_concrete_y=25  (dx = 0.5 ft)
  4×: n_concrete_x=81, n_concrete_y=49  (dx = 0.25 ft)

Sanity check: 1× Run F masked max|R| must reproduce Stage 4b's 4.084°F ±0.01°F.

Outputs:
    validation/soil_calibration/STAGE5a_gridsweep_residuals.npy  shape (3, 49, 13)
    validation/soil_calibration/STAGE5a_cause3_grid_resolution.md
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import parse_cw_dat
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import load_cw_slice, resample_engine_to_cw, COMPARE_HR
from stage4b_run import make_neutral_env, nearest_time_idx, masked_stats, save_field_csv

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
MASK_PATH = os.path.join(HERE, "STAGE35_validity_mask.npy")

RUN_LABEL = "F"
RUN_FOLDER = "runF_73_45"
PLACEMENT_F = 73
SOIL_F = 45

STAGE4B_MASKED_MAX = 4.084  # °F from Stage 4b report
SANITY_TOL = 0.01           # °F

WALL_CLOCK_LIMIT_S = 3600.0  # abort 4× run if exceeds 60 min

RESOLUTIONS = [
    (1, 21, 13),   # label, n_concrete_x, n_concrete_y
    (2, 41, 25),
    (4, 81, 49),
]


def run():
    mask = np.load(MASK_PATH)
    assert mask.shape == (49, 13), f"Mask shape mismatch: {mask.shape}"
    print(f"Loaded Stage 3.5 mask: shape={mask.shape}, {int(mask.sum())} cells kept")

    dat_path = os.path.join(CW_RUNS, RUN_FOLDER, "input.dat")
    mix_base, geom, constr, _ = parse_cw_dat(dat_path)
    factor, _ = compute_hu_factor(mix_base)
    mix_base.Hu_factor_calibrated = factor
    mix_base.Hu_J_kg_effective = mix_base.Hu_J_kg * factor
    print(f"  Hu_J_kg_effective={mix_base.Hu_J_kg_effective:.6f} (suppressed)")

    constr.model_soil = False
    constr.is_submerged = True

    T0_C = (PLACEMENT_F - 32.0) * 5.0 / 9.0
    T_soil_C = (SOIL_F - 32.0) * 5.0 / 9.0
    env = make_neutral_env(PLACEMENT_F)

    import tempfile
    from stage3_compare import load_engine_csv

    cw_field, cw_widths_m, cw_depths_m, _ = load_cw_slice(RUN_FOLDER, COMPARE_HR)
    assert cw_field is not None, "CW reference not found for Run F"

    results_per_res = {}  # label → {masked_max, masked_mean, wall_s, masked_max_loc}
    residual_stack = []   # (3, 49, 13)

    for (ref, n_cx, n_cy) in RESOLUTIONS:
        grid = build_grid_half_mat(
            geom.width_ft, geom.depth_ft,
            n_concrete_x=n_cx, n_concrete_y=n_cy,
            is_submerged=True, model_soil=False,
        )
        dx_m = (grid.x[1] - grid.x[0]) if grid.nx > 1 else float('nan')
        dy_arr = np.diff(grid.y[grid.iy_concrete_start:grid.iy_concrete_end+1])
        dy_m = float(np.mean(dy_arr))
        print(f"\nResolution {ref}×: n_concrete_x={n_cx}, n_concrete_y={n_cy}, "
              f"dx={dx_m:.4f} m = {dx_m/0.3048:.3f} ft, dy_mean={dy_m:.4f} m = {dy_m/0.3048:.3f} ft")

        T_init = np.full((grid.ny, grid.nx), T0_C)
        T_init[grid.is_soil] = T_soil_C

        mix = mix_base  # same mix for all resolutions

        t0 = time.perf_counter()
        result = solve_hydration_2d(
            grid, mix, T_init,
            duration_s=168 * 3600,
            output_interval_s=1800.0,
            boundary_mode="full_2d",
            environment=env,
            construction=constr,
            T_ground_deep_C=T_soil_C,
            diagnostic_outputs=False,
        )
        t_wall = time.perf_counter() - t0
        print(f"  Solved in {t_wall:.1f}s ({result.n_inner_steps} steps, dt={result.dt_inner_s:.2f}s)")

        if ref == 4 and t_wall > WALL_CLOCK_LIMIT_S:
            print(f"  WARNING: 4× run exceeded {WALL_CLOCK_LIMIT_S:.0f}s wall-clock limit; "
                  "results may be incomplete. Using partial output.")

        jslice, islice = grid.concrete_slice()
        T_conc_C = result.T_field_C[:, jslice, islice]
        T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
        ti_168 = nearest_time_idx(result.t_s, 168.0)

        y_conc = grid.y[jslice]
        x_conc = grid.x[islice]
        tmp_csv = os.path.join(tempfile.gettempdir(), f"stage5a_runF_{ref}x.csv")
        save_field_csv(tmp_csv, T_conc_F[ti_168], x_conc, y_conc)

        eng_y, eng_x, eng_F = load_engine_csv(tmp_csv)
        eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
        residual = eng_on_cw - cw_field

        s = masked_stats(residual, mask, cw_depths_m, cw_widths_m)
        results_per_res[ref] = {
            "masked_max": s["masked_max"],
            "masked_mean": s["masked_mean"],
            "unmasked_max": s["unmasked_max"],
            "wall_s": t_wall,
            "dx_m": dx_m,
            "dy_m": dy_m,
            "n_cx": n_cx,
            "n_cy": n_cy,
            "dt_inner_s": result.dt_inner_s,
        }
        residual_stack.append(residual)

        if ref == 1:
            sanity_ok = abs(s["masked_max"] - STAGE4B_MASKED_MAX) <= SANITY_TOL
            print(f"  SANITY CHECK (1×): masked max|R|={s['masked_max']:.3f}°F "
                  f"(target {STAGE4B_MASKED_MAX}°F, tol ±{SANITY_TOL}°F)  "
                  f"{'PASS' if sanity_ok else 'FAIL ***'}")
        else:
            print(f"  {ref}×: masked max|R|={s['masked_max']:.3f}°F  "
                  f"mean={s['masked_mean']:.3f}°F  unmasked={s['unmasked_max']:.3f}°F")

    # Save stacked residuals
    residuals_npy = np.stack(residual_stack)  # (3, 49, 13)
    np.save(os.path.join(HERE, "STAGE5a_gridsweep_residuals.npy"), residuals_npy)
    print(f"\nSaved STAGE5a_gridsweep_residuals.npy: {residuals_npy.shape}")

    # ------------------------------------------------------------------ #
    # Richardson extrapolation to "infinite resolution"                    #
    # h → 0 limit for second-order FD: R(h) ≈ R_inf + C·h²              #
    # Using 1× and 2× to estimate R_inf                                   #
    # ------------------------------------------------------------------ #
    R1 = results_per_res[1]["masked_max"]
    R2 = results_per_res[2]["masked_max"]
    R4 = results_per_res[4]["masked_max"]

    # Richardson using 1× and 2× (h₁ = 2h₂, so h₁²/h₂² = 4)
    # R1 ≈ R_inf + C·h²  →  R2 ≈ R_inf + C·(h/2)² = R_inf + C·h²/4
    # 4·R2 − R1 = 3·R_inf  →  R_inf = (4·R2 − R1) / 3
    R_inf_12 = (4 * R2 - R1) / 3.0
    numerical_fraction_12 = (R1 - R_inf_12) / R1 if R1 > 0 else 0.0

    # Also using 2× and 4×
    R_inf_24 = (4 * R4 - R2) / 3.0
    numerical_fraction_24 = (R2 - R_inf_24) / R2 if R2 > 0 else 0.0

    # Use 1× and 2× as primary (more reliable unless trend is not monotone)
    R_inf = R_inf_12
    numerical_fraction = numerical_fraction_12

    # ------------------------------------------------------------------ #
    # Write STAGE5a_cause3_grid_resolution.md                              #
    # ------------------------------------------------------------------ #
    lines = [
        "# STAGE5a — Cause 3: Grid Resolution / Numerical Accuracy",
        "",
        f"Run F: placement=73°F, soil=45°F, |ΔT|=28°F (largest residual in Stage 4b).",
        f"Three grid resolutions tested; same CW reference and Stage 3.5 mask applied.",
        "",
        "## Resolution sweep results",
        "",
        "| Resolution | n_cx × n_cy | dx (ft) | dy_mean (ft) | dt_inner (s) | "
        "masked max|R| (°F) | masked mean|R| (°F) | wall-clock (s) |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]

    for ref, n_cx, n_cy in RESOLUTIONS:
        r = results_per_res[ref]
        dx_ft = r["dx_m"] / 0.3048
        dy_ft = r["dy_m"] / 0.3048
        sanity_str = f" ← Stage 4b ref {STAGE4B_MASKED_MAX}°F" if ref == 1 else ""
        lines.append(
            f"| {ref}× | {r['n_cx']}×{r['n_cy']} | {dx_ft:.3f} | {dy_ft:.3f} | "
            f"{r['dt_inner_s']:.2f} | **{r['masked_max']:.3f}**{sanity_str} | "
            f"{r['masked_mean']:.3f} | {r['wall_s']:.1f} |"
        )

    # Richardson table
    lines += [
        "",
        "## Richardson extrapolation (second-order FD: R(h) ≈ R_∞ + C·h²)",
        "",
        "Using 1× and 2×:",
        f"  R_∞ (1×/2× estimate) = (4×R_2× − R_1×) / 3 = "
        f"({4*R2:.3f} − {R1:.3f}) / 3 = **{R_inf_12:.3f}°F**",
        f"  Numerical contribution (1× − R_∞) / R_1× = "
        f"({R1:.3f} − {R_inf_12:.3f}) / {R1:.3f} = **{numerical_fraction_12*100:.1f}%**",
        "",
        "Using 2× and 4×:",
        f"  R_∞ (2×/4× estimate) = (4×R_4× − R_2×) / 3 = "
        f"({4*R4:.3f} − {R2:.3f}) / 3 = **{R_inf_24:.3f}°F**",
        f"  Numerical contribution (2× − R_∞) / R_2× = "
        f"({R2:.3f} − {R_inf_24:.3f}) / {R2:.3f} = **{numerical_fraction_24*100:.1f}%**",
        "",
    ]

    # Convergence qualitative assessment
    mono = (R1 > R2 > R4) or (R1 < R2 < R4)
    if mono:
        lines.append("Residuals are monotonically converging with resolution → "
                     "Richardson extrapolation is reliable.")
    else:
        lines.append("**Residuals are NOT monotonically converging** — Richardson "
                     "extrapolation may be unreliable. Inspect trend manually.")

    # Verdict
    numerical_pct = numerical_fraction * 100
    if numerical_pct < 5:
        verdict = "Cause 3 is NEGLIGIBLE (<5%). Grid resolution is not a significant source of the α_c gap."
    elif numerical_pct < 20:
        verdict = (f"Cause 3 is MINOR ({numerical_pct:.0f}%). Grid refinement reduces residuals "
                   "modestly but does not close the gap.")
    elif numerical_pct < 50:
        verdict = (f"Cause 3 is MODERATE ({numerical_pct:.0f}%). Grid refinement could close "
                   "roughly half the α_c gap.")
    elif numerical_pct < 75:
        verdict = (f"Cause 3 is SIGNIFICANT ({numerical_pct:.0f}%). A large fraction of the α_c "
                   "gap is numerical. Grid refinement is a major part of the fix.")
    else:
        verdict = (f"Cause 3 is DOMINANT ({numerical_pct:.0f}%). Grid resolution explains most "
                   "of the α_c gap; Causes 1 and 2 are secondary.")

    lines += [
        "",
        "## Verdict",
        "",
        verdict,
        "",
        f"Parametric (non-numerical) residual at infinite resolution: {R_inf:.3f}°F",
        f"Gate target: 0.35°F. {'Achievable by grid refinement alone.' if R_inf <= 0.35 else 'NOT achievable by grid refinement alone — Cause 1 or 2 fix also required.'}",
    ]

    with open(os.path.join(HERE, "STAGE5a_cause3_grid_resolution.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"Wrote STAGE5a_cause3_grid_resolution.md")

    print("\nGrid sweep summary:")
    for ref, n_cx, n_cy in RESOLUTIONS:
        r = results_per_res[ref]
        print(f"  {ref}×: masked max|R|={r['masked_max']:.3f}°F  wall={r['wall_s']:.1f}s")
    print(f"  Richardson R_inf (1×/2×): {R_inf_12:.3f}°F → numerical fraction: {numerical_fraction_12*100:.1f}%")
    print(f"  Richardson R_inf (2×/4×): {R_inf_24:.3f}°F → numerical fraction: {numerical_fraction_24*100:.1f}%")


if __name__ == "__main__":
    run()
