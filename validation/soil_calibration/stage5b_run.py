#!/usr/bin/env python3
"""Stage 5b — 9-run validation at 6× grid refinement (grid_refinement=6 default).

Stage 5a diagnosed the 2–4°F residuals in Stage 4b as ~90–93% numerical
discretization error from the coarse 1.0 × 6.67 ft native grid.  Stage 5b
applies the fix: the engine default is now grid_refinement=6 (121 × 73 nodes,
~0.17 × 1.10 ft), predicted to close the 0.35°F gate.

Gate: all 9 runs (A–I) masked max|R| ≤ 0.35°F.

Outputs:
    validation/soil_calibration/stage5b_runs/run{A..I}_t168.csv
    validation/soil_calibration/plots/STAGE5b_run{A..I}_engine_vs_cw.png
    validation/soil_calibration/STAGE5b_validation.md
"""
import os
import sys
import json
import time
import traceback
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import parse_cw_dat
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat

from stage3_compare import (
    load_engine_csv,
    load_cw_slice,
    resample_engine_to_cw,
    COMPARE_HR,
    VMIN,
    VMAX,
)
from stage4b_run import (
    make_neutral_env,
    nearest_time_idx,
    save_field_csv,
    masked_stats,
    RUNS,
)

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
ENG_DIR = os.path.join(HERE, "stage5b_runs")
PLOTS = os.path.join(HERE, "plots")
os.makedirs(ENG_DIR, exist_ok=True)
os.makedirs(PLOTS, exist_ok=True)

MASK_PATH = os.path.join(HERE, "STAGE35_validity_mask.npy")

# Stage 4b masked max|R| values for comparison table.
STAGE4B_MASKED_MAX = {
    "A": 0.146, "B": 1.940, "C": 2.391, "D": 1.940,
    "E": 2.391, "F": 4.084, "G": 3.887, "H": 4.084, "I": 3.887,
}

GATE_ALL = 0.35  # °F — 1% of 35°F ΔT_max spec


def plot_comparison_5b(label, placement, soil, cw_field, eng_field, residual,
                       cw_depths_m, cw_widths_m, mask, grid_info_str=""):
    delta_T = soil - placement
    w_min = float(cw_widths_m[-1])
    w_max = float(cw_widths_m[0])
    d_max = float(cw_depths_m[-1])
    extent = [w_min, w_max, d_max, 0.0]

    cw_plot  = cw_field[:, ::-1]
    eng_plot = eng_field[:, ::-1]
    res_plot = residual[:, ::-1]
    sym = max(0.01, float(np.abs(residual).max()))
    masked_sym = max(0.01, float(np.abs(residual[mask]).max()))

    fig, axes = plt.subplots(3, 1, figsize=(6, 14))
    fig.suptitle(
        f"Stage 5b — Run {label}: placement={placement}°F, soil={soil}°F, "
        f"ΔT={delta_T:+d}°F\n"
        f"model_soil=False, is_submerged=True, grid_refinement=6  |  t={COMPARE_HR:.0f} hr\n"
        f"unmasked max|R|={sym:.2f}°F  masked max|R|={masked_sym:.3f}°F  "
        f"gate≤{GATE_ALL}°F  {'PASS' if masked_sym <= GATE_ALL else 'FAIL'}",
        fontsize=8, fontweight="bold",
    )
    panel_data = [
        (cw_plot,  "CW (reference)",             "inferno", (VMIN, VMAX)),
        (eng_plot, "Engine (resampled to CW grid)", "inferno", (VMIN, VMAX)),
        (res_plot, f"Engine−CW (sym={sym:.3f}°F)", "RdBu_r", (-sym, sym)),
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
    fig.tight_layout(rect=[0, 0, 1, 0.91])

    out = os.path.join(PLOTS, f"STAGE5b_run{label}_engine_vs_cw.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def run():
    if not os.path.exists(MASK_PATH):
        raise FileNotFoundError(f"Stage 3.5 mask not found: {MASK_PATH}")

    mask = np.load(MASK_PATH)
    print(f"Loaded Stage 3.5 mask: shape={mask.shape}, {int(mask.sum())} cells kept")

    manifest = {}
    grid_info = None
    wall_clocks = {}

    for label, folder, placement, soil in RUNS:
        print(f"\nRun {label} ({folder}): placement={placement}°F, soil={soil}°F")
        dat_path = os.path.join(CW_RUNS, folder, "input.dat")
        if not os.path.isfile(dat_path):
            print(f"  SKIP: input.dat not found")
            manifest[label] = {"status": "skip"}
            continue

        try:
            from kinetics_correction import compute_hu_factor
            mix, geom, constr, _ = parse_cw_dat(dat_path)
            factor, _ = compute_hu_factor(mix)
            mix.Hu_factor_calibrated = factor
            mix.Hu_J_kg_effective = mix.Hu_J_kg * factor
            print(f"  Hu_J_kg_effective={mix.Hu_J_kg_effective:.4f}")

            constr.model_soil = False
            constr.is_submerged = True

            # Use 6× default (no grid kwargs passed).
            grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                                       is_submerged=True, model_soil=False)
            if grid_info is None:
                jslice, islice = grid.concrete_slice()
                grid_info = {
                    "grid_refinement": 6,
                    "nx_full": grid.nx, "ny_full": grid.ny,
                    "n_concrete_x": grid.n_concrete_x, "n_concrete_y": grid.n_concrete_y,
                    "nx_concrete": len(grid.x[islice]), "ny_concrete": len(grid.y[jslice]),
                    "dx_m": float(grid.dx),
                    "dx_ft": float(grid.dx) / 0.3048,
                    "model_soil": False, "is_submerged": True,
                }
                print(f"  Grid: {grid.nx}×{grid.ny} (n_cx={grid.n_concrete_x}, "
                      f"n_cy={grid.n_concrete_y}, dx={grid.dx/0.3048:.3f} ft)")

            T0_C = (constr.placement_temp_F - 32.0) * 5.0 / 9.0
            T_soil_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            T_initial = np.full((grid.ny, grid.nx), T0_C)
            T_initial[grid.is_soil] = T_soil_C

            env = make_neutral_env(constr.placement_temp_F)

            t0 = time.perf_counter()
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
            t_wall = time.perf_counter() - t0
            wall_clocks[label] = t_wall
            print(f"  Solved in {t_wall:.1f}s ({result.n_inner_steps} steps, "
                  f"dt={result.dt_inner_s:.1f}s)")

            jslice, islice = grid.concrete_slice()
            T_conc_C = result.T_field_C[:, jslice, islice]
            T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
            ti = nearest_time_idx(result.t_s, 168.0)

            csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
            save_field_csv(csv_path, T_conc_F[ti], grid.x[islice], grid.y[jslice])
            manifest[label] = {"status": "ok", "t_wall_s": float(t_wall)}

        except Exception as e:
            print(f"  ERROR: {e}\n{traceback.format_exc()}")
            manifest[label] = {"status": "error", "error": str(e)}

    with open(os.path.join(ENG_DIR, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    if grid_info:
        with open(os.path.join(ENG_DIR, "grid_info.json"), "w") as f:
            json.dump(grid_info, f, indent=2)

    print("\n--- Computing masked residuals and generating plots ---")
    stats = {}
    for label, folder, placement, soil in RUNS:
        csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
        if not os.path.isfile(csv_path):
            print(f"  SKIP Run {label}: CSV not found")
            continue

        eng_y, eng_x, eng_F = load_engine_csv(csv_path)
        cw_field, cw_widths_m, cw_depths_m, _ = load_cw_slice(folder, COMPARE_HR)
        if cw_field is None:
            print(f"  SKIP Run {label}: CW data missing")
            continue

        eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
        residual = eng_on_cw - cw_field
        s = masked_stats(residual, mask, cw_depths_m, cw_widths_m)
        stats[label] = s

        plot_comparison_5b(label, placement, soil, cw_field, eng_on_cw, residual,
                           cw_depths_m, cw_widths_m, mask)
        status = "PASS" if s["masked_max"] <= GATE_ALL else "FAIL ***"
        print(f"  Run {label}: masked={s['masked_max']:.3f}°F  "
              f"mean={s['masked_mean']:.3f}°F  {status}")

    print("\n--- Gate evaluation (all runs ≤ 0.35°F) ---")
    gate_pass = True
    report_lines = [
        "# STAGE5b — Validation Report (6× Grid Refinement)",
        "",
        "**Engine change:** `build_grid_half_mat` default `grid_refinement=6`  ",
        "**Grid:** `n_concrete_x=121, n_concrete_y=73` (≈0.17 ft × 1.10 ft cells)  ",
        "**Config:** `model_soil=False`, `is_submerged=True`  ",
        "**Comparison time:** t=168 hr  ",
        "**Mask:** Stage 3.5 validity mask (|residual_A| < 0.3°F gate region)  ",
        f"**Gate:** All 9 runs ≤ {GATE_ALL}°F masked max|R| (1% of 35°F ΔT_max spec)",
        "",
        "## Residual Table",
        "",
        "| Run | |ΔT| (°F) | Stage 4b max|R| (°F) | Stage 5b max|R| (°F) | "
        "Improvement (°F) | Gate | Wall-clock (s) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]

    for label, folder, placement, soil in RUNS:
        delta_T = abs(soil - placement)
        s4b = STAGE4B_MASKED_MAX.get(label, float("nan"))
        if label not in stats:
            report_lines.append(
                f"| {label} | {delta_T} | {s4b:.3f} | SKIP | — | — | — |"
            )
            continue
        s = stats[label]
        m5b = s["masked_max"]
        ok = m5b <= GATE_ALL
        if not ok:
            gate_pass = False
        improvement = s4b - m5b
        verdict = "PASS ✓" if ok else "FAIL ✗"
        t_wall = wall_clocks.get(label, float("nan"))
        report_lines.append(
            f"| {label} | {delta_T} | {s4b:.3f} | **{m5b:.3f}** | "
            f"{improvement:+.3f} | **{verdict}** | {t_wall:.1f} |"
        )
        print(f"  Run {label} (|ΔT|={delta_T}°F): "
              f"5b={m5b:.3f}°F  {'PASS' if ok else 'FAIL ***'}")

    overall = "**PASS**" if gate_pass else "**FAIL**"
    report_lines += [
        "",
        f"## Gate Verdict: {overall}",
        "",
    ]
    if gate_pass:
        report_lines += [
            "All 9 runs pass the 0.35°F gate at 6× grid resolution.",
            "Stage 5 (suppressed-hydration calibration) is **CLOSED**.",
            "",
            "Next: 14-mix full-stack revalidation (separate session).",
        ]
    else:
        failing = [
            label for label, *_ in RUNS
            if label in stats and stats[label]["masked_max"] > GATE_ALL
        ]
        report_lines += [
            f"Failing runs: {', '.join(failing)}",
            "",
            "Investigate root cause before proceeding. Possible causes:",
            "- Predicted post-fix residual too optimistic; try 8× refinement",
            "- Mask edge effects at higher resolution",
            "- Stage 5a-undetected parametric issue in specific runs",
        ]

    if grid_info:
        report_lines += [
            "",
            "## Grid Info",
            "",
            f"- `grid_refinement`: {grid_info.get('grid_refinement', 6)}",
            f"- `n_concrete_x`: {grid_info.get('n_concrete_x')} "
            f"(`nx_full`: {grid_info.get('nx_full')})",
            f"- `n_concrete_y`: {grid_info.get('n_concrete_y')} "
            f"(`ny_full`: {grid_info.get('ny_full')})",
            f"- `dx`: {grid_info.get('dx_ft', float('nan')):.4f} ft "
            f"({grid_info.get('dx_m', float('nan')):.4f} m)",
        ]

    report_path = os.path.join(HERE, "STAGE5b_validation.md")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"\nReport saved: {report_path}")
    print(f"Overall gate: {overall}")


if __name__ == "__main__":
    run()
