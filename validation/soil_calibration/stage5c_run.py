#!/usr/bin/env python3
"""Stage 5c — CW-grid alignment via blanket_thickness_m=0.0.

Stage 5b failed for runs F, G, H, I (masked max|R| = 0.387–0.423°F). A 4×→10×
refinement plateau confirmed the residuals are structural (comparison-grid
alignment), not numerical.  Root cause: the engine adds a 0.02 m blanket row
above the concrete, placing the concrete bottom node at y ≈ 24.404 m, while
CW's deepest comparison node sits at 24.380 m (24 mm offset).

Stage 5c passes blanket_thickness_m=0.0 at the CW-comparison call site.  With
zero blanket thickness the blanket row still exists at j=0 for index alignment
(pure-R architecture), but y[1]==y[0]==0, so the concrete top aligns with
CW di=0 and the concrete bottom is at y=depth_m (4 mm from CW di=48 — CW's
own grid precision; residual ≈ 0.067°F, well within gate).

Engine source is not modified beyond the CFL guard added in the separate
prep commit (6a7586b).

Gate: all 9 runs (A–I) masked max|R| ≤ 0.35°F.

Outputs:
    validation/soil_calibration/stage5c_runs/run{A..I}_t168.csv
    validation/soil_calibration/plots/STAGE5c_run{A..I}_engine_vs_cw.png
    validation/soil_calibration/STAGE5c_validation.md
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
ENG_DIR = os.path.join(HERE, "stage5c_runs")
PLOTS = os.path.join(HERE, "plots")
os.makedirs(ENG_DIR, exist_ok=True)
os.makedirs(PLOTS, exist_ok=True)

MASK_PATH = os.path.join(HERE, "STAGE35_validity_mask.npy")

# Stage 5b masked max|R| values for comparison table.
STAGE5B_MASKED_MAX = {
    "A": 0.131, "B": 0.230, "C": 0.250, "D": 0.188,
    "E": 0.302, "F": 0.401, "G": 0.387, "H": 0.401, "I": 0.423,
}

GATE_ALL = 0.35  # °F — 1% of 35°F ΔT_max spec


def plot_comparison_5c(label, placement, soil, cw_field, eng_field, residual,
                       cw_depths_m, cw_widths_m, mask):
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
        f"Stage 5c — Run {label}: placement={placement}°F, soil={soil}°F, "
        f"ΔT={delta_T:+d}°F\n"
        f"model_soil=False, is_submerged=True, blanket=0.0, grid_refinement=6  |  t={COMPARE_HR:.0f} hr\n"
        f"unmasked max|R|={sym:.2f}°F  masked max|R|={masked_sym:.3f}°F  "
        f"gate≤{GATE_ALL}°F  {'PASS' if masked_sym <= GATE_ALL else 'FAIL'}",
        fontsize=8, fontweight="bold",
    )
    panel_data = [
        (cw_plot,  "CW (reference)",               "inferno", (VMIN, VMAX)),
        (eng_plot, "Engine (resampled to CW grid)", "inferno", (VMIN, VMAX)),
        (res_plot, f"Engine−CW (sym={sym:.3f}°F)",  "RdBu_r",  (-sym, sym)),
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

    out = os.path.join(PLOTS, f"STAGE5c_run{label}_engine_vs_cw.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def run():
    if not os.path.exists(MASK_PATH):
        raise FileNotFoundError(f"Stage 3.5 mask not found: {MASK_PATH}")

    mask = np.load(MASK_PATH)
    print(f"Loaded Stage 3.5 mask: shape={mask.shape}, {int(mask.sum())} cells kept")
    print("Stage 5c uses Stage 3.5 base mask only — no extensions.")

    manifest = {}
    grid_info = None
    geom_info = None
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

            # blanket_thickness_m=0.0: concrete bottom aligns with CW di=48.
            grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                                       is_submerged=True, model_soil=False,
                                       blanket_thickness_m=0.0)
            if grid_info is None:
                jslice, islice = grid.concrete_slice()
                grid_info = {
                    "grid_refinement": 6,
                    "blanket_thickness_m": 0.0,
                    "nx_full": grid.nx, "ny_full": grid.ny,
                    "n_concrete_x": grid.n_concrete_x, "n_concrete_y": grid.n_concrete_y,
                    "nx_concrete": len(grid.x[islice]), "ny_concrete": len(grid.y[jslice]),
                    "dx_m": float(grid.dx),
                    "dx_ft": float(grid.dx) / 0.3048,
                    "y0": float(grid.y[0]),
                    "y1": float(grid.y[1]),
                    "y_concrete_end": float(grid.y[grid.iy_concrete_end]),
                    "model_soil": False, "is_submerged": True,
                }
                print(f"  Grid: {grid.nx}×{grid.ny} (n_cx={grid.n_concrete_x}, "
                      f"n_cy={grid.n_concrete_y}, dx={grid.dx/0.3048:.3f} ft)")
                print(f"  §5.2 Geometry: y[0]={grid.y[0]:.6f} y[1]={grid.y[1]:.6f} "
                      f"y[iy_concrete_end]={grid.y[grid.iy_concrete_end]:.6f}")

                # Load CW depths for alignment check
                _, cw_widths_m_ref, cw_depths_m_ref, _ = load_cw_slice(folder, COMPARE_HR)
                geom_info = {
                    "y_concrete_end_m":   float(grid.y[grid.iy_concrete_end]),
                    "cw_depths_m_last":   float(cw_depths_m_ref[-1]),
                    "bottom_offset_m":    float(grid.y[grid.iy_concrete_end] - cw_depths_m_ref[-1]),
                    "y_concrete_top_m":   float(grid.y[1]),
                    "cw_depths_m_first":  float(cw_depths_m_ref[0]),
                    "top_offset_m":       float(grid.y[1] - cw_depths_m_ref[0]),
                }
                print(f"  §5.2 Bottom offset: {geom_info['bottom_offset_m']*1000:.1f} mm "
                      f"(engine {geom_info['y_concrete_end_m']:.4f} m vs CW {geom_info['cw_depths_m_last']:.4f} m)")
                print(f"  §5.2 Top offset:    {geom_info['top_offset_m']*1000:.1f} mm")
                bottom_ok = abs(geom_info["bottom_offset_m"]) < 0.01  # < 10 mm
                top_ok    = abs(geom_info["top_offset_m"])    < 1e-6
                print(f"  §5.2 Alignment: bottom {'OK (<10mm)' if bottom_ok else 'WARN (≥10mm)'}, "
                      f"top {'OK' if top_ok else 'WARN'}")
                if not top_ok:
                    raise RuntimeError("§5.2 STOP: top alignment failed — concrete top y[1] != CW di=0")

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
    lateral_stats = {}  # §5.3 lateral check
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

        # §5.3 lateral alignment check: wi columns 8, 9, 10
        for wi in [8, 9, 10]:
            if wi < residual.shape[1]:
                col_mask = mask[:, wi]
                col_max = float(np.abs(residual[col_mask, wi]).max()) if col_mask.any() else 0.0
                lateral_stats.setdefault(label, {})[f"wi{wi}_masked_maxR"] = col_max

        plot_comparison_5c(label, placement, soil, cw_field, eng_on_cw, residual,
                           cw_depths_m, cw_widths_m, mask)
        status = "PASS" if s["masked_max"] <= GATE_ALL else "FAIL ***"
        print(f"  Run {label}: masked={s['masked_max']:.3f}°F  "
              f"mean={s['masked_mean']:.3f}°F  {status}")

    print("\n--- §5.3 Lateral alignment check (wi=8,9,10 masked max|R|) ---")
    for label in ["F", "G", "H", "I"]:
        if label in lateral_stats:
            parts = [f"wi{wi}={lateral_stats[label].get(f'wi{wi}_masked_maxR', float('nan')):.3f}°F"
                     for wi in [8, 9, 10]]
            print(f"  Run {label}: {', '.join(parts)}")

    print("\n--- Gate evaluation (all runs ≤ 0.35°F) ---")
    gate_pass = True
    report_lines = [
        "# STAGE5c — Validation Report (CW-Grid Alignment)",
        "",
        "**Engine change:** `build_grid_half_mat` called with `blanket_thickness_m=0.0`  ",
        "**Grid:** `n_concrete_x=121, n_concrete_y=73` (≈0.17 ft × 1.10 ft cells)  ",
        "**Config:** `model_soil=False`, `is_submerged=True`, `blanket_thickness_m=0.0`  ",
        "**Comparison time:** t=168 hr  ",
        "**Mask:** Stage 3.5 validity mask (base, no extensions)  ",
        f"**Gate:** All 9 runs ≤ {GATE_ALL}°F masked max|R| (1% of 35°F ΔT_max spec)",
        "",
        "## Residual Table",
        "",
        "| Run | |ΔT| (°F) | Stage 5b max|R| (°F) | Stage 5c max|R| (°F) | "
        "Improvement (°F) | Gate | Wall-clock (s) |",
        "| --- | --- | --- | --- | --- | --- | --- |",
    ]

    for label, folder, placement, soil in RUNS:
        delta_T = abs(soil - placement)
        s5b = STAGE5B_MASKED_MAX.get(label, float("nan"))
        if label not in stats:
            report_lines.append(
                f"| {label} | {delta_T} | {s5b:.3f} | SKIP | — | — | — |"
            )
            continue
        s = stats[label]
        m5c = s["masked_max"]
        ok = m5c <= GATE_ALL
        if not ok:
            gate_pass = False
        improvement = s5b - m5c
        verdict = "PASS ✓" if ok else "FAIL ✗"
        t_wall = wall_clocks.get(label, float("nan"))
        report_lines.append(
            f"| {label} | {delta_T} | {s5b:.3f} | **{m5c:.3f}** | "
            f"{improvement:+.3f} | **{verdict}** | {t_wall:.1f} |"
        )
        print(f"  Run {label} (|ΔT|={delta_T}°F): "
              f"5c={m5c:.3f}°F  {'PASS' if ok else 'FAIL ***'}")

    overall = "**PASS**" if gate_pass else "**FAIL**"
    report_lines += [
        "",
        f"## Gate Verdict: {overall}",
        "",
    ]
    if gate_pass:
        report_lines += [
            "All 9 runs pass the 0.35°F gate with CW-grid alignment at 6× grid resolution.",
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
            "Investigate root cause before proceeding.",
        ]

    if geom_info:
        report_lines += [
            "",
            "## §5.2 Geometry Verification",
            "",
            f"- Engine y[1] (concrete top): {geom_info['y_concrete_top_m']:.6f} m  "
            f"vs CW di=0: {geom_info['cw_depths_m_first']:.6f} m  "
            f"(offset: {geom_info['top_offset_m']*1000:.2f} mm)",
            f"- Engine y[iy_concrete_end]: {geom_info['y_concrete_end_m']:.6f} m  "
            f"vs CW di=48: {geom_info['cw_depths_m_last']:.6f} m  "
            f"(offset: {geom_info['bottom_offset_m']*1000:.1f} mm — CW grid precision)",
            f"- Reduction from Stage 5b: 24 mm → {geom_info['bottom_offset_m']*1000:.1f} mm",
        ]

    if lateral_stats:
        report_lines += [
            "",
            "## §5.3 Lateral Alignment Check (wi=8,9,10)",
            "",
            "| Run | wi=8 masked max|R| (°F) | wi=9 masked max|R| (°F) | wi=10 masked max|R| (°F) |",
            "| --- | --- | --- | --- |",
        ]
        for label in ["F", "G", "H", "I"]:
            if label in lateral_stats:
                parts = [f"{lateral_stats[label].get(f'wi{wi}_masked_maxR', float('nan')):.3f}"
                         for wi in [8, 9, 10]]
                report_lines.append(f"| {label} | {parts[0]} | {parts[1]} | {parts[2]} |")

    if grid_info:
        report_lines += [
            "",
            "## Grid Info",
            "",
            f"- `blanket_thickness_m`: {grid_info.get('blanket_thickness_m')}",
            f"- `grid_refinement`: 6",
            f"- `n_concrete_x`: {grid_info.get('n_concrete_x')} "
            f"(`nx_full`: {grid_info.get('nx_full')})",
            f"- `n_concrete_y`: {grid_info.get('n_concrete_y')} "
            f"(`ny_full`: {grid_info.get('ny_full')})",
            f"- `dx`: {grid_info.get('dx_ft', float('nan')):.4f} ft "
            f"({grid_info.get('dx_m', float('nan')):.4f} m)",
        ]

    report_path = os.path.join(HERE, "STAGE5c_validation.md")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"\nReport saved: {report_path}")
    print(f"Overall gate: {overall}")


if __name__ == "__main__":
    run()
