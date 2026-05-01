#!/usr/bin/env python3
"""Stage 4b — Run 9 engine simulations with model_soil=False, is_submerged=True.

This is the CW-matching configuration: Dirichlet T_soil applied directly at
the concrete face, no soil mesh. Validates that the structural fix (Stage 4b)
closes the 9–19°F bottom-CL residuals found in Stage 3.5.

Validation target:
  Run A masked max|R| ≤ 0.5°F
  Runs B–I masked max|R| ≤ 2.0°F

Outputs:
    validation/soil_calibration/stage4b_runs/run{A..I}_t168.csv
    validation/soil_calibration/stage4b_runs/manifest.json
    validation/soil_calibration/stage4b_runs/grid_info.json
    validation/soil_calibration/plots/STAGE4b_run{A..I}_engine_vs_cw.png
    validation/soil_calibration/STAGE4b_validation.md
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

from cw_scenario_loader import parse_cw_dat, CWMixDesign, CWGeometry, CWConstruction, CWEnvironment
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat

# Reuse IO + resampling helpers from Stage 3.
from stage3_compare import (
    load_engine_csv,
    load_cw_slice,
    resample_engine_to_cw,
    COMPARE_HR,
    VMIN,
    VMAX,
)

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
ENG_DIR = os.path.join(HERE, "stage4b_runs")
PLOTS = os.path.join(HERE, "plots")
os.makedirs(ENG_DIR, exist_ok=True)
os.makedirs(PLOTS, exist_ok=True)

MASK_PATH = os.path.join(HERE, "STAGE35_validity_mask.npy")

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

# Stage 3.5 residuals for comparison table.
STAGE35_MASKED_MAX = {
    "A": 0.235, "B": 9.13, "C": 11.71, "D": 9.13,
    "E": 11.71, "F": 18.97, "G": 17.42, "H": 19.566, "I": 17.42,
}

GATE_A = 0.5    # °F
GATE_BI = 2.0   # °F


def make_neutral_env(placement_temp_F):
    n_hrs = 169
    hours = np.arange(n_hrs, dtype=float)
    T_air = np.full(n_hrs, placement_temp_F)
    T_air_C = (T_air - 32.0) * 5.0 / 9.0
    return CWEnvironment(
        hours=hours,
        T_air_F=T_air,
        RH_pct=np.full(n_hrs, 60.0),
        solar_W_m2=np.zeros(n_hrs),
        wind_m_s=np.zeros(n_hrs),
        cloud_cover=np.ones(n_hrs),
        pressure_hPa=np.full(n_hrs, 1013.25),
        T_dp_C=T_air_C - 10.0,
        T_sky_C=T_air_C - 5.0,
        daily_max_F=[placement_temp_F] * 8,
        daily_min_F=[placement_temp_F] * 8,
        cw_ave_max_daily_temp_F=placement_temp_F,
        cw_ave_min_daily_temp_F=placement_temp_F,
        cw_ave_max_solar_W_m2=0.0,
        cw_ave_max_wind_m_s=0.0,
        cw_ave_max_RH_pct=60.0,
        cw_ave_min_RH_pct=60.0,
    )


def nearest_time_idx(t_s, target_hr):
    return int(np.abs(np.asarray(t_s) - target_hr * 3600.0).argmin())


def save_field_csv(path, field_2d_F, x_m, y_m):
    ny, nx = field_2d_F.shape
    with open(path, "w") as f:
        f.write("depth_m\\width_m," + ",".join(f"{x:.4f}" for x in x_m) + "\n")
        for iy in range(ny):
            row = f"{y_m[iy]:.4f}," + ",".join(f"{field_2d_F[iy, ix]:.4f}" for ix in range(nx))
            f.write(row + "\n")


def plot_comparison(label, placement, soil, cw_field, eng_field, residual,
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
        f"Stage 4b — Run {label}: placement={placement}°F, soil={soil}°F, "
        f"ΔT={delta_T:+d}°F\n"
        f"model_soil=False, is_submerged=True  |  t={COMPARE_HR:.0f} hr\n"
        f"unmasked max|R|={sym:.2f}°F  masked max|R|={masked_sym:.2f}°F",
        fontsize=9, fontweight="bold",
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
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    out = os.path.join(PLOTS, f"STAGE4b_run{label}_engine_vs_cw.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def masked_stats(residual, mask, cw_depths_m, cw_widths_m):
    abs_R = np.abs(residual)
    masked_abs = abs_R[mask]
    masked_max = float(masked_abs.max())
    masked_mean = float(masked_abs.mean())
    masked_full = np.where(mask, abs_R, -np.inf)
    flat_idx = int(np.argmax(masked_full))
    di, wi = np.unravel_index(flat_idx, residual.shape)
    return {
        "unmasked_max": float(abs_R.max()),
        "masked_max": masked_max,
        "masked_mean": masked_mean,
        "max_depth_ft": float(cw_depths_m[di]) * 3.28084,
        "max_width_ft": float(cw_widths_m[wi]) * 3.28084,
    }


def run():
    if not os.path.exists(MASK_PATH):
        raise FileNotFoundError(f"Stage 3.5 mask not found: {MASK_PATH}. Run stage35_revalidate.py first.")

    mask = np.load(MASK_PATH)
    print(f"Loaded Stage 3.5 mask: shape={mask.shape}, {int(mask.sum())} cells kept")

    manifest = {}
    grid_info = None

    for label, folder, placement, soil in RUNS:
        print(f"\nRun {label} ({folder}): placement={placement}°F, soil={soil}°F")
        dat_path = os.path.join(CW_RUNS, folder, "input.dat")
        if not os.path.isfile(dat_path):
            print(f"  SKIP: input.dat not found")
            manifest[label] = {"status": "skip"}
            continue

        try:
            mix, geom, constr, raw = parse_cw_dat(dat_path)

            from kinetics_correction import compute_hu_factor
            factor, note = compute_hu_factor(mix)
            mix.Hu_factor_calibrated = factor
            mix.Hu_J_kg_effective = mix.Hu_J_kg * factor

            print(f"  Hu_J_kg_effective={mix.Hu_J_kg_effective:.4f} "
                  f"({'suppressed' if mix.Hu_J_kg_effective < 10 else 'ACTIVE!'})")

            # Stage 4b config: model_soil=False, is_submerged=True
            constr.model_soil = False
            constr.is_submerged = True

            grid = build_grid_half_mat(
                geom.width_ft, geom.depth_ft,
                is_submerged=True,
                model_soil=False,
            )
            if grid_info is None:
                jslice, islice = grid.concrete_slice()
                grid_info = {
                    "nx_full": grid.nx, "ny_full": grid.ny,
                    "nx_concrete": len(grid.x[islice]),
                    "ny_concrete": len(grid.y[jslice]),
                    "x_concrete_m": grid.x[islice].tolist(),
                    "y_concrete_m": grid.y[jslice].tolist(),
                    "model_soil": False, "is_submerged": True,
                }
                print(f"  Grid: ny={grid.ny} (no soil rows), n_soil_y={grid.n_soil_y}")

            T0_C = (constr.placement_temp_F - 32.0) * 5.0 / 9.0
            T_soil_C = (constr.soil_temp_F - 32.0) * 5.0 / 9.0
            T_initial = np.full((grid.ny, grid.nx), T0_C)
            # No soil cells: T_initial[grid.is_soil] is a no-op but harmless.
            T_initial[grid.is_soil] = T_soil_C

            env = make_neutral_env(constr.placement_temp_F)
            T_ground_deep_C = T_soil_C

            t0 = time.perf_counter()
            result = solve_hydration_2d(
                grid, mix, T_initial,
                duration_s=168 * 3600,
                output_interval_s=1800.0,
                boundary_mode="full_2d",
                environment=env,
                construction=constr,
                T_ground_deep_C=T_ground_deep_C,
                diagnostic_outputs=False,
            )
            t_wall = time.perf_counter() - t0
            print(f"  Solved in {t_wall:.1f}s")

            jslice, islice = grid.concrete_slice()
            T_conc_C = result.T_field_C[:, jslice, islice]
            T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
            y_conc = grid.y[jslice]
            x_conc = grid.x[islice]

            ti = nearest_time_idx(result.t_s, 168.0)
            csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
            save_field_csv(csv_path, T_conc_F[ti], x_conc, y_conc)
            print(f"  Saved: {csv_path} (T range: {T_conc_F[ti].min():.1f}–{T_conc_F[ti].max():.1f}°F)")

            manifest[label] = {"status": "ok", "t_wall_s": float(t_wall)}

        except Exception as e:
            tb = traceback.format_exc()
            print(f"  ERROR: {e}\n{tb}")
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

        plot_comparison(label, placement, soil, cw_field, eng_on_cw, residual,
                        cw_depths_m, cw_widths_m, mask)
        print(f"  Run {label}: unmasked={s['unmasked_max']:.2f}°F  "
              f"masked={s['masked_max']:.3f}°F  mean={s['masked_mean']:.3f}°F")

    print("\n--- Gate evaluation ---")
    gate_pass = True
    report_lines = [
        "# STAGE4b — Validation Report (model_soil=False, is_submerged=True)",
        "",
        "**Config:** `model_soil=False`, `is_submerged=True`  ",
        "**Comparison time:** t=168 hr  ",
        "**Mask:** Stage 3.5 Run-A validity mask (|residual_A| < 0.3°F)  ",
        "**Gate:** Run A ≤ 0.5°F, Runs B–I ≤ 2.0°F (masked max|R|)",
        "",
        "## Residual Table",
        "",
        "| Run | Placement/Soil (°F) | Stage 3.5 masked max|R| (°F) | "
        "Stage 4b masked max|R| (°F) | Improvement | Gate |",
        "| --- | --- | --- | --- | --- | --- |",
    ]

    for label, folder, placement, soil in RUNS:
        if label not in stats:
            report_lines.append(f"| {label} | {placement}/{soil} | — | SKIP | — | — |")
            continue
        s = stats[label]
        s35 = STAGE35_MASKED_MAX.get(label, float("nan"))
        m4b = s["masked_max"]
        gate = GATE_A if label == "A" else GATE_BI
        ok = m4b <= gate
        if not ok:
            gate_pass = False
        improvement = s35 - m4b if not np.isnan(s35) else float("nan")
        verdict = "PASS" if ok else "FAIL"
        report_lines.append(
            f"| {label} | {placement}/{soil} | {s35:.3f} | **{m4b:.3f}** | "
            f"{improvement:+.3f} | **{verdict}** |"
        )
        print(f"  Run {label}: masked={m4b:.3f}°F  gate≤{gate}°F  {verdict}")

    overall = "**PASS**" if gate_pass else "**FAIL**"
    report_lines += [
        "",
        f"## Gate Verdict: {overall}",
        "",
    ]
    if gate_pass:
        report_lines.append(
            "All runs pass the Stage 4b validation gate. "
            "The model_soil=False, is_submerged=True configuration matches CW "
            "within the target tolerance."
        )
    else:
        report_lines.append(
            "One or more runs fail the gate. See masked max|R| values above "
            "for the location of remaining residuals."
        )

    report_path = os.path.join(HERE, "STAGE4b_validation.md")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines) + "\n")
    print(f"\nReport saved: {report_path}")
    print(f"Overall gate: {overall}")


if __name__ == "__main__":
    run()
