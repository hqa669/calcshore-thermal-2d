#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-pre — 2D temperature distribution maps.

For each of the 4 α_u targets (F scenario: T_pl=73°F, T_soil=45°F), plots a
2×3 figure: rows = CW / engine, columns = t=0 / 84 / 168 hr.

Outputs (figures/):
  temp_map_alpha020_F_scenario.png
  temp_map_alpha040_F_scenario.png
  temp_map_alpha060_F_scenario.png
  temp_map_alpha080_F_scenario.png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
SC2  = (HERE / "../stage2_calibration").resolve()
ROOT = (HERE / "../../..").resolve()
SC   = (HERE / "../../soil_calibration").resolve()

sys.path.insert(0, str(SC2))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from kinetics_correction import compute_hu_factor
from stage4b_run import make_neutral_env, nearest_time_idx

CW_DATA = SC2 / "cw_data"
FIG_DIR = HERE / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RUNS = [
    ("alpha02_F", 0.20, CW_DATA / "thermal_alpha02_F_73_45"),
    ("alpha04_F", 0.40, CW_DATA / "thermal_alpha04_F_73_45"),
    ("alpha06_F", 0.60, CW_DATA / "thermal_alpha06_F_73_45"),
    ("alpha08_F", 0.80, CW_DATA / "thermal_alpha08_F_73_45"),
]

T_PL_F  = 73.0
T_SO_F  = 45.0
FACTOR  = 0.96
T_TARGETS_HR = [0.0, 84.0, 168.0]


def solve_engine(folder, T_pl_F, T_so_F, factor=0.96):
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))
    fac, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = fac
    mix.Hu_J_kg_effective    = mix.Hu_J_kg * fac
    constr.model_soil   = False
    constr.is_submerged = True
    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0 = (T_pl_F - 32.0) * 5.0 / 9.0
    Ts = (T_so_F - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0)
    Ti[grid.is_soil] = Ts
    original = te2d.K_UC_CALIBRATION_FACTOR_SPRINT7
    try:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = factor
        res = solve_hydration_2d(
            grid, mix, Ti,
            duration_s=168.0 * 3600.0,
            output_interval_s=1800.0,
            boundary_mode="full_2d",
            environment=make_neutral_env(T_pl_F),
            construction=constr,
            T_ground_deep_C=Ts,
            diagnostic_outputs=False,
        )
    finally:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = original
    jsl, isl = grid.concrete_slice()
    eng_x = grid.x[isl]   # widths  (m)
    eng_y = grid.y[jsl]   # depths  (m)
    return res, eng_x, eng_y, jsl, isl


def eng_field_at(res, jsl, isl, target_hr):
    ti = nearest_time_idx(res.t_s, target_hr)
    return res.T_field_C[ti, jsl, isl] * 9.0 / 5.0 + 32.0


def cw_field_at(v, target_hr):
    ti = int(np.abs(v.time_hrs - target_hr).argmin())
    return v.T_field_F[ti]


def make_figure(label, alpha_u, folder):
    print(f"  Solving engine for {label} (α_u={alpha_u:.2f})...", flush=True)
    res, eng_x, eng_y, jsl, isl = solve_engine(folder, T_PL_F, T_SO_F, FACTOR)
    v = parse_cw_temp_output(str(folder / "output.txt"))

    eng_fields = [eng_field_at(res, jsl, isl, t) for t in T_TARGETS_HR]
    cw_fields  = [cw_field_at(v, t) for t in T_TARGETS_HR]

    # shared color limits across all 6 panels
    all_vals = np.concatenate([f.ravel() for f in eng_fields + cw_fields])
    vmin, vmax = float(np.min(all_vals)), float(np.max(all_vals))

    cmap = "RdYlBu_r"
    fig, axes = plt.subplots(2, 3, figsize=(13, 9))

    # Both grids: flip width axis so side face is at x=0 (left),
    # centerline at x=max (right) — consistent orientation for both sources.
    # CW: widths_m[0]=6.1 (center) → flip [:, ::-1] → col 0 = side (x=0)
    # Engine: eng_x[0]=0 (center) → flip [:, ::-1] → col 0 = side (x=max)
    cw_w_max  = float(v.widths_m[0])    # = 6.1 m (half-width)
    eng_w_max = float(eng_x[-1])        # = half-width (same geometry)
    cw_ext    = [0.0, cw_w_max,  float(v.depths_m[-1]), 0.0]
    eng_ext   = [0.0, eng_w_max, float(eng_y[-1]),      0.0]

    cw_t0_label = float(v.time_hrs[0])  # CW has no exact t=0; first output ~5min

    for col, t_hr in enumerate(T_TARGETS_HR):
        # CW
        ax = axes[0, col]
        cw_t_actual = float(v.time_hrs[int(np.abs(v.time_hrs - t_hr).argmin())])
        cw_label = (f"CW  t={cw_t0_label:.2f} hr (first output)"
                    if t_hr == 0.0 else f"CW  t={t_hr:.0f} hr")
        im = ax.imshow(cw_fields[col][:, ::-1], aspect="auto", origin="upper",
                       extent=cw_ext, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_title(cw_label, fontsize=9)
        ax.set_xlabel("Width from side (m)", fontsize=8)
        ax.set_ylabel("Depth (m)", fontsize=8)
        ax.tick_params(labelsize=7)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=6)

        # engine
        ax = axes[1, col]
        im = ax.imshow(eng_fields[col][:, ::-1], aspect="auto", origin="upper",
                       extent=eng_ext, cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        ax.set_title(f"engine  t={t_hr:.0f} hr", fontsize=9)
        ax.set_xlabel("Width from side (m)", fontsize=8)
        ax.set_ylabel("Depth (m)", fontsize=8)
        ax.tick_params(labelsize=7)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=6)

    fig.suptitle(
        f"Temperature distribution — F scenario (T_pl=73°F, T_soil=45°F)  α_u={alpha_u:.2f}\n"
        f"Color scale: {vmin:.1f}–{vmax:.1f}°F  (shared across all panels)",
        fontsize=10, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])

    fname = f"temp_map_alpha{int(alpha_u*100):03d}_F_scenario.png"
    path  = FIG_DIR / fname
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Wrote: {path}")
    return path


def main():
    print(f"\nGenerating 2D temperature maps (F scenario, t=0/84/168 hr)...")
    print(f"4 engine solves × ~60s each ≈ 4 min\n")
    for label, alpha_u, folder in RUNS:
        make_figure(label, alpha_u, folder)
    print("\nDone.")


if __name__ == "__main__":
    main()
