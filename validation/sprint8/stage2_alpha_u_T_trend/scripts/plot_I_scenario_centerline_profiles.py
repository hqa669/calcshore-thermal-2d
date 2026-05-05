#!/usr/bin/env python3
"""§4.11.8 diagnostic — Centerline T(z) profiles, CW vs engine(c(T_pl)), I-scenarios.

For each of the 4 I-scenario datasets (T_pl=100°F, T_soil=73°F at α_u ∈ {0.20, 0.40,
0.60, 0.80}), plots the centerline (wi=0) temperature profile from bottom surface up
through the half-mat at t = 0, 84, 168 hr. Engine is run with the §4.11.8 chosen
quadratic c(T_pl) correction (c(100°F) = 1.3623), Hu_eff = 12,936.96 J/kg,
k_uc factor = 0.96.

Output:
  figures/I_centerline_profiles_v3.png   — 4×3 grid (4 α_u rows × 3 time columns)
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE  = Path(__file__).resolve().parents[1]
DATA  = HERE / "data"
FIGS  = HERE / "figures"
REPO  = Path(__file__).resolve().parents[4]
STAGE = REPO / "validation" / "sprint8" / "stage2_calibration"
SC    = REPO / "validation" / "soil_calibration"
DIAG  = REPO / "validation" / "sprint8" / "stage2_diag_Hu"

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STAGE))
sys.path.insert(0, str(SC))
sys.path.insert(0, str(DIAG / "scripts"))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from stage4b_run import make_neutral_env, nearest_time_idx

CW_DATA = STAGE / "cw_data"

HU_EFF_J_KG     = 12936.96
K_UC_FACTOR     = 0.96
DURATION_HR     = 168.0
OUTPUT_INTERVAL = 1800.0
SNAPSHOT_HRS    = [0.0, 84.0, 168.0]

I_RUNS = [
    ("α_u=0.20",  "thermal_alpha02_I_100_73", 0.20),
    ("α_u=0.40",  "thermal_alpha04_I_100_73", 0.40),
    ("α_u=0.60",  "thermal_alpha06_I_100_73", 0.60),
    ("α_u=0.80",  "thermal_alpha08_I_100_73", 0.80),
]
T_PL = 100.0
T_SO = 73.0


def load_chosen_fit():
    with open(DATA / "c_T_pc_fit_chosen.csv") as f:
        rows = list(csv.reader(f))
    info = {r[0]: r[1] for r in rows if len(r) >= 2 and r[0]}
    a0, a1, a2 = json.loads(info["params_json"])
    def c_of_T(T_F):
        T = np.asarray(T_F, dtype=float)
        return a0 + a1 * T + a2 * T ** 2
    return float(c_of_T(T_PL)), c_of_T


def run_engine(folder, alpha_u_nom, c_factor):
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))
    mix.Hu_factor_calibrated = 1.0
    mix.Hu_J_kg_effective    = HU_EFF_J_KG
    mix.alpha_u              = alpha_u_nom * c_factor

    constr.model_soil   = False
    constr.is_submerged = True

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0 = (T_PL - 32.0) * 5.0 / 9.0
    Ts = (T_SO - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0)
    Ti[grid.is_soil] = Ts

    original = te2d.K_UC_CALIBRATION_FACTOR_SPRINT7
    try:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = K_UC_FACTOR
        res = solve_hydration_2d(
            grid, mix, Ti,
            duration_s=DURATION_HR * 3600.0,
            output_interval_s=OUTPUT_INTERVAL,
            boundary_mode="full_2d",
            environment=make_neutral_env(T_PL),
            construction=constr,
            T_ground_deep_C=Ts,
            diagnostic_outputs=False,
        )
    finally:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = original

    return res, grid


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    c_at_100, _ = load_chosen_fit()
    print(f"§4.11.8 diagnostic — I-scenario centerline T(z) profiles")
    print(f"  c(T_pl=100°F) = {c_at_100:.4f}, Hu_eff = {HU_EFF_J_KG} J/kg, k_uc = {K_UC_FACTOR}")
    print(f"  Snapshots at t = {SNAPSHOT_HRS} hr")
    print("=" * 100)

    fig, axes = plt.subplots(4, 3, figsize=(13, 12), sharey=True)
    color_cw  = "#2166ac"
    color_eng = "#d62728"

    for row, (au_label, folder_name, alpha_u_nom) in enumerate(I_RUNS):
        folder = CW_DATA / folder_name
        print(f"\n  {au_label}  ({folder_name})  α_u_nom={alpha_u_nom}, "
              f"α_u_eng={alpha_u_nom * c_at_100:.4f}")
        res, grid = run_engine(folder, alpha_u_nom, c_at_100)
        jsl, isl = grid.concrete_slice()
        # Engine concrete grid
        eng_y_m = grid.y[jsl]                # depth coords (m), top→bottom
        eng_x_m = grid.x[isl]                # width coords (m)
        # Centerline in engine: x=max (largest x = symmetric center of mat)
        eng_wi = len(eng_x_m) - 1

        # CW data
        v = parse_cw_temp_output(str(folder / "output.txt"))
        cw_depths_m = v.depths_m            # 49 depths, 0..z_max
        cw_widths_m = v.widths_m
        cw_wi = int(np.argmax(cw_widths_m))            # centerline (CW x=max → mat center)

        for col, t_hr in enumerate(SNAPSHOT_HRS):
            ax = axes[row, col]
            # Engine snapshot
            ti_e = nearest_time_idx(res.t_s, t_hr)
            T_eng_C = res.T_field_C[ti_e, jsl, isl]   # (n_depth, n_width)
            T_eng_F = T_eng_C * 9.0 / 5.0 + 32.0
            T_eng_centerline = T_eng_F[:, eng_wi]      # (n_depth,)

            # CW snapshot
            ti_c = int(np.abs(v.time_hrs - t_hr).argmin())
            T_cw_F = v.T_field_F[ti_c]                 # (n_depth, n_width)
            T_cw_centerline = T_cw_F[:, cw_wi]         # (49,)

            # Plot: depth on y (downward positive), T on x.
            # CW: depths_m go from 0 (top) to bottom. We want bottom at top of plot? No —
            # convention: depth on y axis (0=top at top of plot, increasing downward).
            ax.plot(T_cw_centerline, cw_depths_m, "o-",
                    color=color_cw, lw=1.5, markersize=4, label="CW")
            ax.plot(T_eng_centerline, eng_y_m, "s-",
                    color=color_eng, lw=1.5, markersize=3, label="Engine c(T_pl)")
            ax.invert_yaxis()                          # depth grows downward
            ax.set_title(f"{au_label}, t = {t_hr:.0f} hr", fontsize=10)
            if col == 0:
                ax.set_ylabel("depth from top (m)", fontsize=10)
            if row == 3:
                ax.set_xlabel("T (°F)", fontsize=10)
            ax.grid(True, ls="--", alpha=0.4)
            if row == 0 and col == 0:
                ax.legend(fontsize=9, loc="lower right")

            # Quick stats: max diff
            T_eng_on_cw = np.interp(cw_depths_m, eng_y_m, T_eng_centerline)
            max_diff = float(np.max(np.abs(T_eng_on_cw - T_cw_centerline)))
            print(f"    t={t_hr:>5.0f} hr  CW range [{T_cw_centerline.min():.2f}, "
                  f"{T_cw_centerline.max():.2f}]°F  "
                  f"max|eng-CW| at centerline = {max_diff:.4f}°F")

    fig.suptitle(
        f"I-scenario centerline T(z) profiles — CW vs engine with c(T_pl=100°F)={c_at_100:.4f}\n"
        f"(Hu_eff = {HU_EFF_J_KG} J/kg, k_uc = {K_UC_FACTOR}, half-mat center x=max)",
        fontsize=11
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIGS / "I_centerline_profiles_v3.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nWrote {out}")
    print("Done.")


if __name__ == "__main__":
    main()
