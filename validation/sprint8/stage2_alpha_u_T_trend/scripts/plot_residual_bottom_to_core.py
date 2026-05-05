#!/usr/bin/env python3
"""§4.11.8 diagnostic — Centerline residual (CW − engine) at t=168 hr.

Bottom surface → center core for A, F, I scenarios.
One figure with 3 panels (one per scenario); each panel has 4 α_u curves.

Output:
  figures/residual_bottom_to_core_v3.png
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
T_FINAL_HR      = 168.0

SCENARIOS = [
    {
        "label":   "A",
        "T_pl_F":  73.0,
        "T_so_F":  73.0,
        "runs": [
            ("α_u=0.20", "thermal_alpha02_A_73_73", 0.20),
            ("α_u=0.40", "thermal_alpha04_A_73_73", 0.40),
            ("α_u=0.60", "thermal_alpha06_A_73_73", 0.60),
            ("α_u=0.80", "thermal_alpha08_A_73_73", 0.80),
        ],
    },
    {
        "label":   "F",
        "T_pl_F":  73.0,
        "T_so_F":  45.0,
        "runs": [
            ("α_u=0.20", "thermal_alpha02_F_73_45", 0.20),
            ("α_u=0.40", "thermal_alpha04_F_73_45", 0.40),
            ("α_u=0.60", "thermal_alpha06_F_73_45", 0.60),
            ("α_u=0.80", "thermal_alpha08_F_73_45", 0.80),
        ],
    },
    {
        "label":   "I",
        "T_pl_F":  100.0,
        "T_so_F":  73.0,
        "runs": [
            ("α_u=0.20", "thermal_alpha02_I_100_73", 0.20),
            ("α_u=0.40", "thermal_alpha04_I_100_73", 0.40),
            ("α_u=0.60", "thermal_alpha06_I_100_73", 0.60),
            ("α_u=0.80", "thermal_alpha08_I_100_73", 0.80),
        ],
    },
]

AU_COLORS  = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]
AU_LABELS  = ["α_u=0.20", "α_u=0.40", "α_u=0.60", "α_u=0.80"]


def load_c_fit():
    with open(DATA / "c_T_pc_fit_chosen.csv") as f:
        rows = list(csv.reader(f))
    info = {r[0]: r[1] for r in rows if len(r) >= 2 and r[0]}
    a0, a1, a2 = json.loads(info["params_json"])
    def c_of_T(T_F):
        T = np.asarray(T_F, dtype=float)
        return float(a0 + a1 * T + a2 * T ** 2)
    return c_of_T


def run_engine(folder, alpha_u_nom, c_factor, T_pl_F, T_so_F):
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))
    mix.Hu_factor_calibrated = 1.0
    mix.Hu_J_kg_effective    = HU_EFF_J_KG
    mix.alpha_u              = alpha_u_nom * c_factor
    constr.model_soil        = False
    constr.is_submerged      = True
    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0 = (T_pl_F - 32.0) * 5.0 / 9.0
    Ts = (T_so_F - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0)
    Ti[grid.is_soil] = Ts
    orig = te2d.K_UC_CALIBRATION_FACTOR_SPRINT7
    try:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = K_UC_FACTOR
        res = solve_hydration_2d(
            grid, mix, Ti,
            duration_s=DURATION_HR * 3600.0,
            output_interval_s=OUTPUT_INTERVAL,
            boundary_mode="full_2d",
            environment=make_neutral_env(T_pl_F),
            construction=constr,
            T_ground_deep_C=Ts,
            diagnostic_outputs=False,
        )
    finally:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = orig
    return res, grid


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    c_of_T = load_c_fit()

    print("§4.11.8 — residual (CW−engine) at t=168 hr, bottom surface → center core")
    print(f"  Hu_eff={HU_EFF_J_KG} J/kg  k_uc={K_UC_FACTOR}")
    print("=" * 100)

    fig, axes = plt.subplots(1, 3, figsize=(14, 7), sharey=False)

    for ax, sc in zip(axes, SCENARIOS):
        T_pl_F = sc["T_pl_F"]
        T_so_F = sc["T_so_F"]
        c_val  = c_of_T(T_pl_F)
        label  = sc["label"]
        print(f"\nScenario {label}  T_pl={T_pl_F}°F  T_so={T_so_F}°F  c={c_val:.4f}")

        ax.axvline(0, color="k", lw=0.8, ls="--", alpha=0.5)

        for (au_label, folder_name, alpha_u_nom), color in zip(sc["runs"], AU_COLORS):
            folder = CW_DATA / folder_name
            res, grid = run_engine(folder, alpha_u_nom, c_val, T_pl_F, T_so_F)
            jsl, isl  = grid.concrete_slice()
            eng_y_m   = grid.y[jsl]               # 0 (top) → max_depth (bottom)
            eng_x_m   = grid.x[isl]
            eng_wi    = len(eng_x_m) - 1           # x=max → symmetric center

            v         = parse_cw_temp_output(str(folder / "output.txt"))
            cw_depths = v.depths_m                 # 0 (top) → max_depth (bottom)
            cw_wi     = int(np.argmax(v.widths_m)) # widths_m[0]=max → mat center

            max_depth = cw_depths[-1]
            half_depth = max_depth / 2.0

            # Engine at t=168 hr — bottom half only
            ti_e      = nearest_time_idx(res.t_s, T_FINAL_HR)
            T_eng_F   = res.T_field_C[ti_e, jsl, isl] * 9.0 / 5.0 + 32.0
            T_eng_cl  = T_eng_F[:, eng_wi]         # shape (n_y,)
            mask_e    = eng_y_m >= half_depth
            eng_y_bot = eng_y_m[mask_e]            # bottom half, depth from top
            T_eng_bot = T_eng_cl[mask_e]

            # CW at t=168 hr — bottom half only
            ti_c      = int(np.abs(v.time_hrs - T_FINAL_HR).argmin())
            T_cw_cl   = v.T_field_F[ti_c][:, cw_wi]
            mask_c    = cw_depths >= half_depth
            cw_y_bot  = cw_depths[mask_c]          # depth from top
            T_cw_bot  = T_cw_cl[mask_c]

            # Interpolate engine onto CW depth grid (both already in bottom half)
            T_eng_on_cw = np.interp(cw_y_bot, eng_y_bot, T_eng_bot)
            residual    = T_cw_bot - T_eng_on_cw  # CW − engine

            # Convert depth-from-top → depth-from-bottom for y-axis
            y_from_bot = max_depth - cw_y_bot      # 0 = bottom surface, ~12 = center

            max_res = float(np.max(np.abs(residual)))
            print(f"  {au_label}  max|CW−eng|={max_res:.4f}°F  "
                  f"res at center={residual[0]:.4f}°F  res at bottom={residual[-1]:.4f}°F")

            ax.plot(residual, y_from_bot, "o-", color=color, lw=1.5, ms=4, label=au_label)

        ax.set_xlabel("CW − engine (°F)", fontsize=11)
        ax.set_ylabel("depth from bottom surface (m)", fontsize=11)
        ax.set_title(
            f"Scenario {label}\n"
            f"T_pl={T_pl_F}°F  T_so={T_so_F}°F  c={c_val:.4f}",
            fontsize=11,
        )
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(True, ls="--", alpha=0.4)
        ax.axhline(half_depth, color="gray", lw=0.7, ls=":", alpha=0.6)  # center core

    fig.suptitle(
        f"Centerline residual (CW − engine) at t=168 hr — bottom surface → center core\n"
        f"(Hu_eff={HU_EFF_J_KG} J/kg, k_uc={K_UC_FACTOR}, quadratic c(T_pl) correction)",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIGS / "residual_bottom_to_core_v3.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\nWrote {out}")
    print("Done.")


if __name__ == "__main__":
    main()
