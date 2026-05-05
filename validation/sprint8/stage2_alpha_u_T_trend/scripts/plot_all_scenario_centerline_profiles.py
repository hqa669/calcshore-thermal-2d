#!/usr/bin/env python3
"""§4.11.8 diagnostic — Centerline T(z) profiles, CW vs engine(c(T_pl)), all scenarios.

Produces one 4×3 figure (4 α_u rows × 3 time columns) per scenario type (A, F, I).
Engine is run with the §4.11.8 chosen quadratic c(T_pl) correction,
Hu_eff = 12,936.96 J/kg, k_uc factor = 0.96.

Output:
  figures/A_centerline_profiles_v3.png  — T_pl=73°F, T_so=73°F
  figures/F_centerline_profiles_v3.png  — T_pl=73°F, T_so=45°F
  figures/I_centerline_profiles_v3.png  — T_pl=100°F, T_so=73°F
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

SCENARIOS = [
    {
        "label":     "A",
        "T_pl_F":    73.0,
        "T_so_F":    73.0,
        "runs": [
            ("α_u=0.20", "thermal_alpha02_A_73_73", 0.20),
            ("α_u=0.40", "thermal_alpha04_A_73_73", 0.40),
            ("α_u=0.60", "thermal_alpha06_A_73_73", 0.60),
            ("α_u=0.80", "thermal_alpha08_A_73_73", 0.80),
        ],
        "out": "A_centerline_profiles_v3.png",
    },
    {
        "label":     "F",
        "T_pl_F":    73.0,
        "T_so_F":    45.0,
        "runs": [
            ("α_u=0.20", "thermal_alpha02_F_73_45", 0.20),
            ("α_u=0.40", "thermal_alpha04_F_73_45", 0.40),
            ("α_u=0.60", "thermal_alpha06_F_73_45", 0.60),
            ("α_u=0.80", "thermal_alpha08_F_73_45", 0.80),
        ],
        "out": "F_centerline_profiles_v3.png",
    },
    {
        "label":     "I",
        "T_pl_F":    100.0,
        "T_so_F":    73.0,
        "runs": [
            ("α_u=0.20", "thermal_alpha02_I_100_73", 0.20),
            ("α_u=0.40", "thermal_alpha04_I_100_73", 0.40),
            ("α_u=0.60", "thermal_alpha06_I_100_73", 0.60),
            ("α_u=0.80", "thermal_alpha08_I_100_73", 0.80),
        ],
        "out": "I_centerline_profiles_v3.png",
    },
]


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
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = original

    return res, grid


def plot_scenario(sc, c_of_T):
    T_pl_F  = sc["T_pl_F"]
    T_so_F  = sc["T_so_F"]
    c_val   = c_of_T(T_pl_F)
    label   = sc["label"]

    print(f"\n{'='*100}")
    print(f"Scenario {label}  T_pl={T_pl_F}°F  T_so={T_so_F}°F  c(T_pl)={c_val:.4f}")
    print(f"{'='*100}")

    fig, axes = plt.subplots(4, 3, figsize=(13, 12), sharey=True)
    color_cw  = "#2166ac"
    color_eng = "#d62728"

    for row, (au_label, folder_name, alpha_u_nom) in enumerate(sc["runs"]):
        folder = CW_DATA / folder_name
        alpha_u_eng = alpha_u_nom * c_val
        print(f"\n  {au_label}  ({folder_name})  α_u_eng={alpha_u_eng:.4f}")
        res, grid = run_engine(folder, alpha_u_nom, c_val, T_pl_F, T_so_F)
        jsl, isl  = grid.concrete_slice()
        eng_y_m   = grid.y[jsl]
        eng_x_m   = grid.x[isl]
        eng_wi    = len(eng_x_m) - 1   # x=max → symmetric center of mat

        v         = parse_cw_temp_output(str(folder / "output.txt"))
        cw_depths = v.depths_m
        cw_wi     = int(np.argmax(v.widths_m))  # widths_m[0]=max → mat center

        for col, t_hr in enumerate(SNAPSHOT_HRS):
            ax = axes[row, col]

            ti_e = nearest_time_idx(res.t_s, t_hr)
            T_eng_F = res.T_field_C[ti_e, jsl, isl] * 9.0 / 5.0 + 32.0
            T_eng_cl = T_eng_F[:, eng_wi]

            ti_c = int(np.abs(v.time_hrs - t_hr).argmin())
            T_cw_cl = v.T_field_F[ti_c][:, cw_wi]

            ax.plot(T_cw_cl,  cw_depths, "o-", color=color_cw,  lw=1.5, ms=4, label="CW")
            ax.plot(T_eng_cl, eng_y_m,   "s-", color=color_eng, lw=1.5, ms=3, label="Engine c(T_pl)")
            ax.invert_yaxis()
            ax.set_title(f"{au_label}, t = {t_hr:.0f} hr", fontsize=10)
            if col == 0:
                ax.set_ylabel("depth from top (m)", fontsize=10)
            if row == 3:
                ax.set_xlabel("T (°F)", fontsize=10)
            ax.grid(True, ls="--", alpha=0.4)
            if row == 0 and col == 0:
                ax.legend(fontsize=9, loc="lower right")

            T_eng_on_cw = np.interp(cw_depths, eng_y_m, T_eng_cl)
            max_diff = float(np.max(np.abs(T_eng_on_cw - T_cw_cl)))
            print(f"    t={t_hr:>5.0f} hr  CW [{T_cw_cl.min():.2f}, {T_cw_cl.max():.2f}]°F  "
                  f"max|eng-CW|={max_diff:.4f}°F")

    fig.suptitle(
        f"Scenario {label} centerline T(z) — CW vs engine  "
        f"T_pl={T_pl_F}°F  T_so={T_so_F}°F  c(T_pl)={c_val:.4f}\n"
        f"(Hu_eff={HU_EFF_J_KG} J/kg, k_uc={K_UC_FACTOR})",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out = FIGS / sc["out"]
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"\n  → Wrote {out}")


def main():
    FIGS.mkdir(parents=True, exist_ok=True)
    c_of_T = load_c_fit()
    print(f"§4.11.8 diagnostic — centerline T(z) profiles, all scenarios")
    print(f"  Hu_eff = {HU_EFF_J_KG} J/kg, k_uc = {K_UC_FACTOR}")
    print(f"  Snapshots at t = {SNAPSHOT_HRS} hr")
    for sc in SCENARIOS:
        plot_scenario(sc, c_of_T)
    print("\nDone.")


if __name__ == "__main__":
    main()
