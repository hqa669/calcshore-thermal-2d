#!/usr/bin/env python3
"""k×0.96 — All 9 runs residual profiles at t=168 hr.

Generates 9×4 figure: for each run, Side T, Side residual, Bottom T,
Bottom residual.  Prints full residual tables for R1 and R2.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SC   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, SC)

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

CW_RUNS = os.path.join(SC, "cw_runs")
OUT_DIR = os.path.join(SC, "diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

K_FACTOR   = 0.96
COMPARE_HR = 168
GATE_F     = 0.35
R1_DI      = 24
R2_WI      = 0

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


def run_engine(folder, pl_F, so_F):
    mix, geom, constr, _ = parse_cw_dat(os.path.join(CW_RUNS, folder, "input.dat"))
    fac, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = fac
    mix.Hu_J_kg_effective = mix.Hu_J_kg * fac
    mix.thermal_conductivity_BTU_hr_ft_F *= K_FACTOR
    constr.model_soil = False; constr.is_submerged = True
    grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                               is_submerged=True, model_soil=False,
                               blanket_thickness_m=0.0)
    T0 = (pl_F-32)*5/9; Ts = (so_F-32)*5/9
    Ti = np.full((grid.ny, grid.nx), T0); Ti[grid.is_soil] = Ts
    res = solve_hydration_2d(grid, mix, Ti, duration_s=COMPARE_HR*3600,
                             output_interval_s=1800., boundary_mode="full_2d",
                             environment=make_neutral_env(pl_F), construction=constr,
                             T_ground_deep_C=Ts, diagnostic_outputs=False)
    jsl, isl = grid.concrete_slice()
    ti = nearest_time_idx(res.t_s, float(COMPARE_HR))
    TF = res.T_field_C[ti, jsl, isl]*9/5+32
    return grid.y[jsl], grid.x[isl], TF


def load_cw(folder):
    v = parse_cw_temp_output(os.path.join(CW_RUNS, folder, "output.txt"))
    ti = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m


def main():
    # CW reference grid
    _, cw_widths_m, cw_depths_m = load_cw("runA_baseline")
    di_r2   = np.arange(24, 49)
    z_ft    = -cw_depths_m[di_r2] / 0.3048          # −40 to −80 ft
    x_r1_ft = (cw_widths_m[0] - cw_widths_m) / 0.3048  # 0→20 ft (CL→face)
    wi9_x   = x_r1_ft[9]
    di45_z  = -cw_depths_m[45] / 0.3048

    print(f"Running all 9 engines with k×{K_FACTOR} …\n")

    eng_data = {}
    cw_data  = {}
    for label, folder, pl, so in RUNS:
        print(f"  Run {label} …", end=" ", flush=True)
        ey, ex, eT  = run_engine(folder, pl, so)
        cw_f, _, _  = load_cw(folder)
        eng_on_cw   = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)
        eng_data[label] = eng_on_cw
        cw_data[label]  = cw_f
        R = eng_on_cw - cw_f
        r1 = float(np.max(np.abs(R[R1_DI, :])))
        r2 = float(np.max(np.abs(R[24:48, R2_WI])))
        gate = "PASS" if (r1 <= GATE_F and r2 <= GATE_F) else "FAIL"
        print(f"R1={r1:.4f}  R2={r2:.4f}  [{gate}]")

    # ── Residual tables ──
    print(f"\n{'='*90}")
    print(f"R1 — Side profile (di=24) residual, Engine−CW (°F)  [k×{K_FACTOR}]")
    print(f"{'='*90}")
    wi_idx = list(range(13))
    print(f"  {'Run':>4} | " + "  ".join(f"wi={w:>2}" for w in wi_idx))
    print("  " + "-"*(7 + 8*13))
    for label, *_ in RUNS:
        R = eng_data[label] - cw_data[label]
        vals = "  ".join(f"{R[R1_DI, w]:>+6.3f}" for w in wi_idx)
        r1_max = np.max(np.abs(R[R1_DI, :]))
        print(f"  {label:>4} | {vals}   max={r1_max:.4f}")

    print(f"\n{'='*90}")
    print(f"R2 — Bottom profile (wi=0) residual, Engine−CW (°F)  [k×{K_FACTOR}]")
    print(f"{'='*90}")
    sample_di = di_r2[::2]   # every other row (13 values)
    print(f"  {'Run':>4} | " + "  ".join(f"di={d:>2}" for d in sample_di))
    print("  " + "-"*(7 + 8*len(sample_di)))
    for label, *_ in RUNS:
        R = eng_data[label] - cw_data[label]
        vals = "  ".join(f"{R[d, R2_WI]:>+6.3f}" for d in sample_di)
        r2_max = np.max(np.abs(R[24:48, R2_WI]))
        print(f"  {label:>4} | {vals}   max={r2_max:.4f}")

    # ── Figure ──
    _make_figure(eng_data, cw_data, cw_depths_m, cw_widths_m,
                 x_r1_ft, z_ft, di_r2, wi9_x, di45_z)


def _make_figure(eng_data, cw_data, cw_depths_m, cw_widths_m,
                 x_r1_ft, z_ft, di_r2, wi9_x, di45_z):
    out_path = os.path.join(OUT_DIR, "stage5d_k096_all_profiles.png")
    n = len(RUNS)

    # shared y-limits across all runs for each column type
    all_r1_T, all_r2_T = [], []
    for label, *_ in RUNS:
        all_r1_T += eng_data[label][R1_DI, :].tolist() + cw_data[label][R1_DI, :].tolist()
        all_r2_T += eng_data[label][di_r2, R2_WI].tolist() + cw_data[label][di_r2, R2_WI].tolist()

    r1_T_lo = min(all_r1_T) - 1.5
    r1_T_hi = max(all_r1_T) + 1.5
    r2_T_lo = min(all_r2_T) - 1.5
    r2_T_hi = max(all_r2_T) + 1.5
    R_LIM   = 0.5   # fixed ±0.5°F for all residual axes

    fig, axes = plt.subplots(n, 4, figsize=(20, 3.2*n))
    fig.suptitle(
        f"k×{K_FACTOR} — All 9 Runs at t=168 hr\n"
        "Col 1: Side T (di=24)  Col 2: Side residual  "
        "Col 3: Bottom T (wi=0)  Col 4: Bottom residual",
        fontsize=11, fontweight="bold",
    )

    for row, (label, folder, pl, so) in enumerate(RUNS):
        dT  = abs(pl - so)
        eng = eng_data[label]
        cw  = cw_data[label]
        R   = eng - cw
        direction = "cool" if pl > so else ("—" if pl==so else "warm")

        r1_max_val = float(R[R1_DI, np.argmax(np.abs(R[R1_DI, :]))])
        r1_max_wi  = int(np.argmax(np.abs(R[R1_DI, :])))
        r2_col     = R[di_r2, R2_WI]
        r2_max_idx = int(np.argmax(np.abs(r2_col)))
        r2_max_val = float(r2_col[r2_max_idx])
        r2_max_di  = int(di_r2[r2_max_idx])

        gate_r1 = abs(r1_max_val) <= GATE_F
        gate_r2 = abs(r2_max_val) <= GATE_F
        run_color = "green" if (gate_r1 and gate_r2) else "red"

        # ── Col 0: Side T ──
        ax = axes[row, 0]
        ax.plot(x_r1_ft, cw[R1_DI, :],  "b-o", ms=4, lw=1.8, label="CW")
        ax.plot(x_r1_ft, eng[R1_DI, :], "r-o", ms=4, lw=1.8, label="Engine")
        ax.axvline(wi9_x, color="gray", ls="--", lw=0.8)
        ax.set_ylim(r1_T_lo, r1_T_hi)
        ax.set_ylabel("T (°F)", fontsize=7)
        ax.set_title(f"Run {label} | |ΔT|={dT}°F {direction} | Side T (di=24)", fontsize=8,
                     color=run_color, fontweight="bold")
        ax.tick_params(labelsize=6)
        if row == 0: ax.legend(fontsize=6, loc="lower right")

        # ── Col 1: Side residual ──
        ax = axes[row, 1]
        ax.plot(x_r1_ft, R[R1_DI, :], "r-o", ms=4, lw=1.8)
        ax.axhline(0, color="black", lw=0.7)
        ax.axhline( GATE_F, color="green", ls="--", lw=1.1)
        ax.axhline(-GATE_F, color="green", ls="--", lw=1.1)
        ax.axhspan(-GATE_F, GATE_F, color="green", alpha=0.07)
        ax.fill_between(x_r1_ft, R[R1_DI, :], 0,
                        where=(np.abs(R[R1_DI, :]) > GATE_F),
                        color="red", alpha=0.20)
        ax.axvline(wi9_x, color="gray", ls="--", lw=0.8)
        ax.set_ylim(-R_LIM, R_LIM)
        ax.set_ylabel("Eng−CW (°F)", fontsize=7)
        ax.set_title(f"Side residual  max|R|={abs(r1_max_val):.3f}°F @ wi={r1_max_wi}  "
                     f"{'✓' if gate_r1 else '✗'}", fontsize=8,
                     color="darkgreen" if gate_r1 else "darkred")
        ax.annotate(f"{r1_max_val:+.3f}",
                    xy=(x_r1_ft[r1_max_wi], r1_max_val),
                    xytext=(x_r1_ft[r1_max_wi]+1.5, r1_max_val+(0.12 if r1_max_val<0 else -0.12)),
                    fontsize=6.5, color="darkred",
                    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.7))
        ax.tick_params(labelsize=6)

        # ── Col 2: Bottom T ──
        ax = axes[row, 2]
        ax.plot(z_ft, cw[di_r2, R2_WI],  "b-o", ms=4, lw=1.8, label="CW")
        ax.plot(z_ft, eng[di_r2, R2_WI], "r-o", ms=4, lw=1.8, label="Engine")
        ax.axvline(di45_z, color="gray", ls="--", lw=0.8)
        ax.set_xlim(-82, -38); ax.set_ylim(r2_T_lo, r2_T_hi)
        ax.set_ylabel("T (°F)", fontsize=7)
        ax.set_title(f"Run {label} | Bottom T (wi=0)", fontsize=8)
        ax.tick_params(labelsize=6)
        if row == 0: ax.legend(fontsize=6, loc="upper left")

        # ── Col 3: Bottom residual ──
        ax = axes[row, 3]
        ax.plot(z_ft, r2_col, "r-o", ms=4, lw=1.8)
        ax.axhline(0, color="black", lw=0.7)
        ax.axhline( GATE_F, color="green", ls="--", lw=1.1)
        ax.axhline(-GATE_F, color="green", ls="--", lw=1.1)
        ax.axhspan(-GATE_F, GATE_F, color="green", alpha=0.07)
        ax.fill_between(z_ft, r2_col, 0,
                        where=(np.abs(r2_col) > GATE_F),
                        color="red", alpha=0.20)
        ax.axvline(di45_z, color="gray", ls="--", lw=0.8)
        ax.set_xlim(-82, -38); ax.set_ylim(-R_LIM, R_LIM)
        ax.set_ylabel("Eng−CW (°F)", fontsize=7)
        ax.set_title(f"Bottom residual  max|R|={abs(r2_max_val):.3f}°F @ di={r2_max_di}  "
                     f"{'✓' if gate_r2 else '✗'}", fontsize=8,
                     color="darkgreen" if gate_r2 else "darkred")
        ax.annotate(f"{r2_max_val:+.3f}",
                    xy=(z_ft[r2_max_idx], r2_max_val),
                    xytext=(z_ft[r2_max_idx]+3.0, r2_max_val+(0.12 if r2_max_val<0 else -0.12)),
                    fontsize=6.5, color="darkred",
                    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.7))
        ax.tick_params(labelsize=6)

        # shared x-labels on bottom row only
        if row == n - 1:
            axes[row, 0].set_xlabel("Dist from CL (ft)", fontsize=7)
            axes[row, 1].set_xlabel("Dist from CL (ft)", fontsize=7)
            axes[row, 2].set_xlabel("Depth z (ft)", fontsize=7)
            axes[row, 3].set_xlabel("Depth z (ft)", fontsize=7)

    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out_path}")


if __name__ == "__main__":
    main()
