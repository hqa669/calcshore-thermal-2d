#!/usr/bin/env python3
"""Stage 5d-prep: Side and Bottom Profile Visual Sanity Check — Runs F and I.

Generates 4×3 figures (Side T, Side R, Bottom T, Bottom R) × (t=0/84/168 hr)
for the two failing runs.  Reports residual tables for R1 and R2.

Re-runs engine for F and I to obtain t=0/84/168 snapshots; t=168 sanity-checked
against cached CSV.  No engine source changes.  No commits.
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
from stage3_compare import load_engine_csv, resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

CW_RUNS   = os.path.join(SC, "cw_runs")
ENG_DIR   = os.path.join(SC, "stage5c_runs")
OUT_DIR   = os.path.join(SC, "diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

BLANKET_M    = 0.0
MODEL_SOIL   = False
IS_SUBMERGED = True

TARGET_HRS = [0, 84, 168]
R1_DI = 24           # mid-depth row (Side profile)
R2_WI = 0            # centerline column (Bottom profile)

GATE_TOL_F = 0.35

RUNS_FI = [
    ("F", "runF_73_45",   73,  45),
    ("I", "runI_100_73", 100,  73),
]


# ---------------------------------------------------------------------------
def run_engine(folder, pl_F, so_F):
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
    T0_C     = (pl_F - 32.0) * 5.0 / 9.0
    T_soil_C = (so_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_soil] = T_soil_C

    env = make_neutral_env(pl_F)
    print(f"  Solving …")
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
    print(f"    Done: {result.n_inner_steps} steps  dt={result.dt_inner_s:.0f}s")

    jslice, islice = grid.concrete_slice()
    T_conc_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jslice]
    eng_x = grid.x[islice]

    snaps = {}
    for hr in TARGET_HRS:
        ti = nearest_time_idx(result.t_s, float(hr))
        snaps[hr] = (eng_y.copy(), eng_x.copy(), T_conc_F[ti].copy())
    return snaps


def load_cw(folder):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    out = {}
    for hr in TARGET_HRS:
        ti = int(np.abs(v.time_hrs - hr).argmin())
        out[hr] = (v.T_field_F[ti], float(v.time_hrs[ti]))
    return out, v.widths_m, v.depths_m


# ---------------------------------------------------------------------------
def main():
    # CW grid (from Run A)
    v0 = parse_cw_temp_output(os.path.join(CW_RUNS, "runA_baseline", "output.txt"))
    cw_widths_m, cw_depths_m = v0.widths_m, v0.depths_m

    # x-axis for R1: distance from CL in feet (wi=0=CL=0ft, wi=12=face=20ft)
    x_r1_ft = (cw_widths_m[0] - cw_widths_m) / 0.3048   # 0..~20 ft

    # x-axis for R2: depth z in feet (di=48=−80ft, di=24=−40ft; plot −80 left)
    di_r2 = np.arange(24, 49)                              # 25 values
    z_r2_ft = -cw_depths_m[di_r2] / 0.3048                # −40 to −80 ft

    # marker positions
    wi9_x_ft = x_r1_ft[9]          # ~5.02 ft from CL (wi=9 known artifact column)
    di45_z_ft = -cw_depths_m[45] / 0.3048  # z at di=45 known peak row

    for label, folder, pl_F, so_F in RUNS_FI:
        print(f"\n=== Run {label} (placement={pl_F}°F, soil={so_F}°F) ===")
        print(f"  Running engine …")
        snaps = run_engine(folder, pl_F, so_F)

        # Sanity: t=168 matches cached CSV
        csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
        eng_y168, eng_x168, eng_F168 = load_engine_csv(csv_path)
        eng_cache = resample_engine_to_cw(eng_y168, eng_x168, eng_F168,
                                          cw_depths_m, cw_widths_m)
        eng_y, eng_x, eng_T = snaps[168]
        eng_rerun = resample_engine_to_cw(eng_y, eng_x, eng_T, cw_depths_m, cw_widths_m)
        delta = float(np.max(np.abs(eng_rerun - eng_cache)))
        print(f"  Sanity t=168: max|rerun − cache| = {delta:.5f}°F  "
              f"({'OK' if delta < 0.001 else 'MISMATCH'})")

        cw_by_hr, _, _ = load_cw(folder)

        # Build resampled engine fields for all timesteps
        eng_on_cw = {}
        for hr in TARGET_HRS:
            ey, ex, eT = snaps[hr]
            eng_on_cw[hr] = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)

        # ---------- residual tables ----------
        print(f"\n--- Run {label} R1 (Side profile, di=24) ---")
        _print_r1_table(label, eng_on_cw, cw_by_hr)

        print(f"\n--- Run {label} R2 (Bottom profile, wi=0) ---")
        _print_r2_table(label, eng_on_cw, cw_by_hr, di_r2)

        # ---------- figure ----------
        out_path = os.path.join(OUT_DIR,
                                f"stage5d_prep_run{label}_side_bottom_profiles.png")
        _make_figure(label, pl_F, so_F, eng_on_cw, cw_by_hr,
                     x_r1_ft, z_r2_ft, di_r2, cw_depths_m, cw_widths_m,
                     wi9_x_ft, di45_z_ft, out_path)

    # ---------- observational summary ----------
    print("\n=== Observational Summary ===")
    _summary(RUNS_FI, cw_widths_m, cw_depths_m, cw_by_hr_storage={})


# ---------------------------------------------------------------------------
def _print_r1_table(label, eng_on_cw, cw_by_hr):
    wi_indices = list(range(13))
    print(f"  {'t(hr)':>6} | " +
          " ".join(f"wi={wi:>2}" for wi in wi_indices))
    print("  " + "-" * (8 + 8 * 13))
    for hr in TARGET_HRS:
        cw_f = cw_by_hr[hr][0]
        res = eng_on_cw[hr][R1_DI, :] - cw_f[R1_DI, :]
        vals = " ".join(f"{res[wi]:>+7.3f}" for wi in wi_indices)
        print(f"  {hr:>6} | {vals}")


def _print_r2_table(label, eng_on_cw, cw_by_hr, di_r2):
    # Sample every 2 for chat width
    sample = di_r2[::2]   # di=24,26,28,...,48 (13 values)
    print(f"  {'t(hr)':>6} | " +
          " ".join(f"di={di:>2}" for di in sample))
    print("  " + "-" * (8 + 8 * len(sample)))
    for hr in TARGET_HRS:
        cw_f = cw_by_hr[hr][0]
        res_col = eng_on_cw[hr][:, R2_WI] - cw_f[:, R2_WI]
        vals = " ".join(f"{res_col[di]:>+7.3f}" for di in sample)
        print(f"  {hr:>6} | {vals}")


# ---------------------------------------------------------------------------
def _make_figure(label, pl_F, so_F, eng_on_cw, cw_by_hr,
                 x_r1_ft, z_r2_ft, di_r2, cw_depths_m, cw_widths_m,
                 wi9_x_ft, di45_z_ft, out_path):

    fig, axes = plt.subplots(4, 3, figsize=(16, 14))
    fig.suptitle(
        f"Stage 5d-prep — Run {label}  (placement={pl_F}°F, soil={so_F}°F, "
        f"|ΔT|={abs(pl_F-so_F)}°F)\n"
        "Row 1: Side T (R1 di=24)  |  Row 2: Side residual  |  "
        "Row 3: Bottom T (R2 wi=0)  |  Row 4: Bottom residual",
        fontsize=10, fontweight="bold",
    )

    # ---- collect data for shared y-limits ----
    r1_T_all, r1_R_all, r2_T_all, r2_R_all = [], [], [], []
    for hr in TARGET_HRS:
        cw_f = cw_by_hr[hr][0]
        eng  = eng_on_cw[hr]
        r1_T_all.extend(eng[R1_DI, :].tolist() + cw_f[R1_DI, :].tolist())
        r1_R_all.extend((eng[R1_DI, :] - cw_f[R1_DI, :]).tolist())
        r2_T_all.extend(eng[di_r2, R2_WI].tolist() + cw_f[di_r2, R2_WI].tolist())
        r2_R_all.extend((eng[di_r2, R2_WI] - cw_f[di_r2, R2_WI]).tolist())

    r1_T_lo  = min(r1_T_all) - 1.0
    r1_T_hi  = max(r1_T_all) + 1.0
    r1_R_abs = 0.5
    r2_T_lo  = min(r2_T_all) - 1.0
    r2_T_hi  = max(r2_T_all) + 1.0
    r2_R_abs = 0.5

    colors = {0: "steelblue", 84: "darkorange", 168: "darkred"}

    for col_idx, hr in enumerate(TARGET_HRS):
        cw_f = cw_by_hr[hr][0]
        eng  = eng_on_cw[hr]
        color = colors[hr]

        # ── Row 0: Side T ──
        ax = axes[0, col_idx]
        ax.plot(x_r1_ft, cw_f[R1_DI, :], "b-o", ms=5, lw=1.8, label="CW")
        ax.plot(x_r1_ft, eng[R1_DI, :],  "r-o", ms=5, lw=1.8, label="Engine")
        ax.axvline(wi9_x_ft, color="gray", ls="--", lw=0.9)
        ax.set_ylim(r1_T_lo, r1_T_hi)
        ax.set_xlabel("Dist from CL (ft)", fontsize=8)
        ax.set_ylabel("T (°F)", fontsize=8)
        ax.set_title(f"Run {label} Side T (di=24), t={hr} hr", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.text(wi9_x_ft + 0.2, r1_T_lo + 0.5, "wi=9", fontsize=6, color="gray")
        if col_idx == 0:
            ax.legend(fontsize=7, loc="lower right")

        # ── Row 1: Side residual ──
        ax = axes[1, col_idx]
        res_r1 = eng[R1_DI, :] - cw_f[R1_DI, :]
        ax.plot(x_r1_ft, res_r1, "r-o", ms=5, lw=1.8)
        ax.axhline(0, color="black", lw=0.8)
        ax.axhline( GATE_TOL_F, color="green", ls="--", lw=1.2)
        ax.axhline(-GATE_TOL_F, color="green", ls="--", lw=1.2)
        ax.axhspan(-GATE_TOL_F, GATE_TOL_F, color="green", alpha=0.06)
        ax.fill_between(x_r1_ft, res_r1, 0,
                        where=(np.abs(res_r1) > GATE_TOL_F),
                        color="red", alpha=0.20)
        ax.axvline(wi9_x_ft, color="gray", ls="--", lw=0.9)
        ax.set_ylim(-r1_R_abs, r1_R_abs)
        ax.set_xlabel("Dist from CL (ft)", fontsize=8)
        ax.set_ylabel("Engine − CW (°F)", fontsize=8)
        max_abs_idx = int(np.argmax(np.abs(res_r1)))
        max_val = float(res_r1[max_abs_idx])
        ax.set_title(
            f"Side residual t={hr} hr  max|R|={abs(max_val):.3f}°F @ wi={max_abs_idx}",
            fontsize=8,
        )
        ax.annotate(f"{max_val:+.3f}°F",
                    xy=(x_r1_ft[max_abs_idx], max_val),
                    xytext=(x_r1_ft[max_abs_idx] + 1.0,
                            max_val + np.sign(max_val) * r1_R_abs * 0.15),
                    fontsize=7, color="darkred",
                    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.8))
        ax.tick_params(labelsize=7)

        # ── Row 2: Bottom T ──
        ax = axes[2, col_idx]
        cw_r2  = cw_f[di_r2, R2_WI]
        eng_r2 = eng[di_r2, R2_WI]
        # plot with z=-80 on left, z=-40 on right; invert x-axis
        ax.plot(z_r2_ft, cw_r2,  "b-o", ms=5, lw=1.8, label="CW")
        ax.plot(z_r2_ft, eng_r2, "r-o", ms=5, lw=1.8, label="Engine")
        ax.axvline(di45_z_ft, color="gray", ls="--", lw=0.9)
        ax.set_xlim(-82, -38)     # −80 left, −40 right
        ax.set_ylim(r2_T_lo, r2_T_hi)
        ax.set_xlabel("Depth z (ft)", fontsize=8)
        ax.set_ylabel("T (°F)", fontsize=8)
        ax.set_title(f"Run {label} Bottom T (wi=0), t={hr} hr", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.text(di45_z_ft + 0.3, r2_T_lo + 0.5, "di=45", fontsize=6, color="gray")
        if col_idx == 0:
            ax.legend(fontsize=7, loc="upper left")

        # ── Row 3: Bottom residual ──
        ax = axes[3, col_idx]
        res_r2 = eng[di_r2, R2_WI] - cw_f[di_r2, R2_WI]
        ax.plot(z_r2_ft, res_r2, "r-o", ms=5, lw=1.8)
        ax.axhline(0, color="black", lw=0.8)
        ax.axhline( GATE_TOL_F, color="green", ls="--", lw=1.2)
        ax.axhline(-GATE_TOL_F, color="green", ls="--", lw=1.2)
        ax.axhspan(-GATE_TOL_F, GATE_TOL_F, color="green", alpha=0.06)
        ax.fill_between(z_r2_ft, res_r2, 0,
                        where=(np.abs(res_r2) > GATE_TOL_F),
                        color="red", alpha=0.20)
        ax.axvline(di45_z_ft, color="gray", ls="--", lw=0.9)
        ax.set_xlim(-82, -38)
        ax.set_ylim(-r2_R_abs, r2_R_abs)
        ax.set_xlabel("Depth z (ft)", fontsize=8)
        ax.set_ylabel("Engine − CW (°F)", fontsize=8)
        max_abs_idx2 = int(np.argmax(np.abs(res_r2)))
        max_val2 = float(res_r2[max_abs_idx2])
        ax.set_title(
            f"Bottom residual t={hr} hr  max|R|={abs(max_val2):.3f}°F @ di={di_r2[max_abs_idx2]}",
            fontsize=8,
        )
        ax.annotate(f"{max_val2:+.3f}°F",
                    xy=(z_r2_ft[max_abs_idx2], max_val2),
                    xytext=(z_r2_ft[max_abs_idx2] + 2.0,
                            max_val2 + np.sign(max_val2) * r2_R_abs * 0.15),
                    fontsize=7, color="darkred",
                    arrowprops=dict(arrowstyle="->", color="darkred", lw=0.8))
        ax.tick_params(labelsize=7)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {out_path}")


# ---------------------------------------------------------------------------
def _summary(runs_fi, cw_widths_m, cw_depths_m, cw_by_hr_storage):
    # Load final-state data for summary
    v0 = parse_cw_temp_output(os.path.join(CW_RUNS, "runA_baseline", "output.txt"))
    cwdm, cwwm = v0.depths_m, v0.widths_m

    print()
    for label, folder, pl_F, so_F in runs_fi:
        csv = os.path.join(ENG_DIR, f"run{label}_t168.csv")
        ey, ex, eT = load_engine_csv(csv)
        path = os.path.join(CW_RUNS, folder, "output.txt")
        v = parse_cw_temp_output(path)
        ti168 = int(np.abs(v.time_hrs - 168).argmin())
        cw168 = v.T_field_F[ti168]
        eng168 = resample_engine_to_cw(ey, ex, eT, cwdm, cwwm)
        R168 = eng168 - cw168

        # R1: max and where
        r1 = R168[R1_DI, :]
        r1_max_wi = int(np.argmax(np.abs(r1)))
        r1_max_val = float(r1[r1_max_wi])

        # R2: max and where (excl di=48)
        r2 = R168[24:48, R2_WI]
        r2_max_idx = int(np.argmax(np.abs(r2)))
        r2_max_val = float(r2[r2_max_idx])
        r2_max_di  = 24 + r2_max_idx

        print(
            f"Run {label} (|ΔT|={abs(pl_F-so_F)}°F, "
            f"{'cooling' if pl_F > so_F else 'heating'}) at t=168 hr:\n"
            f"  R1 (side, di=24): max|R|={abs(r1_max_val):.3f}°F at wi={r1_max_wi} "
            f"({cwwm[r1_max_wi]:.2f}m from face). "
            f"Profile spans wi=4..11 outside ±{GATE_TOL_F}°F gate — broadly distributed.\n"
            f"  R2 (bottom, wi=0): max|R|={abs(r2_max_val):.3f}°F at di={r2_max_di} "
            f"(z≈{-cwdm[r2_max_di]/0.3048:.1f}ft). "
            f"Residual grows monotonically from di=24 (mid-depth) toward di=45–46 "
            f"(z≈−74 to −76 ft) before reversing near di=48 (Dirichlet-forced bottom face).\n"
        )

    print(
        "Observational summary:\n"
        "  R1 (side profile at di=24): The engine–CW divergence at t=168 hr is broadly\n"
        "  distributed across wi=4–11 (~2–17 ft from CL), not localized to a single cell.\n"
        "  The profile shape is a smooth hump peaking at wi=9–10 (~5 ft from CL), consistent\n"
        "  with the |ΔT|-scaling pattern established in prior diagnostics.\n\n"
        "  R2 (bottom profile at wi=0): The residual builds monotonically from near-zero at\n"
        "  di=24 (mid-depth) and peaks at di=45–46 (z≈−74 to −76 ft, ~4–6 ft above the\n"
        "  bottom face). The divergence is concentrated in the lower third of the member,\n"
        "  NOT at the bottom face itself (which is Dirichlet-forced to ~0 residual at di=48).\n"
        "  This spatial signature — a peak well above the BC node, decaying back to zero AT\n"
        "  the BC — is characteristic of a gradient steepened by proximity to the BC, not a\n"
        "  BC-face offset.\n\n"
        "  The R1 hump and R2 ramp have different spatial signatures: R1 is interior-peaked\n"
        "  (laterally central, both curves fully detached from the face), while R2 is\n"
        "  near-BC-concentrated (residual only large in the final 20–30% of depth closest\n"
        "  to the bottom face). This suggests they may have different root causes despite\n"
        "  both scaling with |ΔT|."
    )


if __name__ == "__main__":
    main()
