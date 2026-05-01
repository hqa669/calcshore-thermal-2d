#!/usr/bin/env python3
"""
CHECK 3 — MIX-01 full-stack engine vs CW centerline trajectory.

Extracts the centerline temperature time series from the engine
(re-running at the same parameters as run_all.py) and from the CW
output.txt (via cw_scenario_loader). Writes a CSV with both
trajectories on the CW time grid, plus the residual = T_engine − T_cw.

Also generates a 2-panel PNG (trajectories on top, residual on bottom).

Read-only on source code; writes only to diagnostics/.
"""

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
os.chdir(REPO_ROOT)

from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d


SCENARIO_DIR = "validation/cw_exports/MIX-01"
OUT_CSV = "validation/post_apr28_full_stack/diagnostics/check3_mix01_trajectories.csv"
OUT_PNG = "validation/post_apr28_full_stack/diagnostics/check3_mix01_full_stack_trajectory.png"


def f_to_c(t_f): return (t_f - 32.0) * 5.0 / 9.0


def main():
    scn = load_cw_scenario(
        os.path.join(SCENARIO_DIR, "input.dat"),
        os.path.join(SCENARIO_DIR, "weather.dat"),
        os.path.join(SCENARIO_DIR, "output.txt"),
    )

    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T_initial = np.full((grid.ny, grid.nx),
                        f_to_c(scn.construction.placement_temp_F))

    res = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=False,
    )

    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_center_F = res.T_field_C[:, iy_mid, grid.nx - 1] * 9.0 / 5.0 + 32.0
    eng_t_s = res.t_s

    val = scn.cw_validation
    n_cw_t, n_cw_d, n_cw_w = val.T_field_F.shape
    cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]
    cw_t_s = val.time_hrs * 3600.0

    eng_on_cw_grid = np.interp(cw_t_s, eng_t_s, eng_center_F)
    residual = eng_on_cw_grid - cw_center_F

    arr = np.column_stack([val.time_hrs, eng_on_cw_grid, cw_center_F, residual])
    np.savetxt(OUT_CSV, arr, delimiter=",", fmt="%.4f",
               header="time_hrs,T_engine_F,T_cw_F,residual_F", comments="")
    print(f"Wrote {OUT_CSV}")

    # Characterize residual signature
    t = val.time_hrs
    abs_resid = np.abs(residual)

    first_exceed_idx = int(np.argmax(abs_resid > 1.0)) if np.any(abs_resid > 1.0) else -1
    first_exceed_hr = float(t[first_exceed_idx]) if first_exceed_idx >= 0 else None

    max_idx = int(np.argmax(abs_resid))
    max_hr = float(t[max_idx])
    max_val = float(residual[max_idx])

    final_resid = float(residual[-1])

    def rms(a): return float(np.sqrt(np.mean(a**2)))
    rms_full = rms(residual)
    mask_late = t >= 48.0
    mask_early = t <= 48.0
    rms_early = rms(residual[mask_early])
    rms_late = rms(residual[mask_late])

    # Diurnal AC component check: look at residual minus a slowly-varying
    # baseline. Compute std of (residual − running_mean_24hr).
    win = max(1, int(round(24.0 / (t[1] - t[0]))))  # samples per 24 hr
    if win < len(residual):
        baseline = np.convolve(residual, np.ones(win) / win, mode="same")
        ac = residual - baseline
        ac_std = float(np.std(ac))
        dc = float(np.mean(residual[len(residual) // 4:]))  # late-time mean
    else:
        ac_std = float("nan"); dc = float("nan")

    sig = (
        f"\nMIX-01 full-stack centerline residual signature\n"
        f"------------------------------------------------\n"
        f"Window: {len(residual)} samples on CW grid, t ∈ [{t[0]:.2f}, "
        f"{t[-1]:.2f}] hr\n"
        f"First |residual| > 1°F: hr {first_exceed_hr}\n"
        f"|residual| max: {max_val:+.2f}°F at hr {max_hr:.1f}\n"
        f"Residual at t=168 hr: {final_resid:+.2f}°F\n"
        f"RMS over full trajectory:       {rms_full:.2f}°F\n"
        f"RMS over [0, 48] hr (early):    {rms_early:.2f}°F\n"
        f"RMS over [48, 168] hr (late):   {rms_late:.2f}°F\n"
        f"DC offset (late-time mean):     {dc:+.2f}°F\n"
        f"AC component (std after 24-hr running-mean detrend): {ac_std:.2f}°F\n"
    )
    print(sig)

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

        ax = axes[0]
        ax.plot(t, eng_on_cw_grid, color="tab:blue", label=f"Engine center (post-apr28)")
        ax.plot(t, cw_center_F, "--", color="tab:orange", label="CW center")
        ax.set_ylabel("T (°F)")
        ax.set_title(f"MIX-01 full-stack centerline — engine vs CW (commit 978f5da)")
        ax.legend()
        ax.grid(alpha=0.3)
        ax.axvline(48.0, color="gray", linestyle=":", alpha=0.5)
        ax.text(48.5, ax.get_ylim()[0] + 5, "RMS window start", color="gray", fontsize=8)

        ax = axes[1]
        ax.plot(t, residual, color="tab:purple", label="T_engine − T_cw")
        ax.axhline(0.0, color="gray", alpha=0.5)
        ax.axhline(-1.0, color="gray", linestyle=":", alpha=0.3)
        ax.axhline(+1.0, color="gray", linestyle=":", alpha=0.3)
        ax.axvline(48.0, color="gray", linestyle=":", alpha=0.5)
        ax.set_xlabel("t (hr)")
        ax.set_ylabel("Residual (°F)")
        ax.set_title("Residual: T_engine − T_cw")
        ax.legend()
        ax.grid(alpha=0.3)

        plt.tight_layout()
        plt.savefig(OUT_PNG, dpi=110)
        plt.close(fig)
        print(f"Wrote {OUT_PNG}")
    except ImportError:
        print("matplotlib not installed; skipping PNG")


if __name__ == "__main__":
    main()
