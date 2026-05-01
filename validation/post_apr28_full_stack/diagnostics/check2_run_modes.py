#!/usr/bin/env python3
"""
CHECK 2 — Compare CW adiabatic-vs-full-stack heat shedding to the
engine's, per mix.

For each library mix (skip MIX-13):
  - Read CW adiabatic peak from validation/adiabatic_references/.
  - Read CW full-stack peak from scn.cw_validation (output.txt).
  - Run engine in adiabatic mode at the full-stack T0 (from input.dat).
  - Read engine full-stack peak from run_all_output.md (re-run if missing).

Compute, per mix:
  - engine_drop_F = engine_adiab_peak − engine_fs_peak  (both at fs T0)
  - cw_drop_F estimated using the T0-shifted CW adiabatic
    (apr28 validation: engine_adiab ≈ CW_adiab at the adiab CSV T0;
     extrapolate by T0 offset to fs T0)
  - drop_diff = engine_drop − cw_drop
    (positive → engine sheds MORE heat than CW)

Read-only: makes no source changes. Writes to diagnostics/.
"""

import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
os.chdir(REPO_ROOT)

from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import (
    build_grid_half_mat, build_grid_rectangular, solve_hydration_2d,
)


MIXES = [f"MIX-{n:02d}" for n in range(1, 16) if n != 13]

ADIAB_REF_DIR = "validation/adiabatic_references"
SCN_ROOT = "validation/cw_exports"


def f_to_c(t_f): return (t_f - 32.0) * 5.0 / 9.0
def c_to_f(t_c): return np.asarray(t_c) * 9.0 / 5.0 + 32.0


def cw_adiab_peak_and_T0(mix_id_lower):
    """Read CW adiabatic peak and starting T0 from the canonical CSV."""
    path = os.path.join(ADIAB_REF_DIR, f"cw_adiabatic_reference_{mix_id_lower}.csv")
    arr = np.genfromtxt(path, delimiter=",", skip_header=1)
    # columns: time_hrs, T_center_F_adiabatic, T_max_xs_F, T_ambient_F
    t0 = float(arr[0, 1])
    peak = float(arr[:, 1].max())
    return t0, peak


def run_engine_adiabatic_at_T0(mix_id, T0_F, duration_hrs=168.0):
    """Run the engine in adiabatic mode at a specified T0 (for the given mix
    composition). Returns peak centerline T in °F."""
    input_dat = os.path.join(SCN_ROOT, mix_id, "input.dat")
    weather   = os.path.join(SCN_ROOT, mix_id, "weather.dat")
    output    = os.path.join(SCN_ROOT, mix_id, "output.txt")
    scn = load_cw_scenario(input_dat, weather, output)  # default calibration on

    # Tiny grid; adiabatic field stays uniform
    grid = build_grid_rectangular(Lx_m=1.0, Ly_m=1.0, nx=3, ny=3)
    T_init = np.full((grid.ny, grid.nx), f_to_c(T0_F), dtype=np.float64)

    res = solve_hydration_2d(
        grid, scn.mix, T_init,
        duration_s=duration_hrs * 3600.0,
        output_interval_s=300.0,
        boundary_mode="adiabatic",
    )
    T_F = c_to_f(res.T_field_C[:, grid.ny // 2, grid.nx // 2])
    return float(T_F.max())


def main():
    rows = []
    for mix_id in MIXES:
        mix_lower = mix_id.replace("-", "").lower()  # MIX-01 → mix01

        # Load full-stack scenario (gives us fs T0, fs CW peak)
        input_dat = os.path.join(SCN_ROOT, mix_id, "input.dat")
        weather   = os.path.join(SCN_ROOT, mix_id, "weather.dat")
        output    = os.path.join(SCN_ROOT, mix_id, "output.txt")
        scn = load_cw_scenario(input_dat, weather, output)
        fs_T0 = scn.construction.placement_temp_F
        cw_fs_peak = float(scn.cw_validation.T_max_xs_F.max())

        # CW adiabatic peak from canonical CSV
        adiab_T0, cw_adiab_peak = cw_adiab_peak_and_T0(mix_lower)

        # Engine peaks
        engine_adiab_peak_at_fs_T0 = run_engine_adiabatic_at_T0(mix_id, fs_T0)
        # Also run at adiab T0 for the apples-to-apples adiab comparison
        engine_adiab_peak_at_adiab_T0 = (
            engine_adiab_peak_at_fs_T0 if abs(fs_T0 - adiab_T0) < 0.5
            else run_engine_adiabatic_at_T0(mix_id, adiab_T0)
        )

        # Engine full-stack peak from run_all_output.md (re-extract via solver
        # so this script is self-contained)
        grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
        T_init = np.full((grid.ny, grid.nx), f_to_c(fs_T0), dtype=np.float64)
        res_fs = solve_hydration_2d(
            grid, scn.mix, T_init,
            duration_s=168 * 3600,
            output_interval_s=1800.0,
            boundary_mode="full_2d",
            environment=scn.environment,
            construction=scn.construction,
            diagnostic_outputs=False,
        )
        jslice, islice = grid.concrete_slice()
        T_conc_F = res_fs.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
        engine_fs_peak = float(T_conc_F.max())

        # Drop computations
        # Engine drop: both runs at fs_T0 — clean.
        engine_drop = engine_adiab_peak_at_fs_T0 - engine_fs_peak

        # CW drop (estimated): CW_adiab is at adiab_T0 not fs_T0. We use the
        # apr28 fact that engine_adiab ≈ CW_adiab at the adiab CSV T0 to
        # extrapolate CW_adiab to fs_T0:
        #   estimated CW_adiab(fs_T0) = engine_adiab(fs_T0)
        # i.e., we trust the apr28 adiabatic validation to hold at any T0.
        cw_adiab_at_fs_T0_estimate = engine_adiab_peak_at_fs_T0
        cw_drop_estimate = cw_adiab_at_fs_T0_estimate - cw_fs_peak

        # Adiabatic-validation residual at adiab T0 (sanity check on apr28)
        adiab_residual = engine_adiab_peak_at_adiab_T0 - cw_adiab_peak

        # Adiabatic rise (for percent normalization, use fs_T0 reference)
        engine_adiab_rise = engine_adiab_peak_at_fs_T0 - fs_T0
        engine_drop_pct = (engine_drop / engine_adiab_rise * 100.0
                           if engine_adiab_rise > 0 else float("nan"))
        cw_drop_pct = (cw_drop_estimate / engine_adiab_rise * 100.0
                       if engine_adiab_rise > 0 else float("nan"))
        drop_diff_pct = engine_drop_pct - cw_drop_pct
        drop_diff_F = engine_drop - cw_drop_estimate

        rows.append(dict(
            mix=mix_id,
            fs_T0=fs_T0, adiab_T0=adiab_T0,
            cw_adiab_peak=cw_adiab_peak,
            cw_fs_peak=cw_fs_peak,
            engine_adiab_at_adiab_T0=engine_adiab_peak_at_adiab_T0,
            engine_adiab_at_fs_T0=engine_adiab_peak_at_fs_T0,
            engine_fs_peak=engine_fs_peak,
            adiab_residual=adiab_residual,
            engine_drop=engine_drop,
            cw_drop_estimate=cw_drop_estimate,
            engine_adiab_rise=engine_adiab_rise,
            engine_drop_pct=engine_drop_pct,
            cw_drop_pct=cw_drop_pct,
            drop_diff_F=drop_diff_F,
            drop_diff_pct=drop_diff_pct,
        ))
        print(f"  {mix_id}: fs_T0={fs_T0:.0f}°F  adiab_T0={adiab_T0:.0f}°F  "
              f"CW_adiab={cw_adiab_peak:.2f}°F  CW_fs={cw_fs_peak:.2f}°F  "
              f"eng_adiab_fsT0={engine_adiab_peak_at_fs_T0:.2f}°F  "
              f"eng_fs={engine_fs_peak:.2f}°F  "
              f"engine_drop={engine_drop:+.2f}°F  cw_drop_est={cw_drop_estimate:+.2f}°F  "
              f"diff={drop_diff_F:+.2f}°F", flush=True)

    # Markdown table
    md = ["| Mix | fs T0 | adiab T0 | CW adiab peak | CW fs peak | "
          "Engine adiab @ adiab T0 | Engine adiab @ fs T0 | Engine fs peak | "
          "Adiab residual (eng−CW @ adiab T0) | Engine drop (adiab−fs @ fs T0) | "
          "CW drop (estimated @ fs T0) | Engine drop − CW drop | "
          "Engine drop %  | CW drop %  | Δ% |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(
            f"| {r['mix']} | {r['fs_T0']:.0f}°F | {r['adiab_T0']:.0f}°F | "
            f"{r['cw_adiab_peak']:.2f}°F | {r['cw_fs_peak']:.2f}°F | "
            f"{r['engine_adiab_at_adiab_T0']:.2f}°F | "
            f"{r['engine_adiab_at_fs_T0']:.2f}°F | "
            f"{r['engine_fs_peak']:.2f}°F | "
            f"{r['adiab_residual']:+.2f}°F | "
            f"{r['engine_drop']:+.2f}°F | "
            f"{r['cw_drop_estimate']:+.2f}°F | "
            f"{r['drop_diff_F']:+.2f}°F | "
            f"{r['engine_drop_pct']:.1f}% | "
            f"{r['cw_drop_pct']:.1f}% | "
            f"{r['drop_diff_pct']:+.1f}% |"
        )
    md_path = os.path.join("validation", "post_apr28_full_stack",
                           "diagnostics", "check2_cw_modes.md")
    with open(md_path, "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"\nWrote {md_path}")

    # Aggregate stats over the 14 mixes
    diffs_pct = np.array([r["drop_diff_pct"] for r in rows], dtype=float)
    diffs_F   = np.array([r["drop_diff_F"]   for r in rows], dtype=float)
    print(f"\nAggregate drop_diff (engine_drop − cw_drop):")
    print(f"  mean = {diffs_F.mean():+.2f}°F  ({diffs_pct.mean():+.2f}% of adiab rise)")
    print(f"  SD   = {diffs_F.std(ddof=1):.2f}°F  ({diffs_pct.std(ddof=1):.2f}%)")
    print(f"  min  = {diffs_F.min():+.2f}°F   max = {diffs_F.max():+.2f}°F")


if __name__ == "__main__":
    main()
