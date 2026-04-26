"""
Shared engine driver for PR 18 reconnaissance diagnostics.

Mirrors compare_to_cw.run_one()'s load+solve+extract sequence, but
returns raw arrays (not just metrics) needed by diagnostics d1-d5.
compare_to_cw.py is not modified.
"""
import os
import sys
import numpy as np
from dataclasses import dataclass

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _REPO)

from cw_scenario_loader import load_cw_scenario, CWScenario
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d, Grid2D, HydrationResult

REFERENCE_MIXES = ("MIX-01", "MIX-02", "MIX-03", "MIX-11", "MIX-12")
SCENARIO_ROOT = os.path.join(_REPO, "validation", "cw_exports")

T_START_RMS_HR = 48.0  # gate window start — matches compare_to_cw.py


@dataclass
class MixCache:
    mix_name: str
    scn: CWScenario
    grid: Grid2D
    result: HydrationResult
    cw_time_hrs: np.ndarray            # (n_cw_t,)
    cw_center_F: np.ndarray            # (n_cw_t,) — mid-depth, width=0
    eng_center_F_at_cw_t: np.ndarray   # (n_cw_t,) — engine interp to CW times
    eng_centerline_col_F_at_cw_t: np.ndarray  # (n_cw_t, n_conc_y), depth-0=top
    cw_centerline_col_F: np.ndarray    # (n_cw_t, n_cw_d), depth-0=top
    gate_mask: np.ndarray              # [48, 168] hr
    mode_a_mask: np.ndarray            # [8, 48] hr
    diurnal_mask: np.ndarray           # [144, 168] hr


def load_mix_cache(mix_name: str, scenario_root: str = SCENARIO_ROOT) -> MixCache:
    scenario_dir = os.path.join(scenario_root, mix_name)

    scn = load_cw_scenario(
        os.path.join(scenario_dir, "input.dat"),
        os.path.join(scenario_dir, "weather.dat"),
        os.path.join(scenario_dir, "output.txt"),
    )

    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)

    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )

    val = scn.cw_validation
    cw_t_s = val.time_hrs * 3600.0
    n_cw_t, n_cw_d, _n_cw_w = val.T_field_F.shape

    # CW centerline mid-depth (compare_to_cw.py:209)
    cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]

    # Engine centerline mid-depth
    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_ctr_C = result.T_field_C[:, iy_mid, grid.nx - 1]
    eng_ctr_F = eng_ctr_C * 9.0 / 5.0 + 32.0
    eng_center_F_at_cw_t = np.interp(cw_t_s, result.t_s, eng_ctr_F)

    # Engine centerline column (all concrete depths) sampled at CW times
    jslice, _islice = grid.concrete_slice()
    n_conc_y = grid.iy_concrete_end - grid.iy_concrete_start + 1
    eng_col = np.empty((n_cw_t, n_conc_y), dtype=float)
    for ti in range(n_cw_t):
        idx = int(np.searchsorted(result.t_s, cw_t_s[ti]))
        idx = min(idx, len(result.t_s) - 1)
        col_C = result.T_field_C[idx, jslice, grid.nx - 1]
        eng_col[ti] = col_C * 9.0 / 5.0 + 32.0

    cw_col_F = val.T_field_F[:, :, 0]   # (n_cw_t, n_cw_d)

    gate_mask    = val.time_hrs >= T_START_RMS_HR
    mode_a_mask  = (val.time_hrs >= 8.0) & (val.time_hrs <= T_START_RMS_HR)
    diurnal_mask = (val.time_hrs >= 144.0) & (val.time_hrs <= 168.0)

    return MixCache(
        mix_name=mix_name,
        scn=scn,
        grid=grid,
        result=result,
        cw_time_hrs=val.time_hrs,
        cw_center_F=cw_center_F,
        eng_center_F_at_cw_t=eng_center_F_at_cw_t,
        eng_centerline_col_F_at_cw_t=eng_col,
        cw_centerline_col_F=cw_col_F,
        gate_mask=gate_mask,
        mode_a_mask=mode_a_mask,
        diurnal_mask=diurnal_mask,
    )
