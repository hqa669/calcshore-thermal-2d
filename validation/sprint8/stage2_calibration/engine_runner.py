"""Sprint 8 Stage 2 — shared engine driver.

Entry point: run_and_compare(folder, T_pl_F, T_so_F, factor).

The k_uc factor is applied by monkey-patching
thermal_engine_2d.K_UC_CALIBRATION_FACTOR_SPRINT7 before each solve and
restoring it afterward.  factor=0.96 exactly replicates Sprint 7's
baked-in calibration; other values are absolute final multipliers.

Region definitions (Structure C, Sprint 7):
  R1 — side profile at di=24 across all widths         gate: max|R| ≤ 0.35°F
  R2 — bottom centerline wi=0, di=24..48               gate: max|R| ≤ 0.35°F
  R3 — bottom-side corner di=43..48, wi=8..12          report RMSE only
"""

import os
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = (HERE / "../../..").resolve()
SC   = (HERE / "../../soil_calibration").resolve()

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from kinetics_correction import compute_hu_factor
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

COMPARE_HR = 168
GATE_F     = 0.35

R1_DI = 24
R2_WI = 0
R2_DI = slice(24, 49)
R3_DI = slice(43, 49)
R3_WI = slice(8, 13)


@dataclass
class ResidualReport:
    r1_max_F:  float
    r2_max_F:  float
    r3_rmse_F: float
    gate_r1:   bool
    gate_r2:   bool
    passes:    bool
    minimax:   float
    R_field:   np.ndarray = field(repr=False)
    eng_F:     np.ndarray = field(repr=False, default=None)
    cw_F:      np.ndarray = field(repr=False, default=None)


def run_and_compare(
    folder: Path,
    T_pl_F: float,
    T_so_F: float,
    factor: float = 0.96,
    compare_hr: float = COMPARE_HR,
) -> ResidualReport:
    """Solve engine at given k_uc factor; return residuals vs CW output.txt."""
    folder = Path(folder)
    original = te2d.K_UC_CALIBRATION_FACTOR_SPRINT7
    try:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = factor
        eng_y, eng_x, eng_T_F = _run_engine(folder, T_pl_F, T_so_F, compare_hr)
    finally:
        te2d.K_UC_CALIBRATION_FACTOR_SPRINT7 = original

    cw_T_F, cw_widths_m, cw_depths_m = _load_cw(folder, compare_hr)
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_T_F, cw_depths_m, cw_widths_m)
    R = eng_on_cw - cw_T_F

    r1   = float(np.max(np.abs(R[R1_DI, :])))
    r2   = float(np.max(np.abs(R[R2_DI, R2_WI])))
    r3   = float(np.sqrt(np.mean(R[R3_DI, R3_WI] ** 2)))
    mini = max(r1, r2)
    return ResidualReport(
        r1_max_F=r1, r2_max_F=r2, r3_rmse_F=r3,
        gate_r1=(r1 <= GATE_F), gate_r2=(r2 <= GATE_F),
        passes=(r1 <= GATE_F and r2 <= GATE_F),
        minimax=mini, R_field=R, eng_F=eng_on_cw, cw_F=cw_T_F,
    )


def _run_engine(folder: Path, T_pl_F: float, T_so_F: float, compare_hr: float):
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))
    fac, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = fac
    mix.Hu_J_kg_effective    = mix.Hu_J_kg * fac
    constr.model_soil    = False
    constr.is_submerged  = True
    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0 = (T_pl_F - 32.0) * 5.0 / 9.0
    Ts = (T_so_F - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0)
    Ti[grid.is_soil] = Ts
    res = solve_hydration_2d(
        grid, mix, Ti,
        duration_s=compare_hr * 3600.0,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=make_neutral_env(T_pl_F),
        construction=constr,
        T_ground_deep_C=Ts,
        diagnostic_outputs=False,
    )
    jsl, isl = grid.concrete_slice()
    ti = nearest_time_idx(res.t_s, float(compare_hr))
    T_F = res.T_field_C[ti, jsl, isl] * 9.0 / 5.0 + 32.0
    return grid.y[jsl], grid.x[isl], T_F


def _load_cw(folder: Path, compare_hr: float):
    v  = parse_cw_temp_output(str(folder / "output.txt"))
    ti = int(np.abs(v.time_hrs - compare_hr).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m
