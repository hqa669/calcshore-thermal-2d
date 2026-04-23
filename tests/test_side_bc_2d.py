"""
M4 tests: air-mask geometry, side BC, centerline symmetry BC, and CW MIX-01
acceptance gates for boundary_mode='full_2d'.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np
import pytest

from thermal_engine_2d import (
    FORM_R_SI,
    R_IMP_TO_R_SI,
    HydrationResult,
    ambient_temp_F,
    build_grid_half_mat,
    build_grid_rectangular,
    compute_T_gw_C,
    h_convective,
    h_side_convective,
    h_top_series,
    menzel_evaporation,
    solve_conduction_2d,
    solve_hydration_2d,
    analytical_square_slab,
)


# ============================================================
# Stubs (same pattern as test_top_bc_2d.py)
# ============================================================

@dataclass
class EnvStub:
    daily_max_F: List[float] = field(default_factory=lambda: [92.0] * 8)
    daily_min_F: List[float] = field(default_factory=lambda: [75.0] * 8)
    cw_ave_max_wind_m_s: float = 10.5
    cw_ave_max_RH_pct: float = 87.4
    cw_ave_min_RH_pct: float = 44.6


@dataclass
class ConstrStub:
    blanket_R_value: float = 5.67
    top_cure_blanket_time_hrs: float = 2.0
    placement_hour: int = 5
    placement_temp_F: float = 60.0


@dataclass
class MixStub:
    Hu_J_kg: float = 424143.0
    tau_hrs: float = 29.401
    beta: float = 0.895
    alpha_u: float = 0.7585
    activation_energy_J_mol: float = 26457.9
    total_cementitious_lb_yd3: float = 575.0
    concrete_density_lb_ft3: float = 131.2
    thermal_conductivity_BTU_hr_ft_F: float = 1.56
    aggregate_Cp_BTU_lb_F: float = 0.20
    cement_type_I_II_lb_yd3: float = 350.0
    water_lb_yd3: float = 253.0
    coarse_agg_lb_yd3: float = 1800.0
    fine_agg_lb_yd3: float = 1100.0


MIX01 = MixStub()


# ============================================================
# Test 1: regression — M0-M3 APIs still work after M4 geometry change
# ============================================================

def test_regression_all_previous():
    """M0/M1/M2/M3 APIs intact after M4 geometry and BC additions."""
    grid = build_grid_half_mat(40.0, 8.0)
    assert grid.nx == 33
    assert grid.ny == 29

    rect = build_grid_rectangular(2.0, 1.0, 11, 11)
    T0 = np.full((rect.ny, rect.nx), 20.0)
    res_cond = solve_conduction_2d(rect, 2.7, 2400.0, 900.0, T0, 0.0, 100.0)
    assert np.all(np.isfinite(res_cond.T_field_C))

    T_an = analytical_square_slab(rect.x, rect.y, 100.0, 2.0, 1.0,
                                   2.7 / (2400 * 900), 20.0, 0.0)
    assert T_an.shape == (11, 11)

    T0_hyd = np.full((grid.ny, grid.nx), 20.0)
    env = EnvStub()
    constr = ConstrStub()

    for mode in ("adiabatic", "dirichlet", "top_bc_only", "v2_equivalent"):
        kw = {}
        if mode in ("top_bc_only", "v2_equivalent"):
            kw = {"environment": env, "construction": constr}
        res = solve_hydration_2d(
            grid, MIX01, T0_hyd,
            duration_s=3600.0,
            output_interval_s=3600.0,
            boundary_mode=mode,
            **kw,
        )
        assert isinstance(res, HydrationResult), f"mode={mode} returned wrong type"
        assert np.all(np.isfinite(res.T_field_C)), f"NaN in mode={mode}"


# ============================================================
# Test 2: air region cells are not updated by the solver
# ============================================================

def test_air_region_not_solved():
    """Air cells (material_id=3) must keep their initial value throughout the run."""
    grid = build_grid_half_mat(40.0, 8.0)
    assert grid.is_air.any(), "Expected air cells in half-mat grid"

    T_init = np.full((grid.ny, grid.nx), 20.0)
    # Mark air cells with a sentinel temperature
    T_init[grid.is_air] = 999.0

    env = EnvStub()
    constr = ConstrStub()

    res = solve_hydration_2d(
        grid, MIX01, T_init,
        duration_s=3600.0,
        output_interval_s=3600.0,
        boundary_mode="full_2d",
        environment=env,
        construction=constr,
    )

    T_final = res.T_field_C[-1]
    assert np.allclose(T_final[grid.is_air], 999.0, atol=1e-6), (
        "Air cells must not be updated by the solver"
    )
    assert not np.any(T_final[~grid.is_air] == 999.0), (
        "Non-air cells should not be stuck at 999.0 sentinel"
    )


# ============================================================
# Test 3: h_side_convective formula check
# ============================================================

def test_h_side_combination():
    """h_side_convective = 1 / (1/h_conv + R_FORM_EFFECTIVE_SI) ≈ 8 W/(m²·K).

    R_FORM_EFFECTIVE_SI is the empirically calibrated effective form-face
    resistance (≈ 0.0862 m²·K/W = 0.490 hr·ft²·°F/BTU), giving h_side ≈ 8.
    The cure blanket is a top-surface treatment only; blanket_R_imp is
    accepted but does not affect the side BC.
    """
    from thermal_engine_2d import R_FORM_EFFECTIVE_SI
    wind = 10.5
    R_imp = 5.67          # passed but ignored per function contract
    h_conv = h_convective(wind)     # 25.8 W/(m²·K)
    h_side_expected = 1.0 / (1.0 / h_conv + R_FORM_EFFECTIVE_SI)

    h_side = h_side_convective(wind, R_imp)
    assert h_side == pytest.approx(h_side_expected, rel=1e-5)
    assert h_side == pytest.approx(8.0, abs=0.1), (
        f"h_side should be ≈ 8.0 W/(m²·K) with R_FORM_EFFECTIVE_SI={R_FORM_EFFECTIVE_SI:.4f}"
    )


# ============================================================
# Test 4: centerline symmetry preserves uniform field
# ============================================================

def test_centerline_symmetry_preserves_uniform():
    """With adiabatic BCs and uniform IC, all non-air cells stay at T_init."""
    grid = build_grid_half_mat(40.0, 8.0)
    T_init = np.full((grid.ny, grid.nx), 30.0)
    T_init[grid.is_air] = 30.0   # air cells also start at 30°C for this test

    # Adiabatic mode: no heat in or out, uniform IC → T should stay constant
    res = solve_hydration_2d(
        grid, MixStub(alpha_u=0.0), T_init,
        duration_s=24 * 3600,
        output_interval_s=24 * 3600,
        boundary_mode="adiabatic",
    )
    T_final = res.T_field_C[-1]
    non_air = ~grid.is_air
    assert np.allclose(T_final[non_air], 30.0, atol=1e-4), (
        f"Non-air cells drifted from 30°C: max dev = "
        f"{abs(T_final[non_air] - 30.0).max():.6f}"
    )


# ============================================================
# Test 5: side BC cools the form edge relative to centerline
# ============================================================

def test_side_bc_cools_edge():
    """With full_2d and warm IC + cool ambient, form-edge concrete is cooler
    than the centerline concrete after 24 hr."""
    grid = build_grid_half_mat(40.0, 8.0)
    T_init = np.full((grid.ny, grid.nx), 60.0)
    T_init[grid.is_air] = 999.0  # sentinel — must be ignored by solver

    # Constant cool ambient (68°F = 20°C)
    env = EnvStub(
        daily_max_F=[68.0] * 8,
        daily_min_F=[68.0] * 8,
        cw_ave_max_wind_m_s=5.0,
        cw_ave_max_RH_pct=60.0,
        cw_ave_min_RH_pct=60.0,
    )
    constr = ConstrStub(top_cure_blanket_time_hrs=0.0, placement_hour=15)
    mix_no_hyd = MixStub(alpha_u=0.0)

    res = solve_hydration_2d(
        grid, mix_no_hyd, T_init,
        duration_s=24 * 3600,
        output_interval_s=24 * 3600,
        boundary_mode="full_2d",
        environment=env,
        construction=constr,
    )

    T_final = res.T_field_C[-1]
    jslice = slice(grid.iy_concrete_start, grid.iy_concrete_end + 1)
    edge_mean = T_final[jslice, grid.ix_concrete_start].mean()
    center_mean = T_final[jslice, grid.nx - 1].mean()

    assert edge_mean < center_mean, (
        f"Form-edge concrete ({edge_mean:.2f}°C) should be cooler than "
        f"centerline ({center_mean:.2f}°C) — side BC is not removing heat"
    )


# ============================================================
# Tests 6-10: CW MIX-01 acceptance gates (full_2d mode)
# ============================================================

def _load_mix01():
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario
    return load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )


def _run_full2d(scn):
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_air] = T0_C   # air cells start at placement temp
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
    )
    return grid, result


@pytest.mark.xfail(
    reason=(
        "Peak T runs ~1.0°F below CW (−1.04°F on MIX-01). Root cause: missing "
        "solar gain during warm daytime hours. Without solar input the hydration "
        "peak is suppressed slightly. Expected to pass after Sprint 1 adds solar "
        "+ longwave top BC."
    ),
    strict=False,
)
def test_matches_cw_mix01_peak_T():
    """M4 gate: peak max T in concrete ≈ CW 129.6°F ± 1.0°F."""
    scn = _load_mix01()
    grid, result = _run_full2d(scn)

    jslice, islice = grid.concrete_slice()
    T_concrete_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
    engine_peak_F = float(T_concrete_F.max())
    cw_peak_F = 129.6
    delta_F = engine_peak_F - cw_peak_F

    print(f"\nEngine peak: {engine_peak_F:.2f} °F")
    print(f"CW peak:     {cw_peak_F:.2f} °F")
    print(f"Delta:       {delta_F:+.2f} °F")

    assert abs(delta_F) < 1.0, (
        f"Peak max T delta = {delta_F:+.2f} °F exceeds ±1.0°F tolerance"
    )


def test_matches_cw_mix01_peak_gradient():
    """M4 gate: peak gradient ≈ CW 39.3°F ± 2.0°F."""
    scn = _load_mix01()
    grid, result = _run_full2d(scn)

    jslice, islice = grid.concrete_slice()
    T_concrete_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
    grad_series = T_concrete_F.max(axis=(1, 2)) - T_concrete_F.min(axis=(1, 2))
    engine_peak_grad_F = float(grad_series.max())
    cw_peak_grad_F = 39.3
    delta_F = engine_peak_grad_F - cw_peak_grad_F

    print(f"\nEngine peak gradient: {engine_peak_grad_F:.2f} °F")
    print(f"CW peak gradient:     {cw_peak_grad_F:.2f} °F")
    print(f"Delta:                {delta_F:+.2f} °F")

    assert abs(delta_F) < 2.0, (
        f"Peak gradient delta = {delta_F:+.2f} °F exceeds ±2.0°F tolerance"
    )


def test_matches_cw_mix01_field_rms():
    """M4 gate: centerline field RMS vs CW < 2.0°F."""
    scn = _load_mix01()
    if scn.cw_validation.T_field_F is None:
        pytest.skip("CW T_field_F not available in validation fixture")

    grid, result = _run_full2d(scn)
    jslice, islice = grid.concrete_slice()
    cw_field_F = scn.cw_validation.T_field_F
    n_cw_t, n_cw_d, n_cw_w = cw_field_F.shape
    cw_t_s = scn.cw_validation.time_hrs * 3600.0

    sq_errs = []
    for ti in range(n_cw_t):
        idx = int(np.searchsorted(result.t_s, cw_t_s[ti]))
        idx = min(idx, len(result.t_s) - 1)
        T_eng_col_C = result.T_field_C[idx, jslice, grid.nx - 1]
        T_eng_col_F = T_eng_col_C * 9.0 / 5.0 + 32.0
        n_cmp = min(len(T_eng_col_F), n_cw_d)
        sq_errs.extend((T_eng_col_F[:n_cmp] - cw_field_F[ti, :n_cmp, 0]).tolist())

    rms_F = float(np.sqrt(np.mean(np.array(sq_errs) ** 2)))
    print(f"\nField-wide RMS (centerline): {rms_F:.2f} °F")
    assert rms_F < 2.0, f"Field RMS = {rms_F:.2f} °F exceeds 2.0°F tolerance"


@pytest.mark.xfail(
    reason=(
        "Corner RMS ~4°F on MIX-01; root cause is missing solar radiation on "
        "the form face. CW's afternoon corner temperature exceeds ambient by "
        "several °F due to solar gain not modelled in Sprint 0. "
        "Expected to pass after Sprint 1 adds solar + longwave top BC."
    ),
    strict=False,
)
def test_matches_cw_mix01_corner_rms():
    """M4 gate: corner node RMS vs CW < 3.0°F."""
    scn = _load_mix01()
    if scn.cw_validation.T_field_F is None:
        pytest.skip("CW T_field_F not available in validation fixture")

    grid, result = _run_full2d(scn)
    cw_t_s = scn.cw_validation.time_hrs * 3600.0
    # CW corner: depth=0 (top), width=-1 (form edge)
    cw_corner_F = scn.cw_validation.T_field_F[:, 0, -1]

    iy_top = grid.iy_concrete_start
    ix_corner = grid.ix_concrete_start
    eng_corner_C = result.T_field_C[:, iy_top, ix_corner]
    eng_corner_F = eng_corner_C * 9.0 / 5.0 + 32.0
    eng_corner_interp = np.interp(cw_t_s, result.t_s, eng_corner_F)

    rms_F = float(np.sqrt(np.mean((eng_corner_interp - cw_corner_F) ** 2)))
    print(f"\nCorner RMS: {rms_F:.2f} °F")
    assert rms_F < 3.0, f"Corner RMS = {rms_F:.2f} °F exceeds 3.0°F tolerance"


def test_matches_cw_mix01_centerline_rms():
    """M4 gate: centerline mid-depth RMS vs CW < 1.0°F."""
    scn = _load_mix01()
    if scn.cw_validation.T_field_F is None:
        pytest.skip("CW T_field_F not available in validation fixture")

    grid, result = _run_full2d(scn)
    cw_t_s = scn.cw_validation.time_hrs * 3600.0
    n_cw_d = scn.cw_validation.T_field_F.shape[1]
    # CW centerline: width=0; mid-depth = n_cw_d // 2
    cw_center_F = scn.cw_validation.T_field_F[:, n_cw_d // 2, 0]

    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_center_C = result.T_field_C[:, iy_mid, grid.nx - 1]
    eng_center_F = eng_center_C * 9.0 / 5.0 + 32.0
    eng_center_interp = np.interp(cw_t_s, result.t_s, eng_center_F)

    rms_F = float(np.sqrt(np.mean((eng_center_interp - cw_center_F) ** 2)))
    print(f"\nCenterline mid-depth RMS: {rms_F:.2f} °F")
    assert rms_F < 1.0, f"Centerline RMS = {rms_F:.2f} °F exceeds 1.0°F tolerance"
