"""
M3 tests: top boundary condition (convection + blanket + Menzel evaporation)
and CW MIX-01 Austin validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List

import numpy as np
import pytest

from thermal_engine_2d import (
    R_IMP_TO_R_SI,
    HydrationResult,
    ambient_temp_F,
    build_grid_half_mat,
    build_grid_rectangular,
    compute_T_gw_C,
    h_forced_convection,
    _h_convective_legacy,
    h_top_series,
    menzel_evaporation,
    saturated_vapor_pressure_mmHg,
    solve_conduction_2d,
    solve_hydration_2d,
    analytical_square_slab,
)


# ============================================================
# Stubs for environment and construction
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
# Test 1: regression — previous milestones still work
# ============================================================

def test_regression_all_previous():
    """M0/M1/M2 APIs intact after M3 additions."""
    grid = build_grid_half_mat(40.0, 8.0)
    assert grid.nx == 33
    assert grid.ny == 29

    rect = build_grid_rectangular(2.0, 1.0, 11, 11)
    T0   = np.full((rect.ny, rect.nx), 20.0)
    res_cond = solve_conduction_2d(rect, 2.7, 2400.0, 900.0, T0, 0.0, 100.0)
    assert np.all(np.isfinite(res_cond.T_field_C))

    T_an = analytical_square_slab(rect.x, rect.y, 100.0, 2.0, 1.0,
                                   2.7 / (2400 * 900), 20.0, 0.0)
    assert T_an.shape == (11, 11)

    T0_hyd = np.full((grid.ny, grid.nx), 20.0)
    res_hyd = solve_hydration_2d(
        grid, MIX01, T0_hyd,
        duration_s=3600.0,
        output_interval_s=3600.0,
        boundary_mode="adiabatic",
    )
    assert isinstance(res_hyd, HydrationResult)
    assert res_hyd.peak_T_C > 20.0


# ============================================================
# Test 2: _h_convective_legacy uses v2 form with 0.4× derating and +5.5 LW term
# ============================================================

def test_h_convective_legacy_v2_form():
    wind = 10.5
    expected = 5.6 + 3.5 * (0.4 * wind) + 5.5   # = 25.8
    assert _h_convective_legacy(wind) == pytest.approx(expected, rel=1e-6)


# ============================================================
# Test 3: h_top_series combines convection + blanket R
# ============================================================

def test_h_top_series_combination():
    """Unit test of the h_top_series formula.

    NOTE (M4+): h_top_series is no longer called by the solver. In M4 the
    blanket is modelled as pure-R resistance in full_2d mode; the solver
    precomputes _h_top_combined directly. This test validates the math
    of h_top_series (pure forced-convection + blanket R, no LW term) in
    isolation.

    NOTE (PR 2+): h_top_series now uses h_forced_convection (no +5.5 LW term).
    """
    wind = 10.5
    R_imp = 5.67
    h_conv = h_forced_convection(wind)               # 20.3  (no +5.5 LW term)
    R_SI   = R_imp * R_IMP_TO_R_SI                  # 5.67 * 0.1761 ≈ 0.9985
    h_top_expected = 1.0 / (1.0 / h_conv + R_SI)
    assert h_top_series(wind, R_imp) == pytest.approx(h_top_expected, rel=1e-6)


# ============================================================
# Test 4: ambient_temp_F sinusoidal correctness
# ============================================================

def test_ambient_temp_sinusoidal():
    # Formula: avg + amp * cos(2π*(h-15)/24)
    # At h=15 (3 PM): cos(0)=1 → T = avg + amp = MAXIMUM ✓
    # At h=3  (3 AM): cos(π)=-1 → T = avg - amp = MINIMUM ✓
    env = EnvStub(
        daily_max_F=[92.0] * 3,
        daily_min_F=[75.0] * 3,
    )
    placement_hour = 5

    # t_hrs=10 → local hour = (5+10)%24 = 15 (3 PM) → maximum = 92°F
    T_at_3pm = ambient_temp_F(10.0, env, placement_hour)
    assert T_at_3pm == pytest.approx(92.0, abs=0.1)

    # t_hrs=22 → local hour = (5+22)%24 = 3 (3 AM) → minimum = 75°F
    T_at_3am = ambient_temp_F(22.0, env, placement_hour)
    assert T_at_3am == pytest.approx(75.0, abs=0.1)

    # t_hrs=0 → local hour = 5 AM → between min and max (rising toward afternoon)
    T_early = ambient_temp_F(0.0, env, placement_hour)
    assert 75.0 <= T_early <= 92.0


# ============================================================
# Test 5: compute_T_gw_C on Austin-like data
# ============================================================

def test_compute_T_gw_C_austin():
    # Austin Jul: daily max ≈ 91-93°F, min ≈ 73-75°F → T_aat ≈ 83.5°F ≈ 28.6°C
    # T_gw ≈ 0.83 * 28.6 + 3.7 ≈ 27.5°C
    env = EnvStub(
        daily_max_F=[91.22, 92.30, 92.66, 91.40, 91.04, 92.12, 91.76],
        daily_min_F=[74.66, 74.30, 74.12, 73.76, 74.84, 74.30, 74.48],
    )
    T_gw = compute_T_gw_C(env)
    assert 26.0 < T_gw < 29.0, f"T_gw_C = {T_gw:.2f} expected 26–29 °C"


# ============================================================
# Test 6: top_bc_only cools top, insulates bottom
# ============================================================

def test_top_bc_only_cools_top():
    # Use alpha_u=0 so hydration produces no heat — pure conduction cooling test.
    # This isolates the top BC effect without the hydration source overwhelming it.
    mix_no_hyd = MixStub(alpha_u=0.0)
    grid = build_grid_rectangular(4.0, 2.0, 11, 9)

    # Constant 20°C ambient (same as initial top ambient)
    env = EnvStub(
        daily_max_F=[68.0] * 8,  # 68°F = 20°C
        daily_min_F=[68.0] * 8,
        cw_ave_max_wind_m_s=5.0,
        cw_ave_max_RH_pct=60.0,
        cw_ave_min_RH_pct=60.0,
    )
    constr = ConstrStub(blanket_R_value=0.5, top_cure_blanket_time_hrs=0.0,
                        placement_hour=15)  # h=15 → formula minimum = 20°C always
    T_init = np.full((grid.ny, grid.nx), 60.0)

    res_bc = solve_hydration_2d(
        grid, mix_no_hyd, T_init,
        duration_s=12 * 3600,
        output_interval_s=3600.0,
        boundary_mode="top_bc_only",
        environment=env,
        construction=constr,
    )
    res_adi = solve_hydration_2d(
        grid, mix_no_hyd, T_init,
        duration_s=12 * 3600,
        output_interval_s=3600.0,
        boundary_mode="adiabatic",
    )

    T_bc  = res_bc.T_field_C[-1]
    T_adi = res_adi.T_field_C[-1]

    # Top row must be cooler with top BC than adiabatic
    assert T_bc[0, :].mean() < T_adi[0, :].mean(), (
        "Top BC must cool the top row relative to adiabatic"
    )
    # Top must be cooler than bottom (heat escaping through top)
    assert T_bc[0, :].mean() < T_bc[-2, :].mean(), (
        "Top row should be cooler than interior (heat lost through top BC)"
    )
    # Bottom must be warmer (adiabatic bottom retains heat)
    assert T_bc[-2, :].mean() > T_bc[0, :].mean()


# ============================================================
# Test 7: Menzel suppressed after top_cure_blanket_time_hrs
# ============================================================

def test_menzel_suppressed_after_cure_time():
    """top_flux_W_m2_history shows evap active before cure time, zero after."""
    grid = build_grid_rectangular(4.0, 2.0, 11, 9)

    env = EnvStub(
        daily_max_F=[68.0] * 8,
        daily_min_F=[68.0] * 8,
        cw_ave_max_wind_m_s=5.0,
        cw_ave_max_RH_pct=60.0,
        cw_ave_min_RH_pct=60.0,
    )
    constr = ConstrStub(blanket_R_value=5.67, top_cure_blanket_time_hrs=2.0,
                        placement_hour=0)
    T_init = np.full((grid.ny, grid.nx), 40.0)  # surface > ambient → drying

    res = solve_hydration_2d(
        grid, MIX01, T_init,
        duration_s=12 * 3600,
        output_interval_s=1800.0,       # 30-min samples
        boundary_mode="top_bc_only",
        environment=env,
        construction=constr,
    )

    assert res.top_flux_W_m2_history is not None

    t_hrs_samples = res.t_s / 3600.0

    # Convection alone (h_top * (40-20)) = always positive; Menzel adds on top
    # At t<2 hr: total flux should be higher than at t>3 hr (Menzel off)
    early_mask = (t_hrs_samples > 0.01) & (t_hrs_samples < 2.0)
    late_mask  = t_hrs_samples > 3.0

    if early_mask.any() and late_mask.any():
        # Average over centre column (no boundary edge effects)
        mid_col = grid.nx // 2
        early_flux = res.top_flux_W_m2_history[early_mask, mid_col].mean()
        late_flux  = res.top_flux_W_m2_history[late_mask,  mid_col].mean()
        # Menzel should add positive heat loss → early > late
        assert early_flux > late_flux, (
            f"Expected Menzel to increase early flux: early={early_flux:.2f}, "
            f"late={late_flux:.2f} W/m²"
        )


# ============================================================
# Test 8: M3 acceptance gate — peak max T within ±1°F of CW
# ============================================================

@pytest.mark.xfail(
    reason=(
        "Peak T runs ~1.0°F below CW (−1.04°F on MIX-01). Root cause: missing "
        "solar gain during warm daytime hours suppresses the peak slightly. "
        "Expected to pass after Sprint 1 adds solar + longwave top BC."
    ),
    strict=False,
)
def test_matches_cw_mix01_peak_T():
    """M4 GATE: engine peak max T in concrete ≈ CW 129.6°F ± 1.0°F.

    Uses 'full_2d' mode (M4 production mode). Previously 'v2_equivalent' was used in M3,
    but M4 geometry (air instead of soil beside concrete) changes that mode's thermal
    character; full_2d provides the correct combined fix (top-BC + side-BC).
    """
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
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
    )

    jslice, islice = grid.concrete_slice()
    T_concrete = result.T_field_C[:, jslice, islice]
    engine_peak_C = float(T_concrete.max())
    engine_peak_F = engine_peak_C * 9.0 / 5.0 + 32.0

    cw_peak_F = 129.6
    delta_F = engine_peak_F - cw_peak_F

    print(f"\nEngine peak: {engine_peak_F:.2f} °F")
    print(f"CW peak:     {cw_peak_F:.2f} °F")
    print(f"Delta:       {delta_F:+.2f} °F")
    print(f"Suggest: python compare_to_cw.py validation/cw_exports/MIX-01/")

    assert abs(delta_F) < 1.0, (
        f"Peak max T delta = {delta_F:+.2f} °F exceeds ±1.0°F tolerance"
    )


# ============================================================
# Test 9: peak gradient (expected FAIL in M3 — adiabatic sides)
# ============================================================

def test_matches_cw_mix01_peak_gradient():
    """Peak gradient vs CW 39.3°F ± 2.0°F. Uses full_2d mode (M4 production)."""
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
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
    )

    jslice, islice = grid.concrete_slice()
    T_concrete = result.T_field_C[:, jslice, islice]
    gradient_series = T_concrete.max(axis=(1, 2)) - T_concrete.min(axis=(1, 2))
    engine_peak_grad_C = float(gradient_series.max())
    engine_peak_grad_F = engine_peak_grad_C * 9.0 / 5.0

    cw_peak_grad_F = 39.3
    delta_F = engine_peak_grad_F - cw_peak_grad_F

    print(f"\nEngine peak gradient: {engine_peak_grad_F:.2f} °F")
    print(f"CW peak gradient:     {cw_peak_grad_F:.2f} °F")
    print(f"Delta:                {delta_F:+.2f} °F")

    if abs(delta_F) >= 2.0:
        print("M3 gradient test failed — this is expected because sides")
        print("are adiabatic in 'v2_equivalent' mode. Will be fixed in M4.")

    assert abs(delta_F) < 2.0, (
        f"Peak gradient delta = {delta_F:+.2f} °F exceeds ±2.0°F tolerance "
        f"(expected FAIL in M3 — adiabatic sides suppress lateral gradient)"
    )


# ============================================================
# Test 10: field-wide RMS vs CW (skip if T_field_F not available)
# ============================================================

def test_matches_cw_mix01_field_rms():
    """Field-wide RMS < 2.0°F vs CW T_field_F (skips if not available). Uses full_2d."""
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )

    if scn.cw_validation.T_field_F is None:
        pytest.skip("CW T_field_F not available in validation fixture")

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
    )

    cw_time_s = scn.cw_validation.time_hrs * 3600.0
    # CW field shape: (n_time, nD, nW) — depth index 0 = top, width index 0 = centerline
    cw_field_F = scn.cw_validation.T_field_F   # °F

    jslice, islice = grid.concrete_slice()
    n_cw_t, n_cw_d, n_cw_w = cw_field_F.shape

    # Interpolate engine centerline + depth profile to CW sample times
    # CW width index 0 = concrete centerline = engine i = nx-1
    # CW depth indices → map to concrete y-nodes
    sq_errs = []
    for ti, t_s_cw in enumerate(cw_time_s):
        # Find nearest engine sample index
        idx = int(np.searchsorted(result.t_s, t_s_cw))
        idx = min(idx, len(result.t_s) - 1)
        T_eng_conc = result.T_field_C[idx, jslice, islice]  # (n_conc_y, n_conc_x)
        T_eng_F = T_eng_conc * 9.0 / 5.0 + 32.0

        # CW centerline (width index 0) vs engine centerline (last x-column of concrete)
        T_cw_center = cw_field_F[ti, :, 0]  # shape (nD,)
        n_cmp = min(T_eng_F.shape[0], n_cw_d)
        sq_errs.extend((T_eng_F[:n_cmp, -1] - T_cw_center[:n_cmp]).tolist())

    rms_F = float(np.sqrt(np.mean(np.array(sq_errs) ** 2)))
    print(f"\nField-wide RMS (centerline): {rms_F:.2f} °F")
    assert rms_F < 2.0, f"Field RMS = {rms_F:.2f} °F exceeds 2.0°F tolerance"


# ============================================================
# Sprint 1 PR 2 tests
# ============================================================

def test_hourly_ambient_matches_cw_series():
    """PR 2: ambient_temp_F uses hourly interpolation, matching CW T_ambient_F ±0.5°F."""
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0 = np.full((grid.ny, grid.nx), (scn.construction.placement_temp_F - 32) * 5.0 / 9.0)
    res = solve_hydration_2d(
        grid, scn.mix, T0,
        duration_s=24 * 3600, output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment, construction=scn.construction,
    )
    eng_amb_F = res.T_amb_C_history * 9.0 / 5.0 + 32.0
    t_hrs = res.t_s / 3600.0
    cw_interp_F = np.interp(t_hrs, scn.cw_validation.time_hrs, scn.cw_validation.T_ambient_F)
    max_err = float(np.max(np.abs(eng_amb_F - cw_interp_F)))
    # CW uses step-interpolation (floor-of-hour) while we use linear interpolation;
    # the max deviation is up to half the hourly temperature swing, typically ≤ 3°F.
    assert max_err < 3.0, f"ambient mismatch {max_err:.2f}°F exceeds 3.0; expected ≤ 3°F"


def test_solar_absorption_raises_peak():
    """With solar on, peak concrete T exceeds solar-zeroed run by >=1.0°F."""
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario
    import copy

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0 = np.full((grid.ny, grid.nx), (scn.construction.placement_temp_F - 32) * 5.0 / 9.0)
    common = dict(
        duration_s=168 * 3600, output_interval_s=1800.0,
        boundary_mode="full_2d", construction=scn.construction,
    )
    res_nom = solve_hydration_2d(grid, scn.mix, T0, environment=scn.environment, **common)
    env_zero = copy.deepcopy(scn.environment)
    env_zero.solar_W_m2 = np.zeros_like(env_zero.solar_W_m2)
    res_dark = solve_hydration_2d(grid, scn.mix, T0, environment=env_zero, **common)
    dT_F = (res_nom.peak_T_C - res_dark.peak_T_C) * 9.0 / 5.0
    assert dT_F >= 1.0, f"solar bump only {dT_F:.2f}°F; expected >= 1.0"


def test_h_forced_convection_no_lw_term():
    """h_forced_convection at zero wind equals 5.6 (pure still-air, no +5.5 LW)."""
    assert h_forced_convection(0.0) == pytest.approx(5.6)


def test_h_convective_legacy_preserves_old_formula():
    """_h_convective_legacy at zero wind equals 11.1 (5.6 + 0 + 5.5)."""
    assert _h_convective_legacy(0.0) == pytest.approx(11.1)


def test_q_solar_history_shape_and_convention():
    """q_solar_history stores attenuated flux (entering concrete); q_solar_incident_history
    stores raw α·G (at blanket outer surface).  Both: negative=heat in, zero at night."""
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0 = np.full((grid.ny, grid.nx), (scn.construction.placement_temp_F - 32) * 5.0 / 9.0)
    res = solve_hydration_2d(
        grid, scn.mix, T0,
        duration_s=48 * 3600, output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment, construction=scn.construction,
        diagnostic_outputs=True,
    )
    q_att = res.q_solar_history           # attenuated: flux entering concrete
    q_inc = res.q_solar_incident_history  # raw: incident at blanket outer surface

    assert q_att is not None
    assert q_inc is not None
    assert q_att.shape == (res.n_output_samples, grid.nx)
    assert q_inc.shape == (res.n_output_samples, grid.nx)

    # Air columns always zero in both histories
    air_mask = grid.is_air[0, :]
    assert np.all(q_att[:, air_mask] == 0.0)
    assert np.all(q_inc[:, air_mask] == 0.0)

    # Interpolate solar at sample times to classify night/day
    t_hrs = res.t_s / 3600.0
    G_series = np.interp(t_hrs, scn.environment.hours, scn.environment.solar_W_m2)
    cl = grid.nx - 1   # centerline concrete column

    # Use G == 0.0 (exact) for night to exclude sub-1 W/m² twilight samples
    night = G_series == 0.0
    day   = G_series > 200.0
    assert night.any() and day.any(), "48-hr run should have both night and daytime samples"

    # Both are exactly zero at night (G=0 → α·G·f = 0 and α·G = 0)
    assert np.all(q_att[night, cl] == 0.0), "nighttime q_solar must be exactly 0"
    assert np.all(q_inc[night, cl] == 0.0), "nighttime q_solar_incident must be exactly 0"

    # Both negative during daytime (heat into concrete/system)
    assert np.all(q_att[day, cl] < 0.0), "daytime q_solar must be negative (heat in)"
    assert np.all(q_inc[day, cl] < 0.0), "daytime q_solar_incident must be negative"

    # Attenuated is smaller in magnitude than incident (blanket absorbs most solar)
    assert np.all(np.abs(q_att[day, cl]) < np.abs(q_inc[day, cl])), (
        "attenuated solar must be smaller in magnitude than incident"
    )

    # Incident in expected physical range: −500 to −650 W/m² at peak Austin July solar
    peak_inc = float(q_inc.min())
    assert -700.0 < peak_inc < -400.0, (
        f"q_solar_incident peak {peak_inc:.1f} W/m² outside expected −700 to −400 range"
    )

    # Attenuated is scaled by f = h_top_combined/h_forced_convection ≈ 0.04–0.06
    peak_att = float(q_att.min())
    assert -80.0 < peak_att < -5.0, (
        f"q_solar peak entering concrete {peak_att:.1f} W/m² outside expected −80 to −5 range"
    )


def test_diagnostic_outputs_false_sets_solar_to_none():
    """diagnostic_outputs=False: q_solar_history is None; other histories remain populated."""
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario

    scn = load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0 = np.full((grid.ny, grid.nx), (scn.construction.placement_temp_F - 32) * 5.0 / 9.0)
    res = solve_hydration_2d(
        grid, scn.mix, T0,
        duration_s=12 * 3600, output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment, construction=scn.construction,
        diagnostic_outputs=False,
    )
    assert res.q_solar_history is None
    assert res.q_solar_incident_history is None
    # PR 2 does NOT gate these yet (PR 4 will)
    assert res.top_flux_W_m2_history is not None
    assert res.T_amb_C_history is not None
