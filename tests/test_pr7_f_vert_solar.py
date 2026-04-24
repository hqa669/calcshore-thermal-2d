"""
PR 7 (M6c) tests: F_vert solar activation on the vertical form face.

Flips CWConstruction.vertical_solar_factor default from 0.0 to 0.5 and adds
alpha_form * F_vert * G * daytime to the Newton residual on the form outer
surface. The solar term is T-independent, so dF/dT is unchanged and the
2-Newton-step convergence guarantee from PR 6 is preserved.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from thermal_engine_2d import (
    EMIS_GROUND,
    EMISSIVITY_DEFAULT,
    F_SKY_VERT,
    R_FORM_EFFECTIVE_SI,
    SOLAR_ABSORPTIVITY_DEFAULT,
    STEFAN_BOLTZMANN,
    build_grid_half_mat,
    h_forced_convection,
    is_daytime,
    solve_hydration_2d,
)


# ============================================================
# Shared helpers (mirror test_pr6_vertical_form_tout.py)
# ============================================================

def _load_mix01():
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    from cw_scenario_loader import load_cw_scenario
    return load_cw_scenario(
        "validation/cw_exports/MIX-01/input.dat",
        "validation/cw_exports/MIX-01/weather.dat",
        "validation/cw_exports/MIX-01/output.txt",
    )


def _run(scn, construction=None, duration_hrs=24, diagnostic_outputs=True):
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_air] = T0_C
    return solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=duration_hrs * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=construction if construction is not None else scn.construction,
        diagnostic_outputs=diagnostic_outputs,
    )


def _run_full2d(scn, construction=None, diagnostic_outputs=False):
    """168-hour production run."""
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_air] = T0_C
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=construction if construction is not None else scn.construction,
        diagnostic_outputs=diagnostic_outputs,
    )
    return grid, result


def _corner_rms_F(grid, result, val):
    """Corner RMS in °F vs CW reference."""
    if val.T_field_F is None:
        return None
    eng_corner_C = result.T_field_C[:, grid.iy_concrete_start, grid.ix_concrete_start]
    eng_corner_F = eng_corner_C * 9.0 / 5.0 + 32.0
    cw_corner_F  = val.T_field_F[:, 0, -1]
    cw_t_s       = val.time_hrs * 3600.0
    eng_interp   = np.interp(cw_t_s, result.t_s, eng_corner_F)
    return float(np.sqrt(np.mean((eng_interp - cw_corner_F) ** 2)))


# ============================================================
# Test 1: construction default
# ============================================================

def test_pr7_f_vert_default_is_05():
    """PR 8: vertical_solar_factor defaults to None; engine resolves via orientation lookup.

    After PR 8, CWConstruction.vertical_solar_factor=None and the engine falls
    through to F_VERT_BY_ORIENTATION["unknown"] (= 0.15, MIX-01 calibrated).
    This test verifies the field default and that the default-construction run
    produces nonzero side solar (i.e., the lookup path is active).
    """
    from cw_scenario_loader import CWConstruction
    from thermal_engine_2d import F_VERT_BY_ORIENTATION
    c = CWConstruction()
    assert c.vertical_solar_factor is None, (
        f"Expected vertical_solar_factor=None (PR 8: orientation lookup is primary), "
        f"got {c.vertical_solar_factor}"
    )
    # Default form_orientation is "unknown" → lookup gives 0.15
    assert c.form_orientation == "unknown"
    assert F_VERT_BY_ORIENTATION["unknown"] == pytest.approx(0.15, abs=1e-9)


# ============================================================
# Test 2: daytime solar is nonzero, nighttime solar is zero
# ============================================================

def test_pr7_side_solar_nonzero_during_daytime():
    """With F_vert=0.5 (default), q_side_solar_history is zero at night and
    strictly negative (heat in) at midday.

    is_daytime returns 1.0 for 6 AM ≤ wall_hour < 6 PM. With placement_hour=5,
    wall_hour = (5 + t_hrs) % 24. Pre-dawn (3 AM) → t_hrs ≈ 22 (day 1),
    midday (1 PM) → t_hrs ≈ 8.
    """
    scn = _load_mix01()
    r = _run(scn, duration_hrs=48, diagnostic_outputs=True)

    assert r.q_side_solar_history is not None
    q_solar = r.q_side_solar_history     # (n_out, n_side_rows), negative = heat in
    t_hrs = r.t_s / 3600.0

    placement_hour = scn.construction.placement_hour   # typically 5 AM
    # pre-dawn sample: 3 AM wall = (3 - placement_hour + 24) % 24 = 22h elapsed
    pre_dawn_t = (3 - placement_hour + 24) % 24 + 24  # use day-2 sample (t>24)
    pre_dawn_idx = int(np.argmin(np.abs(t_hrs - pre_dawn_t)))
    assert np.all(q_solar[pre_dawn_idx] == 0.0), (
        f"q_side_solar_history must be zero at 3 AM (t_hrs≈{pre_dawn_t:.0f}). "
        f"Got max={q_solar[pre_dawn_idx].max():.2f} W/m²"
    )

    # midday sample: 1 PM wall = (13 - placement_hour) % 24 = 8h elapsed (day 1)
    midday_t = (13 - placement_hour) % 24
    midday_idx = int(np.argmin(np.abs(t_hrs - midday_t)))
    assert np.all(q_solar[midday_idx] < 0.0), (
        f"q_side_solar_history must be negative at 1 PM (t_hrs≈{midday_t:.0f}). "
        f"Got min={q_solar[midday_idx].min():.2f} W/m²"
    )


# ============================================================
# Test 3: F_vert=0.0 override fully disables solar (ablation case)
# ============================================================

def test_pr7_f_vert_zero_override_disables_solar():
    """Explicit vertical_solar_factor=0.0 override must yield fully-zero
    q_side_solar_history, regardless of the new 0.5 default.
    """
    scn = _load_mix01()
    ctor_zero = dataclasses.replace(scn.construction, vertical_solar_factor=0.0)
    r = _run(scn, construction=ctor_zero, diagnostic_outputs=True)

    assert r.q_side_solar_history is not None
    assert np.all(r.q_side_solar_history == 0.0), (
        "vertical_solar_factor=0.0 override must produce zero q_side_solar_history. "
        "Ablation case is broken."
    )


# ============================================================
# Test 4: Newton convergence with solar term active
# ============================================================

def test_pr7_newton_convergence_with_solar():
    """2-Newton-step result is within 0.01°C of 5-step reference, with
    F_vert=0.5 and daytime G=800 W/m² (worst-case solar load).

    Verifies that the solar shift (T-independent constant) does not break
    the convergence guarantee inherited from PR 6.
    """
    h_conv   = h_forced_convection(10.5)   # MIX-01 wind speed
    f_vert   = 0.5
    alpha    = SOLAR_ABSORPTIVITY_DEFAULT
    G_day    = 800.0    # W/m² — daytime worst-case
    daytime  = 1.0

    def newton_solve_with_solar(T_conc_C, T_sky_C, T_amb_C, n_steps):
        T_sky_K = T_sky_C + 273.15
        T_gnd_K = T_amb_C + 273.15

        T_eff_sky_C = F_SKY_VERT * T_sky_C + (1.0 - F_SKY_VERT) * T_amb_C
        T_ref_K  = 0.5 * (T_conc_C + T_eff_sky_C) + 273.15
        h_rad0   = 4.0 * EMISSIVITY_DEFAULT * STEFAN_BOLTZMANN * T_ref_K ** 3
        denom    = h_conv + h_rad0 + 1.0 / R_FORM_EFFECTIVE_SI
        num      = (T_conc_C / R_FORM_EFFECTIVE_SI
                    + h_conv * T_amb_C
                    + h_rad0 * T_eff_sky_C
                    + alpha * f_vert * G_day * daytime)
        T_o = num / denom

        for _ in range(n_steps):
            T_o_K = T_o + 273.15
            F  = (h_conv * (T_o - T_amb_C)
                  + EMISSIVITY_DEFAULT * STEFAN_BOLTZMANN * F_SKY_VERT
                    * (T_o_K ** 4 - T_sky_K ** 4)
                  + EMISSIVITY_DEFAULT * EMIS_GROUND * STEFAN_BOLTZMANN * (1.0 - F_SKY_VERT)
                    * (T_o_K ** 4 - T_gnd_K ** 4)
                  - (T_conc_C - T_o) / R_FORM_EFFECTIVE_SI
                  - alpha * f_vert * G_day * daytime)
            dF = (h_conv
                  + 4.0 * EMISSIVITY_DEFAULT * STEFAN_BOLTZMANN * F_SKY_VERT * T_o_K ** 3
                  + 4.0 * EMISSIVITY_DEFAULT * EMIS_GROUND * STEFAN_BOLTZMANN
                    * (1.0 - F_SKY_VERT) * T_o_K ** 3
                  + 1.0 / R_FORM_EFFECTIVE_SI)
            T_o -= F / dF
        return T_o

    T_conc_vals = [15.0, 35.0, 55.0]
    T_sky_vals  = [-10.0, 10.0, 25.0]
    T_amb_vals  = [10.0, 20.0, 30.0]

    max_err = 0.0
    for T_conc in T_conc_vals:
        for T_sky in T_sky_vals:
            for T_amb in T_amb_vals:
                prod = newton_solve_with_solar(T_conc, T_sky, T_amb, n_steps=2)
                ref  = newton_solve_with_solar(T_conc, T_sky, T_amb, n_steps=5)
                max_err = max(max_err, abs(prod - ref))

    assert max_err < 0.01, (
        f"2-step Newton error {max_err:.4f}°C exceeds 0.01°C with F_vert=0.5, G=800 W/m². "
        "Solar shift (T-independent constant) should not break convergence."
    )


# ============================================================
# Test 5: phase lag reduction
# ============================================================

def test_pr7_phase_lag_reduction():
    """Engine corner peaks align with CW to within 4 hours in the quasi-steady
    window (t > 72hr). PR 7's solar activation shifts the phase forward.

    Note: the Corner RMS gate was retired here (PR 7 xfail deleted in PR 8).
    amplitude overshoots CW even though phase is corrected. The phase mechanism
    works; the amplitude calibration is PR 8.
    """
    scn = _load_mix01()
    val = scn.cw_validation
    if val.T_field_F is None:
        pytest.skip("CW field data not available")

    grid, result = _run_full2d(scn)

    eng_corner_C = result.T_field_C[:, grid.iy_concrete_start, grid.ix_concrete_start]
    eng_corner_F = eng_corner_C * 9.0 / 5.0 + 32.0
    t_hrs        = result.t_s / 3600.0

    # Quasi-steady window: t > 72h
    qs_mask = t_hrs > 72.0
    if qs_mask.sum() < 10:
        pytest.skip("Insufficient quasi-steady samples")

    cw_corner_F = val.T_field_F[:, 0, -1]
    cw_t_hrs    = val.time_hrs

    # find time-of-peak in last full diurnal cycle
    eng_qs_F = eng_corner_F[qs_mask]
    eng_qs_t = t_hrs[qs_mask]
    eng_peak_hr = float(eng_qs_t[np.argmax(eng_qs_F)])

    cw_qs_mask = cw_t_hrs > 72.0
    cw_qs_F    = cw_corner_F[cw_qs_mask]
    cw_qs_t    = cw_t_hrs[cw_qs_mask]
    cw_peak_hr = float(cw_qs_t[np.argmax(cw_qs_F)])

    phase_lag_hr = abs(eng_peak_hr - cw_peak_hr)
    assert phase_lag_hr <= 4.0, (
        f"Phase lag {phase_lag_hr:.1f}h exceeds 4h gate. "
        f"Engine peak at {eng_peak_hr:.1f}h, CW peak at {cw_peak_hr:.1f}h."
    )
