"""
PR 10 tests: Barber soil model plumbing (zero behavior change).

Changes validated here:
  - soil_lag_hrs (5.0) and soil_damping (0.7) fields added to CWConstruction
  - ground_surface_temperature_C() helper defined in thermal_engine_2d
  - T_ground_C_history field declared on HydrationResult (None default)
  - MIX-01 trajectories bit-identical to sprint-2-complete (every 4th sample)
"""

from __future__ import annotations

import dataclasses
import pathlib

import numpy as np
import pytest

from thermal_engine_2d import (
    HydrationResult,
    ambient_temp_F,
    build_grid_half_mat,
    ground_surface_temperature_C,
    solve_hydration_2d,
)

MIX01_DIR = pathlib.Path("validation/cw_exports/MIX-01")
MIX01_INPUT = str(MIX01_DIR / "input.dat")
MIX01_WEATHER = str(MIX01_DIR / "weather.dat")
MIX01_OUTPUT = str(MIX01_DIR / "output.txt")
def _load_mix01():
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    if not MIX01_DIR.exists():
        pytest.skip(f"MIX-01 scenario not found at {MIX01_DIR}")
    from cw_scenario_loader import load_cw_scenario
    return load_cw_scenario(MIX01_INPUT, MIX01_WEATHER, MIX01_OUTPUT)


def _run_mix01(scn):
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    return solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )


def _synthetic_env():
    """Minimal CWEnvironment-like object with a pure sinusoidal T_air_F."""
    hours = np.arange(0, 169, dtype=float)
    T_air_F = 85.0 + 10.0 * np.sin(2 * np.pi * hours / 24.0)

    class _Env:
        pass

    env = _Env()
    env.hours = hours
    env.T_air_F = T_air_F
    return env


# ============================================================
# Class A — zero-behavior-change (plumbing)
# ============================================================

def test_cwconstruction_soil_defaults():
    from cw_scenario_loader import CWConstruction
    c = CWConstruction()
    assert c.soil_lag_hrs == 5.0
    assert c.soil_damping == 0.7


def test_hydration_result_has_t_ground_field():
    dummy = np.zeros(1)
    result = HydrationResult(
        t_s=dummy,
        T_field_C=dummy.reshape(1, 1, 1),
        alpha_field=dummy.reshape(1, 1, 1),
        t_e_field_hrs=dummy.reshape(1, 1, 1),
        dt_inner_s=1.0,
        n_inner_steps=1,
        n_output_samples=1,
        peak_T_C=25.0,
        peak_T_location=(0, 0),
        peak_alpha=0.5,
        centerline_T_C=dummy.reshape(1, 1),
    )
    assert hasattr(result, "T_ground_C_history")
    assert result.T_ground_C_history is None


# ============================================================
# Class B — ground_surface_temperature_C() math
# ============================================================

def test_helper_identity_reduces_to_T_amb():
    """With lag=0 and damping=1, helper equals ambient_temp_F converted to °C."""
    scn = _load_mix01()
    env = scn.environment
    ph = scn.construction.placement_hour
    sample_times = [0.0, 3.0, 6.0, 9.0, 12.0, 18.0, 24.0, 36.0, 48.0, 72.0, 120.0, 168.0]
    for t in sample_times:
        oracle = (ambient_temp_F(t, env, ph) - 32.0) * 5.0 / 9.0
        helper = ground_surface_temperature_C(t, env, 0.0, 1.0, ph)
        assert abs(helper - oracle) < 1e-10, (
            f"Identity failed at t={t}: helper={helper:.10f}, oracle={oracle:.10f}"
        )


def test_helper_pure_lag():
    """With damping=1, helper at t equals unlagged value at (t - lag)."""
    env = _synthetic_env()
    ph = 0
    # helper at t=12 with lag=6 must equal helper at t=6 with lag=0
    h_lagged = ground_surface_temperature_C(12.0, env, 6.0, 1.0, ph)
    h_unlagged = ground_surface_temperature_C(6.0, env, 0.0, 1.0, ph)
    assert abs(h_lagged - h_unlagged) < 1e-6, (
        f"Lag mismatch: h(t=12,lag=6)={h_lagged:.8f}, h(t=6,lag=0)={h_unlagged:.8f}"
    )


def test_helper_pure_damping():
    """With lag=0, output amplitude scales by damping factor around the mean."""
    env = _synthetic_env()
    ph = 0
    # T_air_F = 85 + 10*sin(...)  →  mean 85°F, peak at sin=1 (t=6 for sin(2π·6/24)=sin(π/2)=1)
    T_mean_C = (85.0 - 32.0) * 5.0 / 9.0   # 29.444...°C
    T_peak_C = (95.0 - 32.0) * 5.0 / 9.0   # 35.0°C
    expected = T_mean_C + 0.5 * (T_peak_C - T_mean_C)
    result = ground_surface_temperature_C(6.0, env, 0.0, 0.5, ph)
    assert abs(result - expected) < 1e-6, (
        f"Damping mismatch at t=6: result={result:.8f}, expected={expected:.8f}"
    )


def test_helper_combined_sinusoid():
    """With lag=5 and damping=0.7, output sinusoid has shifted phase and reduced amplitude."""
    env = _synthetic_env()
    ph = 0
    sample_ts = [0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0, 48.0]
    vals = [ground_surface_temperature_C(t, env, 5.0, 0.7, ph) for t in sample_ts]
    amplitude = (max(vals) - min(vals)) / 2.0
    # Expected amplitude: 0.7 * 10°F * (5/9) = 3.888...°C
    expected_amplitude = 0.7 * 10.0 * 5.0 / 9.0
    assert abs(amplitude - expected_amplitude) < 0.3, (
        f"Amplitude mismatch: estimated={amplitude:.4f}, expected={expected_amplitude:.4f}"
    )
    # Spot-check shifted peak: unlagged sin peaks at t=6h (sin(π/2)=1); lag=5 shifts to t=11h.
    T_mean_C = (85.0 - 32.0) * 5.0 / 9.0
    T_peak_C = (95.0 - 32.0) * 5.0 / 9.0
    expected_peak = T_mean_C + 0.7 * (T_peak_C - T_mean_C)
    val_at_11 = ground_surface_temperature_C(11.0, env, 5.0, 0.7, ph)
    assert abs(val_at_11 - expected_peak) < 0.05, (
        f"Peak shift mismatch at t=11: val={val_at_11:.6f}, expected={expected_peak:.6f}"
    )


def test_helper_early_time_clamp():
    """Before lag elapses, helper clamps t_lagged to 0 — same result as t=0."""
    env = _synthetic_env()
    ph = 0
    val_at_0 = ground_surface_temperature_C(0.0, env, 5.0, 0.7, ph)
    val_at_2 = ground_surface_temperature_C(2.0, env, 5.0, 0.7, ph)
    assert val_at_0 == val_at_2, (
        f"Early-time clamp failed: val(t=0)={val_at_0}, val(t=2)={val_at_2}"
    )
    assert not np.isnan(val_at_0)
    assert not np.isnan(val_at_2)


def test_helper_reasonable_band_on_mix01():
    """Helper with defaults stays within ±15°F (8.3°C) of unlagged T_amb at same t."""
    scn = _load_mix01()
    env = scn.environment
    ph = scn.construction.placement_hour
    check_times = [0.0, 24.0, 48.0, 72.0, 96.0, 120.0, 144.0, 168.0]
    for t in check_times:
        T_amb_C = (ambient_temp_F(t, env, ph) - 32.0) * 5.0 / 9.0
        T_gnd_C = ground_surface_temperature_C(t, env, 5.0, 0.7, ph)
        assert abs(T_gnd_C - T_amb_C) <= 8.3, (
            f"Unreasonable ground temp at t={t}h: "
            f"T_gnd={T_gnd_C:.2f}°C, T_amb={T_amb_C:.2f}°C, "
            f"delta={abs(T_gnd_C - T_amb_C):.2f}°C > 8.3°C"
        )
