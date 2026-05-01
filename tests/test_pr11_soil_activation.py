"""
PR 11 tests: Barber soil model activation.

Changes validated here:
  - ground_surface_temperature_C() wired into both form-face Newton sites
  - T_ground_C_history populated at sample times (full_2d only; None otherwise)
  - Ablation (damping=1.0, lag=0.0) reproduces sprint-2-complete bit-identically
  - Defaults (damping=0.7, lag=5.0) produce a measurably different trajectory
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from thermal_engine_2d import (
    build_grid_half_mat,
    solve_hydration_2d,
)

MIX01_DIR = pathlib.Path("validation/cw_exports/MIX-01")
MIX01_INPUT = str(MIX01_DIR / "input.dat")
MIX01_WEATHER = str(MIX01_DIR / "weather.dat")
MIX01_OUTPUT = str(MIX01_DIR / "output.txt")
FIXTURE_PATH = pathlib.Path("tests/fixtures/sprint_2_complete_mix01_T_field.npz")


def _load_mix01():
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    if not MIX01_DIR.exists():
        pytest.skip(f"MIX-01 scenario not found at {MIX01_DIR}")
    from cw_scenario_loader import load_cw_scenario
    return load_cw_scenario(MIX01_INPUT, MIX01_WEATHER, MIX01_OUTPUT)


def _run_mix01(scn, boundary_mode="full_2d", duration_hrs=168):
    # model_soil=True: Barber soil tests require the transient soil mesh.
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft, model_soil=True)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    return solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=duration_hrs * 3600,
        output_interval_s=1800.0,
        boundary_mode=boundary_mode,
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )


# ============================================================
# Class A — ablation: degenerate params must reproduce sprint-2-complete
# ============================================================

def test_ablation_reproduces_sprint2():
    """damping=1.0, lag=0.0 are identity params — must reproduce sprint-2-complete."""
    if not FIXTURE_PATH.exists():
        pytest.skip(f"Fixture not found at {FIXTURE_PATH} — run fixture-gen script first")
    scn = _load_mix01()
    scn.construction.soil_damping = 1.0
    scn.construction.soil_lag_hrs = 0.0
    result = _run_mix01(scn)
    fixture = np.load(FIXTURE_PATH)
    np.testing.assert_allclose(
        result.T_field_C[::4], fixture["T_field_C"],
        rtol=0, atol=1e-12,
        err_msg=(
            "Ablation failed: damping=1.0/lag=0.0 did not reproduce sprint-2-complete. "
            "Barber wiring has a math bug."
        ),
    )


# ============================================================
# Class B — activation: defaults must change the trajectory
# ============================================================

def test_nondefault_params_differ_from_sprint2():
    """PR 11 Barber params (damping=0.7, lag=5.0) must diverge from sprint-2-complete.
    Params set explicitly — PR 17 changed defaults to no-op (0.0, 1.0)."""
    if not FIXTURE_PATH.exists():
        pytest.skip(f"Fixture not found at {FIXTURE_PATH} — run fixture-gen script first")
    scn = _load_mix01()
    scn.construction.soil_lag_hrs = 5.0
    scn.construction.soil_damping = 0.7
    result = _run_mix01(scn)
    fixture = np.load(FIXTURE_PATH)
    max_diff = float(np.max(np.abs(result.T_field_C[::4] - fixture["T_field_C"])))
    assert max_diff > 0.01, (
        f"Barber model appears inactive — max |ΔT_field_C| across all nodes/times "
        f"is only {max_diff:.6f}°C (expected > 0.01°C with defaults)."
    )


def test_t_ground_history_populated():
    """T_ground_C_history must be populated with correct shape and sane values."""
    scn = _load_mix01()
    result = _run_mix01(scn)
    assert result.T_ground_C_history is not None, "T_ground_C_history is None in full_2d mode"
    assert result.T_ground_C_history.shape == (result.n_output_samples,), (
        f"Shape mismatch: {result.T_ground_C_history.shape} vs ({result.n_output_samples},)"
    )
    diff = np.abs(result.T_ground_C_history - result.T_amb_C_history)
    max_diff_C = float(diff.max())
    limit_C = 15.0 * 5.0 / 9.0   # 15°F → °C
    assert max_diff_C <= limit_C, (
        f"T_ground_C_history is unreasonably far from T_amb_C_history: "
        f"max |ΔT| = {max_diff_C:.2f}°C (limit {limit_C:.2f}°C / 15°F)"
    )


def test_t_ground_lagged_damped_with_pr11_params():
    """End-to-end: T_ground amplitude ≈ 0.7 · T_amb amplitude in [48, 168]hr window.
    Params set explicitly (lag=5.0, damping=0.7) — PR 17 changed defaults to no-op."""
    scn = _load_mix01()
    scn.construction.soil_lag_hrs = 5.0
    scn.construction.soil_damping = 0.7
    result = _run_mix01(scn)
    t_hrs = result.t_s / 3600.0
    mask = (t_hrs >= 48.0) & (t_hrs <= 168.0)
    T_amb = result.T_amb_C_history[mask]
    T_gnd = result.T_ground_C_history[mask]
    amp_amb = (float(T_amb.max()) - float(T_amb.min())) / 2.0
    amp_gnd = (float(T_gnd.max()) - float(T_gnd.min())) / 2.0
    expected = 0.7 * amp_amb
    rel_err = abs(amp_gnd - expected) / expected if expected > 0 else float("inf")
    assert rel_err < 0.05, (
        f"Amplitude damping end-to-end: amp_gnd={amp_gnd:.4f}°C, "
        f"expected 0.7 * amp_amb={expected:.4f}°C, rel_err={rel_err:.3f} (limit 0.05)"
    )


# ============================================================
# Class C — mode gating
# ============================================================

def test_t_ground_history_none_in_adiabatic():
    """T_ground_C_history must be None when boundary_mode='adiabatic'."""
    scn = _load_mix01()
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=24 * 3600,
        output_interval_s=3600.0,
        boundary_mode="adiabatic",
    )
    assert result.T_ground_C_history is None, (
        "T_ground_C_history should be None in adiabatic mode"
    )
