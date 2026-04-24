"""
PR 5 (M6a) plumbing tests: verify new constants and HydrationResult field
are inert — no solver behavior change vs Sprint 1 baseline.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d


# ============================================================
# Shared helpers
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


# ============================================================
# Test 1: constants importable with exact values
# ============================================================

def test_pr5_constants_declared():
    from thermal_engine_2d import F_SKY_VERT, EMIS_GROUND
    assert F_SKY_VERT == 0.5
    assert EMIS_GROUND == 0.95


# ============================================================
# Test 2: T_outer_form_C_history populated by diagnostic mode
# (PR 5 declared the field; PR 6 activates it)
# ============================================================

def test_pr5_hydration_result_new_field_populated_by_diag():
    scn = _load_mix01()
    result_diag = _run(scn, diagnostic_outputs=True)
    assert result_diag.T_outer_form_C_history is not None, (
        "T_outer_form_C_history should be populated when diagnostic_outputs=True (PR 6)"
    )
    result_no_diag = _run(scn, diagnostic_outputs=False)
    assert result_no_diag.T_outer_form_C_history is None, (
        "T_outer_form_C_history should be None when diagnostic_outputs=False"
    )


# ============================================================
# Test 3: sentinel bit-identity (load-bearing PR 5 test)
# ============================================================

def test_pr5_sentinel_solar_absorptivity_inert_under_zero_fvert():
    """Solar absorptivity is still inert when F_vert=0.0 after PR 6.

    PR 6 activates LW on the form face (emissivity-dependent), so emissivity
    is no longer sentinel-eligible. But solar absorptivity is still gated by
    F_vert: alpha_sol * F_vert = 0 whenever F_vert=0.0, regardless of alpha.
    This test guards that invariant — a regression here means PR 7 or later
    accidentally activated solar without enabling F_vert.
    """
    scn = _load_mix01()

    assert scn.construction.vertical_solar_factor == 0.0, (
        "Sentinel test assumes F_vert=0.0 default — if this changes, the "
        "invariant (alpha_side*F_vert=0 regardless of absorptivity) breaks."
    )

    baseline = _run(scn)

    # Only vary solar absorptivity — emissivity now affects LW and is intentionally
    # physics-coupled in PR 6; changing it SHOULD change T_field_C.
    sentinel_ctor = dataclasses.replace(
        scn.construction,
        solar_absorptivity_side=0.7,    # nonstandard; default is 0.65
    )
    sentinel = _run(scn, construction=sentinel_ctor)

    # Bit-identical for solar: alpha * F_vert = 0 so the path is dead.
    np.testing.assert_array_equal(
        baseline.T_field_C, sentinel.T_field_C,
        err_msg="T_field_C differs when only solar_absorptivity_side changes under "
                "F_vert=0.0 — the solar path is being activated unintentionally.",
    )
    np.testing.assert_array_equal(
        baseline.alpha_field, sentinel.alpha_field,
        err_msg="alpha_field differs — same bug class.",
    )
    assert baseline.peak_T_C == sentinel.peak_T_C
    assert baseline.peak_T_location == sentinel.peak_T_location
