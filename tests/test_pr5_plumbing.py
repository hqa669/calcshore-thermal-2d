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

def test_pr7_solar_absorptivity_now_affects_output():
    """After PR 7 (F_vert=0.5 default), solar_absorptivity_side must
    change side-face solar flux. Inverse of the PR 5/6 sentinel contract.

    PR 5/6 proved the solar path was inert under F_vert=0.0. PR 7 flips the
    default to 0.5, so the solar term is now live and alpha must matter.
    """
    scn = _load_mix01()
    r_base = _run(scn)
    c2 = dataclasses.replace(scn.construction, solar_absorptivity_side=0.7)
    r_mod = _run(scn, construction=c2)
    assert not np.array_equal(
        r_base.q_side_solar_history, r_mod.q_side_solar_history
    ), "solar_absorptivity_side must affect q_side_solar_history at F_vert=0.5"
