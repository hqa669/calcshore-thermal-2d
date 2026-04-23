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
# Test 2: new HydrationResult field stays None in all code paths
# ============================================================

def test_pr5_hydration_result_new_field_default_none():
    scn = _load_mix01()
    for diag in (True, False):
        result = _run(scn, diagnostic_outputs=diag)
        assert result.T_outer_form_C_history is None, (
            f"T_outer_form_C_history should be None in PR 5 "
            f"(diagnostic_outputs={diag})"
        )


# ============================================================
# Test 3: sentinel bit-identity (load-bearing PR 5 test)
# ============================================================

def test_pr5_sentinel_bit_identical():
    """Prove Sprint-1 side-face fields are inert under default F_vert=0.0.

    This is PR 5's load-bearing test. If it fails, the claim 'PR 5 changes
    nothing' is false, and PR 6's isolation of LW physics from plumbing
    bugs is compromised. Debug before proceeding to PR 6.
    """
    scn = _load_mix01()

    assert scn.construction.vertical_solar_factor == 0.0, (
        "Sentinel test assumes F_vert=0.0 default — if this changes, the "
        "invariant (alpha_side*F_vert=0 regardless of absorptivity) breaks."
    )

    baseline = _run(scn)

    sentinel_ctor = dataclasses.replace(
        scn.construction,
        emissivity_side=0.9,            # nonstandard; default is 0.88
        solar_absorptivity_side=0.7,    # nonstandard; default is 0.65
    )
    sentinel = _run(scn, construction=sentinel_ctor)

    # Bit-identical: no tolerance. If this fails, investigate floating-point
    # side channels — do NOT paper over with assert_allclose.
    np.testing.assert_array_equal(
        baseline.T_field_C, sentinel.T_field_C,
        err_msg="T_field_C differs under sentinel side-face fields — "
                "PR 5 plumbing has accidental coupling under F_vert=0.0.",
    )
    np.testing.assert_array_equal(
        baseline.alpha_field, sentinel.alpha_field,
        err_msg="alpha_field differs under sentinel — same bug class.",
    )
    assert baseline.peak_T_C == sentinel.peak_T_C
    assert baseline.peak_T_location == sentinel.peak_T_location
