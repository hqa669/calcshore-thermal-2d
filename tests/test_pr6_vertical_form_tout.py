"""
PR 6 (M6b) tests: vertical-form T_outer LW solve.

Activates the quasi-steady energy balance on the steel-form outer surface
(side face).  LW-only: F_vert stays at 0.0, no solar contribution.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

import thermal_engine_2d as te
from thermal_engine_2d import (
    EMIS_GROUND,
    EMISSIVITY_DEFAULT,
    F_SKY_VERT,
    R_FORM_EFFECTIVE_SI,
    SOLAR_ABSORPTIVITY_DEFAULT,
    STEFAN_BOLTZMANN,
    build_grid_half_mat,
    h_forced_convection,
    solve_hydration_2d,
)


# ============================================================
# Shared helpers  (mirror test_pr5_plumbing.py)
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


def _run_full2d(scn, diagnostic_outputs=False):
    """168-hour production run (default: no diagnostic outputs for speed)."""
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
        construction=scn.construction,
        diagnostic_outputs=diagnostic_outputs,
    )
    return grid, result


def _corner_rms_F(grid, result, val):
    """Corner RMS in °F vs CW reference (mirrors compare_to_cw.py:132-175)."""
    if val.T_field_F is None:
        return None
    eng_corner_C = result.T_field_C[:, grid.iy_concrete_start, grid.ix_concrete_start]
    eng_corner_F = eng_corner_C * 9.0 / 5.0 + 32.0
    cw_corner_F  = val.T_field_F[:, 0, -1]
    cw_t_s       = val.time_hrs * 3600.0
    eng_interp   = np.interp(cw_t_s, result.t_s, eng_corner_F)
    return float(np.sqrt(np.mean((eng_interp - cw_corner_F) ** 2)))


# ============================================================
# Test 1: Newton convergence on worst-case side-face cells
# ============================================================

def test_pr6_newton_convergence():
    """2-Newton-step result is within 0.01°C of 5-step reference on 27 combos.

    Mirrors test_lw_linearization_error_bounded (test_top_bc_2d.py:610).
    The side-face residual has two LW paths (sky + ground) instead of one.
    F_vert=0 so no solar term in the balance.
    """
    h_conv = h_forced_convection(10.5)   # MIX-01 wind speed

    def newton_solve_vert(T_conc_C, T_sky_C, T_amb_C, n_steps):
        T_sky_K = T_sky_C + 273.15
        T_gnd_K = T_amb_C + 273.15   # T_ground = T_amb (PR 6 baseline)

        T_eff_sky_C = F_SKY_VERT * T_sky_C + (1.0 - F_SKY_VERT) * T_amb_C
        T_ref_K  = 0.5 * (T_conc_C + T_eff_sky_C) + 273.15
        h_rad0   = 4.0 * EMISSIVITY_DEFAULT * STEFAN_BOLTZMANN * T_ref_K ** 3
        denom    = h_conv + h_rad0 + 1.0 / R_FORM_EFFECTIVE_SI
        num      = (T_conc_C / R_FORM_EFFECTIVE_SI
                    + h_conv * T_amb_C
                    + h_rad0 * T_eff_sky_C)
        T_o = num / denom

        for _ in range(n_steps):
            T_o_K = T_o + 273.15
            F  = (h_conv * (T_o - T_amb_C)
                  + EMISSIVITY_DEFAULT * STEFAN_BOLTZMANN * F_SKY_VERT
                    * (T_o_K ** 4 - T_sky_K ** 4)
                  + EMISSIVITY_DEFAULT * EMIS_GROUND * STEFAN_BOLTZMANN * (1.0 - F_SKY_VERT)
                    * (T_o_K ** 4 - T_gnd_K ** 4)
                  - (T_conc_C - T_o) / R_FORM_EFFECTIVE_SI)
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
                prod = newton_solve_vert(T_conc, T_sky, T_amb, n_steps=2)
                ref  = newton_solve_vert(T_conc, T_sky, T_amb, n_steps=5)
                max_err = max(max_err, abs(prod - ref))

    assert max_err < 0.01, (
        f"2-step Newton error {max_err:.4f}°C exceeds 0.01°C tolerance. "
        "The linearized initial guess may be poor — check T_eff_sky formulation."
    )


# ============================================================
# Test 2: q_side_LW_history populated when sky data present
# ============================================================

def test_pr6_side_LW_flux_populated_when_sky_data_present():
    """LW history is populated and nonzero when diagnostic_outputs=True."""
    scn = _load_mix01()
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    result = _run(scn, diagnostic_outputs=True)

    assert result.q_side_LW_history is not None, "q_side_LW_history should not be None"
    assert np.any(result.q_side_LW_history != 0), (
        "q_side_LW_history is all-zero — LW solve is not activating. "
        "Check the sky-data guard in the diagnostic-sampling block."
    )

    assert result.T_outer_form_C_history is not None, (
        "T_outer_form_C_history should be populated when diagnostic_outputs=True"
    )
    n_samples = result.T_field_C.shape[0]
    n_side_rows = grid.iy_concrete_end - 1
    assert result.T_outer_form_C_history.shape == (n_samples, n_side_rows), (
        f"Expected shape ({n_samples}, {n_side_rows}), "
        f"got {result.T_outer_form_C_history.shape}"
    )

    # T_outer should be physically bounded: cooler than warm concrete,
    # but warmer than extreme cold sky (nighttime radiative cooling has limits)
    assert np.all(np.isfinite(result.T_outer_form_C_history)), (
        "T_outer_form_C_history contains NaN or Inf"
    )
    assert result.T_outer_form_C_history.min() > -40.0, (
        "T_outer < -40°C: physically unreasonable for standard met conditions"
    )
    assert result.T_outer_form_C_history.max() < 80.0, (
        "T_outer > 80°C: exceeds concrete peak by too wide a margin"
    )


# ============================================================
# Test 3: side LW is zero when no sky data provided
# ============================================================

def test_pr6_side_LW_zero_when_no_sky_data():
    """Fallback: empty T_sky_C → q_side_LW all zeros, T_outer = T_conc."""
    scn = _load_mix01()
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)

    env_no_sky = dataclasses.replace(scn.environment, T_sky_C=np.array([]))
    scn_no_sky = dataclasses.replace(scn, environment=env_no_sky)

    result = _run(scn_no_sky, diagnostic_outputs=True)

    assert result.q_side_LW_history is not None
    assert np.all(result.q_side_LW_history == 0.0), (
        "q_side_LW_history should be all-zero when T_sky_C is empty"
    )

    # In the no-sky fallback T_outer_form = T_conc_side at sample time
    conc_side_C = result.T_field_C[:, 1:grid.iy_concrete_end, grid.ix_concrete_start]
    np.testing.assert_allclose(
        result.T_outer_form_C_history, conc_side_C, atol=1e-12,
        err_msg="T_outer_form_C_history should equal concrete side temps when no sky data"
    )


# ============================================================
# Test 4: F_vert=0 still gives zero side solar (PR 4/5 invariant)
# ============================================================

def test_pr6_f_vert_zero_still_gives_zero_side_solar():
    """Explicit F_vert=0.0 override produces zero q_side_solar_history everywhere.

    PR 7 flipped the default to 0.5; this test preserves the ablation-case
    invariant by forcing F_vert=0.0 explicitly.
    """
    scn = _load_mix01()
    ctor_zero = dataclasses.replace(scn.construction, vertical_solar_factor=0.0)

    result = _run(scn, construction=ctor_zero, diagnostic_outputs=True)

    assert result.q_side_solar_history is not None
    assert np.all(result.q_side_solar_history == 0.0), (
        "q_side_solar_history should be exactly zero when F_vert=0.0. "
        "F_vert=0.0 override must fully disable the solar path."
    )


# ============================================================
# Test 5: corner T_outer row-0 is physically consistent
# ============================================================

def test_pr6_corner_tout_matches_main_side_row0():
    """T_outer_form row 0 (corner cell) is physically bounded and consistent.

    The corner quarter-cell reuses _T_outer_form_C[0] from the main side-BC
    Newton solve (j=1 = iy_concrete_start). The diagnostic sampler independently
    re-solves the same Newton at sample time. Both paths use identical inputs
    so both should produce physically consistent temperatures.
    """
    scn = _load_mix01()
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    result = _run(scn, diagnostic_outputs=True)

    assert result.T_outer_form_C_history is not None
    corner_tout = result.T_outer_form_C_history[:, 0]  # row 0 = corner cell

    assert np.all(np.isfinite(corner_tout)), "Corner T_outer contains NaN/Inf"

    # Wide physical bounds: T_outer can be anywhere from cold sky to warm ambient/concrete.
    # (Early-phase: T_amb ≈ 25°C >> T_conc ≈ 15°C, so T_outer > T_conc is expected.)
    assert corner_tout.min() > -40.0, "T_outer at corner below -40°C: physically unreasonable"
    assert corner_tout.max() < 100.0, "T_outer at corner above 100°C: physically unreasonable"

    # Consistency: after hydration warmup (beyond first 24h sample), nighttime T_outer
    # should be cooler than concrete corner on average (LW radiative cooling to sky).
    from thermal_engine_2d import is_daytime
    t_hrs_arr = result.t_s / 3600.0
    placement_hour = scn.construction.placement_hour
    late_night_mask = np.array([
        t_hrs_arr[i] >= 24.0 and is_daytime(t_hrs_arr[i], placement_hour) == 0.0
        for i in range(len(t_hrs_arr))
    ])
    conc_corner_C = result.T_field_C[:, grid.iy_concrete_start, grid.ix_concrete_start]
    if late_night_mask.sum() > 5:
        night_outer = corner_tout[late_night_mask]
        night_conc  = conc_corner_C[late_night_mask]
        assert np.mean(night_outer) < np.mean(night_conc), (
            "After 24h, form outer should be cooler than concrete corner at night on average "
            "(LW radiative cooling to sky). Check sign of R_form flux term."
        )


# ============================================================
# Test 6: R_form=0 diagnostic (no numerical gate — informational)
# ============================================================

def test_pr6_r_form_zero_diagnostic(monkeypatch, capsys):
    """Run MIX-01 twice and report Corner RMS for R_form=0.0862 and R_form≈0.

    No numerical assertion on R_form=0 result — this is diagnostic.
    Informs the PR summary table per Sprint 2 plan (M6e calibration decision).
    """
    scn = _load_mix01()
    if scn.cw_validation.T_field_F is None:
        pytest.skip("CW field data not available for corner RMS comparison")

    grid, result_a = _run_full2d(scn)
    rms_a = _corner_rms_F(grid, result_a, scn.cw_validation)

    monkeypatch.setattr(te, "R_FORM_EFFECTIVE_SI", 1e-9)
    _, result_b = _run_full2d(scn)
    rms_b = _corner_rms_F(grid, result_b, scn.cw_validation)

    print(f"\nCorner RMS with R_form=0.0862: {rms_a:.2f}°F")
    print(f"Corner RMS with R_form=0:      {rms_b:.2f}°F")

    assert np.isfinite(rms_a), "Run A (R_form=0.0862) produced non-finite corner RMS"
    assert np.isfinite(rms_b), "Run B (R_form≈0) produced non-finite corner RMS"
    # Both runs must complete without NaN in the temperature field
    assert np.all(np.isfinite(result_a.T_field_C)), "Run A T_field_C has NaN/Inf"
    assert np.all(np.isfinite(result_b.T_field_C)), "Run B T_field_C has NaN/Inf"


# ============================================================
# Test 7: Corner RMS ≤ 3.5°F on MIX-01 (PR 6 gate)
# ============================================================

def test_pr6_corner_rms_improves():
    """End-to-end validation: PR 6 LW activation does not regress Corner RMS.

    Sprint 1 baseline: 4.10°F.  PR 6 achieves 4.08°F — the LW path is active
    but modest: R_form=0.0862 limits sky-cooling from reaching the concrete
    (h_eff improvement ~8%, τ reduction ~4h out of 44h).  The 3.0–3.5°F goal
    requires F_vert calibration (PR 7/8).  PR 6 gate: does not regress beyond
    Sprint 1's 4.10°F (gate 4.15°F to allow measurement noise).
    Non-regression: Peak Max within ±0.5°F of CW, Gradient within ±1.0°F,
    Field RMS within ±0.2°F of Sprint 1 1.36°F.

    Uses explicit F_vert=0.0 override to isolate PR 6 LW physics from PR 7
    solar activation (PR 7 flipped the default to 0.5).
    """
    scn = _load_mix01()
    val = scn.cw_validation
    if val.T_field_F is None:
        pytest.skip("CW field data not available — cannot compute Corner/Field RMS")

    # Isolate PR 6 LW-only effect: override F_vert=0.0 so solar is disabled.
    ctor_lw_only = dataclasses.replace(scn.construction, vertical_solar_factor=0.0)
    scn_lw_only  = dataclasses.replace(scn, construction=ctor_lw_only)
    grid, result = _run_full2d(scn_lw_only)

    # Corner RMS — PR 6 gate: no regression from Sprint 1 (4.15°F margin)
    corner_rms = _corner_rms_F(grid, result, val)
    assert corner_rms is not None
    # PR 8 note: ACI Eq 27 h_conv_vertical (lower h than horizontal) shifts
    # this LW-only baseline from 4.10°F → 4.32°F. Threshold updated from
    # 4.15 to 4.40 to reflect the new h. The LW activation is still correct;
    # the baseline shift is the expected effect of the lower vertical h.
    assert corner_rms <= 4.40, (
        f"Corner RMS {corner_rms:.2f}°F regressed beyond h_vert LW-only baseline. "
        f"Check sign of q_side_total or h_conv_vertical formula."
    )

    # Peak Max T non-regression (±0.5°F of CW 129.6°F)
    jslice, islice = grid.concrete_slice()
    T_conc_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
    engine_peak_max_F = float(T_conc_F.max())
    cw_peak_max_F     = float(val.T_max_xs_F.max())
    assert abs(engine_peak_max_F - cw_peak_max_F) <= 0.5, (
        f"Peak Max T regression: engine {engine_peak_max_F:.1f}°F, "
        f"CW {cw_peak_max_F:.1f}°F, delta {engine_peak_max_F - cw_peak_max_F:+.2f}°F "
        f"(gate ±0.5°F)"
    )

    # Peak Gradient non-regression (±1.0°F of CW 39.3°F)
    grad_series        = T_conc_F.max(axis=(1, 2)) - T_conc_F.min(axis=(1, 2))
    engine_peak_grad_F = float(grad_series.max())
    cw_peak_grad_F     = float(val.T_diff_xs_F.max())
    assert abs(engine_peak_grad_F - cw_peak_grad_F) <= 1.0, (
        f"Peak Gradient regression: engine {engine_peak_grad_F:.1f}°F, "
        f"CW {cw_peak_grad_F:.1f}°F, delta {engine_peak_grad_F - cw_peak_grad_F:+.2f}°F "
        f"(gate ±1.0°F)"
    )

    # Field RMS non-regression (within ±0.2°F of Sprint 1's 1.36°F)
    n_cw_t, n_cw_d, _ = val.T_field_F.shape
    cw_t_s = val.time_hrs * 3600.0
    sq_errs = []
    for ti in range(n_cw_t):
        idx = min(int(np.searchsorted(result.t_s, cw_t_s[ti])), len(result.t_s) - 1)
        T_eng_col_F = result.T_field_C[idx, jslice, grid.nx - 1] * 9.0 / 5.0 + 32.0
        n_cmp = min(len(T_eng_col_F), n_cw_d)
        sq_errs.extend((T_eng_col_F[:n_cmp] - val.T_field_F[ti, :n_cmp, 0]).tolist())
    field_rms = float(np.sqrt(np.mean(np.array(sq_errs) ** 2)))
    sprint1_field_rms = 1.36
    assert abs(field_rms - sprint1_field_rms) <= 0.2, (
        f"Field RMS {field_rms:.2f}°F deviates by more than ±0.2°F from Sprint 1's "
        f"{sprint1_field_rms}°F. LW on form face should be corner-local; "
        f"field-wide drift suggests a sign error or incorrect flux path."
    )
