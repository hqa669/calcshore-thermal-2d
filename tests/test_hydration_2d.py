"""
M2 tests: hydration heat generation + variable properties.
All 9 tests listed in the M2 spec must pass.  The MIX-01 fixture is defined
as a local dataclass stub — no import from cw_scenario_loader.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pytest

from thermal_engine_2d import (
    BLANKET_PROPERTIES_2D,
    LB_FT3_TO_KG_M3,
    LB_YD3_TO_KG_M3,
    BTU_HR_FT_F_TO_W_M_K,
    BTU_LB_F_TO_J_KG_K,
    R_IMP_TO_R_SI,
    HydrationResult,
    arrhenius_vec,
    build_grid_half_mat,
    build_grid_rectangular,
    hydration_alpha_vec,
    hydration_rate_vec,
    solve_conduction_2d,
    solve_hydration_2d,
    specific_heat_variable,
    thermal_conductivity_variable,
    analytical_square_slab,
)


# ============================================================
# MIX-01 fixture stub
# ============================================================

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
# In-test 1D reference solver (used only by test_centerline_matches_1d_v2)
# ============================================================
# Implements the same physics as solve_hydration_2d on a single vertical
# column, with adiabatic (zero-flux) BCs on top and bottom.  This is the
# "fallback" oracle described in the M2 spec: if the 2D stencil is correct,
# its centerline column must match this 1D solver.

def _solve_1d_reference(
    y_m: np.ndarray,
    material_id_1d: np.ndarray,
    mix: MixStub,
    T0_C: np.ndarray,
    duration_s: float,
    output_interval_s: float,
    blanket_R_value: float = 5.67,
    blanket_thickness_m: float = 0.02,
) -> np.ndarray:
    """1D adiabatic solver with hydration + variable k/Cp.

    Returns T_out of shape (n_samples, ny).  Uses the same physics
    functions (arrhenius_vec, hydration_alpha_vec, etc.) as the 2D solver.

    Parameters
    ----------
    y_m : (ny,) node positions in metres.
    material_id_1d : (ny,) integer material labels (0=blanket,1=concrete,2=soil).
    mix : MixStub  hydration/thermal parameters.
    T0_C : (ny,) initial temperature °C.
    duration_s : float  total duration seconds.
    output_interval_s : float  sample interval seconds.
    blanket_R_value : float  hr·ft²·°F/BTU.
    blanket_thickness_m : float  metres.

    Returns
    -------
    np.ndarray  shape (n_samples, ny)  temperature °C.
    """
    from thermal_engine_2d import (
        SOIL_PROPERTIES_2D, R_GAS,
    )

    ny = len(y_m)
    is_concrete = material_id_1d == 1
    is_soil     = material_id_1d == 2
    is_blanket  = material_id_1d == 0

    # Unit conversions (same as production solver)
    Hu    = mix.Hu_J_kg
    tau   = mix.tau_hrs
    beta_h = mix.beta
    au    = mix.alpha_u
    Ea    = mix.activation_energy_J_mol
    rho_c = mix.concrete_density_lb_ft3 * LB_FT3_TO_KG_M3
    Cc    = mix.total_cementitious_lb_yd3 * LB_YD3_TO_KG_M3
    k_uc  = mix.thermal_conductivity_BTU_hr_ft_F * BTU_HR_FT_F_TO_W_M_K
    Ca    = mix.aggregate_Cp_BTU_lb_F * BTU_LB_F_TO_J_KG_K
    Wc    = mix.cement_type_I_II_lb_yd3 * LB_YD3_TO_KG_M3
    Ww    = mix.water_lb_yd3 * LB_YD3_TO_KG_M3
    Wa    = (mix.coarse_agg_lb_yd3 + mix.fine_agg_lb_yd3) * LB_YD3_TO_KG_M3

    R_SI = blanket_R_value * R_IMP_TO_R_SI
    k_blanket = blanket_thickness_m / R_SI

    soil = SOIL_PROPERTIES_2D

    # Per-node material properties
    rho = np.empty(ny)
    rho[is_concrete] = rho_c
    rho[is_soil]     = soil['rho']
    rho[is_blanket]  = BLANKET_PROPERTIES_2D['rho']

    Cp = np.empty(ny)
    Cp[is_soil]    = soil['Cp']
    Cp[is_blanket] = BLANKET_PROPERTIES_2D['Cp']
    Cp[is_concrete] = specific_heat_variable(
        Wc, Wa, Ww, Ca, np.zeros(is_concrete.sum()), T0_C[is_concrete], rho_c
    )

    k_arr = np.empty(ny)
    k_arr[is_concrete] = thermal_conductivity_variable(k_uc, 0.0)
    k_arr[is_soil]     = soil['k']
    k_arr[is_blanket]  = k_blanket

    # State
    T      = T0_C.astype(np.float64, copy=True)
    te     = np.full(ny, 0.01)
    alpha  = np.zeros(ny)

    # CFL
    dy_all = np.diff(y_m)
    dy_min = float(dy_all.min())
    alpha_diff_max = max(
        1.33 * k_uc / (rho_c * 800.0),
        soil['k'] / (soil['rho'] * soil['Cp']),
        k_blanket / (BLANKET_PROPERTIES_2D['rho'] * BLANKET_PROPERTIES_2D['Cp']),
    )
    dt_inner = 0.4 * 0.5 * dy_min**2 / alpha_diff_max

    # Output schedule
    n_samples = int(round(duration_s / output_interval_s)) + 1
    t_samples = np.linspace(0.0, duration_s, n_samples)
    t_samples[-1] = duration_s

    T_out = np.empty((n_samples, ny))
    T_out[0] = T
    next_idx = 1
    t = 0.0
    T_new = np.empty_like(T)

    while next_idx < n_samples:
        dt_step = dt_inner
        if t + dt_step > t_samples[next_idx] + 1e-10:
            dt_step = t_samples[next_idx] - t

        # Hydration on concrete
        T_K = T + 273.15
        af = arrhenius_vec(T_K[is_concrete], Ea)
        te[is_concrete] += (dt_step / 3600.0) * af
        alpha[is_concrete] = hydration_alpha_vec(te[is_concrete], tau, beta_h, au)
        da_dte = hydration_rate_vec(te[is_concrete], tau, beta_h, au)
        da_dt_per_s = da_dte * af / 3600.0

        Q = np.zeros(ny)
        Q[is_concrete] = Hu * Cc * da_dt_per_s

        # Variable k, Cp on concrete
        k_arr[is_concrete] = thermal_conductivity_variable(k_uc, alpha[is_concrete])
        Cp[is_concrete]    = specific_heat_variable(
            Wc, Wa, Ww, Ca, alpha[is_concrete], T[is_concrete], rho_c
        )

        # Harmonic-mean face k
        k_face = 2.0 * k_arr[:-1] * k_arr[1:] / (k_arr[:-1] + k_arr[1:])

        # Interior update (non-uniform dy)
        for j in range(1, ny - 1):
            dym = y_m[j]   - y_m[j - 1]
            dyp = y_m[j + 1] - y_m[j]
            flux_plus  = k_face[j]   * (T[j + 1] - T[j]) / dyp
            flux_minus = k_face[j-1] * (T[j] - T[j - 1]) / dym
            dflux = (flux_plus - flux_minus) * 2.0 / (dym + dyp)
            T_new[j] = T[j] + dt_step * (dflux + Q[j]) / (rho[j] * Cp[j])

        # Adiabatic BCs (ghost-node mirror)
        T_new[0]  = T_new[1]
        T_new[-1] = T_new[-2]

        T, T_new = T_new, T
        t += dt_step

        if abs(t - t_samples[next_idx]) < 1e-6:
            T_out[next_idx] = T
            next_idx += 1

    return T_out


# ============================================================
# Tests
# ============================================================

def test_regression_m0_m1():
    """M0/M1 API still importable and functional after M2 additions."""
    grid = build_grid_half_mat(40.0, 8.0)
    assert grid.nx == 33
    assert grid.ny == 29

    rect = build_grid_rectangular(2.0, 1.0, 11, 11)
    T0   = np.full((rect.ny, rect.nx), 20.0)
    res  = solve_conduction_2d(rect, 2.7, 2400.0, 900.0, T0, 0.0, 100.0)
    assert np.all(np.isfinite(res.T_field_C))

    # analytical_square_slab still callable
    T_an = analytical_square_slab(rect.x, rect.y, 100.0, 2.0, 1.0,
                                   2.7 / (2400 * 900), 20.0, 0.0)
    assert T_an.shape == (11, 11)


def test_hydration_rate_positive_early():
    """Hydration progresses, temperature rises, alpha monotone."""
    grid = build_grid_rectangular(1.0, 1.0, 7, 7)
    T0   = np.full((grid.ny, grid.nx), 20.0)
    res  = solve_hydration_2d(
        grid, MIX01, T0, 20.0,
        duration_s=48 * 3600,
        output_interval_s=3600.0,
        boundary_mode="adiabatic",
    )
    alpha_final = res.alpha_field[-1]
    assert alpha_final.max() > 0.4,  "hydration should exceed 0.4 by 48 hr"
    assert alpha_final.max() < MIX01.alpha_u + 1e-6, "alpha must not exceed alpha_u"
    # Monotone in time (allow tiny float rounding)
    dalpha = np.diff(res.alpha_field, axis=0)
    assert np.all(dalpha >= -1e-10), "alpha must be non-decreasing"
    assert res.peak_T_C > 40.0, "temperature should rise above 40 °C"


def test_hydration_adiabatic_no_heat_loss():
    """Energy added by hydration equals change in stored thermal energy (5% tol)."""
    grid  = build_grid_rectangular(1.0, 1.0, 7, 7)
    T0_C  = 20.0
    T0    = np.full((grid.ny, grid.nx), T0_C)
    res   = solve_hydration_2d(
        grid, MIX01, T0, T0_C,
        duration_s=168 * 3600,
        output_interval_s=6 * 3600,
        boundary_mode="adiabatic",
    )

    rho_c = MIX01.concrete_density_lb_ft3 * LB_FT3_TO_KG_M3
    Cc    = MIX01.total_cementitious_lb_yd3 * LB_YD3_TO_KG_M3
    Ca    = MIX01.aggregate_Cp_BTU_lb_F * BTU_LB_F_TO_J_KG_K
    k_uc  = MIX01.thermal_conductivity_BTU_hr_ft_F * BTU_HR_FT_F_TO_W_M_K
    Wc    = MIX01.cement_type_I_II_lb_yd3 * LB_YD3_TO_KG_M3
    Ww    = MIX01.water_lb_yd3 * LB_YD3_TO_KG_M3
    Wa    = (MIX01.coarse_agg_lb_yd3 + MIX01.fine_agg_lb_yd3) * LB_YD3_TO_KG_M3

    # Cell areas for finite-volume energy sum (uniform dx/dy here)
    dx = grid.dx
    dy = np.diff(grid.y)
    # Cell height = midpoint rule: for uniform grid, all = dy
    # For a uniform grid: cell_area = dx * dy for interior; edge cells half-width
    # Since it's vertex-centred, approximate with trapezoid weights:
    wy = np.concatenate([[dy[0] / 2],
                         (dy[:-1] + dy[1:]) / 2,
                         [dy[-1] / 2]])
    wx = np.full(grid.nx, dx)
    wx[0]  = dx / 2
    wx[-1] = dx / 2
    cell_area = np.outer(wy, wx)  # (ny, nx)

    # Initial and final Cp (average across field for energy calc)
    alpha_final = res.alpha_field[-1]
    T_final = res.T_field_C[-1]
    Cp_final = specific_heat_variable(Wc, Wa, Ww, Ca, alpha_final, T_final, rho_c)
    Cp_init  = specific_heat_variable(Wc, Wa, Ww, Ca,
                                       np.zeros_like(alpha_final),
                                       np.full_like(T_final, T0_C), rho_c)

    # Use mean Cp for the energy change estimate
    Cp_mean = (Cp_init + Cp_final) / 2.0

    dT = T_final - T0_C
    dE_stored = float(np.sum(rho_c * Cp_mean * dT * cell_area))

    # Expected heat release: Hu * Cc (kg/m³) * alpha_final * total_volume
    total_vol = float(np.sum(cell_area))
    E_hyd = MIX01.Hu_J_kg * Cc * float(alpha_final.mean()) * total_vol
    assert E_hyd > 0
    rel_err = abs(dE_stored - E_hyd) / E_hyd
    assert rel_err < 0.05, f"Energy balance error {rel_err:.3f} > 5%"


def test_centerline_matches_1d_v2():
    """2D centerline matches an in-test 1D reference with same physics."""
    grid = build_grid_half_mat(40.0, 8.0)
    T0_C = 15.56   # 60°F
    T0   = np.full((grid.ny, grid.nx), T0_C)

    res = solve_hydration_2d(
        grid, MIX01, T0, T0_C,
        duration_s=168 * 3600,
        output_interval_s=3600.0,
        boundary_mode="adiabatic",
    )

    # Build 1D reference from the centerline column's material profile
    y_1d   = grid.y
    mat_1d = grid.material_id[:, -1]
    T0_1d  = np.full(grid.ny, T0_C)

    T1d = _solve_1d_reference(
        y_1d, mat_1d, MIX01, T0_1d,
        duration_s=168 * 3600,
        output_interval_s=3600.0,
    )

    # Compare: result.centerline_T_C has shape (n_out, ny)
    err = res.centerline_T_C - T1d
    max_abs = float(np.abs(err).max())
    rms     = float(np.sqrt(np.mean(err**2)))
    assert max_abs < 2.0, f"max |err| = {max_abs:.3f} °C > 2.0 °C"
    assert rms     < 0.5, f"RMS err = {rms:.3f} °C > 0.5 °C"


def test_variable_k_stencil_reduces_to_constant():
    """When alpha_u=0 (no hydration, Q=0), hydration solver matches analytical solution.

    With alpha_u=0: hydration_alpha_vec → 0, Q=0,
    k_variable(k_uc, alpha=0) = 1.33*k_uc (constant across domain).
    The result must match the analytical Fourier-series solution for 2D diffusion
    with those constant properties — proving the stencil uses the correct uniform k.
    """
    mix_zero = MixStub(alpha_u=0.0)
    k_uc   = MIX01.thermal_conductivity_BTU_hr_ft_F * BTU_HR_FT_F_TO_W_M_K
    k_eff  = 1.33 * k_uc    # k_variable(k_uc, alpha=0)
    rho_c  = MIX01.concrete_density_lb_ft3 * LB_FT3_TO_KG_M3
    Ca     = MIX01.aggregate_Cp_BTU_lb_F * BTU_LB_F_TO_J_KG_K
    Wc     = MIX01.cement_type_I_II_lb_yd3 * LB_YD3_TO_KG_M3
    Ww     = MIX01.water_lb_yd3 * LB_YD3_TO_KG_M3
    Wa     = (MIX01.coarse_agg_lb_yd3 + MIX01.fine_agg_lb_yd3) * LB_YD3_TO_KG_M3
    Cp_init = float(specific_heat_variable(Wc, Wa, Ww, Ca, 0.0, 20.0, rho_c))
    alpha_diff = k_eff / (rho_c * Cp_init)  # m²/s

    Lx, Ly = 2.0, 1.0
    nx, ny = 21, 11
    grid = build_grid_rectangular(Lx, Ly, nx, ny)
    T0   = np.full((ny, nx), 100.0)
    T_bc = 0.0

    # Run to mid-diffusion time (same benchmark used in M1 analytical test)
    t_eval = 0.5 * min(Lx, Ly)**2 / alpha_diff
    res = solve_hydration_2d(
        grid, mix_zero, T0, T_bc,
        duration_s=t_eval,
        output_interval_s=t_eval,
        boundary_mode="dirichlet",
    )

    T_analytical = analytical_square_slab(
        grid.x, grid.y, t_eval, Lx, Ly, alpha_diff, 100.0, T_bc, n_terms=21
    )
    # Interior only (boundary is Dirichlet so edges trivially match)
    err = res.T_field_C[-1][1:-1, 1:-1] - T_analytical[1:-1, 1:-1]
    max_err = float(np.abs(err).max())
    rms_err = float(np.sqrt(np.mean(err**2)))
    assert max_err < 1.0, f"max|err| vs analytical = {max_err:.3f} °C > 1.0 °C"
    assert rms_err < 0.3, f"RMS err vs analytical = {rms_err:.3f} °C > 0.3 °C"


def test_soil_does_not_hydrate():
    """Soil cells never develop any hydration; concrete does."""
    grid = build_grid_half_mat(40.0, 8.0)
    T0   = np.full((grid.ny, grid.nx), 20.0)
    res  = solve_hydration_2d(
        grid, MIX01, T0, 20.0,
        duration_s=24 * 3600,
        output_interval_s=6 * 3600,
        boundary_mode="adiabatic",
    )
    assert np.all(res.alpha_field[:, grid.is_soil] == 0.0), \
        "soil alpha must be zero at all times"
    assert res.alpha_field[-1][grid.is_concrete].min() > 0.0, \
        "concrete alpha must be positive after 24 hr"


def test_blanket_does_not_hydrate():
    """Blanket cells never develop any hydration."""
    grid = build_grid_half_mat(40.0, 8.0)
    T0   = np.full((grid.ny, grid.nx), 20.0)
    res  = solve_hydration_2d(
        grid, MIX01, T0, 20.0,
        duration_s=24 * 3600,
        output_interval_s=6 * 3600,
        boundary_mode="adiabatic",
    )
    assert np.all(res.alpha_field[:, grid.is_blanket] == 0.0), \
        "blanket alpha must be zero at all times"


def test_output_sampling_exact_times():
    """Output samples land on exact multiples of output_interval_s."""
    grid = build_grid_rectangular(1.0, 1.0, 5, 5)
    T0   = np.full((grid.ny, grid.nx), 20.0)
    res  = solve_hydration_2d(
        grid, MIX01, T0, 20.0,
        duration_s=1500.0,
        output_interval_s=300.0,
        boundary_mode="adiabatic",
    )
    assert len(res.t_s) == 6, f"expected 6 samples, got {len(res.t_s)}"
    expected = np.array([0.0, 300.0, 600.0, 900.0, 1200.0, 1500.0])
    assert np.allclose(res.t_s, expected, atol=1e-6), \
        f"sample times not exact: {res.t_s}"


def test_cfl_with_variable_properties():
    """Solver remains stable (no NaN/Inf, peak T < 200°C) with variable properties."""
    grid = build_grid_half_mat(40.0, 8.0)
    T0   = np.full((grid.ny, grid.nx), 20.0)
    res  = solve_hydration_2d(
        grid, MIX01, T0, 20.0,
        duration_s=72 * 3600,
        output_interval_s=3600.0,
        boundary_mode="adiabatic",
    )
    assert np.all(np.isfinite(res.T_field_C)), "NaN or Inf in temperature field"
    assert res.peak_T_C < 200.0, f"peak T = {res.peak_T_C:.1f} °C (sanity exceeded)"
