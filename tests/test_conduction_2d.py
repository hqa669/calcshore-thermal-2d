"""Unit tests for 2D conduction solver and analytical validation (M1)."""
from __future__ import annotations

import numpy as np
import pytest

from thermal_engine_2d import (
    analytical_square_slab,
    build_grid_half_mat,
    build_grid_rectangular,
    solve_conduction_2d,
)

# ---------------------------------------------------------------------------
# Constants matching M1 spec (typical concrete)
# ---------------------------------------------------------------------------
K_W_M_K = 2.7
RHO_KG_M3 = 2400.0
CP_J_KG_K = 900.0
ALPHA = K_W_M_K / (RHO_KG_M3 * CP_J_KG_K)   # ≈ 1.25e-6 m²/s


# ---------------------------------------------------------------------------
# 1. build_grid_rectangular basic checks
# ---------------------------------------------------------------------------

def test_build_grid_rectangular_basic():
    Lx, Ly = 2.0, 1.0
    nx, ny = 21, 11
    g = build_grid_rectangular(Lx, Ly, nx, ny)

    assert g.nx == nx and g.ny == ny
    assert g.x.shape == (nx,) and g.y.shape == (ny,)
    assert pytest.approx(g.x[0], abs=1e-12) == 0.0
    assert pytest.approx(g.x[-1], abs=1e-12) == Lx
    assert pytest.approx(g.y[0], abs=1e-12) == 0.0
    assert pytest.approx(g.y[-1], abs=1e-12) == Ly

    dx_expected = Lx / (nx - 1)
    assert np.allclose(np.diff(g.x), dx_expected, atol=1e-12)
    assert np.allclose(np.diff(g.y), Ly / (ny - 1), atol=1e-12)
    assert pytest.approx(g.dx, abs=1e-12) == dx_expected

    assert g.material_id.shape == (ny, nx)
    assert np.all(g.material_id == 1)
    assert np.all(g.is_concrete)
    assert not np.any(g.is_blanket)
    assert not np.any(g.is_soil)


# ---------------------------------------------------------------------------
# 2. Steady state: uniform T_init == T_bc stays flat
# ---------------------------------------------------------------------------

def test_conduction_steady_state_matches_bc():
    g = build_grid_rectangular(1.0, 1.0, 11, 11)
    T0 = np.full((g.ny, g.nx), 20.0)
    res = solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                               T0, 20.0, 1000.0, output_interval_s=500.0)
    assert np.allclose(res.T_field_C, 20.0, atol=1e-6)


# ---------------------------------------------------------------------------
# 3. Cold slab heats toward BC
# ---------------------------------------------------------------------------

def test_conduction_cold_slab_heats_toward_bc():
    Lx, Ly = 1.0, 1.0
    nx, ny = 21, 21
    g = build_grid_rectangular(Lx, Ly, nx, ny)
    T0 = np.zeros((ny, nx))
    T_bc = 100.0
    # short time: diffusion length ~2√(αt) << 0.5*Ly so centre stays cold
    duration = 0.005 * Ly**2 / ALPHA

    res = solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                               T0, T_bc, duration,
                               output_interval_s=duration / 2)
    T_final = res.T_field_C[-1]

    # Corners are on the Dirichlet boundary — exactly 100 °C
    assert pytest.approx(T_final[0, 0], abs=1e-9) == T_bc

    # Interior centre still cold
    jc, ic = ny // 2, nx // 2
    assert T_final[jc, ic] < 10.0, "centre should still be well below BC"

    # Net heat influx: mean > 0
    assert T_final.mean() > 0.0

    # Physically bounded by BCs
    assert T_final.min() >= -1e-10
    assert T_final.max() <= T_bc + 1e-10


# ---------------------------------------------------------------------------
# 4. Core acceptance gate: match analytical solution at mid-diffusion time
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def mid_diffusion_result():
    Lx, Ly = 2.0, 1.0
    nx, ny = 41, 21
    g = build_grid_rectangular(Lx, Ly, nx, ny)
    T0 = np.full((ny, nx), 100.0)
    L_min = min(Lx, Ly)
    duration = 0.5 * L_min**2 / ALPHA
    res = solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                               T0, 0.0, duration, output_interval_s=300.0)
    return g, Lx, Ly, res


def test_conduction_matches_analytical_mid_diffusion(mid_diffusion_result):
    g, Lx, Ly, res = mid_diffusion_result
    t_final = res.t_s[-1]
    T_num = res.T_field_C[-1]

    T_a = analytical_square_slab(g.x, g.y, t_final, Lx, Ly, ALPHA,
                                  100.0, 0.0, n_terms=11)

    # Compare interior nodes only (boundary is exact by construction)
    err = T_num[1:-1, 1:-1] - T_a[1:-1, 1:-1]
    max_err = float(np.abs(err).max())
    rms_err = float(np.sqrt((err**2).mean()))

    assert max_err <= 0.5, f"Interior max|err| = {max_err:.4f} °C (limit 0.5)"
    assert rms_err <= 0.2, f"Interior RMS err  = {rms_err:.4f} °C (limit 0.2)"


# ---------------------------------------------------------------------------
# 5. Symmetry preserved on a square grid
# ---------------------------------------------------------------------------

def test_conduction_symmetry_preserved():
    n = 21
    g = build_grid_rectangular(1.0, 1.0, n, n)
    T0 = np.full((n, n), 50.0)
    duration = 0.2 * 1.0**2 / ALPHA
    res = solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                               T0, 0.0, duration, output_interval_s=duration)
    T_final = res.T_field_C[-1]
    # On a square domain with uniform IC and BC, T[j,i] == T[i,j]
    assert np.allclose(T_final, T_final.T, atol=1e-8), \
        "Field is not symmetric under transpose (symmetry broken)"


# ---------------------------------------------------------------------------
# 6. CFL stability: safe and unstable cases
# ---------------------------------------------------------------------------

def test_conduction_CFL_stability():
    g = build_grid_rectangular(1.0, 1.0, 11, 11)
    T0 = np.full((g.ny, g.nx), 50.0)
    duration = 3000.0

    # Near stability limit — should complete without NaN/inf
    res = solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                               T0, 0.0, duration,
                               output_interval_s=1000.0, cfl_safety=0.95)
    assert not np.any(np.isnan(res.T_field_C)), "NaN detected at cfl_safety=0.95"
    assert not np.any(np.isinf(res.T_field_C)), "Inf detected at cfl_safety=0.95"
    assert res.T_field_C.min() >= -1e-6
    assert res.T_field_C.max() <= 50.0 + 1e-6

    # Violating the stability condition must raise ValueError
    with pytest.raises(ValueError, match="cfl_safety"):
        solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                             T0, 0.0, duration,
                             output_interval_s=1000.0, cfl_safety=2.0)


# ---------------------------------------------------------------------------
# 7. Output sampling cadence
# ---------------------------------------------------------------------------

def test_output_sampling_cadence():
    g = build_grid_rectangular(1.0, 1.0, 11, 11)
    T0 = np.full((g.ny, g.nx), 50.0)
    res = solve_conduction_2d(g, K_W_M_K, RHO_KG_M3, CP_J_KG_K,
                               T0, 0.0, 900.0, output_interval_s=300.0)

    assert res.n_output_samples == 4
    assert len(res.t_s) == 4
    assert res.T_field_C.shape[0] == 4
    assert res.t_s[0] == pytest.approx(0.0, abs=1e-9)
    assert res.t_s[-1] == pytest.approx(900.0, abs=res.dt_inner_s)


# ---------------------------------------------------------------------------
# 8. Analytical helper sanity checks
# ---------------------------------------------------------------------------

def test_analytical_solver_sanity():
    Lx, Ly = 1.0, 1.0
    nx, ny = 21, 21
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, Ly, ny)
    T_init, T_bc = 100.0, 0.0

    # At t ≈ 0 (very small), interior centre should be near T_init
    # Use n_terms=21 to reduce Gibbs overshoot at boundary
    T_early = analytical_square_slab(x, y, 1.0, Lx, Ly, ALPHA,
                                      T_init, T_bc, n_terms=21)
    cx, cy = nx // 2, ny // 2
    assert abs(T_early[cy, cx] - T_init) / T_init < 0.05, \
        "At t≈0, centre should be within 5% of T_init"

    # At very large t, should decay essentially to T_bc everywhere
    t_large = 100.0 * Ly**2 / ALPHA
    T_late = analytical_square_slab(x, y, t_large, Lx, Ly, ALPHA,
                                     T_init, T_bc, n_terms=11)
    assert np.abs(T_late[1:-1, 1:-1]).max() < 1e-3, \
        "At large t, interior should decay to T_bc within 1e-3"


# ---------------------------------------------------------------------------
# 9. M0 regression: build_grid_half_mat still works correctly
# ---------------------------------------------------------------------------

def test_regression_m0_grid_builder():
    # Pin to native 1× grid; these counts document native resolution.
    g = build_grid_half_mat(40.0, 8.0, grid_refinement=1)
    assert g.nx == 33
    assert g.ny == 14   # n_blanket(1) + n_concrete_y(13); no soil rows
    assert int(g.is_concrete.sum()) == 273
    assert g.n_soil_y == 0
    assert int(g.is_soil.sum()) == 0

    # model_soil=True preserves the original 15-row soil buffer.
    g2 = build_grid_half_mat(40.0, 8.0, grid_refinement=1, model_soil=True)
    assert g2.nx == 33
    assert g2.ny == 29   # n_blanket(1) + n_concrete_y(13) + n_soil_y(15)
    assert int(g2.is_concrete.sum()) == 273
    assert g2.n_soil_y == 15
