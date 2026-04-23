"""Unit tests for Grid2D and build_grid_half_mat (M0)."""

import warnings

import numpy as np
import pytest

from thermal_engine_2d import Grid2D, build_grid_half_mat, build_grid_rectangular


@pytest.fixture(scope="module")
def grid() -> Grid2D:
    return build_grid_half_mat(40.0, 8.0)


# ---------------------------------------------------------------------------
# 1. Shape and extents
# ---------------------------------------------------------------------------

def test_default_grid_shape_and_extents(grid):
    assert grid.nx == 33
    assert grid.ny == 29
    assert grid.x.shape == (33,)
    assert grid.y.shape == (29,)
    assert abs(grid.x[0] - (-3.6576)) < 1e-6
    assert abs(grid.x[-1] - 6.0960) < 1e-6
    assert abs(grid.y[0]) < 1e-6
    assert abs(grid.y[-1] - 5.4584) < 1e-6


# ---------------------------------------------------------------------------
# 2. Uniform dx
# ---------------------------------------------------------------------------

def test_uniform_dx(grid):
    diffs = np.diff(grid.x)
    assert np.allclose(diffs, grid.dx, atol=1e-9)
    assert abs(grid.dx - 0.3048) < 1e-9


# ---------------------------------------------------------------------------
# 3. Non-uniform dy
# ---------------------------------------------------------------------------

def test_nonuniform_dy(grid):
    assert grid.y[1] - grid.y[0] == pytest.approx(0.02, abs=1e-9)     # blanket
    assert grid.y[2] - grid.y[1] == pytest.approx(0.2032, abs=1e-4)   # first concrete dy
    assert np.diff(grid.y)[-1] == pytest.approx(0.2, abs=1e-9)         # last soil dy


# ---------------------------------------------------------------------------
# 4. Material ID counts
# ---------------------------------------------------------------------------

def test_material_id_counts(grid):
    # M4 geometry: soil-extension column above concrete-soil interface is now AIR (id=3)
    # blanket = 1 row × 21 cols = 21
    # concrete = 13 rows × 21 cols = 273
    # soil = 15 rows × 33 cols (below concrete) = 495
    # air = 14 rows × 12 cols (soil-ext column above concrete-soil interface) = 168
    assert int(grid.is_blanket.sum()) == 21
    assert int(grid.is_concrete.sum()) == 273
    assert int(grid.is_soil.sum()) == 495
    assert int(grid.is_air.sum()) == 168
    assert grid.is_blanket.sum() + grid.is_concrete.sum() + grid.is_soil.sum() + grid.is_air.sum() == 957
    assert set(np.unique(grid.material_id)) == {0, 1, 2, 3}
    # Masks are mutually exclusive and exhaustive
    assert not np.any(grid.is_blanket & grid.is_concrete)
    assert not np.any(grid.is_blanket & grid.is_soil)
    assert not np.any(grid.is_blanket & grid.is_air)
    assert not np.any(grid.is_concrete & grid.is_soil)
    assert not np.any(grid.is_concrete & grid.is_air)
    assert not np.any(grid.is_soil & grid.is_air)
    assert np.all(grid.is_blanket | grid.is_concrete | grid.is_soil | grid.is_air)


# ---------------------------------------------------------------------------
# 5. Material topology
# ---------------------------------------------------------------------------

def test_material_topology(grid):
    ix0 = grid.ix_concrete_start
    iy_ce = grid.iy_concrete_end

    # Blanket row over concrete is blanket
    assert np.all(grid.material_id[0, ix0:] == 0)

    # Blanket row over soil-extension is AIR (M4: not soil — concrete sides are air)
    assert np.all(grid.material_id[0, :ix0] == 3)

    # Concrete block
    assert np.all(
        grid.material_id[grid.iy_concrete_start : iy_ce + 1, ix0:] == 1
    )

    # Soil-extension column (x < 0) ABOVE concrete-soil interface is AIR (material_id=3)
    assert np.all(grid.material_id[: iy_ce + 1, :ix0] == 3)

    # Soil-extension column (x < 0) BELOW concrete-soil interface is soil
    assert np.all(grid.material_id[iy_ce + 1 :, :ix0] == 2)

    # Soil below concrete
    assert np.all(grid.material_id[iy_ce + 1 :, ix0:] == 2)


# ---------------------------------------------------------------------------
# 6. Centerline and edge indices
# ---------------------------------------------------------------------------

def test_centerline_and_edge_indices(grid):
    assert grid.iy_centerline == grid.nx - 1
    assert grid.x[grid.ix_concrete_start] == pytest.approx(0.0, abs=1e-6)
    assert grid.x[grid.ix_concrete_end] == pytest.approx(6.0960, abs=1e-6)
    assert grid.y[grid.iy_blanket_end] == pytest.approx(0.02, abs=1e-6)
    assert grid.y[grid.iy_concrete_end] == pytest.approx(2.4584, abs=1e-4)


# ---------------------------------------------------------------------------
# 7. concrete_slice extracts correct subdomain
# ---------------------------------------------------------------------------

def test_concrete_slice_extracts_correct_subdomain(grid):
    ny, nx = grid.ny, grid.nx
    field = np.arange(ny * nx).reshape(ny, nx)
    sub = field[grid.concrete_slice()]
    assert sub.shape == (13, 21)
    assert np.all(grid.material_id[grid.concrete_slice()] == 1)


# ---------------------------------------------------------------------------
# 8. CW temp.txt widths recoverable
# ---------------------------------------------------------------------------

def test_cw_temp_txt_widths_recoverable(grid):
    cw_widths = np.array([
        0.00, 0.30, 0.61, 0.91, 1.22, 1.52, 1.83, 2.13,
        2.44, 2.74, 3.05, 3.35, 3.66, 3.96, 4.27, 4.57,
        4.88, 5.18, 5.49, 5.79, 6.10,
    ])
    ix0 = grid.ix_concrete_start
    concrete_x = grid.x[ix0:]
    for w_cw in cw_widths:
        assert np.any(np.abs(concrete_x - w_cw) < 0.01), (
            f"CW width {w_cw} m not found in grid.x[{ix0}:]; "
            f"closest = {concrete_x[np.argmin(np.abs(concrete_x - w_cw))]:.4f}"
        )


# ---------------------------------------------------------------------------
# 9. Warning on mismatched soil_ext_lateral_m
# ---------------------------------------------------------------------------

def test_soil_ext_lateral_override_emits_warning_on_mismatch():
    with pytest.warns(UserWarning):
        g = build_grid_half_mat(40.0, 8.0, soil_ext_lateral_m=5.0, n_soil_x_ext=12)
    assert g.x[0] == pytest.approx(-3.6576, abs=1e-6)


# ---------------------------------------------------------------------------
# 10. No warning when soil_ext_lateral_m matches exactly
# ---------------------------------------------------------------------------

def test_soil_ext_lateral_override_accepts_matching_value():
    matching = 12 * (40.0 * 0.3048 / 2.0) / (21 - 1)  # = 12 * dx = 3.6576
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        g = build_grid_half_mat(40.0, 8.0, soil_ext_lateral_m=matching, n_soil_x_ext=12)
    assert g.nx == 33
    assert g.ny == 29


# ---------------------------------------------------------------------------
# 11. Smaller mat geometry
# ---------------------------------------------------------------------------

def test_smaller_mat_geometry():
    g = build_grid_half_mat(20.0, 4.0)
    assert g.nx == 33
    assert g.ny == 29
    diffs = np.diff(g.x)
    assert np.allclose(diffs, g.dx, atol=1e-9)


# ---------------------------------------------------------------------------
# 12. Air mask (M4)
# ---------------------------------------------------------------------------

def test_air_mask_exists_and_correct(grid):
    assert hasattr(grid, 'is_air'), "Grid2D must have is_air attribute (M4)"
    assert grid.is_air.shape == (grid.ny, grid.nx)
    assert grid.is_air.dtype == bool
    assert int(grid.is_air.sum()) == 168
    assert np.array_equal(grid.is_air, grid.material_id == 3)


def test_rectangular_grid_has_no_air():
    g = build_grid_rectangular(2.0, 1.0, 11, 9)
    assert hasattr(g, 'is_air')
    assert not g.is_air.any()
