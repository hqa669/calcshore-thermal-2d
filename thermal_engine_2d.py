"""
CalcShore Thermal Engine v3 — 2D Finite-Difference Solver
==========================================================
Milestone M0: Grid builder and domain composition.
No physics, no time stepping, no BC application.

Coordinate convention (half-mat, vertex-centered):
  x = 0      : concrete/form edge (soil-extension at x < 0)
  x = W/2    : centerline (symmetry BC, applied in M4)
  y = 0      : top of blanket (top surface exposed to ambient)
  y > 0      : downward into the structure

Material IDs: 0 = blanket, 1 = concrete, 2 = soil
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np

FT_TO_M = 0.3048


@dataclass
class Grid2D:
    """Half-mat 2D grid for the thermal engine.

    Vertex-centered: nodes sit at material interfaces.
    Node [j, i] is at (x[i], y[j]).  Shape convention: (ny, nx).
    """

    # Coordinate arrays (meters)
    x: np.ndarray  # shape (nx,), monotonically increasing
    y: np.ndarray  # shape (ny,), monotonically increasing (downward)

    # Material ID per cell. Shape (ny, nx). Values: 0=blanket, 1=concrete, 2=soil
    material_id: np.ndarray

    # Boolean convenience masks, all shape (ny, nx)
    is_blanket: np.ndarray
    is_concrete: np.ndarray
    is_soil: np.ndarray

    # Index metadata for later milestones
    ix_concrete_start: int  # first x-index where material is concrete/blanket (= n_soil_x_ext)
    ix_concrete_end: int    # last x-index of concrete/blanket (= nx - 1, inclusive)
    iy_blanket_end: int     # y-index of blanket/concrete interface (= n_blanket; node belongs to concrete)
    iy_concrete_start: int  # first y-index of concrete (= n_blanket)
    iy_concrete_end: int    # last y-index of concrete (= n_blanket + n_concrete_y - 1, inclusive)
    iy_centerline: int      # last x-index (= nx - 1), where symmetry BC is applied in M4

    # Raw build params (for debugging / regeneration)
    nx: int
    ny: int
    dx: float   # uniform across x
    n_blanket: int
    n_concrete_x: int
    n_concrete_y: int
    n_soil_x_ext: int
    n_soil_y: int

    def concrete_slice(self) -> tuple[slice, slice]:
        """Returns (y_slice, x_slice) to extract concrete subdomain from a full field."""
        return (
            slice(self.iy_concrete_start, self.iy_concrete_end + 1),
            slice(self.ix_concrete_start, self.ix_concrete_end + 1),
        )


def build_grid_half_mat(
    width_ft: float,
    depth_ft: float,
    n_concrete_x: int = 21,
    n_concrete_y: int = 13,
    blanket_thickness_m: float = 0.02,
    n_blanket: int = 1,
    soil_depth_below_m: float = 3.0,
    n_soil_y: int = 15,
    n_soil_x_ext: int = 12,
    soil_ext_lateral_m: float | None = None,
) -> Grid2D:
    """Build a half-mat 2D grid for the thermal solver.

    The grid is vertex-centered.  Material boundaries align with node rows /
    columns so interface nodes are unambiguously on the boundary.

    Parameters
    ----------
    width_ft : float
        Full mat width in feet.  Half-width (centerline) is modelled.
    depth_ft : float
        Concrete depth in feet.
    n_concrete_x : int
        Number of nodes across the concrete half-width (vertex-centered,
        so n_concrete_x - 1 intervals span exactly W/2).
    n_concrete_y : int
        Number of nodes through the concrete depth.
    blanket_thickness_m : float
        Physical blanket thickness in metres.  Blanket gets n_blanket nodes.
    n_blanket : int
        Number of y-nodes allocated to the blanket layer.
    soil_depth_below_m : float
        Depth of soil domain below the concrete bottom.
    n_soil_y : int
        Number of y-nodes in the soil below concrete.
    n_soil_x_ext : int
        Number of soil-extension cells to the left of x=0.  Cell edges align
        with the uniform dx so the full x-axis has a single dx.
    soil_ext_lateral_m : float | None
        Optional target lateral extent.  If it disagrees with
        n_soil_x_ext * dx by more than 1 cm, a UserWarning is emitted and the
        exact-alignment value is used instead (dx stays uniform).

    Returns
    -------
    Grid2D
    """
    half_width_m = width_ft * FT_TO_M / 2.0
    depth_m = depth_ft * FT_TO_M

    # Uniform dx derived from concrete node spacing (vertex-centered: n-1 intervals)
    dx = half_width_m / (n_concrete_x - 1)

    # Soil-extension lateral extent aligned to exact integer multiples of dx
    L_ext_actual = n_soil_x_ext * dx

    if soil_ext_lateral_m is not None and abs(soil_ext_lateral_m - L_ext_actual) > 0.01:
        warnings.warn(
            f"soil_ext_lateral_m={soil_ext_lateral_m:.4f} m disagrees with "
            f"n_soil_x_ext * dx = {L_ext_actual:.4f} m by more than 1 cm. "
            f"Using exact-alignment value {L_ext_actual:.4f} m to keep dx uniform.",
            UserWarning,
            stacklevel=2,
        )

    nx = n_soil_x_ext + n_concrete_x
    ny = n_blanket + n_concrete_y + n_soil_y

    # --- x coordinates (uniform, monotonically increasing) ---
    x = np.linspace(-L_ext_actual, half_width_m, nx)

    # --- y coordinates (non-uniform: blanket << concrete/soil dy) ---
    # Blanket: n_blanket nodes from y=0 to y=blanket_thickness_m.
    y_blanket = np.linspace(0.0, blanket_thickness_m, n_blanket + 1)[:n_blanket]

    # Concrete: n_concrete_y nodes from blanket_thickness_m to blanket_thickness_m + depth_m.
    y_concrete = np.linspace(blanket_thickness_m, blanket_thickness_m + depth_m, n_concrete_y)

    # Soil below concrete: n_soil_y nodes, uniform dy_soil, starting one dy_soil below concrete bottom.
    dy_soil = soil_depth_below_m / n_soil_y
    y_concrete_bottom = blanket_thickness_m + depth_m
    y_soil = y_concrete_bottom + dy_soil * np.arange(1, n_soil_y + 1)

    y = np.concatenate([y_blanket, y_concrete, y_soil])

    # --- Index metadata ---
    ix_concrete_start = n_soil_x_ext
    ix_concrete_end = nx - 1
    iy_blanket_end = n_blanket  # interface row; y[n_blanket] = blanket_thickness_m
    iy_concrete_start = n_blanket
    iy_concrete_end = n_blanket + n_concrete_y - 1
    iy_centerline = nx - 1  # last x-index is the symmetry centerline

    # --- Material ID array: default fill = 2 (soil) ---
    material_id = np.full((ny, nx), 2, dtype=np.int8)

    # Blanket: rows 0..n_blanket-1, cols ix_concrete_start..nx-1
    material_id[:n_blanket, ix_concrete_start:] = 0

    # Concrete: rows iy_concrete_start..iy_concrete_end, cols ix_concrete_start..nx-1
    material_id[iy_concrete_start : iy_concrete_end + 1, ix_concrete_start:] = 1

    # The soil-extension strip (cols 0..ix_concrete_start-1) keeps 2 for all rows,
    # including y=0 (bare ground beside mat exposed to ambient — design Risk 1).

    # --- Boolean masks ---
    is_blanket = material_id == 0
    is_concrete = material_id == 1
    is_soil = material_id == 2

    return Grid2D(
        x=x,
        y=y,
        material_id=material_id,
        is_blanket=is_blanket,
        is_concrete=is_concrete,
        is_soil=is_soil,
        ix_concrete_start=ix_concrete_start,
        ix_concrete_end=ix_concrete_end,
        iy_blanket_end=iy_blanket_end,
        iy_concrete_start=iy_concrete_start,
        iy_concrete_end=iy_concrete_end,
        iy_centerline=iy_centerline,
        nx=nx,
        ny=ny,
        dx=dx,
        n_blanket=n_blanket,
        n_concrete_x=n_concrete_x,
        n_concrete_y=n_concrete_y,
        n_soil_x_ext=n_soil_x_ext,
        n_soil_y=n_soil_y,
    )
