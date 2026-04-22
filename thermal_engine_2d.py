"""
CalcShore Thermal Engine v3 — 2D Finite-Difference Solver
==========================================================
Milestone M0: Grid builder and domain composition.
Milestone M1: Explicit 2D conduction solver, constant properties,
              Dirichlet BC, analytical Fourier-series validation helper.

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
class ConductionResult:
    """Output from solve_conduction_2d.

    All temperatures in °C, times in seconds.
    """

    t_s: np.ndarray           # shape (n_out,) — sampled time points
    T_field_C: np.ndarray     # shape (n_out, ny, nx)
    dt_inner_s: float         # actual inner time step used
    n_inner_steps: int        # total inner time steps executed
    n_output_samples: int     # len(t_s)
    peak_T_C: float           # max over all time and space
    min_T_C: float            # min over all time and space


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


def build_grid_rectangular(
    Lx_m: float,
    Ly_m: float,
    nx: int,
    ny: int,
) -> Grid2D:
    """Build a uniform rectangular grid with all nodes labelled as concrete.

    Used for M1 analytical validation and future pavement scenarios.
    x in [0, Lx_m], y in [0, Ly_m], uniform spacing on both axes.

    Parameters
    ----------
    Lx_m : float
        Domain width in metres (x-extent).
    Ly_m : float
        Domain depth in metres (y-extent).
    nx : int
        Number of nodes in x (>= 3).
    ny : int
        Number of nodes in y (>= 3).

    Returns
    -------
    Grid2D
    """
    if Lx_m <= 0.0 or Ly_m <= 0.0:
        raise ValueError(f"Lx_m and Ly_m must be positive, got {Lx_m}, {Ly_m}")
    if nx < 3 or ny < 3:
        raise ValueError(f"nx and ny must be >= 3, got {nx}, {ny}")

    x = np.linspace(0.0, Lx_m, nx)
    y = np.linspace(0.0, Ly_m, ny)
    dx = Lx_m / (nx - 1)

    material_id = np.ones((ny, nx), dtype=np.int8)
    is_concrete = np.ones((ny, nx), dtype=bool)
    is_blanket = np.zeros((ny, nx), dtype=bool)
    is_soil = np.zeros((ny, nx), dtype=bool)

    return Grid2D(
        x=x,
        y=y,
        material_id=material_id,
        is_blanket=is_blanket,
        is_concrete=is_concrete,
        is_soil=is_soil,
        ix_concrete_start=0,
        ix_concrete_end=nx - 1,
        iy_blanket_end=-1,
        iy_concrete_start=0,
        iy_concrete_end=ny - 1,
        iy_centerline=ny // 2,
        nx=nx,
        ny=ny,
        dx=dx,
        n_blanket=0,
        n_concrete_x=nx,
        n_concrete_y=ny,
        n_soil_x_ext=0,
        n_soil_y=0,
    )


def solve_conduction_2d(
    grid: Grid2D,
    k_W_m_K: float,
    rho_kg_m3: float,
    Cp_J_kg_K: float,
    T_initial_C: np.ndarray,
    T_boundary_C: float,
    duration_s: float,
    output_interval_s: float = 300.0,
    cfl_safety: float = 0.4,
) -> ConductionResult:
    """Explicit 5-point 2D conduction solver with constant properties.

    Solves  ρ Cp ∂T/∂t = k (∂²T/∂x² + ∂²T/∂y²)  on the Grid2D domain
    with fixed Dirichlet BCs on all four sides.

    Stencil (D4, non-uniform dy):
        T_new[j,i] = T[j,i] + dt·α·(
            (T[j,i+1] - 2T[j,i] + T[j,i-1]) / dx²
          + (2/(dy+ + dy-)) · ((T[j+1,i]-T[j,i])/dy+ - (T[j,i]-T[j-1,i])/dy-)
        )

    Parameters
    ----------
    grid : Grid2D
    k_W_m_K : float
        Thermal conductivity W/(m·K).
    rho_kg_m3 : float
        Density kg/m³.
    Cp_J_kg_K : float
        Specific heat J/(kg·K).
    T_initial_C : np.ndarray
        Shape (ny, nx), initial temperature field in °C.
    T_boundary_C : float
        Dirichlet BC value applied to all four edges every step.
    duration_s : float
        Total simulation duration in seconds.
    output_interval_s : float
        Time between output samples in seconds.
    cfl_safety : float
        Multiplier on CFL-limited dt; must be <= 1.0.

    Returns
    -------
    ConductionResult
    """
    if cfl_safety > 1.0:
        raise ValueError(
            f"cfl_safety={cfl_safety} > 1.0 violates explicit stability; "
            "use cfl_safety <= 1.0"
        )

    alpha = k_W_m_K / (rho_kg_m3 * Cp_J_kg_K)
    dx = grid.dx
    dy_all = np.diff(grid.y)          # (ny-1,)
    dy_min = float(dy_all.min())

    dt_cfl = 0.5 / (alpha * (1.0 / dx**2 + 1.0 / dy_min**2))
    dt = cfl_safety * dt_cfl

    # per-interior-row spacings broadcast over x automatically
    dy_plus  = dy_all[1:].reshape(-1, 1)    # (ny-2, 1): y[j+1]-y[j] for j=1..ny-2
    dy_minus = dy_all[:-1].reshape(-1, 1)   # (ny-2, 1): y[j]-y[j-1]
    inv_dxsq = 1.0 / dx**2
    y_coef   = 2.0 / (dy_plus + dy_minus)   # (ny-2, 1)

    # For M1, k is uniform so all four face-k values equal k_W_m_K.
    # These scalars make the stencil M2-ready: swap to arrays for harmonic mean.
    k_face = k_W_m_K  # k_x_plus = k_x_minus = k_y_plus = k_y_minus

    # Determine output sample times (include t=0)
    n_samples = int(round(duration_s / output_interval_s)) + 1
    t_samples = np.linspace(0.0, duration_s, n_samples)

    T_out = np.empty((n_samples, grid.ny, grid.nx), dtype=np.float64)
    t_out = np.empty(n_samples, dtype=np.float64)

    T = T_initial_C.astype(np.float64, copy=True)
    T_new = np.empty_like(T)

    # Apply BC to initial field so boundary is consistent from the start
    T[0, :]  = T_boundary_C
    T[-1, :] = T_boundary_C
    T[:, 0]  = T_boundary_C
    T[:, -1] = T_boundary_C

    # Store t=0 sample
    T_out[0] = T
    t_out[0] = 0.0
    next_sample = 1

    n_steps = 0
    sim_time = 0.0

    while sim_time < duration_s - 1e-12:
        # Don't overshoot duration
        dt_step = min(dt, duration_s - sim_time)

        Tc  = T[1:-1, 1:-1]
        Tip = T[1:-1, 2:  ]
        Tim = T[1:-1,  :-2]
        Tjp = T[2:  , 1:-1]
        Tjm = T[ :-2, 1:-1]

        lap_x = (Tip - 2.0 * Tc + Tim) * inv_dxsq * k_face
        lap_y = k_face * y_coef * ((Tjp - Tc) / dy_plus - (Tc - Tjm) / dy_minus)

        T_new[1:-1, 1:-1] = Tc + (dt_step / (rho_kg_m3 * Cp_J_kg_K)) * (lap_x + lap_y)

        T_new[0, :]  = T_boundary_C
        T_new[-1, :] = T_boundary_C
        T_new[:, 0]  = T_boundary_C
        T_new[:, -1] = T_boundary_C

        T, T_new = T_new, T
        sim_time += dt_step
        n_steps += 1

        # Capture output samples whose scheduled time has been reached
        while next_sample < n_samples and sim_time >= t_samples[next_sample] - dt * 0.5:
            T_out[next_sample] = T
            t_out[next_sample] = sim_time
            next_sample += 1

    # Ensure all samples filled (handles edge case where duration is exact)
    while next_sample < n_samples:
        T_out[next_sample] = T
        t_out[next_sample] = sim_time
        next_sample += 1

    return ConductionResult(
        t_s=t_out,
        T_field_C=T_out,
        dt_inner_s=dt,
        n_inner_steps=n_steps,
        n_output_samples=n_samples,
        peak_T_C=float(T_out.max()),
        min_T_C=float(T_out.min()),
    )


def analytical_square_slab(
    x: np.ndarray,
    y: np.ndarray,
    t: float,
    Lx: float,
    Ly: float,
    alpha: float,
    T_init: float,
    T_bc: float,
    n_terms: int = 11,
) -> np.ndarray:
    """Analytical solution for 2D diffusion in a rectangle with Dirichlet BCs.

    Solves ∂T/∂t = α(∂²T/∂x² + ∂²T/∂y²) on [0,Lx]×[0,Ly] with
    T=T_bc on all four sides and T(x,y,0)=T_init (uniform).

    θ(x,y,t) = θ₀·(16/π²)·Σ_m Σ_n sin(mπx/Lx)·sin(nπy/Ly)/(m·n)
                           ·exp(-α·π²·t·(m²/Lx²+n²/Ly²))
    where m, n = 1, 3, 5, … (odd integers only), θ = T - T_bc.

    Parameters
    ----------
    x : np.ndarray
        Shape (nx,), node x-coordinates (domain should start at 0).
    y : np.ndarray
        Shape (ny,), node y-coordinates (domain should start at 0).
    t : float
        Evaluation time in seconds.
    Lx, Ly : float
        Domain extents in metres.
    alpha : float
        Thermal diffusivity m²/s.
    T_init : float
        Uniform initial temperature °C.
    T_bc : float
        Dirichlet boundary temperature °C.
    n_terms : int
        Number of odd Fourier modes (1, 3, …, 2*n_terms-1).

    Returns
    -------
    np.ndarray
        Shape (ny, nx), temperature field in °C.
    """
    odd = 2 * np.arange(n_terms) + 1   # [1, 3, 5, ...]
    theta0 = T_init - T_bc
    prefactor = theta0 * (16.0 / np.pi**2)

    X = x[np.newaxis, :]    # (1, nx)
    Y = y[:, np.newaxis]    # (ny, 1)

    theta = np.zeros((y.size, x.size), dtype=np.float64)
    for m in odd:
        sx = np.sin(m * np.pi * X / Lx)         # (1, nx)
        for n in odd:
            sy = np.sin(n * np.pi * Y / Ly)     # (ny, 1)
            decay = np.exp(-alpha * np.pi**2 * t * (m**2 / Lx**2 + n**2 / Ly**2))
            theta += sx * sy * (decay / (m * n))

    return T_bc + prefactor * theta
