# STAGE4a — Bottom-Soil Engine Audit

## S4a.1 Grid allocation

File: `thermal_engine_2d.py`, function `build_grid_half_mat` (L896–1043).

Bottom soil is allocated uniformly at **L976–981**:

```python
dy_soil = soil_depth_below_m / n_soil_y          # L977
y_concrete_bottom = blanket_thickness_m + depth_m # L978
y_soil = y_concrete_bottom + dy_soil * np.arange(1, n_soil_y + 1)  # L979
```

Defaults (L903–904): `soil_depth_below_m=3.0 m`, `n_soil_y=15` → `dy_soil=0.20 m`.

Grid dimensions for the 40ft wide × 70ft deep mat:
- `ny = n_blanket + n_concrete_y + n_soil_y = 1 + 13 + 15 = 29` rows
- Bottom soil occupies rows 14–28 (0-indexed), last row (j=28) is the Dirichlet boundary
- All cells in rows 14–28 have `material_id=2` (soil) by default fill (L992)
- The bottom soil spans **all x columns**, including the concrete-side columns
  (`ix >= ix_concrete_start`): these are also soil because the default fill is 2 and no concrete
  assignment extends below row `iy_concrete_end` = 13.

## S4a.1 Deep BC application

File: `thermal_engine_2d.py`, inner time-step loop.

### Far-left soil column (L2176–2181)

```python
# (d) Far-left column (i = 0, all soil rows): Dirichlet T_gw.
T_new[grid.is_soil[:, 0], 0] = _T_gw_C
```

Applied only to cells where `material_id==2` in column 0. Uses `is_soil` mask to
automatically include the submerged-case rows (soil extends alongside the concrete when
`is_submerged=True`).

### Bottom row (L2183–2184)

```python
# (e) Bottom row: Dirichlet T_gw
T_new[-1, :] = _T_gw_C
```

**Unconditional full-row overwrite.** Applies to every column in the last row — soil,
concrete-below, and the centerline (CL) column (ix=nx-1). No material mask. No
half-cell correction.

Both BCs share the same source `_T_gw_C` (L1483–1486):

```python
_T_gw_C = (
    T_ground_deep_C if T_ground_deep_C is not None
    else compute_T_gw_C(environment)
)
```

## S4a.1 Side vs. bottom symmetry after Stage 3

Stage 3 (commit f5842af) added a **half-cell stencil correction** for the side
(`is_submerged=True`) at L2085–2098:

```python
if grid.is_soil[grid.iy_concrete_start, _ics - 1]:
    # half-cell correction for side soil contact
    _left_extra = _kxm_side * (T[_js_side, _ics-1] - T[_js_side, _ics]) * inv_dxsq
    T_new[_js_side, _ics] += dt_step * (_lat_extra + _left_extra) / rho_cp[_js_side, _ics]
```

This correction accounts for the half-cell mass at the concrete-soil interface: the
form-face node occupies half a cell, so the stencil double-counts the lateral flux.

**The bottom row has no analogous correction.** The bottom-row Dirichlet
(`T_new[-1, :] = _T_gw_C`) is applied *after* all stencil updates and overwrites
whatever the stencil computed. This means:

1. The stencil computes a physically correct temperature for row -1 (deepest soil row)
2. The Dirichlet assignment overwrites it with `_T_gw_C`

For the **bottom-CL corner** (j=-1, i=nx-1), two conflicting BCs apply:
- The CL symmetry stencil (L2100–2118) computes a zero-flux result
- The bottom-row Dirichlet (L2184) then overwrites to `_T_gw_C`

The Dirichlet wins because it is applied last. The result: **the bottom-CL corner
cell is pinned to `T_gw` regardless of the actual temperature the zero-flux CL BC
would produce.**

## S4a.1 Geometry differences between side and bottom

| Feature | Side (form face, x-direction) | Bottom (y-direction) |
|---|---|---|
| Soil Dirichlet | Far-left column, `is_soil[:, 0]` masked | Full last row, no mask |
| Half-cell correction | Yes, at `ix_concrete_start` when `is_submerged` | None |
| CL corner | Not affected by side Dirichlet | **Overwritten** by bottom Dirichlet |
| BC type | Dirichlet via `is_soil` mask | Dirichlet unconditional full-row |
| `_T_gw_C` source | Same `compute_T_gw_C(environment)` | Same |

The side Dirichlet at column 0 uses an `is_soil` mask, so it **never** touches concrete
cells — it correctly restricts to soil-only. The bottom Dirichlet at row -1 applies to
**all** columns without a material check. In the CW validation geometry, this means the
bottom-CL corner (which is soil) is pinned to `T_gw` while the CL centerline stencil
(zero-flux) has already computed a different temperature.

## S4a.1 Structural hypothesis

The **bottom-row Dirichlet lacks the half-cell correction** that the Stage 3 side fix
added. Furthermore, the unconditional full-row overwrite pins the bottom-CL corner to
`T_gw` unconditionally, conflicting with the zero-flux CL symmetry BC.

If CW handles the bottom BC with a proper corner treatment (Robin or a different
discretization at the CL corner), the engine would produce a spurious temperature
fixation at (w=6.10m, depth=24.38m) — the exact cell where the Stage 3.5 masked
residual peaks (9.1–19.6°F, scaling with |ΔT|).

**Hypothesis (not yet confirmed):** the bottom-CL corner residual is structural, not
parametric. The engine's unconditional Dirichlet overwrite at row -1 conflicts with
CW's likely use of a corner BC that blends the symmetry and soil conditions.

This is documented only. No fix is proposed here — that is Stage 4b's scope.
