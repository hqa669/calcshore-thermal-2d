# STAGE2 — Engine Soil Model Audit

## Summary

**Engine treats soil as a 2D transient region, NOT as a boundary condition.**
**Side and bottom soil are treated ASYMMETRICALLY** — only the bottom has a soil region;
the side of the concrete is against air (form-face atmospheric BC).

---

## 1. How Soil is Modeled in the Engine

### Grid structure (`thermal_engine_2d.py:896–1032, :43`)

The engine builds a half-mat 2D grid:

```
material_id:
  0 = blanket (top row)
  1 = concrete (main domain)
  2 = soil     (below concrete, and left-of-concrete extension below grade)
  3 = air      (left-of-concrete extension ABOVE grade = inactive)
```

Key: the lateral soil extension (`n_soil_x_ext=12` cells to the left of x=0) covers only the
**below-grade region** (`y > iy_concrete_end`). For rows **above** the concrete bottom (i.e.,
the side of the concrete wall), those cells are marked `material_id=3` (air/inactive):

```python
# thermal_engine_2d.py:999-1002
material_id[: iy_concrete_end + 1, :ix_concrete_start] = 3
```

So the concrete's SIDE FACE (x=0, all depths) is adjacent to AIR cells, not soil cells.

### Soil-region PDE (`thermal_engine_2d.py:1504, :1576`)

The soil region (material_id=2) participates in the full 2D explicit FD update.
Properties: `SOIL_PROPERTIES_2D = {'k': 1.50, 'rho': 2100.0, 'Cp': 900.0}` (`:229`).
These are **hardcoded** — no lookup from soil-type string.

### Boundary conditions on the soil region (`thermal_engine_2d.py:2151–2163`)

- **(d)** Far-left soil column (`i=0`, soil rows): Dirichlet at `T_gw`
- **(e)** Bottom row entirely: Dirichlet at `T_gw`

`T_gw` source: computed from `T_ground_deep_C` argument if provided, else from
`compute_T_gw_C(environment)` via CW Eq 44 (`T_gw = 0.83 * T_aat_C + 3.7`).

### Side BC (concrete form face) (`thermal_engine_2d.py:1498–1499, :2040`)

The concrete's SIDE is governed by the form-face atmospheric BC (convection + ground
surface temperature via Barber model — NOT soil Dirichlet). In `full_2d` mode:

```python
_h_side = h_side_convective(environment.cw_ave_max_wind_m_s, _r_blanket, _r_form_contact)
# ground_surface_temperature_C(t_hrs, env, ...) = Barber-style ambient tracking
```

The form face sees `T_gnd_C` = lagged/damped ambient, not `soil_temp_F`.

---

## 2. What `soil_temp_F` Gets vs. Doesn't Get

### What the parser extracts (`cw_scenario_loader.py:88–91`)

```python
'placement_temp_F':    438,   # → CWConstruction.placement_temp_F  ← used
'soil_temp_F':         439,   # → CWConstruction.soil_temp_F       ← NOT used by engine
'footing_subbase':     446,   # → CWConstruction.footing_subbase   ← NOT used by engine
```

### What the engine consumes

- `solve_hydration_2d` never reads `construction.soil_temp_F`.
- `T_ground_deep_C` arg is `None` by default → `compute_T_gw_C(env)` (CW Eq 44).
- For Austin summer: `T_aat ≈ 83.5°F → T_gw ≈ 81.3°F` — fixed at env-mean, regardless of the
  9 different soil temperatures in this dataset.
- `footing_subbase` string ("Clay" at idx 446, "Limestone" in practice due to a CW labeling
  mismatch) is stored on `CWConstruction` but never consulted by the engine. No soil-type →
  k/ρ/Cp lookup table exists anywhere in the codebase.

### Soil temperature profile (input.dat lines 487–495)

The input.dat has 9 additional soil temperature values (one per depth row?) at indices 486–494.
These are NOT read by `parse_cw_dat`. They mirror the `soil_temp_F` value across all depths.

---

## 3. Key Structural Gaps

| Gap | Location | Impact |
| --- | -------- | ------ |
| Side of concrete = AIR, not soil | `material_id[:iy_concrete_end+1, :ix_concrete_start] = 3` | Engine ignores soil on the SIDE of the concrete; CW applies soil_temp_F as Dirichlet on both side and bottom |
| `soil_temp_F` not plumbed to T_gw | `compare_to_cw.py:178-186` (no T_ground_deep_C passed) | Default T_gw ≈ 81.3°F (Austin env mean) regardless of input.dat soil_temp |
| No soil-type property lookup | Entire codebase | `footing_subbase="Clay"` has no effect on k_soil, ρ_soil, Cp_soil; always uses hardcoded defaults |
| Soil region initialized at placement_temp | Caller in `run_all.py` / `compare_to_cw.py` | If soil is not initialized at soil_temp, the 3m-deep soil region acts as a 3000hr thermal buffer, effectively decoupling the soil BC from the concrete for the 168hr window |

---

## 4. Engine-vs-CW Protocol for This Study

For the Stage 1+2 comparison, `run_engine_all.py` applies three adjustments:

1. **Explicit T_ground_deep_C = soil_temp_F** — bypasses the CW Eq 44 default; gives engine
   the correct soil Dirichlet BC value.
2. **Soil region initialized at soil_temp_F** — `T_initial[grid.is_soil] = T_soil_C`; removes the
   3m buffer artifact.
3. **Neutral flat-ambient environment** — constant temperature = placement_temp_F, no solar, no
   wind; isolates soil-coupling physics from top-BC variability.

Adjustment 1 and 2 are **not changes to engine source code** — they are initial conditions and
arguments passed to `solve_hydration_2d`. Adjustment 3 affects the TOP BC only (not soil).

The residuals observed in S2.3 AFTER these adjustments isolate **model-form differences** in
how soil is spatially represented:

- Side-face BC: engine uses atmospheric form-face, CW uses soil Dirichlet
- Spatial discretization: engine 13y×21x concrete grid, CW 49y×13x
- Soil thermal properties: engine hardcoded k=1.5, ρ=2100, Cp=900 vs CW's Clay values (unknown)

---

## 5. File References

| File | Lines | Topic |
| ---- | ----- | ----- |
| `thermal_engine_2d.py` | :43, :229 | Material IDs and soil properties |
| `thermal_engine_2d.py` | :896–1032 | `build_grid_half_mat` — soil region geometry |
| `thermal_engine_2d.py` | :999–1002 | Air cells on concrete SIDE |
| `thermal_engine_2d.py` | :2151–2163 | Soil Dirichlet BCs (bottom + far-left) |
| `thermal_engine_2d.py` | :734–750 | `compute_T_gw_C` (CW Eq 44) |
| `cw_scenario_loader.py` | :42–99 | `CW_DAT_INDEX` — field index map |
| `cw_scenario_loader.py` | :262–274 | `CWConstruction.soil_temp_F` field (stored, not consumed) |
| `compare_to_cw.py` | :178–186 | Engine call — `T_ground_deep_C` not passed |
