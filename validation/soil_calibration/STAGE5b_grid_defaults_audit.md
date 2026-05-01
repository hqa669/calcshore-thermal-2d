# STAGE5b — Grid Resolution Defaults Audit (S5b.1 + S5b.2)

## 1. Where the defaults live

**File:** `thermal_engine_2d.py`

Two module-level constants (added in this commit, lines ~901-906):
```python
NATIVE_N_CONCRETE_X = 21   # vertex-centered: n-1=20 intervals across half-width (1.0 ft @ W=40 ft)
NATIVE_N_CONCRETE_Y = 13   # vertex-centered: n-1=12 intervals across depth      (6.67 ft @ D=8 ft)
```

The function `build_grid_half_mat` (line ~907) now has these parameters:
```python
def build_grid_half_mat(
    width_ft: float,
    depth_ft: float,
    grid_refinement: int = 6,        # NEW — default 6× refinement
    n_concrete_x: int | None = None, # CHANGED — was int = 21; None → derive from grid_refinement
    n_concrete_y: int | None = None, # CHANGED — was int = 13; None → derive from grid_refinement
    ...
```

Derivation at function body (inserted before line 951):
```python
if n_concrete_x is None:
    n_concrete_x = (NATIVE_N_CONCRETE_X - 1) * grid_refinement + 1
if n_concrete_y is None:
    n_concrete_y = (NATIVE_N_CONCRETE_Y - 1) * grid_refinement + 1
```

Vertex-centered formula `(N_native - 1) * k + 1` ensures native node positions are preserved
at all refinement levels (the 6× grid nests the native grid exactly).

## 2. Before / After

| Parameter | Before (Stage 5a) | After (Stage 5b) |
|---|---|---|
| `n_concrete_x` default | `21` (hardcoded) | `None` → `(21-1)*6+1 = 121` via `grid_refinement=6` |
| `n_concrete_y` default | `13` (hardcoded) | `None` → `(13-1)*6+1 = 73` via `grid_refinement=6` |
| `grid_refinement` param | does not exist | `int = 6` |
| Cell size dx (W=40 ft) | `40*0.3048/2 / 20 = 0.3048 m (1.000 ft)` | `40*0.3048/2 / 120 = 0.0508 m (0.167 ft)` |
| Cell size dy (D=8 ft)  | `8*0.3048 / 12 = 0.2032 m (6.667 ft)` | `8*0.3048 / 72 = 0.03387 m (1.111 ft)` |
| Full grid nx (default, model_soil=False) | `33` (12+21) | `133` (12+121) |
| Full grid ny (default, model_soil=False) | `14` (1+13) | `74` (1+73) |
| Full grid ny (model_soil=True)            | `29` (1+13+15) | `89` (1+73+15) |

## 3. Design choice: Option B (grid_refinement kwarg)

Option A: change `n_concrete_x = 21` → `n_concrete_x = 121` directly.
Option B (chosen): add `grid_refinement: int = 6` kwarg, derive n_concrete_x/y from it.

**Why Option B:**
- Self-documenting: `grid_refinement=6` reads as a deliberate physical decision (6× refinement
  to close the 0.35°F gate). Bare `n_concrete_x=121` is a magic number.
- Escape hatch: `grid_refinement=1` instantly restores native resolution without re-computing
  node counts — useful for CI speed, fast iteration, and legacy comparisons.
- Named constants: `NATIVE_N_CONCRETE_X / _Y` make the baseline explicit for future stages.
- Backward-compatible: the only existing caller that overrides n_concrete_x/y explicitly
  (`stage5a_grid_sweep.py:80-85`) still works unchanged since explicit kwargs override
  grid_refinement.

## 4. Test changes

All existing tests that asserted grid shape at native resolution were updated to pass
`grid_refinement=1` explicitly. No test assertions were changed — only the call site was
pinned to native so the test continues to document native behavior.

| File | Change |
|---|---|
| `tests/test_grid_2d.py` | Module fixture → `grid_refinement=1`; 4 individual calls pinned |
| `tests/test_conduction_2d.py` | `test_regression_m0_grid_builder` — 2 calls pinned |
| `tests/test_side_bc_2d.py` | `test_regression_all_previous` — 1 call pinned |
| `tests/test_top_bc_2d.py` | `test_regression_all_previous` — 1 call pinned |
| `tests/test_hydration_2d.py` | `test_regression_m0_m1` — 1 call pinned |
| `tests/test_pr8_calibration.py` | `_run_full2d` pinned; gates calibrated at native resolution |
| `tests/test_pr11_soil_activation.py` | `_run_mix01` pinned; sprint-2 fixture at native resolution |

**Rationale for pinning vs. updating:** regression tests that check specific cell counts,
shape constants, or compare against stored fixtures compute the right answer at native
resolution. Changing the assertions to 6× values would convert a "native grid regression"
into a "6× grid regression" — confusing and incorrect for their intent.

Tests that are purely behavioral (BC modes differ, BCs cool the slab, etc.) continue to run
at the new 6× default without pinning.

**Pre-existing failures (27) are unchanged by this commit.** Zero new failures introduced.
Stage 5b validation grid tests use the 6× default implicitly by calling `build_grid_half_mat`
with no grid resolution arguments.
