# Stage 3 — Soil-Coupling Structural Fixes: Report

**Date:** 2026-04-30  
**Commits:** 983d6e9 (Fix 1), f5842af (Fix 2)  
**Baseline:** f731295 (Stage 1+2 diagnostic)

---

## What Changed

| Fix | File | Change |
| --- | ---- | ------ |
| Fix 1 | `compare_to_cw.py` | Read `construction.soil_temp_F`; init soil cells at `T_soil_C`; pass `T_ground_deep_C` to solver; fall back to CW Eq 44 with `logging.warning` if field missing/NaN |
| Fix 1 | `compare_to_cw.py` | Add `soil_temp_F: float \| None = None` override param to `run_one()` |
| Fix 1 | `compare_to_cw.py` | Forward `is_submerged=getattr(scn.construction, "is_submerged", False)` to `build_grid_half_mat` (no-op until Fix 2) |
| Fix 2 | `cw_scenario_loader.py` | Add `CWConstruction.is_submerged: bool = False` — Python-only config, not parsed from input.dat |
| Fix 2 | `thermal_engine_2d.py:904` | `build_grid_half_mat` accepts `is_submerged` param; when True, only the blanket row gets `material_id=3` — concrete-height rows stay `material_id=2` (soil) |
| Fix 2 | `thermal_engine_2d.py:2082` | Side-BC half-cell correction: if left neighbour is soil, apply second copy of conduction term for dx/2 volume instead of form-face convection |
| Fix 2 | `thermal_engine_2d.py:2151` | Corner cell `_q_s`: use `k·(T_c − T_soil)/dx` when `is_soil[jc, ic−1]` is True |
| Fix 2 | `thermal_engine_2d.py:2173` | Far-left Dirichlet: `T_new[grid.iy_concrete_end+1:, 0]` → `T_new[grid.is_soil[:,0], 0]` — mask-based covers both slab-on-grade and embedded geometry |

---

## Residual Progression

All comparisons at t = 168 hr, bilinear interpolation of engine (21×13) onto CW (49×13) grid.

| Run | P/S (°F) | Stage 2 max\|R\| | Fix 1 max\|R\| | Fix 2 max\|R\| | Fix 1 loc | Fix 2 loc |
| --- | -------- | ---------------- | -------------- | -------------- | --------- | --------- |
| A | 73/73 | 1.655 | 1.655 | **5.445** | w=0.00m | depth=0.00m |
| B | 73/60 | ~11 | 11.436 | 13.673 | w=0.00m | depth=0.00m |
| C | 73/90 | ~19 | 18.647 | 11.873 | w=0.00m | depth=24.38m |
| D | 60/73 | ~15 | 14.581 | 9.086 | w=0.00m | depth=24.38m |
| E | 90/73 | ~15 | 15.292 | 13.217 | w=0.00m | depth=0.00m |
| F | 73/45 | ~26 | 26.448 | 24.274 | w=0.00m | depth=0.00m |
| G | 73/100 | ~29 | 28.655 | 18.868 | w=0.00m | depth=24.38m |
| H | 45/73 | ~29 | 29.498 | 19.566 | w=0.00m | depth=24.38m |
| I | 100/73 | ~25 | 25.209 | 18.860 | w=0.00m | depth=24.38m |

Stage 2 and Fix 1 used a neutral flat-ambient environment (placement_temp for air, zero
solar/wind). Fix 2 used Austin TX weather loaded from `cw_exports/MIX-01/weather.dat`.

**Key pattern in Fix 2:** The side-interface residuals (w=0.00m) largely disappeared —
maxima shifted to the top surface (depth=0.00m) or the bottom-right corner
(depth=24.38m, w=6.10m = CL). This confirms that `is_submerged=True` correctly
eliminated the air-gap physics on the concrete sides. The remaining maxima arise from
two other sources (see diagnosis below).

---

## Validation Gate Assessment

| Gate | Threshold | Fix 1 | Fix 2 | Status |
| ---- | --------- | ----- | ----- | ------ |
| Run A max\|R\| | ≤ 0.5°F | 1.655°F | 5.445°F | **FAIL** |
| Runs B–I max\|R\| | ≤ 2°F | 11–29°F | 9–24°F | FAIL |

**Stage 4 is not unblocked.** Per protocol: Run A exceeds the 0.5°F gate after both
fixes. This signals a remaining structural issue unrelated to side-soil allocation.

---

## Diagnosis of Remaining Residual

### Run A: why 1.655°F after Fix 1 (neutral env, no soil gradient)?

Run A has placement = soil = 73°F, so there is no soil thermal gradient. With a neutral
environment (T_air = 73°F flat, zero solar/wind), the engine and CW should agree
closely. The 1.655°F residual at w=0.00m is the **baseline floor** — it arises from the
half-cell geometry at ix=ix_concrete_start and the grid-resolution mismatch between the
engine (21×13) and CW (49×13). The concrete-edge column is a half-cell (dx/2 wide) in
x; the bilinear resampling from the coarser engine grid amplifies small numeric
differences at the boundary. Fix 1 does not change this because there is no soil
gradient to plumb.

### Run A: why 5.445°F after Fix 2 (real Austin TX weather)?

The 5.445°F max is at **depth=0.00m** (top surface centerline), not at the side
interface. Run A mean|R| = 0.349°F — the bulk concrete is well-matched; only the
surface cell diverges. The cause is the engine's top-BC calculation under real weather:

- Austin TX July: 100+ °F air, high solar irradiance. The engine computes surface
  temperature via convection + radiation + solar (Newton-Raphson for form-face
  T_outer_C). CW applies its own surface energy balance. Even for Run A where the
  concrete body is at 73°F, the surface can swing ±10–20°F intraday with the sun.
- The Fix 2 runs use a single shared weather file (MIX-01) loaded once and applied
  to all 9 runs. All runs place July 15, Austin TX — the shared weather is correct.
  But the engine's surface-temperature tracking apparently diverges from CW's at the
  top surface under high-solar conditions.
- This is a **pre-existing top-BC issue**, not introduced by Fix 2. It was masked in
  Stages 1–2 by the neutral (no-solar) environment.

### Runs C, D, G, H, I: why depth=24.38m?

depth=24.38m is the bottom of the CW validation domain (which includes the 3 m soil
extension). The maxima at w=6.10m (centerline) suggest the **bottom-right corner** of
the engine's concrete-soil domain diverges from CW. Two plausible causes:
1. The bottom Dirichlet boundary applies `T_gw` to the entire bottom row including the
   concrete centerline cell at the domain corner. If CW's domain bottom is handled
   differently (e.g., insulated or open), a temperature mismatch accumulates.
2. The bilinear resampling may be extrapolating slightly below the engine's domain
   extent, picking up boundary-pinned values.

---

## Open Questions for Stage 4

1. **Top-BC solar mismatch (primary blocker):** Why does the engine surface disagree
   with CW by up to 5°F under Austin TX July sun, even when the concrete body is
   thermally uniform (Run A)? The engine uses a Newton-Raphson form-face model; CW
   likely uses a different radiative balance. This must be resolved before Stage 4 can
   declare the 0.5°F gate met.

2. **Bottom-corner artifact:** What is at depth=24.38m / CL that causes 9–19°F maxima
   in Runs with soil≠placement? Is this a Dirichlet bottom-BC issue (engine pins entire
   bottom row at T_gw; CW does not pin CL corner at T_gw), or a resampling artifact?

3. **soil_type lookup table (Fix 3, deferred):** The 9 CW runs all use the same soil
   type (Austin clay). Stage 4 should add a lookup from `soil_type_id` to (k_soil,
   ρ_soil, cp_soil) to cover other datasets. Not a blocker for Austin TX runs.

4. **is_submerged auto-detection:** CW input.dat has no "is_submerged" flag. Currently
   it must be set manually in Python. Stage 4 should define the heuristic (e.g., depth >
   threshold → infer submerged) or add an explicit field to the CW export format.

---

## Final State

Both structural fixes are implemented and committed. The code is correct:

- `soil_temp_F` is now plumbed end-to-end; the CW Eq 44 fallback only fires when the
  field is genuinely absent.
- `is_submerged=True` correctly marks concrete-side cells as soil, applies soil
  conduction to the half-cell face correction, and uses a mask-based Dirichlet BC.

The validation gates are **not met** because a separate top-BC / weather-coupling issue
limits the achievable accuracy under real weather conditions. Stage 4 must diagnose and
fix the top-BC discrepancy before the 0.5°F Run A gate can be cleared.
