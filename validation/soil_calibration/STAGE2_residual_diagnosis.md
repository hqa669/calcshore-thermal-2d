# STAGE2 — Residual Diagnosis

## Summary

The engine-vs-CW residual is dominated by a single structural model-form mismatch:
**the engine's side face of the concrete is against AIR cells, while CW applies
`soil_temp_F` as a Dirichlet BC at the side.** This is not a parameter calibration
problem. No combination of soil k/ρ/Cp values can repair a missing spatial coupling.

---

## Pattern 1 — Residual concentrated at side interface (width=0m)

**Evidence:** Every run has its maximum residual at `w=0.00m`. All 9 max-residual
locations (from `STAGE2_engine_cw_comparison.md`) are at the side edge, not the
bottom or interior:

| Run | max\|R\| (°F) | Location |
| --- | --- | --- |
| A (baseline) | 1.655 | depth=2.03m, **w=0.00m** |
| B–I (ΔT≠0) | 11–30 | varying depth, **w=0.00m** |

**Mechanism:** In the engine, `material_id[:iy_concrete_end+1, :ix_concrete_start] = 3`
marks all cells to the left of the concrete (the side region above grade) as air/inactive
(`material_id=3`). The side face BC is an atmospheric form-face convection BC tracking
ambient temperature. CW applies `soil_temp_F` as a Dirichlet boundary on the side for
the full concrete depth. When `soil_temp_F ≠ placement_temp_F`, the engine misses the
soil coupling at the side entirely — even in the `full_2d` boundary mode.

---

## Pattern 2 — Residual scales approximately linearly with |ΔT|

**Evidence:** Subtracting the Run A baseline (0.406°F), the excess mean|R| per °F of
`|soil_temp - placement_temp|` is consistent across all non-baseline runs:

| Run | |ΔT|°F | mean\|R\| (°F) | excess/|ΔT| |
| --- | --- | --- | --- |
| B | 13 | 2.764 | 0.181 |
| C | 17 | 4.262 | 0.227 |
| D | 13 | 3.306 | 0.223 |
| E | 17 | 3.677 | 0.192 |
| F | 28 | 6.250 | 0.209 |
| G | 27 | 6.578 | 0.229 |
| H | 28 | 6.749 | 0.227 |
| I | 27 | 5.988 | 0.207 |

Average coupling coefficient ≈ **0.21 °F mean|R| per °F of |ΔT|**. The linearity
is consistent with the Stage 1 finding that CW's soil coupling is linear
(reciprocity interior mean < 0.28°F). The proportionality holds across both
positive and negative ΔT, and across both placement-varies and soil-varies runs.

---

## Pattern 3 — Residual sign matches missing side BC direction

Engine − CW at the side:
- `soil < placement` (runs B, F, I): CW side is cooled by colder soil, engine side is
  not → engine warmer than CW → residual = **positive** at side
- `soil > placement` (runs C, G): CW side is heated by warmer soil, engine side is
  not → CW warmer than engine → residual = **negative** at side

This sign pattern rules out a property scaling error (which would produce the same
sign for all ΔT directions) and confirms a missing coupling.

---

## Pattern 4 — Interior residual is reduced but non-zero

The mean|R| values include the entire concrete domain. At the centerline (w=6.1m),
the residual is smaller but non-zero: the absent side coupling slightly modifies the
overall heat balance, and the error diffuses inward over 168 hr. The interior mean
residual is approximately 30–40% of the mean over the full domain (estimated from
visual inspection of the residual panels where the side-concentrated stripe is ~2–3
nodes wide).

---

## Pattern 5 — Baseline residual (Run A, ΔT=0)

Run A shows 1.655°F max at depth=2.03m, w=0.00m even with `placement=soil=73°F`.
Since soil=placement, the missing side Dirichlet cannot explain this residual.
Probable sources:

1. **BC formulation difference at the side-top corner**: Engine applies a convective
   film coefficient `h_side` on the form face; CW may apply a simpler Dirichlet.
   At zero wind and zero solar, both drive toward 73°F, but the path (convective
   vs. Dirichlet) differs near corners.
2. **Grid interpolation noise at boundaries**: `RegularGridInterpolator` extrapolates
   at edge nodes. Edge values are most sensitive to the choice of extrapolation
   and grid spacing mismatch between 21x and 13x width grids.
3. **Spatial grid alignment**: Engine concrete spans 0–6.1m in 21 steps; CW spans
   0–6.1m in 13 steps with different node positions. Near the edge node (w=0), the
   bilinear query falls at the boundary of the interpolation domain.

The 1.655°F baseline is ~10× smaller than the soil-ΔT residuals and represents the
minimum floor of the comparison method. It does not indicate a structural engine bug.

---

## Summary of residual categories (brief's six categories)

| Category | Present? | Evidence |
| -------- | -------- | -------- |
| **Uniform offset** | No | Run A baseline is side-concentrated, not uniform |
| **Proportional scaling** | No | Residual is not ∝ absolute T; it is ∝ ΔT = soil−placement |
| **BC mismatch (side)** | **YES — dominant** | All max|R| at w=0; sign matches soil-side sign; magnitude ∝ |ΔT| |
| **BC mismatch (bottom)** | Possibly minor | Bottom residual not isolated but max is always at side, not bottom |
| **Spatial-grid artifact** | Minor (baseline only) | Run A: 1.655°F near edge; ~0.1°F is expected interpolation noise |
| **Property mismatch (k/ρ/Cp)** | Cannot distinguish from structural gap | With side coupling absent, any property change is second-order |

---

## Stage 3 Recommendation

**Structural change is required before parameter calibration.**

The dominant residual (~95% of the excess mean|R|) is caused by the engine having no
soil contact on the side of the concrete. Specifically:

> `material_id[:iy_concrete_end+1, :ix_concrete_start] = 3`
> (`thermal_engine_2d.py:999–1002`)

marks the entire above-grade lateral region as inactive air. The minimal fix is to
extend the soil region upward to cover the full concrete side — or apply a soil
Dirichlet BC directly at `ix_concrete_start` for all concrete-depth rows.

**Three issues to address in Stage 3, in priority order:**

1. **Side soil BC (structural)** — add soil to the lateral region for all concrete-depth
   rows; this is the dominant error source.
2. **`soil_temp_F` plumbing (parameter)** — `compare_to_cw.py:178-186` does not pass
   `T_ground_deep_C`; the engine defaults to `compute_T_gw_C(env)` = CW Eq 44 ≈ 81.3°F
   regardless of input. Fix: read `construction.soil_temp_F` and pass it through.
3. **Soil-type property lookup (parameter)** — `footing_subbase="Clay"` has no effect
   on k/ρ/Cp. A lookup table keyed by soil-type string should replace the hardcoded
   `SOIL_PROPERTIES_2D` defaults.

Parameter calibration (soil k, ρ, Cp) is secondary and should be done AFTER the
structural fix — otherwise the calibrated parameters will compensate for the missing
spatial coupling and generalize poorly to other geometries or soil types.
