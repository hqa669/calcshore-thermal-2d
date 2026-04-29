# Soil Calibration — Stage 1+2 Diagnostic Report

**Date:** 2026-04-29  
**Dataset:** 9 CW runs, placement × soil ΔT matrix, hydration suppressed  
**Scope:** Characterize CW behavior (Stage 1); audit engine soil model and compare (Stage 2)  
**No engine source modifications made in this study.**

---

## 1. Dataset Summary

Nine ConcreteWorks runs covering a 5×5 ΔT matrix at placement and soil temperatures
from 45–100°F. All runs share the same geometry, mix, and hydration parameters.

| Label | Placement (°F) | Soil (°F) | ΔT (soil−placement) |
| ----- | -------------- | --------- | -------------------- |
| A | 73 | 73 | 0 (baseline) |
| B | 73 | 60 | −13 |
| C | 73 | 90 | +17 |
| D | 60 | 73 | +13 |
| E | 90 | 73 | −17 |
| F | 73 | 45 | −28 |
| G | 73 | 100 | +27 |
| H | 45 | 73 | +28 |
| I | 100 | 73 | −27 |

**Controlled constants across all 9 runs:**
- Geometry: 40×80×40 ft (width × depth × length)
- Mix: 575 lb/yd³ cement, Type I/II, no SCMs
- Hydration: `Hu_J_kg=1` (effectively suppressed); Ea=50000, τ=200, β=0.1, αu=0.1
- Form: Steel, removed at t=168 hr, R=5.67
- Soil type: Clay/Clay (side and bottom)
- Weather: Austin TX, July 15 placement (used only for top-BC; isolated from soil
  comparison via neutral flat-ambient environment in Stage 2)
- Analysis window: 0–168 hr

**Mirror pairs** (clean ΔT sign-reversal): B↔D, C↔E, F↔H, G↔I.
Note: the brief's warning that Run I had soil=60°F (ΔT=−40) was incorrect.
Actual `input.dat` inspection confirms Run I has placement=100, soil=73, ΔT=−27,
making G↔I a clean mirror at |ΔT|=27°F.

Full field map and consistency table: `STAGE1_input_field_map.md`.
All 9 runs pass geometry, hydration, form, and soil-string checks.

---

## 2. CW Behavior Characterization

### 2.1 Linearity of Soil–Concrete Coupling (S1.3)

The corrected reciprocity residual `R_corr = (T_pos − T_place_pos) + (T_neg − T_place_neg)`
measures whether each pair's deviation from its own placement temperature is equal and opposite.

| Pair | |ΔT|°F | max\|R_corr\| interior°F | mean\|R_corr\| interior°F |
| ---- | ------- | ------------------------ | ------------------------- |
| B↔D | 13 | — | 0.195 |
| C↔E | 17 | — | 0.254 |
| F↔H | 28 | — | 0.169 |
| G↔I | 27 | — | 0.275 |

**CW's soil–concrete coupling is linear to within 0.28°F (interior mean).** Max residuals
of 5–7°F appear at depth=0m (top surface) at t≈2 hr — this is CW's non-linear atmospheric
top BC (solar + evaporation), not soil-coupling non-linearity.

The linearity of the soil coupling means the 9-run dataset has 4 effective independent
degrees of freedom (Run A + one per |ΔT| level), and the engine comparison residuals are
unambiguously attributable to model-form differences, not CW non-linearities.

### 2.2 Cross-Run Temperature Trajectories (S1.4)

Three diagnostic points: side-near (w=0.51m, d=12.19m), bottom-near (CL, d=23.88m),
deep-interior (CL, d=12.19m). Key observations from `STAGE1_diagnostic_trajectories.md`:

- **Side-near and bottom-near**: temperature drifts toward soil temperature over 168 hr.
  At t=168 hr, T_final ≈ placement + 0.70×ΔT (approximately two-thirds of the way toward
  soil temperature). This is consistent with diffusion penetration (~2m from both interfaces
  at t=168 hr) that has not yet reached the interior.
- **Deep interior (CL, d=12m)**: unaffected by soil temperature variations — reaches
  placement_temp (±0.2°F) at t=168 hr regardless of soil_temp. Bottom soil does not penetrate
  to mid-depth in 168 hr. Side soil (at w=0.51m from side, d=12m) does partially penetrate.
- **Symmetry**: the deviation from placement temperature at each diagnostic point responds
  linearly and symmetrically to ΔT (consistent with S1.3 linearity verdict).

### 2.3 Thermal Front Penetration (S1.5)

1°F front from soil interface vs √t, fit over full 0–168 hr window.

| Statistic | Side | Bottom |
| --------- | ---- | ------ |
| Slope range (m/√hr) | 0.167–0.206 | 0.169–0.207 |
| R² range | 0.91–0.95 | 0.91–0.95 |
| Penetration at t=168 hr | 2.03–2.54 m | 2.03–2.54 m |

**Key findings:**
- CW treats side and bottom interfaces with essentially the same effective diffusivity:
  side slope ≈ bottom slope (within ±0.002 m/√hr) for every run. This is the CW
  behavior the engine should match.
- Slopes are mildly |ΔT|-dependent (0.167 at |ΔT|=13 vs 0.206 at |ΔT|=28), which is
  consistent with the 0.51m grid quantization — larger ΔT causes the front to cross the
  next 0.51m node sooner, inflating the slope estimate. Not a material non-linearity.
- R² of 0.91–0.95 reflects grid quantization (0.51m step function), not model deviation
  from √t behavior. Values are consistent with a single effective diffusivity across the
  full run.
- Estimated effective diffusivity: α_eff ≈ (slope/2)² ≈ 0.005–0.008 m²/hr (middle of
  the range for Clay at k=1.0–1.5 W/m·K, Cp=900 J/kg·K, ρ=1500–2000 kg/m³).

Full data in `STAGE1_penetration_analysis.md`; figures in `plots/penetration_{side|bottom}.png`.

---

## 3. Engine Soil Model Audit

Full audit: `STAGE2_engine_soil_audit.md`. Key findings:

### Grid structure (`thermal_engine_2d.py:896–1032`)

The engine builds a half-mat 2D grid with `material_id`:
- `0` = blanket (top row)
- `1` = concrete (main domain)
- `2` = soil (below concrete, and left-of-concrete BELOW grade only)
- `3` = air/inactive (left-of-concrete ABOVE grade — the side of the concrete)

**Critical:** `material_id[:iy_concrete_end+1, :ix_concrete_start] = 3` marks the entire
region to the left of the concrete, from the surface to the bottom of the concrete, as
inactive air. The side of the concrete has no soil contact in the engine.

### Soil Dirichlet BCs (`thermal_engine_2d.py:2151–2163`)

Only two Dirichlet BCs exist for soil: far-left column and bottom row, both set to `T_gw`.
`T_gw` defaults to `compute_T_gw_C(env)` (CW Eq 44: T_gw = 0.83×T_aat + 3.7) if
`T_ground_deep_C` is not passed explicitly.

### `soil_temp_F` not consumed by engine

`CWConstruction.soil_temp_F` is parsed at `cw_scenario_loader.py:439` but never passed
to `solve_hydration_2d`. In the production code (`compare_to_cw.py:178-186`), the engine
call does not include `T_ground_deep_C`, so the engine always uses the CW Eq 44 ambient
mean regardless of input.dat soil temperature.

### Structural gaps table

| Gap | Code location | Impact |
| --- | ------------- | ------ |
| Side = AIR, not soil | `thermal_engine_2d.py:999–1002` | Entire side face lacks soil coupling |
| `soil_temp_F` not plumbed | `compare_to_cw.py:178–186` | T_gw ≈ 81.3°F (Austin env mean) regardless of input |
| No soil-type property lookup | Entire codebase | `footing_subbase="Clay"` has no effect on k/ρ/Cp |
| Soil region initialized at placement_temp | `run_engine_all.py` (workaround applied) | Without fix, soil acts as 3150hr thermal buffer |

---

## 4. Engine-vs-CW Residual Summary

Engine runs used three adjustments (described in `STAGE2_engine_soil_audit.md §4`):
1. `T_ground_deep_C = soil_temp_F` passed explicitly (bypasses CW Eq 44 default)
2. Soil region initialized at `soil_temp_F` (removes 3m buffer artifact)
3. Neutral flat-ambient environment (isolates soil-coupling physics)

These adjustments fix the two BC plumbing gaps; residuals reflect model-form differences only.

| Run | Placement/Soil | max\|R\| (°F) | mean\|R\| (°F) | RMS R (°F) | Location of max |
| --- | -------------- | -------------- | --------------- | ------------ | --------------- |
| A | 73/73 | 1.655 | 0.406 | 0.617 | depth=2.03m, **w=0.00m** |
| B | 73/60 | 11.436 | 2.764 | 4.511 | depth=0.00m, **w=0.00m** |
| C | 73/90 | 18.647 | 4.262 | 7.064 | depth=2.03m, **w=0.00m** |
| D | 60/73 | 14.581 | 3.306 | 5.502 | depth=2.54m, **w=0.00m** |
| E | 90/73 | 15.292 | 3.677 | 6.021 | depth=0.00m, **w=0.00m** |
| F | 73/45 | 26.448 | 6.250 | 10.292 | depth=0.00m, **w=0.00m** |
| G | 73/100 | 28.655 | 6.578 | 10.912 | depth=2.03m, **w=0.00m** |
| H | 45/73 | 29.498 | 6.749 | 11.246 | depth=4.06m, **w=0.00m** |
| I | 100/73 | 25.209 | 5.988 | 9.853 | depth=13.21m, **w=0.00m** |

**All 9 maximum residuals are at w=0.00m (the side interface).**

For the non-baseline runs, after subtracting the 0.406°F baseline, excess mean|R| scales
approximately as **0.21 °F per °F of |ΔT|**.

Figures: `plots/run{A–I}_engine_vs_cw.png`. Summary: `STAGE2_engine_cw_comparison.md`.

---

## 5. Residual Diagnosis

Full diagnosis: `STAGE2_residual_diagnosis.md`.

### Dominant pattern: missing side soil BC (structural)

Every residual maximum falls at the side interface (w=0m). The sign of the residual
(engine warmer when soil < placement; engine cooler when soil > placement) confirms this
is not a property scaling error but an absent coupling. No parameter change to soil k, ρ,
or Cp can fix a completely absent spatial connection.

### Secondary pattern: proportional scaling with |ΔT|

mean|R| − baseline ≈ 0.21 × |ΔT|. This proportionality holds across positive and
negative ΔT and across placement-varies vs. soil-varies runs, consistent with the linearity
of CW's coupling (S1.3) and with the interpretation that the missing side BC error is
linear in the applied temperature difference.

### Baseline residual (Run A)

1.655°F at w=0m with placement=soil=73°F is residual method noise: bilinear interpolation
at the boundary between the 21x-wide engine grid and the 13x-wide CW grid, plus BC
formulation differences at the corner. This is the floor of the comparison method.

### Brief's six categories: assessment

| Category | Finding |
| -------- | ------- |
| Uniform offset | Not present |
| Proportional scaling | Not present (residual ∝ ΔT, not ∝ absolute T) |
| Side BC mismatch | **Dominant — confirmed** |
| Bottom BC mismatch | Minor or absent (max always at side, not bottom) |
| Spatial-grid artifact | Minor (≤1.7°F; baseline only) |
| Property mismatch (k/ρ/Cp) | Cannot isolate with side BC absent; second-order |

---

## 6. Stage 3 Recommendation

**Structural change first, then property calibration.**

### Priority 1 — Add soil to the concrete side (structural fix)

Modify `thermal_engine_2d.py` to allow soil cells on the lateral face of the concrete.
Specifically, remove or narrow the air region defined by:

```python
material_id[:iy_concrete_end + 1, :ix_concrete_start] = 3
```

The minimal change: for concrete-depth rows (0 through `iy_concrete_end`), set lateral
cells to `material_id=2` (soil) instead of `material_id=3` (air). This extends the soil
region to wrap around the concrete side, consistent with CW's Dirichlet BC there.

This is the dominant residual source. Without this fix, all downstream calibration is
compensating for missing geometry rather than tuning physical properties.

### Priority 2 — Plumb `soil_temp_F` through production code

In `compare_to_cw.py:178-186` (and any other caller of `solve_hydration_2d`), read
`construction.soil_temp_F` and pass it as `T_ground_deep_C`. This fix is already
implemented in `run_engine_all.py` as a workaround; it needs to become the default.

### Priority 3 — Add soil-type property lookup

`footing_subbase="Clay"` (or the raw soil type at input.dat line 466/467) should drive
a lookup of k, ρ, Cp for the soil region. The hardcoded `SOIL_PROPERTIES_2D` at
`thermal_engine_2d.py:229` should be replaced with a table keyed by soil type.
Clay values (typical: k=1.0–1.5 W/m·K, ρ=1500–2000 kg/m³, Cp=900–1200 J/kg·K) may
differ from the current defaults (k=1.50, ρ=2100, Cp=900).

### Priority 4 — Parameter calibration (after priorities 1–3)

Once the side BC is structural-corrected and soil_temp_F is plumbed:
- Re-run the 9-run matrix with the updated engine
- Residuals will now reflect soil property mismatches (not structural gaps)
- Fit k_soil, ρ_soil, Cp_soil to minimize residuals; use pairs B↔D and C↔E (|ΔT|=13/17)
  for calibration and F↔H and G↔I (|ΔT|=28/27) for hold-out validation

### Expected residual after Priority 1

The mean|R| should drop from ~4–7°F (non-baseline runs) to near the Run A baseline
floor (~0.4–0.6°F plus any remaining property mismatch). The max|R| at w=0m should
drop from 11–30°F to near the baseline 1.7°F. If significant residual remains at
w=0m after the structural fix, the soil k/ρ/Cp values are the next variable to tune.

---

## File Inventory

### Stage 1

| File | Contents |
| ---- | -------- |
| `check_inputs.py` | S1.1 input.dat consistency check |
| `STAGE1_input_field_map.md` | Field map + consistency table (all 9 pass) |
| `plot_cw_distributions.py` | S1.2 per-run 3-panel T field plots |
| `plots/run{A-I}_distribution.png` | CW T fields at t=0/84/168 hr (9 files) |
| `plots/all_runs_t168.png` | 3×3 composite at t=168 hr |
| `reciprocity_check.py` | S1.3 linearity check |
| `STAGE1_reciprocity_check.md` | Reciprocity results — interior mean < 0.28°F (linear) |
| `plots/reciprocity_{BD\|CE\|FH\|GI}.png` | Residual fields (4 files) |
| `diagnostic_trajectories.py` | S1.4 cross-run cooling profiles |
| `STAGE1_diagnostic_trajectories.md` | Symmetry + saturation + linearity tables |
| `plots/trajectory_{side_near\|bottom_near\|deep_interior}.png` | T(t) curves (3 files) |
| `penetration_analysis.py` | S1.5 √t front analysis |
| `STAGE1_penetration_analysis.md` | Slope / R² table — symmetric side/bottom coupling |
| `plots/penetration_{side\|bottom}.png` | Penetration vs √t (2 files) |

### Stage 2

| File | Contents |
| ---- | -------- |
| `STAGE2_engine_soil_audit.md` | Engine soil model audit — code locations, gaps |
| `run_engine_all.py` | S2.2 engine runs on all 9 input.dat files |
| `engine_runs/run{A-I}_t{0\|84\|168}.csv` | Engine T field slices (gitignored) |
| `engine_runs/manifest.json` | Run status + metadata (gitignored) |
| `engine_runs/grid_info.json` | Concrete subgrid dimensions (gitignored) |
| `compare_engine_cw.py` | S2.3 engine-vs-CW residual figures |
| `STAGE2_engine_cw_comparison.md` | Residual summary table |
| `plots/run{A-I}_engine_vs_cw.png` | CW \| engine \| residual (9 files) |
| `STAGE2_residual_diagnosis.md` | S2.4 root cause diagnosis |
| `STAGE1_2_REPORT.md` | This file |
