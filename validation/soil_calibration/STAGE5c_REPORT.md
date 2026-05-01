# Stage 5c — Final Report

**Date:** 2026-05-01  
**Engine prep commits:** `6a7586b` (CFL guard), `478ab92` (dy_minus guard)  
**Stage 5c work:** compare_to_cw.py + stage5c_run.py call-site changes

---

## 1. Approach

Stage 5b identified the root cause of runs F/G/H/I failures as a 24 mm vertical
offset between the engine's concrete bottom node (at y = depth_m + blanket_m =
24.384 + 0.020 = 24.404 m) and CW's deepest comparison node (di=48 at 24.380 m).
A 4×–10× refinement plateau confirmed the residual was structural, not numerical.

Stage 5c addressed this by passing `blanket_thickness_m=0.0` to `build_grid_half_mat`
at CW-comparison call sites (`compare_to_cw.py:172`, `stage5c_run.py`). With zero
blanket thickness, the blanket row j=0 still exists at y=0 for index alignment (pure-R
architecture), but y[1]==y[0]==0, so the concrete bottom sits at y=depth_m = 24.384 m —
4 mm from CW di=48 (CW's own grid precision; not further reducible without altering CW data).

### Engine prep commits required

Two engine guards were needed before the call-site change could run without NaN:

**Commit 1 — `6a7586b`**: CFL dy_min guard
- `if _use_pure_r_blanket and dy_all[0] == 0.0: dy_min = dy_all[1:].min()`
- Prevents dt_cfl = 0.0 / inf → 0 when blanket has zero thickness
- Production unaffected (dy_all[0]=0.02 for default blanket)

**Commit 2 — `478ab92`**: dy_minus[0] guard
- `if _use_pure_r_blanket and dy_minus[0, 0] == 0.0: dy_minus[0, 0] = 1.0`
- Prevents `kym * dT / 0 = 0 * dT / 0 = NaN` in the FD stencil and centerline BC
- Setting to 1.0 is safe: k_y_face[j=0] = harmonic_mean(k_blanket=0, k_concrete) = 0
- Production unaffected (dy_minus[0,0]=0.02 for default blanket)

---

## 2. Geometry Verification (§5.2)

| Coordinate | Engine (blanket=0.0) | CW | Offset |
|---|---|---|---|
| Concrete top y[1] | 0.000000 m | cw_d[0] = 0.000000 m | **0.0 mm** (perfect) |
| Concrete bottom y[iy_concrete_end] | 24.384000 m | cw_d[48] = 24.380000 m | **4.0 mm** (CW grid precision) |

The 24 mm vertical offset from Stage 5b was reduced by 83% to 4 mm. The remaining 4 mm
is the discrepancy between `depth_ft * 0.3048 = 24.384 m` and CW's 24.38 m display precision.

---

## 3. Nine-Run Gate Results

**Config:** `model_soil=False`, `is_submerged=True`, `blanket_thickness_m=0.0`, `grid_refinement=6`  
**Gate:** all 9 runs masked max|R| ≤ 0.35°F

| Run | \|ΔT\| (°F) | Stage 5b max\|R\| (°F) | Stage 5c max\|R\| (°F) | Δ | Gate |
|---|---|---|---|---|---|
| A | 0 | 0.131 | **0.131** | 0.000 | **PASS ✓** |
| B | 13 | 0.230 | **0.289** | +0.059 | **PASS ✓** |
| C | 17 | 0.250 | **0.248** | -0.002 | **PASS ✓** |
| D | 13 | 0.188 | **0.173** | -0.015 | **PASS ✓** |
| E | 17 | 0.302 | **0.366** | +0.064 | **FAIL ✗** |
| F | 28 | 0.401 | **0.524** | +0.123 | **FAIL ✗** |
| G | 27 | 0.387 | **0.360** | -0.027 | **FAIL ✗** |
| H | 28 | 0.401 | **0.416** | +0.015 | **FAIL ✗** |
| I | 27 | 0.423 | **0.546** | +0.123 | **FAIL ✗** |

### Gate Verdict: **FAIL** (runs E, F, G, H, I)

---

## 4. Root Cause Diagnosis

### 4.1 Vertical fix: confirmed effective

At di=48, wi=6 (bottom row, mid-width — away from the lateral artifact):

| Run | Stage 5b di=48 | Stage 5c di=48 | Reduction |
|---|---|---|---|
| B | +0.172°F | +0.022°F | 87% |
| E | +0.231°F | +0.035°F | 85% |
| F | +0.391°F | +0.069°F | 82% |
| G | -0.378°F | -0.066°F | 83% |
| H | -0.392°F | -0.069°F | 82% |
| I | +0.370°F | +0.058°F | 84% |

The 24 mm → 4 mm vertical alignment fix works as designed. The di=48 artifact is no longer
the dominant residual source.

### 4.2 Lateral artifact: pre-existing, now dominant

The gate failures in Stage 5c are entirely from a flat-with-depth residual at the wi=9
comparison column (~4.987 ft from the concrete edge). This artifact was present in Stage 5b
but hidden by the larger di=48 residual.

At wi=9, di=20 (mid-depth, away from top/bottom BCs):

| Run | Stage 5b | Stage 5c | Change |
|---|---|---|---|
| A | -0.094°F | -0.094°F | 0.000 |
| B | -0.210°F | -0.210°F | 0.000 |
| E | -0.291°F | -0.290°F | 0.001 |
| F | -0.357°F | -0.355°F | 0.002 |
| G | +0.157°F | +0.155°F | 0.002 |
| H | +0.232°F | +0.231°F | 0.001 |
| I | -0.397°F | -0.395°F | 0.002 |

The wi=9 residual is essentially unchanged between Stage 5b and Stage 5c. The blanket
thickness change (y-axis fix) had no effect on the x-axis lateral artifact, as expected.

**Artifact characteristics:**
- ~0.35–0.40°F constant offset across di=9 to di=40 (flat with depth → structural, not numerical)
- Magnitude scales with |ΔT| (larger lateral gradient → larger offset)
- Present at wi=8/9/10 (column spread); worst at wi=9

**Why Stage 5c FAILS worse than Stage 5b at some runs:**
The Stage 5c gate max is now at the near-bottom wi=9 cells (di=43–47). With `blanket=0.0`,
the entire y-grid shifted 20 mm upward relative to Stage 5b (concrete top at y=0 vs y=0.02).
The near-bottom temperature gradient steepens closer to the soil BC. The pre-existing lateral
x-offset at wi=9, interacting with the steeper near-bottom gradient, produces larger residuals
at di=43–47 in Stage 5c than in Stage 5b. These near-bottom cells were already in the mask.

---

## 5. §5.3 Lateral Alignment Check

| Run | wi=8 masked max\|R\| (°F) | wi=9 masked max\|R\| (°F) | wi=10 masked max\|R\| (°F) |
|---|---|---|---|
| E | 0.353 | 0.366 | 0.291 |
| F | 0.484 | 0.524 | 0.410 |
| G | 0.290 | 0.360 | 0.285 |
| H | 0.360 | 0.416 | 0.336 |
| I | 0.511 | 0.546 | 0.438 |

The lateral artifact is concentrated at wi=9 with spread to wi=8/10. It does NOT improve
with grid refinement (Stage 5b plateau analysis) and does NOT improve with vertical
alignment fixes (Stage 5c). It is a comparison-grid x-alignment structural artifact.

---

## 6. Stage 5 Status

**Stage 5 is NOT closed.** The blanket=0.0 vertical alignment fix is confirmed effective
and correct. The gate is now blocked by a pre-existing lateral comparison artifact at wi=9
(~4.987 ft from concrete edge) that was always present but previously dominated by the
di=48 vertical artifact.

---

## 7. Open Question for Stage 5d

**Root cause of wi=9 lateral artifact:**

CW comparison widths: [6.1, 5.59, ..., 1.52, ..., 0.0] m — 13 columns spaced 0.508 m apart,
spanning 6.1 m. Engine concrete half-width: 20 ft = 6.096 m. The 4 mm difference (6.100 vs
6.096 m) means CW's physical centerline and the engine's physical centerline are offset by 4 mm.

The flat-depth constant offset at wi=9 is consistent with:
1. A 4 mm x-coordinate misalignment between CW's column grid and the engine's vertex grid
   (CW half-width 6.100 m vs engine 6.096 m → 4 mm offset at CL, scaling to ~3.7 mm at wi=9)
2. A different treatment of the submerged side-BC geometry (CW may apply T_soil at a
   slightly different effective x than the engine's vertex at x=0.0 m)
3. CW grid being cell-centered vs engine vertex-centered near the lateral boundary

Stage 5d should: (a) confirm whether width_ft in geom matches CW's comparison grid span,
(b) check whether adjusting build_grid_half_mat to match CW's 6.1 m span (or adjusting
the resampling x-mapping) eliminates the wi=9 artifact for a single run, then re-run gate.

---

## 8. Summary

| Item | Status |
|---|---|
| Engine prep Commit 1 (CFL guard) | **Done** — `6a7586b` |
| Engine prep Commit 2 (dy_minus guard) | **Done** — `478ab92` |
| compare_to_cw.py call-site change | **Done** (uncommitted) |
| stage5c_run.py harness | **Done** (uncommitted) |
| Vertical alignment fix (blanket=0.0) | **Effective** — di=48 reduced 82–87% |
| 9-run gate | **FAIL** (E, F, G, H, I) |
| Root cause of gate failure | Lateral x-alignment artifact at wi=9 |
| Recommended next step | Stage 5d: investigate wi=9 x-coordinate origin |
