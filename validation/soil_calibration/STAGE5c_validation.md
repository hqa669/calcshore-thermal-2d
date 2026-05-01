# STAGE5c — Validation Report (CW-Grid Alignment)

**Engine change:** `build_grid_half_mat` called with `blanket_thickness_m=0.0`  
**Grid:** `n_concrete_x=121, n_concrete_y=73` (≈0.17 ft × 1.10 ft cells)  
**Config:** `model_soil=False`, `is_submerged=True`, `blanket_thickness_m=0.0`  
**Comparison time:** t=168 hr  
**Mask:** Stage 3.5 validity mask (base, no extensions)  
**Gate:** All 9 runs ≤ 0.35°F masked max|R| (1% of 35°F ΔT_max spec)

## Residual Table

| Run | |ΔT| (°F) | Stage 5b max|R| (°F) | Stage 5c max|R| (°F) | Improvement (°F) | Gate | Wall-clock (s) |
| --- | --- | --- | --- | --- | --- | --- |
| A | 0 | 0.131 | **0.131** | +0.000 | **PASS ✓** | 1.2 |
| B | 13 | 0.230 | **0.289** | -0.059 | **PASS ✓** | 1.2 |
| C | 17 | 0.250 | **0.248** | +0.002 | **PASS ✓** | 1.2 |
| D | 13 | 0.188 | **0.173** | +0.015 | **PASS ✓** | 1.2 |
| E | 17 | 0.302 | **0.366** | -0.064 | **FAIL ✗** | 1.2 |
| F | 28 | 0.401 | **0.524** | -0.123 | **FAIL ✗** | 1.2 |
| G | 27 | 0.387 | **0.360** | +0.027 | **FAIL ✗** | 1.2 |
| H | 28 | 0.401 | **0.416** | -0.015 | **FAIL ✗** | 1.2 |
| I | 27 | 0.423 | **0.546** | -0.123 | **FAIL ✗** | 1.2 |

## Gate Verdict: **FAIL**

Failing runs: E, F, G, H, I

Investigate root cause before proceeding.

## §5.2 Geometry Verification

- Engine y[1] (concrete top): 0.000000 m  vs CW di=0: 0.000000 m  (offset: 0.00 mm)
- Engine y[iy_concrete_end]: 24.384000 m  vs CW di=48: 24.380000 m  (offset: 4.0 mm — CW grid precision)
- Reduction from Stage 5b: 24 mm → 4.0 mm

## §5.3 Lateral Alignment Check (wi=8,9,10)

| Run | wi=8 masked max|R| (°F) | wi=9 masked max|R| (°F) | wi=10 masked max|R| (°F) |
| --- | --- | --- | --- |
| F | 0.484 | 0.524 | 0.410 |
| G | 0.290 | 0.360 | 0.285 |
| H | 0.360 | 0.416 | 0.336 |
| I | 0.511 | 0.546 | 0.438 |

## Grid Info

- `blanket_thickness_m`: 0.0
- `grid_refinement`: 6
- `n_concrete_x`: 121 (`nx_full`: 133)
- `n_concrete_y`: 73 (`ny_full`: 74)
- `dx`: 0.1667 ft (0.0508 m)
