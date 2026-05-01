# STAGE5b — Validation Report (6× Grid Refinement)

**Engine change:** `build_grid_half_mat` default `grid_refinement=6`  
**Grid:** `n_concrete_x=121, n_concrete_y=73` (≈0.17 ft × 1.10 ft cells)  
**Config:** `model_soil=False`, `is_submerged=True`  
**Comparison time:** t=168 hr  
**Mask:** Stage 3.5 validity mask (|residual_A| < 0.3°F gate region)  
**Gate:** All 9 runs ≤ 0.35°F masked max|R| (1% of 35°F ΔT_max spec)

## Residual Table

| Run | |ΔT| (°F) | Stage 4b max|R| (°F) | Stage 5b max|R| (°F) | Improvement (°F) | Gate | Wall-clock (s) |
| --- | --- | --- | --- | --- | --- | --- |
| A | 0 | 0.146 | **0.131** | +0.015 | **PASS ✓** | 8.4 |
| B | 13 | 1.940 | **0.230** | +1.710 | **PASS ✓** | 8.5 |
| C | 17 | 2.391 | **0.250** | +2.141 | **PASS ✓** | 8.6 |
| D | 13 | 1.940 | **0.188** | +1.752 | **PASS ✓** | 8.7 |
| E | 17 | 2.391 | **0.302** | +2.089 | **PASS ✓** | 8.7 |
| F | 28 | 4.084 | **0.401** | +3.683 | **FAIL ✗** | 8.7 |
| G | 27 | 3.887 | **0.387** | +3.500 | **FAIL ✗** | 8.7 |
| H | 28 | 4.084 | **0.401** | +3.683 | **FAIL ✗** | 8.7 |
| I | 27 | 3.887 | **0.423** | +3.464 | **FAIL ✗** | 8.7 |

## Gate Verdict: **FAIL**

Failing runs: F, G, H, I

Investigate root cause before proceeding. Possible causes:
- Predicted post-fix residual too optimistic; try 8× refinement
- Mask edge effects at higher resolution
- Stage 5a-undetected parametric issue in specific runs

## Grid Info

- `grid_refinement`: 6
- `n_concrete_x`: 121 (`nx_full`: 133)
- `n_concrete_y`: 73 (`ny_full`: 74)
- `dx`: 0.1667 ft (0.0508 m)
