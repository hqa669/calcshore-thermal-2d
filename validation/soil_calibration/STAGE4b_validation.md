# STAGE4b — Validation Report (model_soil=False, is_submerged=True)

**Config:** `model_soil=False`, `is_submerged=True`  
**Comparison time:** t=168 hr  
**Mask:** Stage 3.5 Run-A validity mask (|residual_A| < 0.3°F)  
**Gate:** Run A ≤ 0.5°F, Runs B–I ≤ 2.0°F (masked max|R|)

## Residual Table

| Run | Placement/Soil (°F) | Stage 3.5 masked max|R| (°F) | Stage 4b masked max|R| (°F) | Improvement | Gate |
| --- | --- | --- | --- | --- | --- |
| A | 73/73 | 0.235 | **0.146** | +0.089 | **PASS** |
| B | 73/60 | 9.130 | **1.935** | +7.195 | **PASS** |
| C | 73/90 | 11.710 | **2.363** | +9.347 | **FAIL** |
| D | 60/73 | 9.130 | **1.811** | +7.319 | **PASS** |
| E | 90/73 | 11.710 | **2.539** | +9.171 | **FAIL** |
| F | 73/45 | 18.970 | **4.084** | +14.886 | **FAIL** |
| G | 73/100 | 17.420 | **3.801** | +13.619 | **FAIL** |
| H | 45/73 | 19.566 | **3.979** | +15.587 | **FAIL** |
| I | 100/73 | 17.420 | **3.977** | +13.443 | **FAIL** |

## Gate Verdict: **FAIL**

One or more runs fail the gate. See masked max|R| values above for the location of remaining residuals.
