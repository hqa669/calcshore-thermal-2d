# Sprint 8 Stage 2-prep — §3.2 + §3.3 Ratio Test Results

## Method

In a linear thermal model with Van Breugel k(α), ΔT ∝ Δk_c = k_uc·0.33·Δα.
Therefore ΔT_A/ΔT_B ≈ Δα_A/Δα_B at each timestep (same BCs, geometry, k_uc).
Discrepancy = (observed_ratio/expected_ratio − 1) × 100%.
Computed for both ratio tests (A/B and A/C) and all three ΔT statistics.

## Results

| Test | Stat | t (hr) | Observed ratio | Expected ratio | Discrepancy% |
|---|---|---|---|---|---|
| A/B | max|ΔT| | 24 | 1.6531 | 1.6105 | +2.64% |
| A/B | mean|ΔT| | 24 | 1.7071 | 1.6105 | +5.99% |
| A/B | RMS ΔT | 24 | 1.6893 | 1.6105 | +4.89% |
| A/B | max|ΔT| | 84 | 2.0364 | 2.0126 | +1.18% |
| A/B | mean|ΔT| | 84 | 2.1382 | 2.0126 | +6.24% |
| A/B | RMS ΔT | 84 | 2.1058 | 2.0126 | +4.63% |
| A/B | max|ΔT| | 168 | 2.1053 | 2.0570 | +2.34% |
| A/B | mean|ΔT| | 168 | 2.1673 | 2.0570 | +5.36% |
| A/B | RMS ΔT | 168 | 2.1400 | 2.0570 | +4.03% |
| A/C | max|ΔT| | 24 | 2.4545 | 2.6379 | -6.95% |
| A/C | mean|ΔT| | 24 | 2.4143 | 2.6379 | -8.48% |
| A/C | RMS ΔT | 24 | 2.4323 | 2.6379 | -7.79% |
| A/C | max|ΔT| | 84 | 1.9310 | 1.9876 | -2.84% |
| A/C | mean|ΔT| | 84 | 1.8783 | 1.9876 | -5.50% |
| A/C | RMS ΔT | 84 | 1.9005 | 1.9876 | -4.38% |
| A/C | max|ΔT| | 168 | 1.9048 | 1.9460 | -2.12% |
| A/C | mean|ΔT| | 168 | 1.8557 | 1.9460 | -4.64% |
| A/C | RMS ΔT | 168 | 1.8753 | 1.9460 | -3.64% |

## §3.3 Magnitude Verdict

**Max |discrepancy| across all stats, timesteps, and ratio tests: 8.48%**

**FAIL (>8.5%) — significant discrepancy; Sprint 8 dataset design needs compensation**

Thresholds: <1% = PASS, 1–5% = CAUTION, >5% = FAIL
