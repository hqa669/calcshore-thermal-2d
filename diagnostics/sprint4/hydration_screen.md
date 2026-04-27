# Sprint 4 — Hydration-Fit Screen (Cluster B)

Early window: [12, 36]hr. Threshold = 1.5× max(Reference) per metric.

## Metrics

| Mix | Set | Centerline RMS [12,36]hr (°F) | Engine rise (°F/hr) | CW rise (°F/hr) | Δ rise (°F/hr) | Engine t½ (hr) | CW t½ (hr) | Δ t½ (hr) | Result |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MIX-01 | Reference | 1.14 | 1.57 | 1.51 | 0.06 | 32.00 | 32.83 | 0.83 | ref |
| MIX-02 | Reference | 0.78 | 1.16 | 1.30 | 0.13 | 35.50 | 36.92 | 1.42 | ref |
| MIX-03 | Reference | 0.87 | 1.53 | 1.51 | 0.02 | 36.50 | 37.58 | 1.08 | ref |
| MIX-11 | Reference | 1.20 | 1.56 | 1.51 | 0.05 | 31.50 | 32.75 | 1.25 | ref |
| MIX-12 | Reference | 1.13 | 1.57 | 1.51 | 0.06 | 32.00 | 32.92 | 0.92 | ref |
| MIX-04 | Cluster B | 0.50 | 2.12 | 2.16 | 0.04 | 30.00 | 30.58 | 0.58 | hydration_pass |
| MIX-08 | Cluster B | 0.74 | 2.19 | 2.16 | 0.03 | 27.50 | 27.92 | 0.42 | hydration_pass |
| MIX-14 | Cluster B | 1.81 | 2.09 | 2.16 | 0.07 | 24.00 | 24.33 | 0.33 | hydration_fail (RMS) |
| MIX-15 | Cluster B | 0.46 | 1.17 | 1.30 | 0.13 | 42.50 | 44.33 | 1.83 | hydration_pass |

## Thresholds derived from Reference

- Centerline RMS [12,36]hr : **1.80 °F**
- Δ peak rise rate         : **0.20 °F/hr**
- Δ time-to-half-peak      : **2.13 hr**

## Cluster B classification

| Mix | Result | Failed metrics |
| --- | --- | --- |
| MIX-04 | hydration_pass | — |
| MIX-08 | hydration_pass | — |
| MIX-14 | hydration_fail | RMS |
| MIX-15 | hydration_pass | — |

**Stage 1 survivors → Stage 2 construction-parameter diff:** MIX-04, MIX-08, MIX-15
