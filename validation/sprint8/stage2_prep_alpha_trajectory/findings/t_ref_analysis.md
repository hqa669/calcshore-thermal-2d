# Sprint 8 Stage 2-prep — §3.4 T_REF_K Back-calculation

## Method

Engine uses T_REF_K = 296.15 K (23°C). CW's reference temperature is undocumented.
At T_pl = 73.0°F (295.93 K), different T_ref values shift the Arrhenius
factor, stretching/compressing the equivalent-age axis nonlinearly across datasets.
Best-fit T_ref minimizes RMS discrepancy across all ratio tests, statistics, and timesteps.

## Arrhenius Factor at T_pl = 73°F

| T_ref | T_ref (K) | Arrhenius factor |
|---|---|---|
| 20.00°C (ACI 308 / Schindler-Folliard 2002 calibration) | 293.15 | 1.212358 |
| 21.00°C | 294.15 | 1.130685 |
| 21.11°C (70°F — Schindler-Folliard paper placement temp) | 294.26 | 1.122077 |
| 22.00°C | 295.15 | 1.055013 |
| 23.00°C (engine default T_REF_K) | 296.15 | 0.984866 ← engine |
| 24.00°C | 297.15 | 0.919809 |
| 25.00°C | 298.15 | 0.859444 |

## Ratio Discrepancy Summary

| T_ref | max\|disc\| | RMS disc | mean disc |
|---|---|---|---|
| 20.00°C (ACI 308 / Schindler-Folliard 2002 calibration) | 5.39% | 3.23% | -0.04% |
| 21.00°C | 5.42% | 3.21% | -0.18% |
| 21.11°C (70°F — Schindler-Folliard paper placement temp) | 5.46% | 3.25% | -0.20% |
| 22.00°C | 5.81% | 3.86% | -0.34% |
| 23.00°C (engine default T_REF_K) | 8.51% | 5.05% | -0.50% ← engine |
| 24.00°C | 12.19% | 6.63% | -0.69% |
| 25.00°C | 16.22% | 8.49% | -0.88% |

## Best-fit T_ref

**21.00°C** — minimizes RMS discrepancy across all 18 ratio test values.
Max |discrepancy| at best-fit: **5.42%**

## Implied |Δα|/α vs Engine (23°C)

| Dataset | t=24hr | t=84hr | t=168hr | max (0–168hr) |
|---|---|---|---|---|
| lowhyd | 1.7% | 1.5% | 1.4% | 2.4% |
| highhyd | 5.5% | 1.9% | 1.0% | 121.3% |
| highhyd_b010 | 1.3% | 1.1% | 1.0% | 1.7% |