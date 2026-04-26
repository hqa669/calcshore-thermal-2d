# s4 — ε_top (Emissivity) Magnitude Sweep

**Knob**: `construction.emissivity_top` override via `dataclasses.replace()`.
Engine reads this via `getattr(construction, 'emissivity_top',
EMISSIVITY_DEFAULT)` at `thermal_engine_2d.py:1359`.

Default: ε_top = 0.88  (ASHRAE blanketed concrete outer surface).
ε_top affects LW radiation balance (nighttime cooling primarily).
Expected signature: one-sided amplitude correction (nighttime trough
shifts, daytime peak less affected). ASHRAE values cluster tightly
near 0.88, so a large sweep range is not physically justified.

## Sweep table

CenterRMS = S0 gate metric, window [48, 168] hr.
Authority threshold: cluster-mean CenterRMS variation ≥ 0.24°F.

| Mix | ε_top | D1 ΔT̄ (°F) | D2 lag (hr) | D3 ratio | CenterRMS (°F) | D4 TopRMS (°F) | D4 spread (°F) | S0 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.50 | +0.523 | +0.00 | 1.7586 | 0.7702 | 2.2589 | 2.0368 | **5/5** |
| MIX-11 | 0.50 | +0.619 | +0.00 | 1.6632 | 0.8364 | 2.4121 | 2.1001 | **5/5** |
| MIX-12 | 0.50 | +0.504 | +0.00 | 1.7134 | 0.7577 | 2.2325 | 2.0236 | **5/5** |
| MIX-03 | 0.50 | +0.662 | +0.00 | 1.7338 | 0.8463 | 2.2994 | 2.0263 | 4/5 |
| MIX-01 | 0.70 | +0.466 | +0.00 | 1.8422 | 0.7506 | 2.0732 | 1.8838 | **5/5** |
| MIX-11 | 0.70 | +0.558 | +0.00 | 1.7365 | 0.8098 | 2.2047 | 1.9572 | **5/5** |
| MIX-12 | 0.70 | +0.448 | +0.00 | 1.7964 | 0.7395 | 2.0509 | 1.8495 | **5/5** |
| MIX-03 | 0.70 | +0.603 | +0.00 | 1.8083 | 0.8165 | 2.1281 | 1.8922 | **5/5** |
| MIX-01 | 0.88 | +0.419 | +0.00 | 1.9102 | 0.7381 | 1.9338 | 1.6966 | **5/5** |
| MIX-11 | 0.88 | +0.508 | +0.00 | 1.7961 | 0.7915 | 2.0449 | 1.8672 | **5/5** |
| MIX-12 | 0.88 | +0.402 | +0.00 | 1.8640 | 0.7281 | 1.9154 | 1.6527 | **5/5** |
| MIX-03 | 0.88 | +0.556 | +0.00 | 1.8689 | 0.7950 | 2.0022 | 1.7338 | **5/5** |
| MIX-01 | 0.95 | +0.402 | +0.00 | 1.9350 | 0.7343 | 1.8860 | 1.6202 | **5/5** |
| MIX-11 | 0.95 | +0.490 | +0.00 | 1.8179 | 0.7856 | 1.9889 | 1.8209 | **5/5** |
| MIX-12 | 0.95 | +0.385 | +0.00 | 1.8887 | 0.7247 | 1.8692 | 1.5754 | **5/5** |
| MIX-03 | 0.95 | +0.538 | +0.00 | 1.8910 | 0.7878 | 1.9598 | 1.6680 | **5/5** |

## Authority assessment

- Cluster-mean CenterRMS variation: **0.0399°F**
- Authority (≥0.24°F): **NO**
- Cluster optimum: ε_top = **0.95**  (cluster-mean CenterRMS = 0.7482°F)
- MIX-03 Δ at optimum: -0.0072°F  (generalizes: YES)
