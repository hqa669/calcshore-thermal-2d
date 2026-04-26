# s3 — α_top (Solar Absorptivity) Magnitude Sweep

**Knob**: `construction.solar_absorptivity_top` override via `dataclasses.replace()`.
Engine reads this via `getattr(construction, 'solar_absorptivity_top',
SOLAR_ABSORPTIVITY_DEFAULT)` at `thermal_engine_2d.py:1358`.

Default: α_top = 0.65  (ASHRAE light-gray concrete/steel).
Expected signature: one-sided amplitude error (daytime peak shifts,
nighttime trough less affected). D3's symmetric over-oscillation makes
α_top a lower-probability dominant knob than h_conv.

## Sweep table

CenterRMS = S0 gate metric, window [48, 168] hr.
Authority threshold: cluster-mean CenterRMS variation ≥ 0.24°F.

| Mix | α_top | D1 ΔT̄ (°F) | D2 lag (hr) | D3 ratio | CenterRMS (°F) | D4 TopRMS (°F) | D4 spread (°F) | S0 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.30 | +0.167 | +0.00 | 2.2479 | 0.7228 | 1.5191 | 1.1259 | **5/5** |
| MIX-11 | 0.30 | +0.240 | +0.00 | 2.0920 | 0.7463 | 1.4763 | 1.0675 | **5/5** |
| MIX-12 | 0.30 | +0.152 | +0.00 | 2.1997 | 0.7187 | 1.5284 | 1.1352 | **5/5** |
| MIX-03 | 0.30 | +0.296 | +0.00 | 2.1678 | 0.7250 | 1.6735 | 1.3449 | **5/5** |
| MIX-01 | 0.45 | +0.276 | +0.00 | 2.1024 | 0.7170 | 1.6239 | 1.2624 | **5/5** |
| MIX-11 | 0.45 | +0.356 | +0.00 | 1.9645 | 0.7528 | 1.6561 | 1.3053 | **5/5** |
| MIX-12 | 0.45 | +0.259 | +0.00 | 2.0551 | 0.7105 | 1.6201 | 1.2587 | **5/5** |
| MIX-03 | 0.45 | +0.408 | +0.00 | 2.0390 | 0.7435 | 1.7410 | 1.4148 | **5/5** |
| MIX-01 | 0.55 | +0.348 | +0.00 | 2.0061 | 0.7235 | 1.7592 | 1.4173 | **5/5** |
| MIX-11 | 0.55 | +0.432 | +0.00 | 1.8801 | 0.7681 | 1.8344 | 1.6264 | **5/5** |
| MIX-12 | 0.55 | +0.331 | +0.00 | 1.9593 | 0.7153 | 1.7475 | 1.4058 | **5/5** |
| MIX-03 | 0.55 | +0.482 | +0.00 | 1.9537 | 0.7657 | 1.8508 | 1.5227 | **5/5** |
| MIX-01 | 0.65 | +0.419 | +0.00 | 1.9102 | 0.7381 | 1.9338 | 1.6966 | **5/5** |
| MIX-11 | 0.65 | +0.508 | +0.00 | 1.7961 | 0.7915 | 2.0449 | 1.8672 | **5/5** |
| MIX-12 | 0.65 | +0.402 | +0.00 | 1.8640 | 0.7281 | 1.9154 | 1.6527 | **5/5** |
| MIX-03 | 0.65 | +0.556 | +0.00 | 1.8689 | 0.7950 | 2.0022 | 1.7338 | **5/5** |
| MIX-01 | 0.75 | +0.491 | +0.00 | 1.8149 | 0.7600 | 2.1374 | 1.9479 | **5/5** |
| MIX-11 | 0.75 | +0.584 | +0.00 | 1.7126 | 0.8221 | 2.2779 | 1.9833 | **5/5** |
| MIX-12 | 0.75 | +0.472 | +0.00 | 1.7692 | 0.7483 | 2.1134 | 1.9225 | **5/5** |
| MIX-03 | 0.75 | +0.629 | +0.00 | 1.7845 | 0.8304 | 2.1859 | 1.9443 | **5/5** |
| MIX-01 | 0.85 | +0.562 | +0.00 | 1.7200 | 0.7886 | 2.3619 | 2.0836 | **5/5** |
| MIX-11 | 0.85 | +0.660 | +0.00 | 1.6294 | 0.8591 | 2.5267 | 2.2213 | **5/5** |
| MIX-12 | 0.85 | +0.542 | +0.00 | 1.6750 | 0.7752 | 2.3333 | 2.0782 | **5/5** |
| MIX-03 | 0.85 | +0.702 | +0.00 | 1.7005 | 0.8712 | 2.3938 | 2.0659 | 4/5 |

## Authority assessment

- Cluster-mean CenterRMS variation: **0.0809°F**
- Authority (≥0.24°F): **NO**
- Cluster optimum: α_top = **0.45**  (cluster-mean CenterRMS = 0.7267°F)
- MIX-03 Δ at optimum: -0.0515°F  (generalizes: YES)
