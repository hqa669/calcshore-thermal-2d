# PR 15 — B1 Diagnostic Sweeps

Phase 1 diagnostic for Sprint 4 PR 15 (B1 calibration: MIX-04, MIX-08).
Mixes: MIX-04 + MIX-08 (B1 targets), MIX-01 (Reference control).
Gates (S0): PeakMax Δ < ±1.0°F, PeakGrad Δ < ±2.0°F, FieldRMS < 2.0°F, CenterRMS < 1.0°F, CornerRMS < 3.0°F.

## Winner Summary

**Q1 — Any F_vert value drives both B1 mixes to S0 PASS without breaking MIX-01?**
NO — No single F_vert value achieves 5/5 on all three mixes. B1 optimum (min mean CornerRMS): F_vert = 0.00 → MIX-04 CornerRMS = 2.16°F, MIX-08 CornerRMS = 2.90°F.

**Q2 — Any R_form value drives both B1 mixes to S0 PASS without breaking MIX-01?**
NO — No single R_form value achieves 5/5 on all three mixes. B1 optimum (min mean CornerRMS): R_form = 0.0400 → MIX-04 CornerRMS = 1.71°F, MIX-08 CornerRMS = 1.26°F.

**Q3 — If neither sweep alone fixes both B1 mixes: do their optima agree (direction + magnitude)?**
F_vert optimum = 0.00 (lower than default 0.15). R_form optimum = 0.0400 (lower than default 0.0862). Direction AGREES — both sweeps push the form-face heat balance in the same direction.

**Q4 — Is residual CornerRMS at the optimum consistent with H3 catch-all?**
PARTIAL — At F_vert = 0.00, one or both B1 mixes cross the S0 threshold (MIX-04=2.16°F, MIX-08=2.90°F). Check Q3 joint fix potential (Decision C).

## Sweep 1 — F_vert

Default: `F_VERT_BY_ORIENTATION["unknown"] = 0.15` (Sprint 2 PR 8 calibration). Override via `construction.vertical_solar_factor` (no source edit required).

| Mix | F_vert | PeakMax Δ (°F) | PeakGrad Δ (°F) | FieldRMS (°F) | CenterRMS (°F) | CornerRMS (°F) | S0 pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.00 | -0.29 | +0.39 | 0.88 | 0.74 | 1.71 | **5/5** |
| MIX-01 | 0.05 | -0.29 | +0.10 | 0.88 | 0.74 | 1.74 | **5/5** |
| MIX-01 | 0.10 | -0.29 | -0.19 | 0.88 | 0.74 | 1.94 | **5/5** |
| MIX-01 | 0.15 | -0.29 | -0.48 | 0.88 | 0.74 | 2.26 | **5/5** |
| MIX-01 | 0.20 | -0.29 | -0.77 | 0.88 | 0.74 | 2.65 | **5/5** |
| MIX-01 | 0.25 | -0.29 | -1.06 | 0.88 | 0.74 | 3.10 | 4/5 |
| MIX-01 | 0.30 | -0.29 | -1.35 | 0.88 | 0.74 | 3.57 | 4/5 |
| MIX-01 | 0.40 | -0.29 | -1.92 | 0.88 | 0.74 | 4.56 | 4/5 |
| MIX-01 | 0.50 | -0.29 | -2.49 | 0.88 | 0.74 | 5.59 | 3/5 |
| MIX-04 | 0.00 | -1.46 | -1.64 | 1.37 | 0.93 | 2.16 | 4/5 |
| MIX-04 | 0.05 | -1.46 | -1.92 | 1.37 | 0.93 | 2.36 | 4/5 |
| MIX-04 | 0.10 | -1.46 | -2.18 | 1.37 | 0.93 | 2.66 | 3/5 |
| MIX-04 | 0.15 | -1.46 | -2.45 | 1.37 | 0.93 | 3.03 | 2/5 |
| MIX-04 | 0.20 | -1.46 | -2.72 | 1.37 | 0.93 | 3.45 | 2/5 |
| MIX-04 | 0.25 | -1.46 | -2.98 | 1.37 | 0.93 | 3.89 | 2/5 |
| MIX-04 | 0.30 | -1.46 | -3.24 | 1.37 | 0.93 | 4.36 | 2/5 |
| MIX-04 | 0.40 | -1.46 | -3.77 | 1.37 | 0.93 | 5.33 | 2/5 |
| MIX-04 | 0.50 | -1.46 | -4.29 | 1.37 | 0.93 | 6.33 | 2/5 |
| MIX-08 | 0.00 | -1.43 | -2.74 | 1.33 | 0.85 | 2.90 | 3/5 |
| MIX-08 | 0.05 | -1.43 | -3.03 | 1.33 | 0.85 | 3.24 | 2/5 |
| MIX-08 | 0.10 | -1.43 | -3.32 | 1.33 | 0.85 | 3.63 | 2/5 |
| MIX-08 | 0.15 | -1.43 | -3.61 | 1.33 | 0.85 | 4.05 | 2/5 |
| MIX-08 | 0.20 | -1.43 | -3.90 | 1.33 | 0.85 | 4.50 | 2/5 |
| MIX-08 | 0.25 | -1.43 | -4.19 | 1.33 | 0.85 | 4.96 | 2/5 |
| MIX-08 | 0.30 | -1.43 | -4.47 | 1.33 | 0.85 | 5.44 | 2/5 |
| MIX-08 | 0.40 | -1.43 | -5.05 | 1.33 | 0.85 | 6.43 | 2/5 |
| MIX-08 | 0.50 | -1.43 | -5.59 | 1.33 | 0.85 | 7.44 | 2/5 |

## Sweep 2 — R_form

Default: `R_FORM_CONTACT_SI = 0.0862` m²·K/W (ADR-04: ACI 347 steel form + wet concrete film). Override via monkey-patch `thermal_engine_2d.R_FORM_CONTACT_SI` (restored after sweep).

| Mix | R_form (m²·K/W) | PeakMax Δ (°F) | PeakGrad Δ (°F) | FieldRMS (°F) | CenterRMS (°F) | CornerRMS (°F) | S0 pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.0400 | -0.29 | +2.91 | 0.88 | 0.74 | 2.03 | 4/5 |
| MIX-01 | 0.0600 | -0.29 | +1.31 | 0.88 | 0.74 | 1.04 | **5/5** |
| MIX-01 | 0.0862 | -0.29 | -0.48 | 0.88 | 0.74 | 2.26 | **5/5** |
| MIX-01 | 0.1000 | -0.29 | -1.31 | 0.88 | 0.74 | 3.11 | 4/5 |
| MIX-01 | 0.1200 | -0.29 | -2.06 | 0.88 | 0.74 | 4.28 | 3/5 |
| MIX-01 | 0.1500 | -0.29 | -3.02 | 0.88 | 0.74 | 5.83 | 3/5 |
| MIX-01 | 0.2000 | -0.29 | -4.34 | 0.88 | 0.74 | 7.98 | 3/5 |
| MIX-04 | 0.0400 | -1.46 | +1.30 | 1.37 | 0.93 | 1.71 | 4/5 |
| MIX-04 | 0.0600 | -1.46 | -0.47 | 1.37 | 0.93 | 1.29 | 4/5 |
| MIX-04 | 0.0862 | -1.46 | -2.45 | 1.37 | 0.93 | 3.03 | 2/5 |
| MIX-04 | 0.1000 | -1.46 | -3.21 | 1.37 | 0.93 | 4.00 | 2/5 |
| MIX-04 | 0.1200 | -1.46 | -4.03 | 1.37 | 0.93 | 5.30 | 2/5 |
| MIX-04 | 0.1500 | -1.46 | -5.10 | 1.37 | 0.93 | 7.02 | 2/5 |
| MIX-04 | 0.2000 | -1.46 | -6.58 | 1.37 | 0.93 | 9.40 | 2/5 |
| MIX-08 | 0.0400 | -1.43 | +0.51 | 1.33 | 0.85 | 1.26 | 4/5 |
| MIX-08 | 0.0600 | -1.43 | -1.49 | 1.33 | 0.85 | 1.79 | 4/5 |
| MIX-08 | 0.0862 | -1.43 | -3.61 | 1.33 | 0.85 | 4.05 | 2/5 |
| MIX-08 | 0.1000 | -1.43 | -4.59 | 1.33 | 0.85 | 5.17 | 2/5 |
| MIX-08 | 0.1200 | -1.43 | -5.73 | 1.33 | 0.85 | 6.66 | 2/5 |
| MIX-08 | 0.1500 | -1.43 | -6.97 | 1.33 | 0.85 | 8.63 | 2/5 |
| MIX-08 | 0.2000 | -1.43 | -8.68 | 1.33 | 0.85 | 11.34 | 2/5 |

