# PR 16 — B2 Diagnostic Ablation (MIX-15)

Phase 1 diagnostic for Sprint 4 PR 16 (B2 calibration: MIX-15 cold placement).
Mixes: MIX-15 (B2 target), MIX-01 (Reference control).
Gates (S0): PeakMax Δ < ±1.0°F, PeakGrad Δ < ±2.0°F, FieldRMS < 2.0°F, CenterRMS < 1.0°F, CornerRMS < 3.0°F.

**Ablation question (binary):** does MIX-15 PeakMax Δ vary across the sweep,
or is it invariant the way B1 (MIX-04, MIX-08) was?

## Phase 1 Verdict

| | F_vert sweep | R_form sweep |
| --- | --- | --- |
| PeakMax Δ at default values | -3.01°F | -3.01°F |
| PeakMax Δ range (max−min) | 0.41°F | 0.08°F |
| PeakGrad Δ range | 2.24°F | 2.92°F |
| FieldRMS range | 0.00°F | 0.00°F |
| CornerRMS range | 1.79°F | 1.80°F |

**MIX-01 holds 5/5 at defaults:** YES

**Authority thresholds:** ≥0.5°F PeakMax variation → boundary-physics-authoritative; <0.1°F → hydration-routed (B1 precedent); 0.1–0.5°F → marginal.

**F_vert verdict:** MARGINAL (variation 0.41°F in [0.1, 0.5)°F — user decision required)
**R_form verdict:** HYDRATION-ROUTED (variation 0.08°F < 0.1°F threshold — B1 precedent)
**Phase 1 verdict (combined):** **MARGINAL**

**Proposed next step:** WAIT — marginal authority. Report to user for decision before Phase 2.

## Sweep 1 — F_vert

Default: `F_VERT_BY_ORIENTATION["unknown"] = 0.15` (Sprint 2 PR 8 calibration).
Override via `construction.vertical_solar_factor` (no source edit required).

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
| MIX-15 | 0.00 | -3.05 | -7.32 | 4.06 | 2.07 | 3.38 | 0/5 |
| MIX-15 | 0.05 | -3.05 | -7.65 | 4.06 | 2.07 | 2.90 | 1/5 |
| MIX-15 | 0.10 | -3.05 | -7.98 | 4.06 | 2.07 | 2.45 | 1/5 |
| MIX-15 | 0.15 | -3.01 | -8.25 | 4.06 | 2.07 | 2.05 | 1/5 |
| MIX-15 | 0.20 | -2.96 | -8.53 | 4.06 | 2.07 | 1.75 | 1/5 |
| MIX-15 | 0.25 | -2.91 | -8.80 | 4.06 | 2.07 | 1.59 | 1/5 |
| MIX-15 | 0.30 | -2.86 | -9.08 | 4.06 | 2.07 | 1.61 | 1/5 |
| MIX-15 | 0.40 | -2.76 | -9.56 | 4.06 | 2.07 | 2.13 | 1/5 |
| MIX-15 | 0.50 | -2.64 | -9.33 | 4.06 | 2.07 | 2.99 | 1/5 |

## Sweep 2 — R_form

Default: `R_FORM_CONTACT_SI = 0.0862` m²·K/W (ADR-04: ACI 347 steel form + wet concrete film).
Override via monkey-patch `thermal_engine_2d.R_FORM_CONTACT_SI` (restored after sweep).

| Mix | R_form (m²·K/W) | PeakMax Δ (°F) | PeakGrad Δ (°F) | FieldRMS (°F) | CenterRMS (°F) | CornerRMS (°F) | S0 pass |
| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.0400 | -0.29 | +2.91 | 0.88 | 0.74 | 2.03 | 4/5 |
| MIX-01 | 0.0600 | -0.29 | +1.31 | 0.88 | 0.74 | 1.04 | **5/5** |
| MIX-01 | 0.0862 | -0.29 | -0.48 | 0.88 | 0.74 | 2.26 | **5/5** |
| MIX-01 | 0.1000 | -0.29 | -1.31 | 0.88 | 0.74 | 3.11 | 4/5 |
| MIX-01 | 0.1200 | -0.29 | -2.06 | 0.88 | 0.74 | 4.28 | 3/5 |
| MIX-01 | 0.1500 | -0.29 | -3.02 | 0.88 | 0.74 | 5.83 | 3/5 |
| MIX-01 | 0.2000 | -0.29 | -4.34 | 0.88 | 0.74 | 7.98 | 3/5 |
| MIX-15 | 0.0400 | -2.98 | -7.11 | 4.06 | 2.07 | 3.68 | 0/5 |
| MIX-15 | 0.0600 | -2.99 | -7.64 | 4.06 | 2.07 | 2.73 | 1/5 |
| MIX-15 | 0.0862 | -3.01 | -8.25 | 4.06 | 2.07 | 2.05 | 1/5 |
| MIX-15 | 0.1000 | -3.02 | -8.54 | 4.06 | 2.07 | 1.99 | 1/5 |
| MIX-15 | 0.1200 | -3.03 | -8.91 | 4.06 | 2.07 | 2.17 | 1/5 |
| MIX-15 | 0.1500 | -3.05 | -9.40 | 4.06 | 2.07 | 2.74 | 1/5 |
| MIX-15 | 0.2000 | -3.05 | -10.03 | 4.06 | 2.07 | 3.78 | 0/5 |

