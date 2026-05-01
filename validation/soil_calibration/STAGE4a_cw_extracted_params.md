# STAGE4a — CW-Extracted Thermal Parameters

## Structural finding (preface to parameter extraction)

**CW's model** for the 40ft×80ft geometry:
- 80ft (24.38m) of concrete, all-concrete domain
- Side BC: Dirichlet T_soil at w=0.00m for ALL depths (0–24.38m)
- Bottom BC: Dirichlet T_soil at depth=24.38m for all widths
- No soil domain below the concrete

**Engine's model** (Stage 3 Fix 2):
- 80ft (24.40m) of concrete
- 3m of soil below the concrete (y=24.40m to 27.40m)
- Dirichlet T_gw at depth=27.40m (3m below the concrete bottom)

The soil buffer prevents the concrete bottom from reaching T_gw as quickly as CW's direct Dirichlet.  This causes the 9–19°F residual.

## Method summary

- **M1** (side penetration slope): 1°F front vs √t from side Dirichlet. Slope → α_c via `slope = 2 × erfinv(1−1/|ΔT|) × √α_c`. Grid resolution 0.508m quantizes the front, degrading precision.
- **M2** (bottom penetration slope): same method from bottom Dirichlet. Both M1 and M2 penetrate CONCRETE (CW has no soil below), so both give α_c.
- **M3** (erfc profile fit): `T(x,t) − T_placement = A × erfc(x/(2√(α_c t)))` at t=24,48,96,168hr with `scipy.curve_fit`. A ≈ ΔT (confirms Dirichlet at x=0). α_c recovered directly — most reliable estimator.
- **M4** (structural proof): compare engine vs CW concrete-bottom temperature at t=168hr, CL column. Gap = soil-buffer artifact.

## Per-run α_c results

| Run | ΔT (°F) | M1 α_c (m²/hr) | M2 α_c (m²/hr) | M3 α_c (m²/hr) |
| --- | --- | --- | --- | --- |
| B | -13 | 0.00447 | 0.00455 | 0.00507 |
| C | +17 | 0.00466 | 0.00473 | 0.00551 |
| D | +13 | 0.00475 | 0.00483 | 0.00542 |
| E | -17 | 0.00411 | 0.00418 | 0.00515 |
| F | -28 | 0.00447 | 0.00452 | 0.00517 |
| G | +27 | 0.00488 | 0.00493 | 0.00544 |
| H | +28 | 0.00475 | 0.00480 | 0.00530 |
| I | -27 | 0.00433 | 0.00439 | 0.00520 |

## M3 α_c by timestamp

| Run | t=24hr | t=48hr | t=96hr | t=168hr |
| --- | --- | --- | --- | --- |
| B | nan | 0.00517 | 0.00505 | 0.00500 |
| C | 0.00575 | 0.00551 | 0.00541 | 0.00537 |
| D | nan | 0.00550 | 0.00539 | 0.00536 |
| E | 0.00541 | 0.00517 | 0.00505 | 0.00497 |
| F | nan | 0.00526 | 0.00515 | 0.00509 |
| G | 0.00568 | 0.00544 | 0.00534 | 0.00530 |
| H | nan | 0.00539 | 0.00527 | 0.00523 |
| I | 0.00547 | 0.00521 | 0.00509 | 0.00503 |

## Aggregate α_c statistics

| Method | Mean (m²/hr) | Std (m²/hr) | Rel std |
| --- | --- | --- | --- |
| M1 side slope | 0.00455 | 0.00024 | 5.2% |
| M2 bot slope | 0.00462 | 0.00024 | 5.2% |
| M3 erfc fit (mean over timestamps) | 0.00528 | 0.00015 | 2.8% |

## M4 — Structural proof: bottom-temperature gap

| Run | ΔT (°F) | CW bot T (°F) | Eng bot T (°F) | Gap (°F) |
| --- | --- | --- | --- | --- |
| B | -13 | 60.01 | 69.04 | +9.03 |
| C | +17 | 90.00 | 78.18 | -11.82 |
| D | +13 | 73.00 | 63.96 | -9.05 |
| E | -17 | 73.00 | 84.82 | +11.82 |
| F | -28 | 45.00 | 64.48 | +19.48 |
| G | +27 | 100.00 | 81.22 | -18.78 |
| H | +28 | 73.00 | 53.52 | -19.48 |
| I | -27 | 73.00 | 91.78 | +18.78 |

The CW concrete-bottom temperature equals T_soil (Dirichlet at depth 24.38m).
The engine concrete-bottom temperature is ~7–10°F warmer, buffered by the
3m soil layer below the concrete.  This gap scales with |ΔT| and matches
the Stage 3.5 masked residual (9.1–19.6°F) at the bottom-CL corner.

## Key interpretation

1. **α_s**: NOT extractable (CW has no soil domain). Irrelevant for the diagnosis.
2. **e_c/e_s admittance**: NOT applicable. CW uses Dirichlet BCs at all
   soil-concrete interfaces; the admittance concept does not apply.
3. **α_c**: Extractable. CW M3 erfc fit gives α_c ≈ 0.0050–0.0056 m²/hr.
   Engine α_c (at α_hyd=0.8) ≈ 0.00449 m²/hr — within ~10–20% of CW.
   This is a secondary factor; see STAGE4a_comparison.md for the full table.
