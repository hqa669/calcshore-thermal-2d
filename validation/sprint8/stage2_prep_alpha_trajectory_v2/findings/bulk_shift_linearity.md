# Sprint 8 Stage 2-prep v2 — §D.1 Bulk-Shift Linearity Sanity

**Gate:** CV(ΔT_bulk / α) < 5% at all 5 timesteps confirms the linear inversion model ΔT_bulk = C(t)·α is valid for Path D.

| t (hr) | R_lowhyd | R_highhyd | R_highhyd_b010 | Mean C | Std | CV% | Gate |
|---|---|---|---|---|---|---|---|
| 24 | 3.187 | 3.402 | 3.447 | 3.345 | 0.113 | 3.4% | **PASS** |
| 48 | 3.500 | 3.464 | 3.483 | 3.482 | 0.015 | 0.4% | **PASS** |
| 84 | 3.289 | 3.523 | 3.441 | 3.418 | 0.097 | 2.8% | **PASS** |
| 120 | 3.167 | 3.533 | 3.457 | 3.386 | 0.158 | 4.7% | **PASS** |
| 168 | 3.552 | 3.551 | 3.479 | 3.527 | 0.034 | 1.0% | **PASS** |

**Overall gate: PASS — linearity assumption valid**