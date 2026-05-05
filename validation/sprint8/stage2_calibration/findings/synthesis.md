# Sprint 8 Stage 2 — Synthesis

## Data summary

| §  | Finding |
|----|---------|
| 3.1 | All 12 new CW datasets copied and present (7.51–7.69 MB each) ✓ |
| 3.2 | All 12 consistent: τ=5hr, β=0.85, Hu=1.0 J/kg; geometry 40×80ft; α_u matches folder name ✓ |
| 3.3 | α(t) at t=168hr reaches 95.0–97.8% of α_u (73°F datasets at 95.0%, 100°F datasets at 97.8%) |
| 3.4 | Baseline residuals at factor=0.96: 0.65–0.89°F (α=0.20), 1.34–1.79°F (α=0.40), 2.02–2.69°F (α=0.60), 2.72–3.59°F (α=0.80) — ALL FAIL Structure C gate |
| 3.5 | k_uc sweep sensitivity near-zero: 0.37°F variation across 12% k change at α=0.20; only 0.06°F variation at α=0.80; optimum at left sweep edge (0.92) for all 4 targets — gate NOT achievable |
| 3.6 | Spread across 5 α points = 0.040 → formally Outcome A-marginal, but see caveat below |
| 3.8 | Final validation at factor=0.96: 9/21 pass (all 9 Sprint 7 datasets), all 12 new datasets fail |
| 3.9 | Sprint 7 regression: ALL 9 PASS ✓ (R1 max=0.229°F, R2 max=0.219°F) |

## Synthesis paragraph

**(a) What is the optimal k_uc calibration?**

Neither Outcome A (constant factor) nor Outcome B (α-dependent function) correctly describes the result. The k_uc sweep reveals that the k_uc factor has near-zero sensitivity to the observed residuals at all four new α targets: moving the factor across the full {0.92..1.04} range changes the worst-case minimax by only 0.37°F at α=0.20, 0.18°F at α=0.40, 0.10°F at α=0.60, and 0.06°F at α=0.80. At α=0.80, a 12% change in k_uc changes the gate metric by less than 0.1°F — the knob is effectively inert. The "optimal factor" of 0.92 reported by the sweep decision script is an artifact of the sweep reaching its left edge, not a physically meaningful calibration point. No k_uc factor in any reasonable range can reduce the residuals to ≤0.35°F at the new α targets.

**(b) Is 0.96 the right constant factor?**

The Sprint 7 calibration (factor = 0.96) does not generalize to higher α values. Residuals grow approximately linearly with α_u: 0.65–0.89°F at α=0.20, growing to 2.72–3.59°F at α=0.80 — roughly 4× the gate threshold. The factor spread across 5 α points (0.036 to 0.80) is formally 0.040 (Outcome A-marginal), but this spread is numerically meaningless because the gate is far from being achieved at any factor.

**(c) What is the dominant source of residual?**

The residuals are driven by the `make_neutral_env` approximation, not by k_uc. The CW simulations used real Austin TX July 2026 weather (diurnal T_air peaking ~100°F, averaging ~80°F over 168hr). The engine uses a flat 73°F neutral environment. At α=0.036, the concrete conductivity is near its maximum (k ≈ 1.32 × k_uc), and the concrete fully equilibrates to its 73°F soil boundaries by t=168hr regardless of surface weather — the weather approximation has negligible impact. At higher α_u, the concrete conductivity drops (k(α_u=0.80) ≈ 1.08 × k_uc), the thermal transient at t=168hr is larger, and the history of boundary conditions matters. The Austin summer weather heats the concrete interior, particularly at the side and top faces. The engine (flat 73°F) does not replicate this heating. Inspection of the residual field at di=24 confirms: the CW temperature at the outer-width cells (wi=0..6, the side-face region) is systematically 0.65–2.72°F above what the engine computes, while the centerline (wi=12) shows near-zero residual. This offset pattern is consistent with weather-driven heating from the surface and cannot be reduced by adjusting k_uc, because the engine's Dirichlet soil BCs pin the outer-face temperature at exactly T_soil while CW's outer-face temperature floats above T_soil due to summer heat.

The residual grows linearly with α_u at approximately +3.4°F per unit of α_u. This linear scaling is consistent with the Van Breugel relationship: Δk(α_u) = k_uc × 0.33 × α_u, meaning the drop in conductivity from α=0 to α_u is proportional to α_u. The weather effect that survives to t=168hr in the concrete interior is therefore proportional to how insulating the concrete has become — i.e., to α_u.

**(d) Do all 21 data points pass Structure C gate after calibration?**

No. 9 pass (the 9 Sprint 7 datasets at α≈0.036), 12 fail (all 12 new datasets at α=0.20 to 0.80). The failures are not marginal: the smallest failure margin is 0.30°F above the gate (alpha02_A: R1=0.652°F vs gate 0.35°F), and the largest is 3.24°F above the gate (alpha08_I: R1=3.59°F). No k_uc factor can bridge these margins.

**(e) Does Sprint 7's 9-run gate still pass?**

Yes. All 9 Sprint 7 datasets pass Structure C at factor=0.96 with substantial margin: R1 max=0.229°F (vs gate 0.35°F), R2 max=0.219°F. Sprint 7 is not regressed. The make_neutral_env approximation works at α≈0.036 because the concrete's near-maximum conductivity ensures the interior equilibrates to the soil boundaries before the 168hr measurement, masking the weather difference.

## Implication

Sprint 8 Stage 2 calibration objective — close Structure C gate at α_u ∈ {0.20, 0.40, 0.60, 0.80} via k_uc calibration — cannot be achieved with the current `make_neutral_env` comparison methodology. The calibrated factor at each α target is NOT a property of the concrete k_uc; it is a property of the mismatch between the CW weather environment and the engine's flat-temperature approximation. Calibrating k_uc against these datasets would produce a k_uc that compensates for a weather effect, not a material conductivity mismatch.

The user should decide whether to:
- Accept that the Sprint 7 k_uc = 0.96 calibration is valid only at α≈0.036 and the validity range is narrower than originally assumed
- Rethink the Stage 2 dataset design to use weather-neutral CW runs (set CW location to a constant T_air environment rather than Austin TX July) or to compare at a different reference time (e.g., t=24hr when weather effects have not yet penetrated into the interior)
- Or accept the make_neutral_env approximation and widen the gate tolerance for higher-α scenarios

This synthesis does not propose which path to take.
