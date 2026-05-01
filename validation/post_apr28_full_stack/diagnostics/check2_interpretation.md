# CHECK 2.3 — Interpretation

**Outcome: 2B — the engine sheds significantly MORE heat than CW does going from adiabatic to full-stack mode. BC over-loss is confirmed.**

## Aggregate statistics

Across the 14 mixes (skip MIX-13):

- **Mean (engine_drop − cw_drop) = +3.95°F**  (SD 0.95°F)
- **Mean (drop_diff as % of adiabatic rise) = +5.33%**  (SD 1.42%)
- Range: +2.19°F (MIX-14, 75°F warm placement) to +6.34°F (MIX-15, 45°F cold placement)
- All 14 mixes have positive drop_diff — the direction is uniform.

The full per-mix table is in `check2_cw_modes.md`. The adiab-residual column (engine_adiab − CW_adiab at the adiabatic T0 of the CSV) ranges from -0.01°F to +1.34°F with a mean near +0.4°F, confirming the apr28 adiabatic validation holds and that the extrapolation `CW_adiab(fs T0) ≈ engine_adiab(fs T0)` used to estimate `cw_drop` is well-grounded (the extrapolation error is bounded by ~1°F and is much smaller than the +3.95°F mean drop_diff).

## Why this rules out the alternatives

**Alternative B (CW mode mismatch) — not supported.** If CW behaved very differently between adiabatic and full-stack (different α_u, different temperature dependence, different effective Hu), we would expect to see large unexplained variation across mixes — composition-dependent or temperature-dependent inconsistencies in the drop comparison. Instead, the drop_diff is tightly clustered (SD 0.95°F) for mid-range placements (60°F mixes range from +3.09°F MIX-04 to +4.88°F MIX-07). The cluster looks like a uniform BC-side phenomenon, not a CW-internal modeling discrepancy. The cold-placement outlier (MIX-15) does suggest the cold regime has additional CW behavior we don't capture, but that's a secondary effect on top of the dominant uniform BC over-loss.

**Alternative A (wiring bug)** was ruled out by CHECK 1.

**The remaining candidates are**:
- BC over-loss in the engine (Hypothesis: engine's top-surface convection / form R-value / soil conduction collectively remove more heat than CW's equivalents).
- A T0-dependent kinetics term (Alternative C). CHECK 3 will discriminate by examining the time-domain residual.

## Detailed observations

**MIX-15 stands out dramatically.** The CW full-stack peak (117.12°F) is essentially identical to the CW adiabatic peak (117.05°F) — CW's full-stack model at 45°F cold placement says the BC layer removes essentially zero heat. The engine, at the same 45°F placement, drops by 6.27°F. This is the largest engine-vs-CW BC mismatch in the library and is consistent with the post-apr28 REPORT.md finding that MIX-15 has the worst residuals. The cold-placement regime appears to be a place where CW and the engine model BC heat loss very differently — a candidate for soil model / cold-air convection investigation in Sprint 7 follow-up.

**MIX-14 (75°F warm placement) has the smallest drop_diff** (+2.19°F). At warmer placement, both CW and engine have larger absolute drops (engine 10.36°F, CW 8.16°F) but they're closer to each other. The warmer concrete loses more heat to BCs in absolute terms, but the engine-CW gap is narrower. This is consistent with BC heat loss scaling roughly with (T_concrete − T_ambient): at 75°F placement the concrete-ambient gradient is similar in early hours, while at 60°F the placement-ambient gradient is smaller and the engine's persistent over-loss becomes more visible relative to the smaller absolute drop.

**The mid-range mixes (60°F placement, all SCM levels) cluster between +3.1°F and +4.9°F drop_diff.** No clear monotonic relationship with SCM percentage. This suggests the BC over-loss is primarily geometric / climate-driven, not composition-driven. The composition does not appear to mediate the engine's BC over-loss in any meaningful way at fixed T0.

## What this implies for Sprint 7

The +3.95°F mean drop_diff IS the regression. Pre-apr28, the engine's over-estimated kinetics (Hu factor ~1.0 instead of ~0.95) generated extra heat that masked this BC over-loss; the engine and CW full-stack peaks agreed on reference mixes by accidental cancellation. With the apr28 calibration removing the kinetics over-supply, the BC over-loss is now visible as a uniform ~3.5°F downward shift in engine PeakMax.

The Sprint 7 primary target should be the BC layer: identify which BC term (top-surface convection, top-surface LW radiation, form-side R-value, or soil-bottom conduction) is over-removing heat. The fact that MIX-14 (warm placement) and MIX-15 (cold placement) show different drop_diff magnitudes (+2.19°F vs +6.34°F) suggests the over-loss is temperature-dependent: the term most likely to scale this way is convective heat transfer (q_conv ∝ T − T_amb), pointing at top-surface convection and/or form-face convection coefficients.

**Proceed to CHECK 3** to confirm the residual signature in the time domain.
