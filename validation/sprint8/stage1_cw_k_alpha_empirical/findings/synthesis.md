# Sprint 8 Stage 1 — Synthesis

## (a) Outcome verdict

**Outcome 1: CW empirically uses Van Breugel α-dependent thermal conductivity.** The data is decisive. Pair A (highhyd − lowhyd), which represents the largest possible α(t) spread under these kinetic parameters, shows max|ΔT| = **1.46°F at t=24 hr**, rising to **2.02°F at t=84 hr**, and **2.16°F at t=168 hr** over the masked domain (di=4–47, all wi). Pair B (highhyd_β=0.10 − lowhyd) yields max|ΔT| = 1.03°F at t=168 hr, also well above the Outcome 2 threshold. Both pairs grow monotonically with time, consistent with k(α) effects accumulating as α diverges across the run. An Outcome 2 signature (constant k) would require all pairs to remain below ~0.1°F; the observed values are 10–20× larger. Sprint 8's calibration design path — multi-α plateau datasets — is confirmed.

## (b) Predicted α and k_c reduction: does Pair A magnitude align?

From the engine isothermal trajectory (which models what CW ran, since Hu=1 J/kg), at **t=168 hr**:
- **lowhyd**: α = 0.0361
- **highhyd**: α = 0.6384

Applying Van Breugel k(α) = k_uc · (1.33 − 0.33·α):
- k_c(lowhyd)  = k_uc · 1.3181
- k_c(highhyd) = k_uc · 1.1193
- Relative reduction: (1.3181 − 1.1193) / 1.3181 = **~15.1%**

A 15% reduction in thermal conductivity over an 80-ft (24.4 m) deep mat being heated from below by 100°F soil (27°F above 73°F placement temperature) is fully consistent with a 2+ °F steady-state temperature field difference. The direction and order of magnitude of Pair A (highhyd warmer than lowhyd near the heated boundary) match the expected outcome: lower k_c in the high-hydration run impedes heat transport from the soil BC, so interior cells near the bottom boundary accumulate slightly more heat. The argmax location (di=8, wi=0 at t=168 — near the top, on the centerline) suggests that after 168 hr the largest T difference is where the competing boundary conditions (cold top, hot bottom) have had the most time to interact through reduced-conductivity concrete. Quantitative alignment is plausible; a detailed model-matching calculation is beyond Stage 1 scope.

## (c) Pair C: what β tells us about CW's α(t) trajectory

Pair C (highhyd − highhyd_β=0.10, same α_u=0.70 and τ=10 hr but β=0.85 vs β=0.10) shows max|ΔT| = **0.59°F at t=24 hr** rising to **1.13°F at t=168 hr**. The engine's isothermal trajectories confirm the mechanism: despite identical α_u and τ, β=0.85 reaches α=0.639 by t=168 while β=0.10 reaches only α=0.329. The high-β run hydrates faster and more completely in the first 24 hr (β controls the curvature of the Schindler function near α=0), so k_c drops earlier and more steeply for highhyd than for highhyd_β=0.10. The 1.13°F difference between these two runs — which share the same τ and α_u — quantifies β's independent contribution to CW's thermal field through k(α). Sprint 8 multi-α calibration will need to distinguish τ-driven from β-driven α(t) trajectories; Pair C confirms that this distinction matters (~1°F at t=168 hr, 40–50% of Pair A's signal).

## (d) Consistency with documented Van Breugel formula

The empirical pattern is entirely consistent with Stage 0's finding that CW V3 uses Van Breugel (Eq. 23 of the CW V3 manual). Three quantitative checks:

1. **Direction**: highhyd runs are warmer than lowhyd in all pairs and at all timesteps. Higher α → lower k_c → slower heat diffusion → warmer interior. This is the expected Van Breugel direction.

2. **Magnitude ordering**: Pair A > Pair B at all timesteps. Pair A has Δα(168) = 0.602 while Pair B has Δα(168) = 0.293 (roughly half). Pair A max|ΔT| = 2.16°F, Pair B = 1.03°F — approximately a 2:1 ratio, matching the ~2:1 ratio in Δα. This linear scaling is exactly what Van Breugel predicts: Δk_c = k_uc · 0.33 · Δα, so doubling Δα doubles Δk_c and (approximately) doubles ΔT.

3. **Temporal growth**: ΔT grows over time in all pairs. Under isothermal kinetics, α increases monotonically toward α_u. The growing α divergence between high- and low-hydration runs means k_c diverges continuously through the run, producing the observed monotonically increasing pairwise ΔT. This is inconsistent with constant-k behavior (which would yield near-zero ΔT at all times) and consistent with Van Breugel.

No evidence of CW deviating from the documented formula is found. The empirical behavior is consistent with a Van Breugel implementation applied to isothermal (Hu=1 J/kg) kinetics over the 168-hour analysis window.
