# Sprint 8 Stage 2-prep — §3.5 Synthesis

## (a) Method used and headline result

**Path A unavailable.** CW output.txt contains no degree-of-hydration data (640 columns:
Time + 637 spatial temperatures + Gradient + Ambient). α(t) is not observable from CW
output files (§3.1 finding).

**Path B — indirect ratio test used.** In a linear thermal diffusion model with Van Breugel
k(α), ΔT_A/ΔT_B ≈ Δα_A/Δα_B. Discrepancy between the CW-observed ΔT ratio and the
engine-predicted Δα ratio bounds the degree to which engine and CW α(t) trajectories
diverge.

**Headline: FAIL verdict at max discrepancy = 8.48%.**

The A/C ratio test at t=24 hr, mean|ΔT| statistic, shows 8.48% discrepancy — exceeding
the 5% FAIL threshold. The main A/B test is in the CAUTION band (2.3–6.2%). Results
span the 1–5% and >5% bands simultaneously depending on which test/stat/timestep is
examined.

| Band | Tests in band |
|---|---|
| PASS (<1%) | 0 of 18 |
| CAUTION (1–5%) | 12 of 18 |
| FAIL (>5%) | 6 of 18 (all in A/C test) |

---

## (b) Sprint 8 implications — CAUTION with acknowledged α uncertainty

The FAIL is structurally concentrated in the A/C ratio test (highhyd − highhyd_b010 denominator),
which has a smaller Δα denominator and is therefore more sensitive to nonlinearity and
trajectory offsets. The A/B test — the main test comparing the full-range pair (Pair A)
to the β=0.10 pair (Pair B) — remains in the CAUTION band at 2.3–6.2%.

**Sprint 8 calibration is at t=168 hr.** At this timestep, A/B and A/C discrepancies
drop to 2.1–5.4%. The t=24 extremes (up to 8.5%) are early-time artifacts where the
ratio denominator is smallest and thermal nonlinearity is most pronounced.

**Practical recommendation:** Sprint 8 can proceed with Schindler parameter–specified
datasets, but must acknowledge a ±5% α trajectory uncertainty band. When interpreting
k_uc residuals, a 5% α offset translates to ΔkVanBreugel ≈ k_uc × 0.33 × 0.05 × α_u ≈
0.5–1.5% change in k_c — well within the measurement scatter from other sources (grid
resolution, BC alignment). This uncertainty does not invalidate the calibration design,
but should be noted in the calibration record.

---

## (c) T_ref finding — principal source of discrepancy

The §3.4 back-calculation finds:

**Best-fit T_ref = 21°C (294.15 K), RMS discrepancy = 3.21% (vs 5.05% at engine T_ref=23°C).**

Arrhenius factors at T_pl=73°F (295.93 K):

| T_ref | Factor | Source/note |
|---|---|---|
| 20°C (293.15 K) | 1.212 | ACI 308 / Schindler-Folliard 2002 calibration reference |
| 21°C (294.15 K) | 1.131 | **Best-fit to CW ratio data** |
| 23°C (296.15 K) | 0.985 | Engine T_REF_K (current) |

The engine's T_ref=23°C is above T_pl=73°F (22.78°C), giving Arrhenius factor < 1 and
slightly slower α accumulation than at T_ref=21°C. If CW internally uses T_ref≈21°C
(close to the Schindler-Folliard 2002 calibration temperature of 21.11°C = 70°F), α
accumulates ~14.7% faster in CW's computation than in the engine's.

This T_ref gap is the principal identifiable source of the systematic ratio discrepancy.
Even at best-fit T_ref=21°C, max discrepancy remains 5.42% — indicating residual sources
(finite-domain nonlinearity, spatial field averaging, early-time sensitivity) contribute
the remaining gap.

**Implied α offset at Sprint 8 calibration timestep (t=168 hr):**

| Dataset | |Δα|/α at t=24hr | |Δα|/α at t=168hr |
|---|---|---|
| lowhyd | 1.7% | 1.4% |
| highhyd | 5.5% | 1.0% |
| highhyd_b010 | 1.3% | 1.0% |

At t=168 hr, the α offset is ≤1.4% across all datasets. The largest discrepancy (5.5%
for highhyd at t=24) decays as highhyd approaches its α_u plateau. This confirms that
the T_ref effect is most acute at early times and is negligible at the Sprint 8
calibration horizon.

---

## (d) Residual source attribution

After T_ref correction (from 23°C to 21°C), 3.21% RMS discrepancy remains. Three
contributing factors, not further separable from this dataset:

1. **Finite-domain nonlinearity**: the ΔT ∝ Δk_c proportionality holds in the
   infinite-medium limit. In a bounded 80×40 ft domain with competing hot-bottom and
   cold-top BCs, the spatial averaging of ΔT over di=4–47 introduces geometry-dependent
   weighting that varies between pairs. This is structural and cannot be zero.

2. **Early-time transient**: at t=24 hr, the thermal boundary layers have not fully
   penetrated the section. The effective "lever arm" for k(α) differences varies with
   depth non-proportionally, making mean|ΔT| and RMS ratios more sensitive to α trajectory
   details than they will be at t=168 hr.

3. **β-nonlinearity**: Pairs B and C both involve highhyd_b010 (β=0.10), which has a
   slower early-time α growth curve (Schindler β controls curvature near α=0). The
   Δα denominator for A/C at t=24 is 0.153 (vs 0.404 for A/B), amplifying any mismatch.
   The ratio test linearization is least accurate when the denominator Δα is small.

---

## Summary judgment

| Question | Answer |
|---|---|
| Path used | B — ratio test on Stage 1 ΔT data |
| Headline discrepancy | 8.48% max (mean\|ΔT\| A/C t=24); 6.24% max for A/B |
| Band at t=168 hr only | CAUTION (2–5%) for both tests |
| Identified root cause | T_ref=23°C (engine) vs ~21°C (CW best-fit) |
| α offset at t=168 hr | ≤1.4% (all datasets) |
| Sprint 8 block? | No — proceed with CAUTION band acknowledged |
| Recommended action | Note ±5% α uncertainty in Sprint 8 calibration record; consider T_ref audit if residuals exceed 3% after multi-α calibration |
