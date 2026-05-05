# compute_hu_factor — Sprint 9 Reference Note

**Date:** 2026-05-04
**Sprint:** 9 (written at sprint close)
**Scope:** Usage, values, bug history, and open questions for Sprint 10.

---

## 1. What It Does

`compute_hu_factor` returns a dimensionless multiplicative scaling factor for the heat-of-hydration
parameter Hu_raw read from CW's `input.dat`. The factor represents the cementitious mixture's
heat-release capacity relative to a pure-cement baseline, weighted by the composition of the mix
(Type I/II cement, fly ash, slag proportions). The output Hu_factor is then applied in the wrapper:

```
Hu_factored = Hu_raw × Hu_factor
```

before the Sprint 8 B1 inverse-compensation step (`Hu_eff = Hu_factored / c`).

The physical interpretation: CW's Hu parameter reflects the total heat per unit mass of
cementitious material. `compute_hu_factor` adjusts for differences between the mix's actual
cementitious blend and the composition assumed by the baseline engine calibration.

---

## 2. Where It Is Called in the Pipeline

The function sits in the Sprint 9 wrapper layer, applied before the engine run:

```python
Hu_factor    = compute_hu_factor(mix)        # Sprint 9 addition
Hu_factored  = Hu_raw * Hu_factor            # multiply before B1 step
c            = ...                           # Sprint 8 c(T_pl) quadratic
Hu_eff       = Hu_factored / c               # B1 inverse-compensation (preserves Q(t))
alpha_u_eff  = c * alpha_u_raw
mix.Hu_J_kg_effective = Hu_eff
mix.alpha_u           = alpha_u_eff
```

The engine (`thermal_engine_2d.py`) reads `mix.Hu_J_kg_effective` at line 1510 and
`mix.alpha_u` at line 1513. It has no knowledge of Hu_factor — the scaling is resolved
entirely in the wrapper before the solve.

The Hu_residual conditional (Sprint 8, Hu < 10,000 → 12,937 J/kg) is logically prior to
Hu_factor; for production mixes with realistic Hu >> 10,000 J/kg, the conditional is a
passthrough and the ordering has no effect in practice.

---

## 3. Observed Values for Pilot Mixes

From `PILOT_REPORT.md` §1.2 (post-fix, Sprint 9):

| Mix | Composition | Hu_raw (J/kg) | Hu_factor | Hu_factored (J/kg) | Hu_eff (J/kg) |
|---|---|---|---|---|---|
| MIX-01 | Type I/II + 22% FA + 17% slag | 424,143 | 0.951403 | 403,531 | 383,178 |
| MIX-07 | Type I/II + 22% FA + 70% slag | 463,076 | 0.894603 | 414,269 | 393,375 |

(`Hu_eff = Hu_factored / c(73°F) = Hu_factored / 1.0532`)

MIX-07's lower Hu_factor (0.8946 vs MIX-01's 0.9514) reflects its higher slag content:
slag contributes less heat per kg than Type I/II cement, so the blended-mix Hu_factor
falls below unity more steeply as slag fraction rises.

---

## 4. Sprint 9 Bug History

**Symptom.** The first pilot run (before the fix) showed MIX-07 with a ~12.5°F core
temperature overshoot relative to CW (engine warmer than CW at di=24, t=168). MIX-01
missed the gate by ~0.5°F at di=47, which is consistent with the structural stencil
mechanism. MIX-07's miss was far outside the expected envelope and spatial diagnosis
pointed to excess bulk heat — a systematic Hu scaling error, not a stencil artifact.

**Root cause.** `compute_hu_factor` used an incorrect denominator in its composition-weighted
Hu calculation. The effect was to inflate Hu_factored for high-slag mixes (where the
denominator error is larger), producing the large core-temperature overshoot in MIX-07.
MIX-01 has a lower slag fraction and was less affected.

**Fix applied mid-sprint.** After correcting the denominator, the output values became
0.951403 (MIX-01) and 0.894603 (MIX-07). MIX-07 core temperature aligned with CW and
max|R| dropped to 1.5843°F — the expected stencil-dominated residual profile, with
R_di36 = +0.7598°F (bulk) and R_di47 = −1.5843°F (stencil dip).

**Implication.** The denominator fix is the only change made to `compute_hu_factor` during
Sprint 9. The corrected function is the version used in all sprint closure artifacts
(PILOT_REPORT.md, SWEEP_V2_REPORT.md).

---

## 5. Open Questions for Sprint 10

**(a) Generalization beyond the pilot kinetic envelope.**
MIX-01 has 17% slag; MIX-07 has 70% slag. These bracket the pilot, but Phase A will introduce
MIX-04, MIX-06, and MIX-13, which may have different FA/slag proportions. The function's
denominator logic should be verified against those compositions before running Phase A sweeps.

**(b) Source location — wrapper vs. engine.**
`compute_hu_factor` currently lives in the Sprint 9 wrapper. If the Hu_factor correction is
considered a permanent part of the calibration (not a run-specific override), it may belong in
the engine source alongside `k_uc × 0.96`. No decision has been made; keeping it wrapper-side
is appropriate until the Phase A audit is complete.

**(c) Specification — production-table coupling.**
The function reads from `mix` object fields that derive from CW's `input.dat`. The exact
mapping from CW's per-ingredient weight fractions to the function's composition arguments
should be documented explicitly. The Sprint 9 pilot only exercised two compositions; edge
cases (e.g., ternary blends with multiple SCMs, or mixes with silica fume) are untested.

**(d) Silica fume compositions.**
Does `compute_hu_factor` handle silica fume content correctly? Sprint 9 did not test SF-bearing
mixes. SF's pozzolanic activity and heat release timing differ from cement, FA, and slag, and
the function's formula may not account for SF-specific heat-fraction behavior. Validation
against MIX-13 (and any other SF mix in the production set) is deferred until formal
specification establishes the formula and its expected SF behavior.

**(e) Interaction with c(T_pl) extrapolation.**
The Sprint 8 c(T_pl) quadratic was calibrated at T_pl ∈ [40, 110]°F using two SCM compositions
(not matched to the production mix set). The cross-product of composition-range uncertainty in
compute_hu_factor and T_pl-range limits in c() has not been characterized. This is a second-order
concern for the pilot but relevant as the validation campaign broadens.

---

## 6. Recommendation

Before extending to Phase A's 5-mix scope:

1. **Audit the denominator logic** against each Phase A mix's production-table composition.
   Confirm the output Hu_factor values are physically plausible (factor should be between the
   heat-fraction contributions of the lightest and heaviest SCM blends in the set).

2. **Lock the specification** — document explicitly which fields from `input.dat` feed the
   function and how the weighting formula is derived from those fields, so the function
   can be reviewed independently of the wrapper code.

3. **Spot-check with a third mix** during Sprint 10's stencil-fix validation pass:
   run one additional production mix (e.g., MIX-04 if available) and confirm the Hu_factor
   output is consistent with expectation before the full 25-run Phase A sweep.
