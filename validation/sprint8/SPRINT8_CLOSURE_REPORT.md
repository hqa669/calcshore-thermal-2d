# Sprint 8 Closure Report — k(α) Curve Calibration

**Date:** 2026-05-03
**Sprint:** 8 (k(α) curve calibration across production hydration range)
**Status:** CLOSED ✓
**Predecessor:** Sprint 7 (soil-concrete BC calibration at α≈0.036, closed 2026-05-01)

---

## 1. Sprint 8 Goal

Determine whether Sprint 7's `k_uc × 0.96` calibration (validated at α ≈ 0.036) holds across the
production-relevant hydration range (α = 0.04 to α ≈ 0.7), or whether the optimal calibration
varies with α. Sprint 7 only validated at the very tip of the production hydration curve;
Sprint 8 covers the rest.

---

## 2. Final Result

Sprint 8 closes with two empirical corrections layered on top of the Sprint 7 calibration:

| Component | Value | Origin |
|---|---|---|
| `k_uc × 0.96` | unchanged | Sprint 7 calibration, carried forward |
| `Hu_residual` | 12,937 J/kg | CW solver minimum-Hu floor compensation |
| `c(T_pl)` | quadratic in T_pl (see §3) | Engine vs CW kinetics shape correction |

Applied together, these close the Sprint 7 Structure C gate (R1, R2 ≤ 0.35°F) for:
- All 32 A-scenario datasets (T_pl = T_soil) at T_pl ∈ {40, 50, ..., 110}°F, α_u ∈ {0.20, 0.40, 0.60, 0.80}
- All Sprint 7 datasets at α_u = 0.10 (no regression)
- 10 of the 12 Sprint 8 multi-α datasets (A, F scenarios across all α_u; I scenarios at α_u ≤ 0.40)

The 2 remaining datasets (I scenario at α_u ∈ {0.60, 0.80}) fail the gate. Spatial diagnostic
isolates the cause as a pre-existing bottom-side BC stencil asymmetry, deferred to Sprint 9.

---

## 3. Calibration Components

### 3.1 Hu_residual = 12,937 J/kg

CW substitutes a hidden minimum heat-release floor when user-set Hu falls below approximately
3,000–10,000 J/kg. Below the floor, CW behaves as if Hu were the floor value; above, CW honors
the input linearly. The engine has no such floor.

The residual was characterized via 14 CW runs at Hu ∈ {1, 10, 100, 1000, 10000, 50000, 100000} J/kg
× α_u ∈ {0.20, 0.80} in the A scenario. Linear fit `ΔT_core = a × Hu + b` per α_u gave:

- α_u=0.20: a = 2.656×10⁻⁵ °C/(J/kg), b = 0.36361 °C → Hu_residual = b/a = 13,691 J/kg
- α_u=0.80: a = 1.106×10⁻⁴ °C/(J/kg), b = 1.50328 °C → Hu_residual = b/a = 13,592 J/kg

Geometric mean: Hu_residual = 13,641.5 J/kg (initial estimate).

Self-calibrated to **12,937 J/kg** via probe solve at the c=1 anchor point (T_pc=70°F, α_u=0.20)
during the §4.11 sweep. The 5% difference reflects refinement to the actual A-scenario behavior;
this is the value used in the closing calibration.

Applied as: `mix.Hu_J_kg_effective = 12,937` in the engine wrapper when CW comparison is the
target. This is wrapper-only; the engine source kinetics path remains unchanged.

### 3.2 c(T_pl) = α_u correction for kinetics shape disagreement

Engine and CW Schindler-style kinetics agree well near T = 70°F but diverge at temperature
extremes — engine is "faster" at cold T, "slower" at hot T relative to CW. The disagreement is
α_u-independent (cross-α_u slope std/mean = 0.68%) and well-modeled as a multiplicative correction
on engine α_u:

```
α_u_engine_corrected = c(T_pl) × α_u_nominal
```

The c factor was characterized via 32 A-scenario CW datasets at
T_pc ∈ {40, 50, 60, 70, 80, 90, 100, 110} °F × α_u ∈ {0.20, 0.40, 0.60, 0.80}. For each (T_pc, α_u)
point, c was tuned to minimize L2 norm of the engine vs CW core temperature trajectory over
t ∈ [0, 168] hr. After dense local refinement (0.01 step in c, ±0.10 around bracket), c values
at all 32 points are interior optima with ±0.005 resolution.

Closed-form fit to α_u-averaged c values:

```
c(T_pl) = −1.3025 + 0.04746 × T_pl − 2.081×10⁻⁴ × T_pl²
```

(T_pl in °F; c dimensionless)

| Fit quality | Value |
|---|---|
| R² | 0.99993 |
| max\|residual\| | 0.0055 |
| Anchored value at T_pl=70°F | 1.000 |
| c(73°F) | 1.0530 |
| c(100°F) | 1.3623 |
| c(110°F) | 1.4008 |

**Important caveat — fit validity range.** This quadratic is fit to data in T_pl ∈ [40, 110] °F.
Extrapolation outside this range is not validated. Production runs at placement temperatures below
40°F or above 110°F should not use this correction without re-extending the calibration sweep.
Future work (deferred to Sprint 9 or beyond) may extend the fitting range to capture broader
operating conditions and potentially replace the quadratic form with a physically motivated
saturating function (e.g., Arrhenius difference) that extrapolates more safely.

The quadratic form is empirical and has no direct physical meaning — it captures the curvature in
the calibrated range but is not a model of the underlying kinetics disagreement. The shape is
consistent with a two-Arrhenius-integral ratio with slightly different reference temperatures, but
that mapping has not been verified.

Applied as: `α_u_engine = c(T_pl) × α_u` in the engine wrapper for CW comparison runs. Wrapper-only;
engine source kinetics path remains unchanged.

### 3.3 Note on Sprint 7 calibration

Sprint 7's `k_uc × 0.96` was derived without the floor or kinetics corrections of this sprint. A
back-decomposition shows the 0.04 reduction had two contributing sources:

- CW floor heat at α(t=168) ≈ 0.036 produced ~0.13°F bulk warming (Sprint 7 §"CW bulk noise floor").
  This was empirically absorbed into the Sprint 7 k_uc factor.
- Real k_uc material disagreement between engine and CW physics (residual after floor accounted for).

The exact split was not re-derived in Sprint 8 — `k_uc × 0.96` is carried forward unchanged as the
calibrated value, with the understanding that it now operates in conjunction with the Hu_residual
floor compensation. Sprint 7's gate at α_u=0.10 still passes with the new corrections layered on
top (verified during §4.11 closure).

---

## 4. Path Traversed

**Stage 0 — CW formula verification.** Confirmed CW V3 Equation 23: `k_c(α) = k_uc × (1.33 − 0.33α)`,
identical to engine. No reverse engineering needed.

**Stage 1 — CW α-dependence empirical confirmation.** 3 CW datasets at varying hydration parameters
confirmed CW responds to α with magnitude matching Van Breugel predictions. Pair-A:B ratio = 2.10
vs Δα ratio 2.06. CW empirically uses Van Breugel.

**Stage 2-prep — Engine vs CW α(t) trajectory.** 5–8% disagreement at intermediate times, narrowing
to ~1% at t=168. Curve-fit suggested T_ref disagreement (engine ~23°C vs CW ~21°C). Decision:
don't change T_ref; absorb into calibration. T_ref alignment deferred. (Later confirmed in §4.10
follow-up that 2°C T_ref shift produces only ~0.005°F residual change — a real but second-order
effect, not the dominant mechanism.)

**Stage 2 — Multi-α calibration (initial attempt, BLOCKED).** 12 CW datasets at
α_u ∈ {0.20, 0.40, 0.60, 0.80} × 3 scenarios (A, F, I) at Hu=1. Engine k_uc sweep across
{0.92, ..., 1.04} produced near-zero leverage (variation 0.06–0.37°F vs gate miss 0.30–3.24°F).
Residuals scaled linearly with α_u at ~4.4°F per unit α. Initial agent attributed this to weather
mismatch (engine `make_neutral_env(73°F)` vs CW Austin TX). Hypothesis was incorrect — ambient was
already overwritten to T_soil in all CW inputs.

**Stage 2-diag-pre — Spatial diagnostic.** Centerline T(z) profile plot at four α_u in F scenario
showed bulk-uniform residual of ~0.7–2.5°F (CW warmer than engine), flat across di=5..40, scaling
linearly with α_u. Pattern inconsistent with BC mismatch (would peak at face) or k(α) shape error
(would bend the profile); consistent with bulk-distributed heat source the engine isn't applying.

**Stage 2-diag-Hu — Hu floor characterization.** 14 CW runs at Hu ∈ {0.01, 1, 10, ..., 100000} J/kg
× α_u ∈ {0.20, 0.80}, A scenario. Results showed:
- CW has a hard floor at Hu ≈ 3,000–10,000 J/kg below which input is ignored.
- Engine has no floor (clean linear response from Hu=0.01 to Hu=100,000).
- Disagreement = CW floor magnitude, near-constant in Hu, scaling linearly with α_u.
- F-scenario cross-check reproduced the bulk gap from Stage 2-diag-pre within 0.2°F.

**Stage 2-floor-test — Apply Hu_residual.** Engine reruns at Hu=13,641.5 J/kg on the original 12
Sprint 8 datasets reduced residuals by 74–99%. 9/12 datasets passed the gate. Remaining failures:
all 4 I-scenarios at α_u ≥ 0.40, residual ≈ 0.20 + 0.205 × α_u °F. T_ref shift to 21°C tested
and rejected (effect ~100× too small).

**Stage 2-alpha-u-T-trend — α_u(T_pl) characterization.** 32 A-scenario CW datasets at 8 T_pl × 4 α_u.
Per-point α_u tuning via L2 trajectory matching across t ∈ [0, 168] hr. Initial sweep at c ∈ [0.80, 1.20]
clamped at edges for 6 of 8 T_pl values. Extended sweep at c ∈ [0.40, 1.80] resolved 7 of 8 T_pl
interior; T_pc=40°F still clamped at 0.40. Final extension at c ∈ [0.20, 0.60] resolved T_pc=40°F
to c ≈ 0.26. Dense local refinement at 0.01 step around all 32 optima confirmed precision.

**Stage 2 §4.11.8 closure — Fit and validate.** Quadratic fit to α_u-averaged c values produced
R² = 0.99993, max residual 0.0055. Applied as α_u_corrected = c(T_pl) × α_u_nominal to the
original 12 Sprint 8 datasets. Result: A and F scenarios pass at all α_u (8/8); I scenarios pass
at α_u ≤ 0.40 but fail at α_u ≥ 0.60.

**Spatial diagnostic at I-scenario failure.** Centerline T(z) and CW−engine residual profile
isolated the failure mechanism: center core residual ≈ 0°F across all α_u (kinetics correctly
matched), with residual concentrated near the bottom soil boundary growing from −0.065°F at
α_u=0.20 to −0.49°F at α_u=0.80. Mechanism is the bottom-side corner BC stencil asymmetry
documented as deferred work in Sprint 7 §"Bottom-side corner stencil correction" — pure
strong Dirichlet write at the bottom corner without the quarter-cell + half-cell treatment
applied at the top-side corner. The asymmetry was small at Sprint 7's α_u=0.10 (0.085°F R3 RMSE);
at Sprint 8's α_u=0.80 it propagates to ~0.49°F at the near-bottom region.

---

## 5. Final Validation

### 5.1 32-point A-scenario sweep at fitted c(T_pl)

| T_pc °F | α_u=0.20 | α_u=0.40 | α_u=0.60 | α_u=0.80 |
|---|---|---|---|---|
| 40 | 0.0102 | 0.0100 | 0.0118 | 0.0144 |
| 50 | 0.0095 | 0.0108 | 0.0137 | 0.0194 |
| 60 | 0.0095 | 0.0118 | 0.0193 | 0.0207 |
| 70 | 0.0102 | 0.0129 | 0.0157 | 0.0202 |
| 80 | 0.0134 | 0.0150 | 0.0188 | 0.0258 |
| 90 | 0.0148 | 0.0257 | 0.0285 | 0.0414 |
| 100 | 0.0281 | 0.0442 | 0.0683 | 0.0887 |
| 110 | 0.0405 | 0.0739 | 0.1079 | 0.1599 |

(Values are max|T_engine − T_CW| over t ∈ [0, 168] hr at the geometric core, in °F. All 32 points
are well under the 0.35°F gate.)

### 5.2 12-dataset Sprint 8 Stage 2 validation at fitted c(T_pl)

| Scenario | T_pl/T_soil | α_u=0.20 | α_u=0.40 | α_u=0.60 | α_u=0.80 |
|---|---|---|---|---|---|
| A | 73/73 | PASS | PASS | PASS | PASS |
| F | 73/45 | PASS | PASS | PASS | PASS |
| I | 100/73 | PASS | PASS | FAIL | FAIL |

10/12 pass. The 2 failures (I scenario at α_u ≥ 0.60) are driven by bottom-side BC stencil
asymmetry, isolated via spatial diagnostic. Center core residual ≈ 0°F at all α_u confirms
kinetics correction is doing its job; the residual is pure BC artifact.

### 5.3 Sprint 7 regression check

Sprint 7's 9-run gate at α_u=0.10 still passes with the Sprint 8 corrections layered on top.
R1 max = 0.229°F, R2 max = 0.219°F (within Sprint 7's reported margins).

---

## 6. Engine State at Sprint 8 Close

| Item | Status | Stage / Origin |
|---|---|---|
| `k_uc × 0.96` (constant) | committed in `thermal_engine_2d.py` | Sprint 7 |
| `model_soil` flag (default False) | committed | Sprint 7 |
| `is_submerged` toggle | committed | Sprint 7 |
| `blanket_thickness_m` parameter | committed | Sprint 7 |
| CFL guard, `dy_minus` guard | committed | Sprint 7 |
| `compute_hu_factor` / kinetics correction | committed | Sprint 7 |
| **Hu_residual = 12,937 J/kg** | wrapper-only | Sprint 8 |
| **c(T_pl) quadratic correction** | wrapper-only | Sprint 8 |
| Bottom-side corner BC stencil | unchanged (Sprint 7 deferred) | — |

The Sprint 8 corrections (Hu_residual and c(T_pl)) are applied as wrapper overrides for CW
comparison runs and have not been committed to the engine source. Whether to integrate them
into the engine source — and in what form (interpolation table, fitted function, or a more
elaborate kinetics path) — is a separate decision deferred from this closure.

---

## 7. Key Findings Outside the Calibration

**CW solver Hu floor.** Unprompted minimum heat-release floor at approximately Hu_min ∈ (3000, 10000) J/kg.
Inputs below this are silently substituted to the floor value. The floor produces ΔT_floor at the
core scaling as `Hu_residual × α(t) × constants` — at α_u=0.10 this is the ~0.13°F bulk warming
documented in Sprint 7 as a CW limitation; at α_u=0.80 it grows to ~1.85°C ≈ 3.33°F. The floor was
not a known feature of CW prior to Sprint 8 characterization.

**Engine vs CW kinetics shape disagreement.** When both solvers are operated at matched effective Hu
(i.e., engine at Hu_residual = 12,937 J/kg), their α(t) trajectories diverge at temperature extremes
in a smooth, α_u-independent pattern. Empirical c(T_pl) quadratic fit captures the divergence to
~1% precision in the [40, 110]°F range. The mechanism is consistent with Arrhenius kinetics
disagreement (the engine and CW likely differ in T_ref or activation energy by a small amount), but
no direct verification of which parameter differs has been performed. The 2°C T_ref disagreement
identified in Stage 2-prep produces only ~0.005°F residual change when corrected, ruling it out as
the dominant mechanism. The actual cause remains an open question.

**Bottom-side corner BC stencil asymmetry (pre-existing, now quantified).** Sprint 7 documented this
as a deferred limitation. Sprint 8 spatial diagnostic at I-scenario α_u=0.80 measured the residual
contribution: ~0.49°F near the bottom soil boundary, decaying to ~0°F at the geometric core. This
is the only remaining residual contributor after Sprint 8 calibrations are applied. Magnitude
scales with α_u (heat flux through the bottom BC).

---

## 8. What Was Deferred

1. **Bottom-side corner BC stencil correction.** Sprint 7's open question, now quantified as the
   sole driver of Sprint 8's 2 remaining gate failures. Implementation: mirror the top-side
   quarter-cell + half-cell stencil to the bottom corner with appropriate sign conventions.
   Estimated effort: 2-4 hours including validation. Recommended as Sprint 9's first task.

2. **c(T_pl) fitting range extension.** The quadratic was fit to T_pl ∈ [40, 110]°F. Production
   runs outside this range are not validated. Future work may extend the calibration sweep to
   T_pl ∈ [30, 120]°F or wider, and may replace the quadratic with an Arrhenius-form fit that
   extrapolates more safely. Effort: ~32 additional CW runs at the new T_pl values plus engine
   re-sweep, ~2-3 hours total.

3. **T_ref alignment.** Engine's T_ref ≈ 23°C, CW's T_ref ≈ 21°C from Stage 2-prep curve-fit.
   2°C shift produces ~0.005°F residual change — small but real. If a physically-motivated
   Arrhenius form is desired for the c(T_pl) correction (replacing the empirical quadratic),
   T_ref alignment would be a natural prerequisite. Not blocking Sprint 8 closure.

4. **ρ_concrete and Cp_concrete cross-validation.** Untouched in Sprint 7 and Sprint 8. Spatial
   diagnostic in Sprint 8 §4.11.8 showed no evidence of ρ/Cp disagreement (core temperature
   matches between engine and CW after corrections), but a parametric sweep would tighten the
   bound. Low priority unless production discrepancies appear.

5. **Engine source integration of Sprint 8 corrections.** Hu_residual and c(T_pl) currently apply
   as wrapper overrides. Whether and how to commit these to `thermal_engine_2d.py` (e.g., as a
   `cw_compatibility_mode` flag, a lookup table, or built-in kinetics modification) is a separate
   design decision. Not blocking production use of the calibration; production wrappers can apply
   the corrections as needed.

6. **CW Hu floor — manual verification.** The floor's value, mechanism, and α_u dependence were
   characterized empirically. Reading the CW V3 manual section on hydration kinetics to confirm
   the floor's origin (e.g., a numerical safeguard, a documented minimum, or an undocumented
   solver behavior) would be valuable for understanding but is not blocking.

7. **Hydration kinetics calibration (τ, β, α_u, E_a).** Sprint 7 carried these from prior sprints;
   Sprint 8 did not touch them. The c(T_pl) correction empirically captures the engine-vs-CW
   kinetics disagreement; whether to refit the underlying kinetics parameters directly is a
   separate, larger effort.

---

## 9. Next Phase

Sprint 9 is recommended to begin with the bottom-side corner BC stencil correction (deferred
work item 1 above). With that closed, the I-scenario α_u ≥ 0.60 residuals should drop into the
gate, completing the calibration coverage across the full Sprint 8 dataset. After that, Sprint 9
is open to take on whatever the next calibration axis requires.

The Sprint 8 calibration record consists of:

```python
# Engine wrapper for CW comparison (or CW-equivalent production runs)
mix.Hu_J_kg_effective = 12_937.0  # Sprint 8 §4.10 — CW floor compensation
T_pl_F = ...  # placement temperature in °F
c = -1.3025 + 0.04746 * T_pl_F - 2.081e-4 * T_pl_F**2  # Sprint 8 §4.11.8
# Validity: T_pl ∈ [40, 110] °F; extrapolation not validated
mix.alpha_u_effective = c * mix.alpha_u_nominal
# k_uc factor: unchanged (Sprint 7 calibration: 0.96)
```

Sprint 8 closes.

---

## Appendix A — File Inventory

Sprint 8 working tree under `validation/sprint8/`:

- `stage2_calibration/` — Stage 2 initial attempt (12 datasets, k_uc sweep)
- `stage2_diag_pre/` — Centerline profile diagnostic
- `stage2_diag_Hu/` — Hu floor characterization (§4.1–§4.9), engine side-by-side, plus follow-up
  stages (§4.10 floor test, §4.11 α_u(T_pl) characterization with v2 extended sweep, v3 refinement,
  fit + validate, spatial diagnostic)
- `stage2_alpha_u_T_trend/` — 32-point A-scenario sweep at engine-side, c_optimal tables and plots

All work uncommitted to main during Sprint 8 per scope discipline. Engine source unchanged.

The unified Sprint 8 report is at
`validation/sprint8/stage2_diag_Hu/STAGE2_DIAG_HU_REPORT.md` (§4.1 through §4.11.8).
This closure document supersedes individual stage reports for the calibration record.

---

**Sprint 8 closed.**
