# Phase A Scoping Brief — Engineering-Tolerance Validation Campaign

**Status:** DEFERRED — pending Sprint 10 stencil fix and compute_hu_factor audit
**Predecessor:** Sprint 9 (realistic mix pilot, 2026-05-04)
**Intended executor:** New Claude Code session after pre-launch checklist (§12) is complete

---

## 1. Goal

Validate the Sprint 9 calibrated engine (k_uc×0.96 + Hu_residual conditional + c(T_pl)
quadratic + Hu_factor composition correction + B1 inverse-comp) against a representative
sweep of production placement scenarios, using 5 production mixes × 5 temperature scenarios =
25 CW datasets.

The core question: after the Sprint 10 stencil correction is committed, does the engine
meet the engineering-tolerance gate (|T_engine − T_CW| < [gate TBD, see §6]) across the
full production operating envelope for both bulk core temperature and maximum temperature
differential?

This phase is **engineering validation**, not calibration. No new corrections are developed
during Phase A. If the gate is not met, the outcome is a diagnosis brief for the next
targeted fix, not an in-session calibration tuning.

---

## 2. Mixes (5)

| Mix ID | Composition | α_u | Hu (J/kg) | τ (hr) | β | Kinetic role |
|---|---|---|---|---|---|---|
| MIX-01 | Type I/II + 22% FA + 17% slag | 0.7585 | 424,143 | 29.4 | 0.895 | Reference moderate kinetics |
| MIX-04 | Type I/II, no slag, high cement (450 lb) | TBD | TBD | TBD | 1.225 | Fast-kinetics, high-heat extreme |
| MIX-06 | Type I/II + 22% FA + 50% slag | TBD | TBD | TBD | 0.637 | Mid slag-extreme, sub-baseline β |
| MIX-07 | Type I/II + 22% FA + 70% slag | 0.8935 | 463,076 | 75.1 | 0.516 | Slow-kinetics, high-slag extreme — pilot-validated |
| MIX-10 | Type I/II + 35% FA + 17% slag | TBD | TBD | TBD | 0.871 | Higher fly ash, moderate kinetics |

MIX-01 and MIX-07 are carried from the Sprint 9 pilot. MIX-04, MIX-06, and MIX-10 complete
the production kinetic envelope; their exact parameters come from the user's mix design table.
The user should provide `input.dat` files for each mix with correct T_pl, T_soil, α_u, Hu, τ,
β, and E_a hand-set to match the production mix table.

**MIX-13 (silica fume) deferred from Phase A.** `compute_hu_factor`'s behavior on SF-bearing
compositions is uncharacterized — SF's pozzolanic activity and heat release timing differ from
cement/FA/slag, and including an SF mix before the function is formally specified would conflate
engine validation with SF-regime uncertainty in `compute_hu_factor`. SF mix coverage deferred
to a later phase after `compute_hu_factor` formal specification (Sprint 10 §5.2).

---

## 3. Scenarios (5)

| Scenario ID | T_pl (°F) | T_soil (°F) | T_ambient (°F) | Description |
|---|---|---|---|---|
| S1 | 50 | 50 | 50 | Cold placement, cold ground |
| S2 | 73 | 73 | 73 | Nominal (Sprint 9 pilot condition) |
| S3 | 73 | 85 | 85 | Nominal placement, warm ground (Sprint 9 pilot condition) |
| S4 | 90 | 60 | 60 | Hot placement, cool ground |
| S5 | 60 | 90 | 90 | Cool placement, hot ground |

T_ambient = T_soil for all scenarios (weather-override discipline preserved from Sprint 7/8/9).
These 5 scenarios provide cold/neutral/hot coverage for both placement temperature (which
drives c(T_pl)) and soil temperature (which drives the bottom BC heat flux direction).

---

## 4. Geometry

Same as Sprint 7/8/9 baseline:
- 40 ft wide × 80 ft deep, half-symmetric
- `is_submerged = True`
- `model_soil = False` (T_soil applied as direct Dirichlet at di=48)
- `blanket_thickness_m = 0.0`
- 6× grid refinement for concrete domain (as established in Sprint 7)

No geometry changes during Phase A.

---

## 5. Metrics

**Primary engineering metrics** (used for the gate; defined formally in §6):
- `T_max` — maximum temperature in the validated region (di ∈ [5, 48], all wi, t ∈ [0, 168 hr])
- `ΔT_max` — maximum spatial differential (max − min over di ∈ [5, 48], all wi) at any
  timestep, then max over t

Both metrics are computed independently for the engine and CW; the gate applies to the
absolute differences `|T_max_engine − T_max_CW|` and `|ΔT_max_engine − ΔT_max_CW|`.

**Diagnostic metrics** (computed for every run; not gated):
```
R(di, t)  = T_engine(di, wi=0, t) − T_CW(di, wi=0, t)  [°F]
max|R|    = max_{di ∈ [24, 48]} |R(di, t=168)|         [centerline t=168 residual]
R(di) profile at t=168 over di ∈ [24, 48]              [for spatial diagnosis]
T_core    = T_engine(di=24, wi=0, t=168)               [hydration-completeness sanity check]
```

The diagnostic metrics support failure-mode isolation when an engineering metric exceeds the
gate, and provide continuity with Sprint 9's stencil/structural-decoupling analysis. They are
not pass/fail criteria.

---

## 6. Gate

**Engineering gate (subject to Tyler Ley confirmation):**
- `|T_max_engine − T_max_CW| < 2°F` across all 25 runs
- `|ΔT_max_engine − ΔT_max_CW| < 2°F` across all 25 runs

**T_max definition:** maximum temperature in the validated region (di ∈ [5, 48], all wi,
t ∈ [0, 168 hr]).

**ΔT_max definition:** maximum spatial differential (max − min over di ∈ [5, 48], all wi) at
any timestep, then max over t.

The validity mask `di ∈ [5, 48]` follows Sprint 7 §3.5 — di=0..4 excluded due to surface-BC
incompatibility, deferred to Phase C top-surface work. If Phase C modifies the validity mask,
Phase A metric definitions must be updated.

`max|R|` (per-node residual at t=168 over di ∈ [24, 48], wi=0) is computed as a diagnostic
metric for understanding scenario behavior, not as a gate. Residual concentrations near di=47
are expected per Sprint 9 §5.1 (deferred stencil work) and do not affect Phase A pass/fail.
If a scenario fails the T_max or ΔT_max gate, max|R| and per-node residuals help diagnose
which physics axis is responsible.

---

## 7. CW Dataset Generation (Out of Scope for Engine Session)

The user generates all 25 CW datasets (5 mixes × 5 scenarios) before the Phase A engine
session begins. The agent does not generate CW datasets.

**CW dataset requirements:**
- One `input.dat` per mix, modified for each scenario's T_pl and T_soil.
- Mix kinetic parameters (α_u, Hu, τ, β, E_a) hand-set to match the production mix table
  at lines 389–390 and the relevant parameter lines (same practice as Sprint 9 pilot).
- T_ambient = T_soil for each scenario (weather-override lines 519–531 set to flat T_soil).
- Geometry unchanged from Sprint 9 baseline (40 × 80 ft submerged).
- Nominal 168 hr simulation window.
- Expected file size per dataset: ~7–8 MB (output.txt).

**Staging:** copy datasets into `validation/sprint10/phase_a/cw_data/{mix01,mix04,mix06,mix07,mix13}/`
with subdirectories per scenario, before handing off to the engine session.

---

## 8. Engine Pipeline

**Per-run wrapper (Sprint 9 closing wrapper + Sprint 10 stencil fix):**

```python
Hu_factor    = compute_hu_factor(mix)         # composition-derived
Hu_factored  = Hu_raw * Hu_factor
c            = -1.3025 + 0.04746*T_pl_F - 2.081e-4*T_pl_F**2  # Sprint 8
alpha_u_eff  = c * alpha_u_raw
Hu_eff       = Hu_factored / c                # B1 inverse-compensation
mix.alpha_u           = alpha_u_eff
mix.Hu_J_kg_effective = Hu_eff
# model_soil=False, is_submerged=True, blanket=0.0, k_uc×0.96
# + Sprint 10 di=47 stencil fix committed in thermal_engine_2d.py
```

**Per-run runtime:** ~1.1 s (from Sprint 9 sweep v2 profiling; may change slightly with
Sprint 10 stencil fix if the fix alters CFL behavior).

**Full sweep compute estimate:** 25 runs × 1.1 s ≈ 28 s total; well within a single session.

**Parallelism:** Phase A runs are independent — a batch loop over all 25 (mix, scenario) pairs
can run sequentially in a single pass. No parallelization required.

---

## 9. Out of Scope for Phase A

The following are explicitly out of scope for the Phase A engine session:

- **No engine-source changes.** The Sprint 10 stencil fix must already be committed before
  Phase A begins. If a new residual mechanism is discovered during Phase A, stop and write
  a diagnosis brief; do not patch the engine in-session.
- **No calibration tuning.** No changes to k_uc, Hu_residual, c(T_pl), or compute_hu_factor
  during Phase A. The calibration is fixed; Phase A is a pass/fail validation.
- **No new CW datasets.** The 25 CW datasets must exist before the session begins. If a
  dataset is missing or corrupt, flag and stop; do not generate CW from within the session.
- **No full-width residual metrics** (di < 24). The gate applies to di ∈ [24, 48] (bottom half)
  only, consistent with Sprint 9 pilot convention.
- **No multi-scenario aggregation.** Report results per (mix, scenario) pair. Do not compute
  cross-scenario averages or aggregate performance statistics — those belong in the Phase B
  production-rollout decision, not Phase A.

---

## 10. Estimated Effort

**Engine sweep and analysis:** ~1 day.
- Wrapper setup and per-mix verification: ~2 hr
- 25-run sweep (compute ≈28 s): ~1 hr including sanity checks and logging
- Residual heatmap generation (5 metrics × 5 scenarios × 2 mixes representative set): ~2 hr
- Table population and gate verdict per cell: ~1 hr

**Phase A report writing:** ~1 day.
- Results tables (max|R|, T_max_core, ΔT_max per cell): ~2 hr
- Spatial diagnostics for any gate-fail cells: ~2 hr
- Synthesis and pass/fail summary: ~2 hr

**Total:** ~2 working days.

---

## 11. Decision Rule

After all 25 runs complete:

**If all 25 cells pass the gate (max|R| < gate threshold):**
→ Phase A passes. Write Phase A closure report. Proceed to Phase B (production rollout
brief). No further calibration work needed at this scope.

**If 1–5 cells fail, residual concentrated at a known mechanism (e.g., di=47 stencil
if Sprint 10 fix is incomplete, or di=24 bulk for one scenario type):**
→ Isolate the failure mechanism via spatial diagnostic. Write a targeted diagnosis brief
for the next sprint (Sprint 11 or equivalent). Do not loop back into Phase A calibration.

**If >5 cells fail, or failures are spatially dispersed (bulk-uniform, not di=47-local):**
→ Likely a new mechanism not covered by Sprint 7–9 calibration. Write a full diagnosis
brief. Escalate to the user before proceeding.

In all cases: do not modify the wrapper or engine source during Phase A. Report outcomes as
a pass/fail grid and surface it for the user's decision on next steps.

---

## 12. Pre-Launch Checklist

Before starting the Phase A engine session, confirm ALL of the following:

- [ ] **Sprint 10 stencil fix committed.** `thermal_engine_2d.py` contains the bottom-side
  di=47 stencil correction. Verify by running a smoke test: MIX-01 at S3 (T_pl=73, T_soil=85)
  should show max|R| substantially below the 1.05°F Sprint 9 baseline.
- [ ] **compute_hu_factor specification audited.** Denominator logic verified for all 5
  Phase A mixes. Output values documented and physically plausible for each mix composition.
  See `COMPUTE_HU_FACTOR_NOTE.md` §6.
- [ ] **All 25 CW datasets present.** Verify `validation/sprint10/phase_a/cw_data/` (or
  equivalent staging path) contains one folder per mix per scenario with non-empty output.txt.
  Check file sizes (expect ~7–8 MB each).
- [ ] **Gate threshold confirmed.** User has explicitly approved the gate value (2°F placeholder,
  or a different value if the engineering tolerance requirement is different).
- [ ] **c(T_pl) validity range checked.** Scenarios S1 (T_pl=50°F) and S4 (T_pl=90°F) are
  within the Sprint 8 calibration range [40, 110]°F. No extrapolation needed for the 5 chosen
  scenarios.
- [ ] **No uncommitted engine source changes.** `git status` shows a clean `thermal_engine_2d.py`
  (the Sprint 10 stencil fix committed, no further modifications).
