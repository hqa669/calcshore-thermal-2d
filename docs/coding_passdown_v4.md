# CalcShore Thermal Engine v3 — Passdown v4 (Sprint 3+)

## How to use this document

This is the active planning + decision document for the CalcShore thermal
engine starting from Sprint 3. Sprint 0/1/2 retrospectives and the
original mission scope are in `coding_passdown_v3.md`; v4 picks up from
`sprint-2-complete` (commit ~dca1dca, tagged 2026-04-24).

Structure:
- **§1 Mission**: what the engine is for, what "done" looks like
- **§2 Active sprint scope**: Sprint 3 concrete scope
- **§3 Sprint roadmap**: Sprints 3–6 outline (mutable)
- **§4 Risk register**: known-unknowns tracked across sprints
- **§5 Architecture decisions (ADRs)**: finalized decisions, findable
  without digging through retrospectives
- **§6 Workflow discipline**: practices that worked in Sprints 0–2,
  updated with Sprint 2 lessons
- **§7 Sprint 3 retrospective**: populated at Sprint 3 close
- **§8 Sprint 4+ retrospectives**: populated as sprints close

---

## §1 Mission

The CalcShore thermal engine generates temperature predictions for mass
concrete placements at timestep resolution, outputting trajectories that
feed the TCP (Temperature Control Plan) document generator and the AI
compliance evaluator.

"Done" for the engine means:

1. **Accurate**: predictions match ConcreteWorks (CW) reference
   trajectories to within S0 tolerances across a 15-mix validation
   library covering diverse climates and mix designs
2. **Fast**: runtime on typical geometry (half-mat, 168hr) under 3
   seconds (currently 1.3s)
3. **Robust**: inputs that deviate from MIX-01 conventions handled
   gracefully; no silent failures on edge-case mixes
4. **Parameterized**: form material, form orientation, site latitude,
   and other customer-varying parameters exposed cleanly

Current state at Sprint 2 close: MIX-01 validation complete at S0
tolerances (5/5 metrics pass on [48, 168]hr window). 92 tests. 15-mix
validation, production robustness, and dual-peak hydration all pending.

---

## §2 Active sprint scope — Sprint 3

**Theme**: Barber soil + S1-aspire recalibration on MIX-01 only.

**Explicit deferrals (not Sprint 3 scope)**:
- Hydration-rise shape divergence — deferred to Sprint 5+ (requires
  dual-peak hydration model, structural change not closable within
  current single-peak α(t_e) formulation)
- 15-mix validation — deferred to Sprint 4 (requires CW export work)
- F_VERT_BY_ORIENTATION stub replacement — deferred to Sprint 4 (depends
  on 15-mix findings)
- Form-material R_form parameterization — deferred to Sprint 6+

### Sprint 3 items

#### Item 3.1: Barber soil model (primary deliverable)

Replace `T_ground = T_amb` in the form-face LW ground term with a
diurnal-lagged ground surface temperature.

**Current state (post-Sprint-2)**: The form-face LW balance at the
vertical form face includes a ground-view term `ε_form·ε_ground·σ·
(1 − F_SKY_VERT)·(T_form⁴ − T_ground⁴)`. `T_ground` is pinned to
`T_amb` because Sprint 2 deferred soil modeling. This is a conservative
baseline — ground surface in summer Austin is actually warmer than
ambient during the day (solar absorption) and cooler at night (sky
radiation), with ~4–6 hour diurnal phase lag.

**Physics approach**: Implement a simplified Barber-style model. Two
options:

**Option 3.1.A — Sinusoidal lag model (simpler, recommended start)**:
```
T_ground(t) = T_amb_daily_mean + (T_amb(t − lag) − T_amb_daily_mean) · damping
```
Two empirical parameters: `lag_hrs` (diurnal phase lag, ~4–6hr typical)
and `damping` (amplitude attenuation, 0.6–0.8 typical for bare soil).
Easy to implement, validates quickly, handles most of the physics.

**Option 3.1.B — Barber 1957 diffusion model (first-principles)**:
Compute ground surface temperature from air temperature forcing via the
semi-infinite soil heat equation with known soil thermal diffusivity.
Five-parameter model (α_soil, ρ_soil, c_soil, skin-layer depth, surface
albedo). Physically principled; may be overkill for the magnitude of
effect expected.

**Recommendation**: Start with Option 3.1.A. If validation shows the
sinusoidal model is insufficient, upgrade to 3.1.B in a follow-up PR.
The effect on Corner RMS is expected to be small (0.05–0.2°F), so a
5-parameter model is unlikely to outperform a 2-parameter model enough
to justify the complexity.

**Implementation scope**:
- New helper function (e.g., `ground_surface_temperature(t_hrs, env)`)
  that takes the environment's hourly T_amb array and returns
  lag-and-damp-modulated T_ground at arbitrary time
- Wire into the form-face Newton residual (currently uses `_T_amb_C`
  for ground term)
- Add parameters to CWConstruction or CWEnvironment (`soil_lag_hrs:
  float = 5.0`, `soil_damping: float = 0.7` as reasonable defaults)
- Populate a new diagnostic field `T_ground_C_history` on HydrationResult
- Update `compare_to_cw.py` if relevant (probably a small annotation
  rather than a new panel)

**Expected PR split**:
- PR 10: Ground surface temperature helper + CWConstruction fields +
  plumbing. Zero behavior change if `soil_lag_hrs=0, soil_damping=1`
  (which reduces to `T_ground = T_amb`). Tests for the helper's lag
  and damping behavior.
- PR 11: Activate non-default soil parameters (`soil_lag_hrs=5,
  soil_damping=0.7`); wire into form-face Newton residual; validate
  against CW.

**Success metric (primary)**: All 5 S0 gates remain PASS on MIX-01.
**Success metric (secondary)**: Corner RMS improves by ≥0.05°F vs
Sprint 2 baseline (2.22°F → ≤2.17°F). Zero regression on other metrics.

**Success metric (aspirational)**: Ground surface temperature lag and
damping produce a visible diurnal-phase improvement in form-face LW
flux trace (panel k of the 4×3 comparison grid).

#### Item 3.2: S1-aspire recalibration (small follow-on)

Update the S1-aspire gates in `compare_to_cw.py` per Sprint 2 close
decisions:

| Metric | Current S1-aspire | New S1-aspire | Sprint 2 result | Status |
|---|---|---|---|---|
| Peak Max T | ±0.5°F | ±0.5°F (unchanged) | −0.3°F | ✓ met |
| Peak Gradient | ±1.0°F | ±1.0°F (unchanged) | −0.3°F | ✓ met |
| Field RMS | ≤1.0°F | ≤1.0°F (unchanged) | 0.88°F | ✓ met |
| Centerline RMS | ≤0.5°F | ≤0.5°F (unchanged) | 0.74°F | ✗ missed (kept to pull) |
| Corner RMS | ≤3.0°F | ≤2.0°F (tightened) | 2.22°F | ✗ missed (creates pull) |

**Expected PR**: PR 12 — small, comment-only or one-constant-change.

### Sprint 3 target

- 3 PRs (PR 10, PR 11, PR 12)
- 2-week duration
- 92 tests → ~98 tests (estimated +5 for Barber helper + integration)
- Tag: `sprint-3-complete`

---

## §3 Sprint roadmap (mutable)

This is the plan of record. Each sprint's scope may change based on
findings from prior sprints.

### Sprint 4 — Multi-mix thermal-physics validation on hydration-clean evaluation set (~6–10 weeks)

Updated 2026-04-25 from "15-mix validation and F_vert stub replacement"
based on §7.5 reconnaissance findings.

**Theme**: validate 2D thermal/boundary physics on the subset of mixes
where hydration tracks CW. Mixes where engine and CW disagree on heat
generation (high-SCM, slag-heavy) are out of Sprint 4 scope and route to
Sprint 5's dual-peak hydration work.

**Evaluation set (8 mixes)**: Reference (MIX-01, 02, 03, 11, 12) +
B1 (MIX-04, 08) + B2 (MIX-15). See §7.5.3 for partitioning rationale.

**Routed to Sprint 5 (6 mixes)**: Cluster A (MIX-05, 06, 07, 09, 10) +
MIX-14 (borderline hydration miss).

**Deferred (1 mix)**: MIX-13 (no CW output; re-export pending).

#### PR sequence

- **PR 13 — `run_one()` refactor + display fixes**. Lift `main()` body
  in `compare_to_cw.py` into `run_one(scenario_dir) -> dict` (Approach A,
  per §7.5 decision). Separate `print_gate_table(result)`. Replace
  hardcoded MIX-01 peak-time constants (lines 41/43) with per-scenario
  extraction. Replace hardcoded PNG path (line 437) with parameterized
  output path. Add `run_all.py` driver with PNG-off-by-default flag.
  MIX-13 sentinel (no `output.txt` → `{"skipped": True, ...}`). MIX-01
  bit-identical regression carries the scaffolding-only contract. ~150
  LOC. Tests for sentinel handling and parameterized PNG path.

- **PR 14 — Baseline data commit**. Run all 14 evaluable mixes via PR 13
  driver. Commit `validation/sprint4_baseline.md` (per-mix gate table,
  evaluation-set partitioning columns). Write §7.6 interim retro
  documenting which mixes hit which gates. No engine changes. Commit IS
  the deliverable.

- **PR 15 — B1 calibration (composition-isolated thermal physics)**.
  Targets MIX-04 + MIX-08. Both under-predict Peak Max + Peak Gradient
  with high Corner RMS at low/mid SCM and 60°F placement. Shared
  signature → likely form-face or top-boundary term with composition
  dependence we treat as constant. Hypotheses to investigate (do
  diagnostics first per §6.3, no pre-committing the fix):
    - F_vert calibration composition-dependent (R1/R2 still in play)
    - R_form composition-dependent (currently constant 0.0862, ADR-04)
    - Some boundary heat-flux term scaling with mix composition
  Stop conditions explicit per §6.6. Ablation of any new physics term
  per §6.4.

- **PR 16 — B2 calibration (placement-temperature-isolated boundary
  physics)**. Targets MIX-15. Cold-placement (45°F vs 60°F Reference)
  with all other parameters identical to Reference. Round-1 PeakGrad Δ
  −8.3°F is off-charts vs Sprint 0–3 experience. Hypotheses:
    - Initial-condition-dependent boundary term (some flux coefficient
      evaluated at placement temp rather than current temp)
    - Cold-placement-specific physics (condensation latent heat?)
    - IC propagation issue when |T₀ − T_amb| is large
  Diagnostics-first per §6.3. Bit-identical regression on Reference set
  required.

- **PR 17 — R4 disposition**. With multi-mix evidence in hand:
  `soil_lag_hrs` default → 0.0; `soil_damping` deprecated (default → 1.0,
  parameter retained for future cold-climate exports). Per §7's Sprint 3
  finding generalized over 8 evaluation-set mixes. ~5-line change +
  comment. Updates ADR-08.

#### Sprint 4 close criteria

- **Floor**: PR 13 + PR 14 land. Reference set (5/5) holds bit-
  identically. R4 decided. PR 15 lands either fixing B1 or
  characterizing B1's failure mode and routing residual to Sprint 6
  with named risk.
- **Target**: B1 (MIX-04, 08) S0 PASS after PR 15. B2 (MIX-15) S0 PASS
  after PR 16. R4 cleaned up in PR 17. Evaluation set 8/8 S0 PASS.
- **Stretch**: B1 + B2 reach S1-aspire on Centerline and Corner RMS.

#### Open risks (see §4)

R1/R2 (F_vert calibration) — owned by PR 15 if B1 hypothesis lands.
R5/R7 (form material, RuntimeWarning cleanup) — Sprint 6 still.
R8 (CW hydration parameter inheritance) — accepted-by-design.

#### Tag plan

- `pr-13-complete`, `pr-14-complete`, `pr-15-complete`, `pr-16-complete`,
  `pr-17-complete`, `sprint-4-complete`. Per §6.8 (new Sprint 3 lesson),
  verify each tag exists locally before trusting Claude Code's report.

### Sprint 5 — Reference-set residual characterization (~3–5 weeks)

Updated 2026-04-25 from "dual-peak hydration model" based on Sprint 4
findings. The dual-peak hydration work has been routed to a separate
hydration sprint series (see "Explicit deferrals" below), and Sprint 5
takes on a narrower, characterization-first scope on the Reference set
that closed Sprint 4 at S0 5/5 but only S1-aspire 3/5.

**Theme**: Reference-set residual characterization — distinguish
thermal-physics terms from inherited kinetics drift before deciding on
a fix. Sprint 5 is characterization-first, fix-conditional.

**Evaluation set (5 mixes)**: Reference (MIX-01, 02, 03, 11, 12). All
39.1% SCM, 60°F placement, all S0 5/5 at sprint-4-complete, all 3/5
S1-aspire (Centerline RMS ✗ on 5/5; Corner RMS ✗ on 4/5).

**Two investigations**:

1. **Centerline RMS systematic offset** (~0.74°F across all 5 Reference
   mixes; MIX-02 outlier at 0.57°F). Systematic across diverse mix
   designs but identical construction → likely a thermal-physics term
   or shared kinetics drift, not mix-specific calibration.
2. **Corner RMS calibration**. PR 15 found R_form has clean authority
   on Corner RMS within Reference's tolerance window; R_form = 0.060
   looked promising. Sprint 5 evaluates this on the full Reference set.

#### PR sequence (provisional, characterization-first)

- **PR 18 — Reconnaissance spike for Investigation 1**. Characterize
  what kind of error the Centerline RMS residual is *before* any code
  change. Diagnostic candidate set:
    - Constant offset (engine vs CW differ by a near-constant ΔT
      throughout the steady-state window)
    - Phase lag (engine and CW have similar amplitude but offset in time)
    - Amplitude scaling (engine's diurnal swing larger or smaller than
      CW's)
    - Depth-localization (residual concentrated at specific depths
      rather than uniform)
    - **Mode-A overlap** (does the Reference centerline residual share
      the same time/temperature signature as PR 16 Phase 2.6 Mode A —
      engine warmer than CW during active hydration, t = 8–48 hr? If
      yes, Sprint 5's centerline work partially characterizes R9's
      Mode A — it strengthens R9 routing, doesn't replace it.)
  Spike outputs go to `validation/diagnostics/pr18_centerline_recon/`.
  Not committed to engine source. Per §6.3 (diagnostic-script-before-
  fix) and §7's "run the parameter sweep before narrating physics"
  lesson generalized.

- **PR 19 — Investigation 1 fix (conditional)**. Land if PR 18 finds
  a thermal-physics term with sufficient authority on Centerline RMS
  (≥0.24°F on Reference if S1-aspire 0.5°F closure is the goal, or a
  proportionate target if PR 18 reframes the goal). Skip and route to
  engine v3 release notes as a known limitation if PR 18 finds the
  residual is dominated by inherited kinetics drift (R9 Mode A
  overlap) — that path keeps Sprint 5 disciplined and routes the
  mechanism cleanly.

- **PR 20 — Investigation 2**. R_form Reference-set evaluation. PR 15
  found R_form=0.060 brings B1 to 4/5 with MIX-01 holding 5/5; PR 20
  evaluates the full Reference set at R_form=0.060 vs the 0.0862
  baseline. Commit conditional on Reference set holding S0 5/5 at
  the new value AND ≥3/5 mixes improving on Corner RMS without
  Centerline degradation. Updates ADR-04 if landed.

- **PR 21 — Sprint 5 close + retrospective**. Sprint roadmap update
  reflecting post-recon Sprint 5 disposition; §8 Sprint 5
  retrospective; tag `sprint-5-complete`. Mirrors PR 17's role in
  Sprint 4.

#### Sprint 5 close criteria

- **Floor**: PR 18 lands characterizing the centerline + corner
  residual mechanisms. R_form recalibration committed if PR 20's
  Reference-set evaluation confirms PR 15's finding (or routed to
  engine v3 release notes if it doesn't). §3 roadmap updated to
  reflect post-recon Sprint 5 scope.
- **Target**: Reference set residuals improve on ≥3/5 mixes for
  Centerline AND Corner. Engine v3 ships with documented residuals.
- **Stretch**: Reference set hits S1-aspire 5/5 on all 5 metrics for
  all 5 mixes.

#### Open risks (see §4)

R5 (R_form form material parameterization) — Sprint 6 still.
R7 (RuntimeWarning cleanup) — Sprint 6 still.
R9 (Arrhenius rate factor divergence above ~30°C) — open; Sprint 5
PR 18 may strengthen routing via Mode-A overlap test but does not
own R9's resolution.

#### Tag plan

`pr-18-complete`, `pr-19-complete` (if landed), `pr-20-complete`
(if landed), `pr-21-complete`, `sprint-5-complete`. Per §6.8 and
§6.9, all tags are local until user pushes; verify locally before
reporting done.

### Sprint 6 — Cleanup + engine v3 release (~3–4 weeks)

Final sprint of the thermal sprint series. Closes the series with
engine v3 as the deliverable.

**Theme**: production-readiness cleanup on the validated Reference-set
baseline.

**Scope**:
- R5: parameterize R_form by form material
  (steel/plywood/plastic-liner) via `CWConstruction.form_type`
  field. The field pre-exists on the dataclass (loader populates from
  CW input.dat index 432); PR 22 wires consumption via
  `resolve_r_form()` resolver replacing the hardcoded
  `R_FORM_CONTACT_SI` constant in the engine. Steel value is 0.0862
  m²·K/W (ADR-04 reinforced by PR 20 evaluation; R_form=0.060 blocked
  by MIX-02 kinetics-driven contradiction). Plywood 0.17 m²·K/W from
  ACI 306R-88 Table 7.3.5 / §7.3 (cw_validated=False,
  OUT-OF-ENVELOPE). Plastic-liner value unpinned (no defensible
  source found; raises NotImplementedError). R5 does NOT attempt to
  fix the MIX-02 cluster contradiction (kinetics-driven; routes to
  hydration series).
- R7: suppress cosmetic `RuntimeWarning: invalid value encountered in
  divide` at harmonic-mean k-divide (`thermal_engine_2d.py:1619, 1624`).
  Currently mitigated via `np.where` guard; cleanup is removing the
  pre-guard divide pattern.
- Engine v3 release: extends `docs/engine_v3_release_notes.md` (PR 20
  started the Centerline residual section). Sprint 6 adds:
  Corner-side residuals (MIX-02 documented as a known mix-regime
  limitation; cluster CornerRMS ~2.0°F is the validated envelope
  statement, not closed); R5 result section; R7 cleanup section;
  final validated envelope (geometry, climate, mix range, S0
  tolerances, known residuals enumerated); out-of-envelope conditions
  (cold placement [R9], high-SCM [Cluster A / hydration series],
  non-half-mat geometry, non-Austin climate, kinetics-divergent mixes
  [MIX-02-class]). API cleanup and customer-facing docstrings.
- Tag `sprint-6-complete` and `engine-v3` on the same commit.

#### Sprint 6 close criteria

- **Floor**: R5 lands; R7 lands; engine v3 release notes
  document the validated envelope and known residuals.
- **Target**: Floor + customer-facing API documented.

### Explicit deferrals (separate sprint series, not bundled into thermal sprints)

These are tracked here so they don't get mistakenly pulled back into
Sprint 5 or Sprint 6 scope. Each is a separate future sprint series,
not a single sprint.

- **Hydration sprint series** — dual-peak hydration model, Cluster A
  recovery (MIX-05, 06, 07, 09, 10), MIX-14 borderline recovery, and
  R9 disambiguation/closure. The deferred first-48hr shape divergence
  identified in Sprint 2 lives here. R3 (accepted, deferred) belongs
  to this series. R8 (CW parameter inheritance) is the inherited-
  calibration adjacent risk. R9 (Arrhenius rate factor divergence
  above ~30°C) routes here if disambiguation lands as kinetics-aligned.
- **Geometry coverage series** — full-mat, slab, column, and other
  non-half-mat geometries. R6 belongs to this series.
- **Multi-climate validation series** — non-Austin sites, cold-climate
  exports, R2 (F_VERT_BY_ORIENTATION non-"unknown" values) testable
  here once non-"unknown" form orientations appear in a CW library.

The thermal sprint series ends at `sprint-6-complete` / `engine-v3`.

---

## §4 Risk register

Known-unknowns tracked explicitly so they don't get lost across sprints.
Each risk has an owner-sprint (which sprint is expected to address it)
and a status (open / mitigated / accepted).

| # | Risk | Owner sprint | Status |
|---|---|---|---|
| R1 | F_VERT_BY_ORIENTATION["unknown"]=0.15 is MIX-01 calibrated only; may not generalize | Sprint 4 PR 15 | Open — diff confirms all evaluation-set mixes have form_orientation="unknown"; if PR 15's B1 fix needs F_vert recalibration, candidate value comes from MIX-04+MIX-08 jointly |
| R2 | F_VERT_BY_ORIENTATION stub values (south=0.35, east/west=0.42, north=0.20) are geometric guesses, unvalidated | Sprint 6+ | Deferred — entire 14-mix library uses form_orientation="unknown"; non-unknown orientations need new CW exports before R2 is testable |
| R3 | First-48hr hydration-rise shape divergence between engine and CW; likely requires dual-peak model | Sprint 5 | Accepted (deferred) |
| R4 | Barber soil parameters (lag, damping) chosen from literature; may not match CW's internal soil model | Sprint 4 PR 17 | **Closed** (tag pr-17-complete). Defaults set to soil_lag_hrs=0.0, soil_damping=1.0 (no-op pair; reduces to T_amb behavior). Helper code retained for future cold-climate data. Reference §7.5.2 (no climate variation in library) and Sprint 3 sweep (damping unidentifiable on MIX-01). ADR-08 updated. |
| R5 | R_FORM_CONTACT_SI=0.0862 hardcoded assuming steel form; customer pilots may use plywood or plastic-lined forms | Sprint 6 PR 22 | **Closed** (Sprint 6 PR 22). Steel value 0.0862 unchanged (ADR-04 reinforced by PR 20). Plywood 0.17 m²·K/W added (ACI 306R-88 Table 7.3.5/§7.3, `cw_validated=False`, OUT-OF-ENVELOPE). Plastic_liner value unpinned (no defensible source found; raises `NotImplementedError`). The MIX-02 cluster contradiction is hydration-series scope (§8.2 finding #1), explicitly not addressed here. |
| R6 | Engine tested on half-mat geometry only; behavior on full-mat, slab, or column geometries unvalidated | Sprint 6+ | Open |
| R7 | Cosmetic RuntimeWarning at harmonic-mean k-divide masks potential future numerical issues | Sprint 6 PR 23 | **Closed** (Sprint 6 PR 23). Engine fix replaced post-hoc np.where mask with np.divide(out=, where=) at the harmonic-mean k-divide sites; run_all.py interim suppression removed; pytest config gained filterwarnings = error::RuntimeWarning to lock the resolution in against future regressions. UserWarning category preserved (PR 22's resolve_r_form cw_validated=False path). |
| R8 | Hydration model uses CW's parameters (τ, β, α_u, Ea, Hu) directly; if those parameters are miscalibrated in CW, our engine inherits the error | Sprint 4–5 | Accepted (inherited by design) |
| R9 | Engine and CW diverge on Arrhenius rate factor at concrete temperatures above ~30°C. Source unknown — possibly different effective Ea, different reference temperature, or a CW-side rate cap not present in engine. Affects cold-placed mixes amplified by sustained high-T residence. MIX-15 −3°F PeakMax deficit traced to this divergence per PR 16 Phase 2.6. Mode A (small early calibration drift, ≤1.7°F, both mixes) is separate and non-blocking. Mode B (late-time divergence, cold-placed mixes only) is the unresolved item. | Sprint 5 if hydration-kinetics-aligned; Sprint 6 if inherited-calibration cleanup | Open |
| R10 | Committed validation baselines (e.g. `validation/sprint4_baseline.md`) and committed test fixtures are downstream artifacts of engine defaults. PRs that change ADR-recorded defaults must regenerate every committed numerical artifact derived from those defaults; otherwise the committed file goes stale and false-positive bit-identity failures surface on later verification runs. | Sprint 5+ (workflow rule) | Mitigated — §6.10 codifies the rule. Precedent: PR 17 changed soil defaults but did not regenerate `validation/sprint4_baseline.md`; staleness surfaced during PR 18's manual verification, fixed in a separate commit on `main`. |

---

## §5 Architecture decisions (ADRs)

Decisions made during Sprints 0–2 that are now load-bearing for future
work. Each entry documents what was decided, when, and why, so future
sprints can find them without digging through retrospectives.

**ADR-01 (Sprint 0)**: 2D half-mat geometry with explicit finite-
difference solver. Decision anchored in v2 validation; not to be
re-litigated without a separate design discussion.

**ADR-02 (Sprint 1 PR 3)**: Quasi-steady T_outer Newton solve on the
top boundary condition (blanket outer surface). 2-step Newton from
linearized initial guess. Convergence verified to <0.01°C.

**ADR-03 (Sprint 2 PR 6)**: Same T_outer pattern applied to the form
outer surface with 4 substitutions: R_form (contact film), F_SKY_VERT
(0.5), ground-view term, α_form. Mathematical twin of ADR-02.

**ADR-04 (Sprint 2 PR 6)**: R_FORM_CONTACT_SI = 0.0862 m²·K/W is
physical contact resistance (steel form + wet concrete film per ACI
347), not legacy calibration. Validated via R_form=0 ablation: Corner
RMS 8.85°F with R_form=0 vs 4.08°F with 0.0862.

*Sprint 6 PR 22 amendment*: The constant is now resolved via
`R_FORM_BY_FORM_TYPE[form_type]` through `resolve_r_form(construction)`
rather than read as a module constant. Steel value 0.0862 unchanged
(PR 20 reinforced). Plywood (0.17 m²·K/W, ACI 306R-88 Table 7.3.5/§7.3)
and plastic_liner (value unpinned, raises `NotImplementedError`) added
with `cw_validated=False`, marked OUT-OF-ENVELOPE in engine v3 release
notes. The `FormTypeRForm` dataclass carries the `cw_validated` flag +
`source` citation as the structural way to hold "documented but
unvalidated" form types forward. R_FORM_CONTACT_SI module attribute
retained as deprecated alias for back-compat with diagnostic harnesses.

**ADR-05 (Sprint 2 PR 8)**: F_VERT_BY_ORIENTATION lookup as primary
path; CWConstruction.vertical_solar_factor is optional override. F_vert
is a low-sensitivity calibration knob (sweep 0.0→0.5 produces only
0.36°F Corner RMS variation); load-bearing physics is LW cooling + h_conv.

**ADR-06 (Sprint 2 PR 8)**: ACI 207.2R Eq 27 orientation-dependent
convection. `h_forced_convection()` for horizontal faces (5.6 + 3.5·0.4·
wind). `h_forced_convection_vertical()` for vertical faces (4.0 +
2.5·0.4·wind).

**ADR-07 (Sprint 2 PR 8.5)**: Validation RMS metrics evaluated on
steady-state window t ∈ [48, 168]hr via `T_START_RMS_HR = 48.0` in
`compare_to_cw.py`. First-48hr transient is hydration-rise scope
(deferred), not boundary physics scope. Window is permanent production
convention going forward; Sprint 5+ dual-peak hydration work may close
the transient but the window stays as a "boundary physics isolated"
gate.

**ADR-08 (Sprint 2 close, updated Sprint 3 + Sprint 4 PR 17)**: Ground temperature in form-face LW term:
- Sprint 2: defaulted to T_amb (no soil model)
- Sprint 3: upgraded to simplified Barber (sinusoidal lag + damping); defaults `soil_lag_hrs=5.0`, `soil_damping=0.7`
- **Sprint 4 PR 17**: defaults reverted to `soil_lag_hrs=0.0`, `soil_damping=1.0` (no-op pair; reduces to T_amb behavior). Helper code retained per R4 disposition. Empirical justification: §7.5.2 confirms 14-mix library has zero climate variation (all TX Austin); Sprint 3 sweep on MIX-01 showed lag dominates Corner RMS sensitivity (~0.03°F/hr) and damping is unidentifiable (<0.01°F authority). Without identifiability data across climates, defaults that introduce phase offset degrade fit monotonically. The helper remains a latent capability for future cold-climate exports where diurnal ground-air differentials become identifiable.
- Full Barber 1957 diffusion model available as a follow-on if future cold-climate data shows simplified Barber is insufficient.

**ADR-09 (Sprint 4 close, established by §7.5/§7.6 work)**: Mixes are
admitted to a thermal-validation evaluation set via an early-window
hydration-fit screen on the [12, 36]hr window — Centerline RMS, Δ peak
rise rate, and Δ time-to-half-peak, with thresholds derived from the
sprint's Reference set at 1.5× worst-case (§7.5.2 records the protocol
and the thresholds used in Sprint 4: 1.80°F, 0.20°F/hr, 2.13 hr). Mixes
that fail the screen are routed out of the active sprint's evaluation
set, not bundled into thermal-physics calibration work. The screen is a
methodology rule about *which mixes are admitted*, not a categorical
claim about *which physics terms are thermal vs hydration* — that
distinction is for diagnostics to determine on a finding-by-finding
basis (PR 16's Phase 2.6 Mode-A drift on a screen-passing mix is the
canonical example of why the rule is admission-shaped, not term-shaped).
Future thermal sprints reuse this screen; the Reference set carried
into a sprint may change, in which case the 1.5× worst-case thresholds
recompute against that sprint's Reference set.

**ADR-10 (Sprint 6 PR 22)**: `form_type` is normalized to lowercase at
the loader trust boundary (`cw_scenario_loader.py`, `parse_cw_dat`).
Engine-side lookup keys in `R_FORM_BY_FORM_TYPE` are lowercase. The
loader calls `.strip().lower()` on the raw CW input.dat string at parse
time, before the value is assigned to `CWConstruction.form_type`. The
dataclass default is also lowercase (`"steel"`).

Rationale: trust-boundary discipline — one canonical form decided once,
at the boundary between CW's text format and the engine's typed Python.
This prevents case-mismatch lookup failures as new form types are added.

Cites ADR-05's `F_VERT_BY_ORIENTATION` lookup as precedent (its keys
are lowercase strings: `"south"`, `"east"`, `"unknown"`, etc.). Note
distinction: ADR-05's lowercase came from the authored dict, not a
normalization step. PR 22 makes the lowercase convention explicit and
enforces it uniformly at the trust boundary so future loaders do not
need to remember to author lowercase keys.

---

## §6 Workflow discipline

Practices that held in Sprints 0–2 and carry forward, updated with
Sprint 2 lessons.

### §6.1 Narrow-PR discipline

Each PR has one job. PR names map to physics changes (or scaffolding),
not to file-edit scope. When a PR's change feels like "while we're in
there, also…", that's a signal to split.

### §6.2 Commit-before-reporting-done

Every PR prompt must include explicit `git add && git commit && git tag
&& git push` steps in the verification section. Claude Code's
"completion" notion does not include git bookkeeping unless it is asked
for. Sprint 2 PR 5's uncommitted-changes surprise was caught only by
manual Check 3; subsequent PRs baked this into the prompt.

After each push, verify clean state:
```
git status    # should show "working tree clean"
git log --oneline <prev-tag>..HEAD    # should show expected commit count
```

### §6.3 Diagnostic scripts before fixes

When a metric moves unexpectedly (either direction), write a 30-line
diagnostic script to decompose the effect BEFORE assuming which
hypothesis is right. Sprint 2 examples:

- `verify_pr6.py` — distinguished "RMS didn't move" from "amplitude
  actually closed, phase lag remaining"
- `verify_pr7.py` — distinguished flat-gate over-delivery (H1) from
  F_vert magnitude overshoot (H4)
- `verify_pr8_floor_v2.py` — distinguished top-BC bias, sampling
  mismatch, and transient contamination

Archived in `diagnostics/sprintN/` for future reuse. Pattern reuses
across sprints (window decomposition, sampling-location sensitivity,
BC-column uniformity).

### §6.4 Ablation tests within physics PRs

When a physics change adds a new term, include an ablation run in the
PR validation (not just the activated run). Examples:

- Sprint 2 PR 6: R_form=0 ablation confirmed R_form=0.0862 is real
  contact physics (Corner RMS 8.85°F vs 4.08°F)
- Sprint 2 PR 7: F_vert=0.0 ablation would have been useful to quantify
  solar contribution; added retroactively in PR 8's sweep
- Sprint 2 PR 8: F_vert sweep across 0.00–0.50 revealed calibration
  insensitivity (sweep range 0.36°F on Corner RMS)

Report ablation numbers in the PR commit message. They inform future
Claude Code sessions about which physics terms are load-bearing vs
marginal.

### §6.5 Verify field names before writing scripts

Sprint 2 had two drafting errors where field names were assumed from
memory rather than verified:
- PR 5 prompt: `T_C_history` / `alpha_history` vs actual `T_field_C` /
  `alpha_field`
- `verify_pr8_floor.py` v1: `cw.t_hrs` / `cw.corner_T_F` vs actual
  `val.time_hrs` / `val.T_field_F[:, 0, -1]`

Before writing any script that reads from a dataclass, `grep` the
dataclass definition for actual field names. 30 seconds of verification
saves an iteration cycle.

### §6.6 Gates as stop conditions, not aspirations

PR prompts specify concrete stop conditions with explicit off-ramps.
Example from Sprint 2 PR 8:
- Corner RMS ≤ 2.5°F: skip M6d (off-ramp)
- 2.5–3.0°F: run M6d as planned polish
- 3.0–3.5°F: run M6d with investigation flag
- > 3.5°F: stop, don't commit, regroup

Claude Code honored the off-ramp at PR 8 (3.96°F → did not commit). Gate
language in prompts must be literal, not aspirational.

### §6.7 Session discipline

- **Narrow PRs stay in one session** until the PR is merged
- **Sprint boundaries start a new session** — context shifts, file-touch
  patterns differ, diagnostic approaches change
- **Debugging sessions can carry context** across multiple exchanges,
  but production PRs fresh-start on the prompt rather than accumulating
- Claude Code session state is ephemeral; `git log` is the durable
  record

### §6.8 Verify tags exist locally before trusting Claude Code's report

Sprint 3 had two cases (PR 11, PR 12) where Claude Code reported a tag
created but `git tag -l` showed it absent locally. Reason traced to
session timing — tag created in the same session that pushed prior
work, but the tag step was skipped silently when the prompt's tag step
ran on already-tagged commits.

After every PR that includes a tag step, run:
git tag -l <expected-tag>          # should print the tag name
git ls-remote --tags origin <tag>  # should show same SHA

If the local check fails, the tag was never created — re-run the tag
step explicitly. Sprint 4 will tag pr-13-complete, pr-14-complete,
pr-15-complete, pr-16-complete, pr-17-complete, sprint-4-complete (six
tag sites; six places to bite).

### §6.9 Manual push and tag discipline (Sprint 5+)

Claude Code commits directly on `main` and creates local tags. Claude
Code does *not* push to remote. The user reviews the local commit
and any test/diff output Claude Code reports, then pushes manually.

Workflow per PR:

1. Claude Code: from `main` at the previous PR's tag, implement
   scope, run gates, commit on `main`, create local tag(s).
2. Claude Code: report commit SHA, tag name, gate results, anything
   notable in the diff. Stop.
3. User: review the commit locally (`git show <SHA>`, `git diff
   <prev-tag>..HEAD --stat`). Verify the change matches the prompt's
   scope.
4. User: if review passes, push: `git push origin main` and (if a
   tag was created) `git push origin <tag>`. If review fails, decide
   whether to amend (`git commit --amend`), reset (`git reset --hard
   <prev-tag>`), or retry; nothing reaches origin until review passes.
5. User: §6.8 verification (`git tag -l <tag>` AND
   `git ls-remote --tags origin <tag>`).

**Why**: the human review checkpoint between commit and push is a
safety gate. Claude Code's gate-passing is necessary but not sufficient
— the diff itself can contain scope creep, missed edge cases, or test
fixtures that look right but were generated against the wrong reference
(caught on PR 13's MIX-01 baseline JSON spot-review and again on PR
18's stale-baseline surfacing).

**Why commit-on-main, not branch+merge**: Sprint 4 used branch+merge
under §6.9's original draft; mid-Sprint-5 the convention shifted to
commit-on-main because the branch overhead added no review value
(reviews happen on the commit, not the branch). Commit-on-main keeps
the rollback path simple (`git reset --hard origin/main` until the
user pushes) without losing the human-review gate.

Claude Code prompts should reflect this:
- Specify "commit on `main` + create local tag, do **NOT** push" explicitly.
- Report SHA + tag name + gate output for user review.
- Skip the §6.8 verification step in the prompt itself (user runs
  it post-push).
- Do **NOT** create a branch.

### §6.10 Regenerate committed numerical artifacts when defaults change

When a PR changes an ADR-recorded default (or any engine parameter
that flows through to a committed validation file, test fixture, or
diagnostic markdown), regenerating every downstream artifact is part
of the PR's scope. The regeneration must be in the same commit as
the default change, not a follow-on PR.

Inventory of committed numerical artifacts that derive from engine
defaults (Sprint 5 snapshot — extend as new files are added):

- `validation/sprint4_baseline.md` — generated by `run_all.py
  --group all`, sensitive to soil_lag_hrs, soil_damping, F_vert*,
  R_FORM_CONTACT_SI, h_conv functional form, α_top, ε_top, and any
  ADR-recorded engine constant.
- `tests/fixtures/mix01_pr17_baseline.json` (and predecessor PR-tagged
  fixtures) — generated for bit-identical regression tests; sensitive
  to any engine default change.
- `validation/diagnostics/*/` — recon spike outputs; tied to the
  engine state at the time of the spike. Listed for awareness, NOT
  required to regenerate when defaults change (recon outputs are
  permanently associated with the spike's source-commit SHA, not
  with current `main`).

Workflow per default-changing PR:

1. Identify which committed numerical artifacts are sensitive to the
   default being changed.
2. Regenerate each. For `validation/sprint4_baseline.md`, that means
   `python run_all.py --group all --output-md
   validation/sprint4_baseline.md --quiet` from the new engine state.
3. Commit the regeneration in the same PR as the default change.
4. The PR's commit message lists which artifacts were regenerated.

If an artifact's derivation chain is unclear, default to regenerating
it. False-positive staleness wastes a verification cycle (PR 18's
case); false-negative staleness (committing a default change without
catching a downstream artifact) wastes a calibration cycle and may
contaminate later analysis.

Precedent: PR 17 (soil_lag_hrs 5.0→0.0, soil_damping 0.7→1.0) did not
regenerate `validation/sprint4_baseline.md`. The staleness surfaced
during PR 18's manual bit-identity verification on the Reference set
— PeakGrad shifted ~0.19°F and CornerRMS shifted 0.04–0.08°F across
all 5 Reference rows, with a form-face-localized signature consistent
with ADR-08's no-op-pair change. Resolved by a separate commit on
`main` regenerating the baseline against current engine state.

**Follow-on lesson (Sprint 5)**: regeneration must verify structural
contract, not just schema. The 9f275e4 baseline regeneration (the §6.10
mitigation for PR 17's stale baseline) dropped the MIX-13 skipped row
that PR 14's manually-curated baseline included; two tests asserting on
the row failed silently until PR 18's verification surfaced the
regression. Fix landed in 1eeb67d (`run_all.py` `GROUPS["all"]` → main()
append pattern).

The lesson: §6.10's "regenerate every committed numerical artifact" rule
is necessary but not sufficient — the regenerated file must also reproduce
the *structural contract* (row count, sentinel handling, column
conventions) of the previous version. A regeneration that produces the
right schema but drops a structurally-required row is a §6.10-class miss
caught only by downstream tests.

Future regeneration commits should: (1) regenerate per §6.10's primary
protocol, (2) `git diff` the regenerated file for structural changes (row
count, schema, sentinel rows), (3) re-run `pytest` to catch test-contract
regressions. The diff + pytest pair catches both the §6.10 primary case
(stale numerical values from default change) and the follow-on case
(structural drift from regeneration mechanism).

---

## §7 Sprint 3 retrospective

Sprint 3 closed 2026-04-24 at tag `sprint-3-complete` (commit 6b8a1a7).
Three PRs: PR 10 (plumbing, 318c400), PR 11 (activation + ablation,
d113e8f), PR 12 (S1-aspire recalibration, 6b8a1a7). All three landed
on scope per §2 plan.

### Deliverables

- Barber soil model activated at both form-face Newton sites (primary
  integration line ~1860, sample-time re-solve line ~2121)
- `ground_surface_temperature_C()` helper added to thermal_engine_2d
- `T_ground_C_history` diagnostic wired into HydrationResult
- `CWConstruction` gained `soil_lag_hrs` (5.0) and `soil_damping`
  (0.7) fields
- S1-aspire Corner RMS tightened from 3.0°F to 2.0°F in
  `compare_to_cw.py`
- Tests: 92 → 105 (13 net new; PR 10's bit-identical regression was
  intentionally retired in PR 11 as its "plumbing-only" contract is
  broken by design once the helper is wired)

### MIX-01 gate numbers at Sprint 3 close

| Metric | Sprint 2 | Sprint 3 | S0 gate | S1-aspire |
|---|---|---|---|---|
| Peak Max T | −0.3°F | −0.3°F | ±1.0°F ✓ | ±0.5°F ✓ |
| Peak Gradient | −0.3°F | −0.5°F | ±2.0°F ✓ | ±1.0°F ✓ |
| Field RMS | 0.88°F | 0.88°F | ≤2.0°F ✓ | ≤1.0°F ✓ |
| Centerline RMS | 0.74°F | 0.74°F | ≤1.0°F ✓ | ≤0.5°F ✗ |
| Corner RMS | 2.22°F | 2.26°F | ≤3.0°F ✓ | ≤2.0°F ✗ |

All five S0 gates remain PASS. S1-aspire: 3/5 met (was 4/5 at Sprint
2 close; Corner RMS moved out of S1-aspire when the threshold
tightened from 3.0°F to 2.0°F).

Corner RMS moved +0.04°F in PR 11 (2.22 → 2.26°F) — the "commit-
conditional" direction per §6.4 prompt gate language. Committed with
documented reasoning rather than reverted.

### Parameter sweep reconnaissance (not committed)

Manual sweep at sprint close, via python-anchor-inject pattern on
`compare_to_cw.py`:

| damping | lag | Corner RMS |
|---|---|---|
| 1.0 | 0.0 | 2.22°F (sprint-2 identity — validates ablation live) |
| 0.9 | 2.0 | 2.22°F |
| 0.7 | 5.0 | 2.26°F (PR 11 defaults) |
| 0.5 | 5.0 | 2.26°F (damping is a near-null knob) |
| 0.7 | 8.0 | 2.31°F |

Key findings:

1. **Lag dominates Corner RMS sensitivity on MIX-01.** Each ~3h of
   lag costs ~0.03°F. Monotonic with lag magnitude.
2. **Damping is effectively unidentifiable on MIX-01.** Sweep from
   0.5 → 1.0 at fixed lag=5 shows no measurable RMS change. On
   Austin-summer geometry, the Barber damping parameter has <0.01°F
   authority.
3. **CW's internal ground treatment appears phase-aligned with T_amb.**
   The sweep shape is consistent with CW using `T_ground ≡ T_amb` (or
   equivalent simple model). Our Barber defaults (5h/0.7) introduce a
   phase offset that degrades fit monotonically.
4. **If Sprint 4 confirms this pattern across climates**, the correct
   MIX-01-only calibration is `soil_lag_hrs = 0.0, soil_damping = 1.0`
   — making the Barber model a no-op by default and the helper a
   latent capability reserved for future cold-climate data where
   diurnal ground-air differentials are large enough to be
   identifiable.

### What worked (carrying to Sprint 4)

- **Bit-identical fixture regression** (PR 10) caught zero regressions,
  which is itself the contract-satisfaction proof. Pattern worth
  reusing whenever a "plumbing-only" claim needs teeth.
- **Ablation-in-activation PR** (PR 11) validated helper math via the
  degenerate-parameter identity (damping=1.0, lag=0.0 → reproduces
  sprint-2-complete bit-identically) before relying on aggregate
  gate numbers for correctness.
- **Python-anchor-inject for parameter sweeps** (grep-unique anchor
  string + `src.replace(anchor, override + anchor, 1)`) was more
  robust than shell `sed` once shell-side fragility was diagnosed.
  Template for future reconnaissance work.

### What didn't work (don't repeat)

- **Initial sweep script used `sed` + `set -e` + pipe-to-grep.** The
  chain was fragile on macOS; argparse errors and `RuntimeWarning`s
  both broke it silently. Manual edit-and-rerun was what actually
  produced data. For ≤6 data points, the manual path is faster than
  scripting.
- **Panel-(k) interpretation in PR 11 result review was directionally
  wrong.** Initial reading suggested Barber damping was over-radiating
  the form face. The sweep showed damping has no effect; the real
  mismatch is phase. Lesson: run the parameter sweep before narrating
  physics from a single parameter point.

---

## §7.5 Sprint 4 reconnaissance findings (pre-PR 13)

Two throwaway diagnostic spikes ran before Sprint 4 PR planning, per
"data-first" principle (see §7's "run the parameter sweep before
narrating physics from a single parameter point" lesson, generalized).
Spike artifacts in `diagnostics/sprint4/` — not committed to main.

### §7.5.1 Round 1: 14-mix S0 baseline

Subprocess-driven loop over `compare_to_cw.py` for MIX-01..15. MIX-13
skipped (no `output.txt`). Climate column degenerate — all 15 mixes are
TX, Austin (designed factorial, not natural sample).

| Mix | SCM% | S0 pass | Notable failures |
|---|---|---|---|
| MIX-01 | 39.1% | 5/5 | — |
| MIX-02 | 39.1% | 5/5 | — |
| MIX-03 | 39.1% | 5/5 | — |
| MIX-04 | 21.7% | 2/5 | PeakMax −1.5, PeakGrad −2.5, Corner 3.03 |
| MIX-05 | 52.2% | 4/5 | CenterRMS 1.49 |
| MIX-06 | 71.3% | 1/5 | PeakMax +2.2, PeakGrad +2.9, Field 2.35, Center 2.81 |
| MIX-07 | 91.3% | 1/5 | PeakMax +4.1, PeakGrad +4.8, Field 3.82, Center 4.45 |
| MIX-08 | 17.4% | 2/5 | PeakGrad −3.6, Corner 4.05 |
| MIX-09 | 41.7% | 1/5 | PeakMax +3.8, Field 4.15, Center 4.84, Corner 4.58 |
| MIX-10 | 52.2% | 4/5 | Center 1.03 |
| MIX-11 | 39.1% | 5/5 | — |
| MIX-12 | 39.1% | 5/5 | — |
| MIX-13 | — | skipped | no CW output |
| MIX-14 | 39.1% | 0/5 | PeakMax +1.5, PeakGrad −3.7, Field 3.20, Center 2.29, Corner 4.86 |
| MIX-15 | 39.1% | 1/5 | PeakMax −3.0, PeakGrad −8.3, Field 4.06, Center 2.07 |

5/14 pass all S0 (MIX-01, 02, 03, 11, 12 — all 39.1% SCM at 60°F
placement temp). Failure clustering by gate:

- **Centerline RMS** fails on 5/5 high-SCM mixes (MIX-05, 06, 07, 09,
  10). Monotone with SCM% on Peak Max (+0.7→+2.2→+4.1°F as SCM goes
  52→71→91%). Hydration mismatch signature.
- **Other failures** (MIX-04, 08, 14, 15) at mid/low SCM, mixed gate
  patterns. Not hydration-cluster-aligned.

### §7.5.2 Round 2: hydration screen + construction-parameter diff

Round 1 surfaced two distinct failure clusters. Round 2 partitioned them
empirically before PR planning.

**Stage 1 — hydration-fit screen.** Computed three early-window
[12, 36]hr metrics (Centerline RMS, Δ peak rise rate, Δ time-to-half-
peak) for the 5 Reference mixes and 4 Cluster B candidates (MIX-04, 08,
14, 15). Threshold = 1.5× Reference worst-case per metric.

| Threshold | Value |
|---|---|
| Centerline RMS [12, 36]hr | 1.80°F |
| Δ peak rise rate | 0.20°F/hr |
| Δ time-to-half-peak | 2.13 hr |

Results:
- MIX-04: hydration_pass (RMS 0.50, Δrise 0.04, Δt½ 0.58)
- MIX-08: hydration_pass (RMS 0.74, Δrise 0.03, Δt½ 0.42)
- MIX-15: hydration_pass (RMS 0.46, Δrise 0.13, Δt½ 1.83)
- MIX-14: hydration_fail (RMS 1.81, **0.01°F over threshold** —
  borderline; Sprint 5 routing per the rule but record borderline status
  for Sprint 5's reference)

**Stage 2 — construction-parameter diff.** Tabulated all 39 non-array
fields across `CWConstruction`, `CWGeometry`, `CWEnvironment` for
Reference (5) + hydration-pass survivors (3). Result: exactly **one
field varies** anywhere in the 8-mix slice.

| Field | Reference | MIX-04 | MIX-08 | MIX-15 |
|---|---|---|---|---|
| `construction.placement_temp_F` | 60 | 60 | 60 | **45** |

All other 38 fields are identical: 40×60×8ft footing, 2026/7/15 @ 5am,
Austin lat/lon, Steel/Red/Limestone construction, 0.65/0.88
absorptivity/emissivity, 5.67 R-blanket, 5.0/0.7 soil, identical hourly
weather arrays. The dataset is a designed factorial isolating mix-design
and (one case) placement temperature.

### §7.5.3 Sprint 4 evaluation-set partitioning (decided pre-PR 13)

Driven by §7.5 findings:

- **Reference set** (5 mixes — Sprint 4 calibration anchor):
  MIX-01, 02, 03, 11, 12. All 39.1% SCM, 60°F placement, S0 5/5.
- **B1 — composition-isolated thermal physics** (2 mixes):
  MIX-04 (21.7% SCM) and MIX-08 (17.4% SCM). Construction byte-identical
  to Reference; only mix-design differs. Hydration tracks in [12, 36]hr.
  Round-1 S0 failures: under-predict Peak Max (−1.5/−1.4°F), under-
  predict Peak Gradient (−2.5/−3.6°F), high Corner RMS (3.03/4.05°F).
  Same direction both mixes — real physics signature.
- **B2 — placement-temperature-isolated boundary physics** (1 mix):
  MIX-15. Same 39.1% SCM and construction as Reference except 45°F
  placement (vs 60°F). Hydration cleanest of any Cluster B candidate
  (RMS 0.46°F). Round-1 S0 failures: PeakMax −3.0, PeakGrad −8.3
  (off-charts vs Sprint 0–3 priors), Field 4.06.
- **Cluster A — Sprint 5 routed** (5 mixes — high-SCM hydration
  mismatch): MIX-05, 06, 07, 09, 10. Monotone Peak Max Δ with SCM%.
  Engine's single-peak α(t_e) cannot capture deferred kinetics of
  slag/high-fly-ash binders. Confirms R3 empirically.
- **Sprint 5 routed (borderline)** (1 mix): MIX-14. Failed hydration
  screen by 0.01°F. Threshold-touching; Sprint 5 should expect easy
  recovery once dual-peak hydration lands.
- **Deferred** (1 mix): MIX-13. No CW output; awaits re-export.

**Sprint 4 evaluation set: 8 mixes** (Reference 5 + B1 2 + B2 1).
Sprint 5 inherits 6 mixes (Cluster A 5 + MIX-14). MIX-13 deferred.

### §7.5.4 Pre-existing issues surfaced (not Sprint 4 scope)

Round 1 noted three quirks in `compare_to_cw.py`. Routed:

- **Hardcoded PNG path** at `compare_to_cw.py:437`
  (`cw_comparison_MIX-01.png`) overwritten on every run. **Fix in
  PR 13** as part of `run_one()` refactor.
- **Hardcoded MIX-01 peak-time constants** at `compare_to_cw.py:41` and
  `:43` (`CW_PEAK_MAX_T_HR`, `CW_PEAK_GRAD_T_HR`). Used as display
  values in print statements at lines 257/259. Causes silently-wrong
  CW peak-time labels on non-MIX-01 mixes (deltas and pass/fail values
  unaffected, computed from runtime CW data). **Fix in PR 13.**
- **`RuntimeWarning: invalid value encountered in divide`** at
  `thermal_engine_2d.py:1619` and `:1624` on first invocation
  (conductivity averaging). Tracked as R7 (§4), Sprint 6 cleanup.
  **NOT fixed in PR 13** — would break PR 13's bit-identical regression
  contract. Suppress at `run_all.py`'s stderr-parse layer if noisy.

---

## §7.6 PR 13 + 14 closeout — population characterization committed

PR 13 (`67a3ebb` → tag `pr-13-complete`) lifted `compare_to_cw.py:main()`
into `run_one()` (Approach A), added `run_all.py` multi-mix driver, and
fixed three §7.5.4 display issues (hardcoded MIX-01 peak-time constants,
hardcoded PNG path, lack of MIX-13 sentinel handling). Bit-identical
regression on MIX-01 PASS, full suite 115/115. PR 14 (this commit) ran
the full 14-mix evaluable set through `run_all.py` and committed the
result to `validation/sprint4_baseline.md`. Population characterization
is now in committed form; §7.5 reconnaissance is durably promoted to
production evidence.

### §7.6.1 Reconnaissance vs production reconciliation

Round-1 spike (subprocess + parser) → PR 14 in-process → numbers should
be identical at the engine level. Confirmed:

| Mix   | Spike S0 | PR 14 S0 | Largest metric drift |
|-------|----------|----------|----------------------|
| MIX-01 | 5/5     | 5/5      | 0.00 (same in-process engine path) |
| MIX-02 | 5/5     | 5/5      | 0.00 |
| MIX-03 | 5/5     | 5/5      | 0.00 |
| MIX-04 | 2/5     | 2/5      | 0.05 (PeakGrad: −2.45 vs spike −2.5; 1-decimal rounding) |
| MIX-05 | 4/5     | 4/5      | 0.00 |
| MIX-06 | 1/5     | 1/5      | 0.05 (PeakGrad: +2.85 vs spike +2.9; 1-decimal rounding) |
| MIX-07 | 1/5     | 1/5      | 0.03 (PeakGrad: +4.77 vs spike +4.8) |
| MIX-08 | 2/5     | 2/5      | 0.03 (PeakMax: −1.43 vs §7.5.3 ref −1.4) |
| MIX-09 | 1/5     | 1/5      | 0.02 (PeakMax: +3.78 vs spike +3.8) |
| MIX-10 | 4/5     | 4/5      | 0.00 |
| MIX-11 | 5/5     | 5/5      | 0.00 |
| MIX-12 | 5/5     | 5/5      | 0.00 |
| MIX-13 | skipped | skipped  | n/a (no CW output)   |
| MIX-14 | 0/5     | 0/5      | 0.01 (PeakMax: +1.49 vs spike +1.5) |
| MIX-15 | 1/5     | 1/5      | 0.05 (PeakGrad: −8.25 vs spike −8.3; 1-decimal rounding) |

(Largest drift is the maximum absolute difference across the 5 metrics —
reported in °F. §7.5.1 reported failure values rounded to 1 decimal
place; all non-zero drift values are attributable to that rounding.
Subprocess and in-process are byte-equivalent paths through the same
`run_one()` body, confirmed.)

### §7.6.2 Sprint 4 evaluation set: production-confirmed

8 mixes, partitioning per §7.5.3:

- **Reference (5/5 each)**: MIX-01, 02, 03, 11, 12. 39.1% SCM, 60°F
  placement.
- **B1 (2/5 each)**: MIX-04 (21.7% SCM), MIX-08 (17.4% SCM). 60°F
  placement, mid/low SCM. Hydration tracks.
- **B2 (1/5)**: MIX-15. 39.1% SCM, 45°F placement (only mix with cold
  placement). Hydration tracks.

5/8 evaluation set passes S0 today. 3/8 are PR 15/16 calibration
targets.

### §7.6.3 PR 15 hypothesis worksheet (B1 — composition-isolated thermal physics)

**PR 15 finding — diagnostic-only close, Decision D.**

PeakMax Δ is boundary-physics-invariant. Across 48 sweep points
(F_vert 0.00–0.50 × 3 mixes; R_form 0.040–0.200 × 3 mixes), MIX-04
PeakMax Δ held at −1.46°F and MIX-08 PeakMax Δ at −1.43°F at every
value without exception. Neither F_vert nor R_form moves the core peak
temperature. B1's dominant failures (PeakMax, PeakGrad at current
defaults) are hydration-heat-generation-driven, not boundary-physics-
driven. **Routed to Sprint 5 (dual-peak hydration work).** Full sweep
data: `validation/diagnostics/pr15_b1_sweeps.md`.

B1 failure is dual-source — core and corner are opposite-signed:

1. **Core too cold** (PeakMax −1.46°F / −1.43°F, PeakGrad −2.45°F /
   −3.61°F at current defaults) — hydration heat generation, Sprint 5
   scope. Immovable by boundary-physics calibration.
2. **Corner too hot** (CornerRMS 3.03°F / 4.05°F at current defaults)
   — engine over-predicts form-face corner temperature. Lower F_vert
   and lower R_form both reduce CornerRMS monotonically. The original
   §7.6.3 framing ("engine under-cools the form face") was incorrect
   at the corner; the corner is over-predicted, not under-predicted.

**R_form=0.060 secondary finding (diagnostic observation, not a
pre-committed Sprint 5 target):** At R_form=0.060, MIX-01 stays 5/5,
MIX-04 reaches 4/5 (PeakGrad=−0.47°F PASS, CornerRMS=1.29°F PASS;
only PeakMax=−1.46°F fails), MIX-08 reaches 4/5 (PeakGrad=−1.49°F
PASS, CornerRMS=1.79°F PASS; only PeakMax=−1.43°F fails). If Sprint 5
fixes the hydration model (PeakMax), an R_form recalibration from
0.0862 → 0.060 may then complete B1. Sprint 5 should evaluate this
jointly — not pre-commit it now.

**ADR-04 load-bearing note:** The R_form sweep shows MIX-01 is 5/5 at
both 0.060 and 0.0862, and drops to 4/5 at 0.100 (CornerRMS=3.11°F
exceeds the 3.0°F S0 threshold). ADR-04's value of 0.0862 sits at the
high end of MIX-01's 5/5 tolerance window — not wrong, but with limited
upward margin. This is load-bearing information for Sprint 5, not a
calibration claim.

Original hypotheses and their dispositions:

**H1 — F_vert calibration is composition-dependent.** FALSIFIED for
PeakMax (zero authority across 0.00–0.50 sweep). Confirmed for
CornerRMS (monotone: lower F_vert → lower CornerRMS for B1). Not a
viable single-knob fix because PeakMax is the binding constraint.

**H2 — R_form is composition-dependent.** FALSIFIED for PeakMax (zero
authority across 0.040–0.200 sweep). R_form has significant authority
over PeakGrad and CornerRMS. At R_form=0.060, B1 PeakGrad and
CornerRMS both pass. Not a viable single-knob fix because PeakMax is
the binding constraint.

**H3 — Some boundary heat-flux term scales with mix composition.**
MOOT at the boundary-physics level. The dominant residual (PeakMax
invariance) is not a catch-all boundary term — it is the hydration
heat generation itself, which is Sprint 5 scope.

Phase 1 gate checks: passed. MIX-01 holds 5/5 at default F_vert=0.15
and default R_form=0.0862. No stop condition triggered.

### §7.6.4 PR 16 hypothesis worksheet (B2 — placement-temperature-isolated boundary physics)

**PR 16 finding — diagnostic-only close, Decision G.**

**V2 (verification) — MIX-15 and MIX-01 share bit-identical mix designs.**
All 10 mix-design fields identical to 4+ decimal places; placement_temp_F
is the only varying field across all 38 construction parameters (§7.5.2
confirmed, now strengthened). The CW dataset is a purer factorial than
originally indicated.

**Phase 1 ablation — MARGINAL boundary-physics authority.** F_vert sweep
(9 values): MIX-15 PeakMax Δ range = 0.41°F (below the 0.50°F authority
threshold). R_form sweep (7 values): PeakMax Δ range = 0.08°F (invariant).
FieldRMS and CenterRMS: 0.00°F variation across both sweeps — completely
boundary-physics-invariant. MIX-01 held 5/5 at defaults (gate passed).
Full sweep data: `validation/diagnostics/pr16_b2_ablation.md`.

**H6 (two-tier) — cold-IC confirmed as trigger; blanket-pin ruled out.**
- H6a (warm IC override: placement_temp_F=60°F, compare to MIX-01 CW):
  S0 5/5, PeakMax Δ −0.29°F, FieldRMS 0.88°F, CenterRMS 0.74°F. Perfect
  Reference-level match. Cold IC confirmed as the causal trigger.
- H6b (cold IC, blanket pin patched to track T_amb, compare to MIX-15 CW):
  Results bit-identical to baseline (S0 1/5 unchanged). The blanket-cell
  pin (lines 2085, 2088) is physically inert — k_cell[grid.is_blanket]=0.0
  at line 1438 thermally decouples blanket cells from the concrete stencil.
  V3 counter confirmed 16,464 patch executions with zero concrete-metric
  change.
- **Blanket-cell pin is NOT a Sprint 6 cleanup item.** The air-cell pin
  at line 2085 is load-bearing (placeholder rho_cp=1; do not change).
  The blanket pin at line 2088 is cosmetically misleading but physically
  harmless. Do not conflate the two.
- H6 artifacts: `validation/diagnostics/pr16_h6_ic_test.md`.

**H4 (grep audit) — code-misuse hypothesis falsified.** Grep for
placement_temp_F / T_initial / T0_C across three production files found
zero non-trivial, non-blanket, non-air candidates in the concrete physics
path. T_initial_C is never re-read mid-simulation for concrete cells.
The only post-initialization use (line 1415, initial Cp seeding at α=0)
is overwritten on the first iteration. H4-class bugs ruled out.
Audit: `validation/diagnostics/pr16_h4_audit.md`.

**Phase 2.6 (equivalent-age trajectory comparison) — two-mode divergence.**
Engine vs CW centerline temperature trajectories (mid-depth) for MIX-01
and MIX-15 compared against CW-implied te/α (engine's Arrhenius model
applied to CW temperature trajectory):

- **Mode A** (both mixes): Engine runs +1.0–1.7°F warmer than CW during
  active hydration (t=8–48 hr). Small calibration drift in early heat
  generation. MIX-01 still passes S0 (PeakMax Δ −0.29°F). Not blocking.
- **Mode B** (MIX-15 only): After a crossover at t≈55–60 hr (concrete
  temperature ~88–93°F ≈ 31–34°C, just above T_REF=23°C), the engine
  progressively under-predicts CW. Δte crosses zero at t≈72 hr and
  reaches −7.6 hr by t=168 hr; ΔT reaches −3.04°F. Engine loses
  equivalent-age vs CW after the crossover, consistent with a divergence
  in the Arrhenius rate factor at concrete temperatures above ~30°C.

Mode B is distinct from Sprint 5's dual-peak hydration scope (which
targets the early-window two-peak structure for high-SCM mixes). Mode B
is a temperature-dependent rate calibration mismatch above ~30°C,
amplified in cold-placed mixes because their temperature curve spends
more time in the high-rate regime during the late-time peak.
**Routed to R9.**
Trajectory data: `validation/diagnostics/pr16_te_alpha_comparison.md`.

Original hypothesis dispositions:
**H4 — IC-dependent boundary term.** FALSIFIED. No placement_temp_F /
T_initial code-misuse found in the concrete physics path. See H4 audit.
**H5 — Cold-placement-specific physics not modeled.** PARTIALLY SUPPORTED.
Phase 2.6 Mode B divergence is consistent with H5's "calibration artifact
for cold-placed concrete" sub-hypothesis. Routed to R9 for Sprint 5/6.
**H6 — IC propagation issue.** CONFIRMED as causal (H6a); MECHANISM
REFINED by H6b (blanket is not the mechanism) and Phase 2.6 (Arrhenius
rate divergence above ~30°C is the root, not a single IC-pinned code term).

### §7.6.5 Sprint 4 close criteria — final disposition (PR 17)

After PR 17, all Sprint 4 PRs have landed. Final state:

- **Floor (achieved)**: PR 13 + 14 + 17 land. Reference 5/5 holds (re-verified post-PR-17). R4 cleaned. PR 15 + 16 land as diagnostic-only with traceable mechanisms (B1 → Sprint 5 dual-peak; B2 → R9 Arrhenius rate divergence above ~30°C).
- **Target (recalibrated and met)**: Originally "B1 + B2 both S0 PASS." Updated to "B1 routed via Sprint 5 hypothesis; B2 routed via R9 with traced mechanism; R4 cleaned." Both sub-targets achieved.
- **Stretch (deferred)**: S1-aspire targets for B1/B2 don't apply since neither was calibrated in Sprint 4. Deferred to Sprint 5/6 alongside the routed mechanisms.

Sprint 4 produced two genuinely new mechanism-level findings (B1 boundary-physics PeakMax invariance; R9 cold-IC Arrhenius rate divergence above ~30°C) that did not exist as hypotheses at sprint-3-complete. The structural shift from "calibration sprint" to "characterization sprint" is recorded in §8 Sprint 4 retrospective.

---

## §8 Sprint 4+ retrospectives

### §8.1 Sprint 4 retrospective

Sprint 4 closed 2026-04-25 at tag `sprint-4-complete` (commit 5596c92).
Five PRs landed across the sprint: PR 13 (run_one() refactor +
display fixes), PR 14 (baseline data commit), PR 15 (B1 ablation
sweeps), PR 16 (B2 ablation + IC test + Phase 2.6 trajectory
comparison), PR 17 (R4 disposition + sprint close).

The sprint planned as a "calibration sprint" (B1 + B2 close to S0)
and shipped as a "characterization sprint" (B1 + B2 mechanism-traced
and routed). The structural shift is the dominant story of the
sprint and the reason the close criteria were recalibrated mid-flight.

#### Deliverables

- `compare_to_cw.py` refactored: `main()` body lifted into
  `run_one(scenario_dir) -> dict`; `print_gate_table(result)` separated;
  hardcoded MIX-01 peak-time constants (lines 41/43) replaced with
  per-scenario extraction; hardcoded PNG path (line 437)
  parameterized.
- `run_all.py` multi-mix driver added (PR 13). Group definitions
  (`reference`, `b1`, `b2`, `cluster_a`, `evaluation_set`, `all`)
  static in `GROUPS`, not in dataclasses. PNG-off-by-default. MIX-13
  sentinel handling (`{"skipped": True, ...}`).
- `validation/sprint4_baseline.md` committed (PR 14): per-mix gate
  table for the full 14-mix evaluable set, partitioned by group.
- `validation/diagnostics/pr15_b1_sweeps.md`: F_vert × R_form sweep
  data for Reference + B1 (48 sweep points across 6 mixes).
- `validation/diagnostics/pr16_b2_ablation.md`: F_vert × R_form sweep
  data for B2.
- `validation/diagnostics/pr16_h6_ic_test.md`: cold-IC override and
  blanket-pin patch tests.
- `validation/diagnostics/pr16_h4_audit.md`: grep audit for IC code-
  misuse hypotheses across `thermal_engine_2d.py`, `cw_scenario_loader.py`,
  `compare_to_cw.py`.
- `validation/diagnostics/pr16_te_alpha_comparison.md`: equivalent-age
  trajectory comparison (engine vs CW) for MIX-01 and MIX-15 across
  168 hr, surfacing Mode A and Mode B divergences.
- ADR-08 updated (PR 17): `soil_lag_hrs=0.0`, `soil_damping=1.0`
  defaults committed. R4 closed.
- Tests: 105 → 119 → 122 across Sprint 4 (17 net new total). PR 13
  added the largest tranche covering sentinel handling and
  parameterized PNG output paths; PR 17 added 3 covering the
  no-op-pair contract (soil_lag_hrs default, soil_damping default,
  no-op pair behavior). The interim "181" figure recorded in earlier
  drafts was a misread; verified pytest baseline at sprint-4-complete
  is 122.

#### Sprint 4 close — gate numbers

Reference set (5/5 mixes hold S0 5/5 at sprint-4-complete):

| Mix | PeakMax | PeakGrad | FieldRMS | CenterRMS | CornerRMS | S0 |
|---|---|---|---|---|---|---|
| MIX-01 | −0.29°F | −0.45°F | 0.88°F | 0.74°F | 2.26°F | 5/5 |
| MIX-02 | (Ref) | (Ref) | (Ref) | 0.57°F | (Ref) | 5/5 |
| MIX-03 | (Ref) | (Ref) | (Ref) | (Ref) | (Ref) | 5/5 |
| MIX-11 | (Ref) | (Ref) | (Ref) | (Ref) | (Ref) | 5/5 |
| MIX-12 | (Ref) | (Ref) | (Ref) | (Ref) | (Ref) | 5/5 |

(Reference-set values matching MIX-01 to first decimal indicated as
"(Ref)"; full per-mix table lives in `validation/sprint4_baseline.md`.)

B1 (MIX-04, MIX-08): S0 2/5 each at sprint-4-complete (unchanged from
PR 14 baseline; PR 15 closed diagnostic-only — no engine change).
B2 (MIX-15): S0 1/5 at sprint-4-complete (unchanged from PR 14
baseline; PR 16 closed diagnostic-only — no engine change).

#### The structural shift: calibration sprint → characterization sprint

PR 15 was scoped as "B1 calibration: pick a knob that closes B1's S0
gates." PR 15's diagnostic Phase (per §6.3) ran F_vert and R_form
sweeps before committing to a fix. The sweep finding — **PeakMax Δ
held at −1.46°F (MIX-04) and −1.43°F (MIX-08) at every sweep point
across F_vert 0.00–0.50 and R_form 0.040–0.200, no exception** —
falsified both H1 (F_vert composition-dependent) and H2 (R_form
composition-dependent) for the binding constraint (PeakMax).

The honest read of that sweep: B1's PeakMax failure is not a
boundary-physics knob away. It is hydration-heat-generation-driven
and not closable inside Sprint 4's thermal-physics scope. PR 15
landed as Decision D (diagnostic-only close) — committed the sweep
data, did not commit a calibration change, routed B1 PeakMax to a
hydration sprint.

PR 16 followed the same pattern. Phase 1 ablation: F_vert sweep
range 0.41°F on PeakMax (below the 0.50°F authority threshold);
R_form sweep range 0.08°F (invariant). H6a (warm IC override at 60°F)
recovered Reference-level S0 (5/5, PeakMax Δ −0.29°F) — confirming
cold IC as the causal trigger but not the *mechanism*. H6b (blanket-pin
patch) was bit-identical to baseline — falsifying the proposed
mechanism. Phase 2.6's te/α trajectory comparison surfaced the actual
mechanism: Mode A small early calibration drift (≤1.7°F, both mixes,
non-blocking) and Mode B late-time Arrhenius rate divergence above
~30°C (cold-placed mixes only, MIX-15 −3.04°F). PR 16 landed as
Decision G (diagnostic-only close), R9 opened.

The sprint shifted, between PR 14 close and PR 15 mid-flight, from
"close 3 of 8 evaluation-set mixes to S0" to "characterize 3 of 8
evaluation-set mixes' failure modes well enough to route them to the
correct downstream sprint with a named mechanism." The recalibrated
close criteria (§7.6.5) reflect that pivot honestly. The pivot was
the right call: a calibration commit on B1 or B2 without first
falsifying the boundary-physics single-knob hypotheses would have
been the same anti-pattern §7's Sprint 3 retrospective warned about
("run the parameter sweep before narrating physics from a single
parameter point"), generalized from MIX-01 to a multi-mix evaluation
set.

#### New mechanism-level findings

Two findings that did not exist as hypotheses at sprint-3-complete:

1. **B1 PeakMax invariance to boundary physics.** Across 48 sweep
   points (F_vert × R_form, Reference + B1), MIX-04 and MIX-08 PeakMax
   Δ held to within ±0.005°F of −1.46/−1.43°F. The signature is
   composition-dependent (Reference 5 don't show it) and
   boundary-physics-invariant. Routed to the hydration sprint series;
   the engine's single-peak α(t_e) Arrhenius-weighted hydration model
   is the suspected source. Documented in §7.6.3.

2. **R9 — Arrhenius rate factor divergence above ~30°C.** Phase 2.6's
   equivalent-age comparison on MIX-15 shows a clean crossover at
   t≈55–60 hr (concrete temperature ≈31–34°C, just above T_REF=23°C),
   after which the engine progressively under-predicts CW. Δte reaches
   −7.6 hr by t=168 hr; ΔT reaches −3.04°F. Mode A (engine warmer than
   CW during active hydration, t=8–48 hr, both mixes) is a separate
   small calibration drift, not the same mechanism. R9 is opened with
   mechanism traced (Phase 2.6 finding), routing deferred pending
   kinetics-vs-inherited-calibration disambiguation. R9 is *not*
   "routed to the hydration series" — that conflates the two mechanism
   classes PR 16 worked to distinguish. Sprint 5 PR 18 D5 confirmed
   the Mode-A signature on 5/5 Reference mixes — positive bias +0.77
   to +1.21°F in [8,48] hr — with peak ΔT shifted from MIX-15's [12,24]
   hr to [39,47] hr, consistent with 60°F-placement hydration timing.
   The strict [12,24] hr time-of-max criterion was MIX-15-specific (cold
   placement) and does not generalize as-is; the placement-temp-adjusted
   signature matches. R9's Mode A is broader than PR 16 thought —
   present on warm-placement Reference too, just kinetics-timing-shifted.

These two findings are the durable scientific output of Sprint 4 —
more durable than the calibration commits the sprint was originally
scoped to deliver.

#### What worked (carrying to Sprint 5)

- **Ablation-first restructure after PR 15.** When PR 15's sweep
  falsified the single-knob hypothesis, the sprint did not chase a
  multi-knob fit. PR 15 closed diagnostic-only; PR 16 inherited the
  same ablation-first discipline (Phase 1 sweeps before Phase 2
  hypotheses). This is §6.3 generalized: not just "diagnostic before
  fix" but "diagnostic-strong-enough-to-route is itself a valid PR
  outcome." Sprint 5 PR 18 is shaped this way by design.

- **Two-phase wait points in PR 15 / PR 16.** Both PRs paused after
  Phase 1 ablation for a scope check before continuing into Phase 2
  hypotheses. PR 15's Phase 1 → Phase 2 wait point was the moment the
  calibration→characterization pivot was decided. PR 16's wait point
  caught H6a's "too clean" warm-IC result (5/5 with −0.29°F PeakMax —
  suspiciously identical to Reference) and triggered the H6b /
  Phase 2.6 follow-on rather than committing on H6a alone. Phase
  boundaries are good think-points; future PRs with diagnostic+fix
  scope should structure them this way.

- **Phase 2.5 verification battery in PR 16.** After H6a's clean
  result, the prompt added an explicit verification phase before
  proceeding: blanket-pin patch (H6b), code-misuse audit (H4), and
  trajectory comparison (Phase 2.6). The battery is what surfaced
  R9 — H6a alone would have routed B2 to hydration with the wrong
  mechanism. Sprint 5 PR 18 should structure a similar verification
  phase if its Phase 1 finding is "too clean."

- **Designed-factorial recognition (PR 14 → PR 15/16 framing).** The
  §7.5.2 finding that the 14-mix library varies on exactly one
  construction field (placement_temp_F: 60 vs 45 in 1/14 mixes)
  reframed the library from "natural sample" to "designed factorial
  isolating mix-design and one cold-placement case." That reframing
  is what made B1 and B2's separation clean enough for ablation
  conclusions to land. Future thermal sprints should run a
  construction-parameter diff before partitioning the evaluation set.

- **Group-static partitioning in `run_all.py`.** GROUPS lives in the
  driver, not in dataclasses; reverse-lookup `_MIX_TO_GROUP` is
  first-match across reference/b1/b2/cluster_a. That keeps Sprint 5's
  group rebalance (Reference-only evaluation set) a 5-line edit
  rather than a refactor.

#### What didn't work (don't repeat)

- **SSH-passphrase blocking push on all 5 PRs.** Every PR in Sprint 4
  hit the same passphrase-prompt-doesn't-render-in-sandbox issue at
  the push step. Manual SSH push (user shell, not Claude Code's
  shell) resolved it each time. The recurrence across 5 PRs without
  documenting the pattern was a workflow miss. **§6.9 codifies the
  fix as a positive practice**: Claude Code commits and tags locally,
  user pushes after review. The discipline is valuable independent
  of SSH and is the default from Sprint 5 onward.

- **PR 15's original "calibration" framing made the diagnostic-only
  close feel like a failure.** The sprint plan in §3 (pre-Sprint-4)
  said "PR 15 — B1 calibration." When PR 15 closed without a
  calibration commit, the prompt had to be re-read as "characterize-
  and-route" mid-PR. Sprint 5's PR 18 / PR 19 structure (recon spike
  PR + conditional-fix PR) is the corrected shape: scope the recon
  PR as recon, scope the fix PR as conditional, don't conflate them.

- **The §7.5.4 "RuntimeWarning suppress at run_all.py stderr-parse
  layer if noisy" guidance under-specified.** PR 13 ended up using
  `warnings.catch_warnings()` filter in `run_all.py` (not a stderr
  parse) because the warning surfaces as a Python `RuntimeWarning`
  before reaching stderr. The §7.5.4 wording set up a wrong mental
  model. R7 cleanup (Sprint 6) is removing the pre-guard divide
  pattern; the suppression in `run_all.py` is the interim.

#### Risk register movement

- **R1 (F_VERT_BY_ORIENTATION["unknown"]=0.15 MIX-01 calibrated only)**:
  Narrowed, not closed. PR 15's sweep showed F_vert has clean Corner
  RMS authority (monotone) but zero PeakMax authority for B1. R1's
  "may not generalize" concern is now scoped: F_vert generalization
  matters for Corner RMS calibration (Sprint 5 PR 20 in scope) and
  doesn't matter for PeakMax (B1's binding constraint, not F_vert's
  domain).
- **R2 (F_VERT_BY_ORIENTATION non-"unknown" stub values)**: Deferred to
  the multi-climate validation sprint series (still requires non-
  "unknown" CW exports). Status unchanged in §4 row.
- **R4 (Barber soil parameter literature picks)**: **Closed** in PR 17.
  Defaults reverted to no-op pair (`soil_lag_hrs=0.0`,
  `soil_damping=1.0`); helper code retained for future cold-climate
  data. ADR-08 updated.
- **R9 (Arrhenius rate factor divergence above ~30°C)**: **Opened** in
  PR 16 Phase 2.6. Mechanism traced. Routing deferred pending
  kinetics-vs-inherited-calibration disambiguation. Sprint 5 PR 18
  may strengthen the Mode-A overlap evidence on the Reference set
  but does not own R9's resolution.

R3, R5, R6, R7, R8 unchanged from sprint-3-complete state.

#### Sprint-series-closure framing (recorded for Sprint 5/6 planning)

Sprint 4's closure changes the shape of the remaining thermal series.
At sprint-3-complete, Sprint 5 was planned as "dual-peak hydration
model" inside the thermal series. Sprint 4's findings — particularly
the B1 PeakMax invariance and R9's rate-divergence-above-30°C signal
— make clear that the hydration work is its own sprint series, not a
single sprint. The thermal series ends at engine v3 (sprint-6-complete)
on a Reference-set baseline; the hydration series, geometry coverage
series, and multi-climate validation series are independent successors.
§3's "Explicit deferrals" subsection records the boundary so future
planning sessions don't pull the deferred work back into thermal-series
scope by default.

### §8.2 Sprint 5 retrospective

Sprint 5 closed 2026-04-25 at tag `sprint-5-complete` (commit [SHA]).
Four PRs landed across the sprint: PR 18 (Centerline recon spike), PR
19 (top-surface BC ablation, negative result), PR 20 (R_form Reference
evaluation + engine v3 release notes draft, conditional-commit deferred),
PR 21 (this commit — close + retrospective). Plus four small support
commits: baseline regeneration (9f275e4, post-PR-17 default-flow fix),
§8 retrospective + ADR-09 + §6.9 passdown amendments (884471f), MIX-02
mini-spike (8dc7131, disposition C), MIX-13 row fix (1eeb67d).

**Theme realization**: Sprint 5 planned as "Reference-set residual
characterization, fix-conditional" and shipped as that exactly, with
both fix-conditional outcomes coming back NEGATIVE. The sprint produced
no engine code commits; what it produced is more valuable — three
structural findings that close hypothesis space for engine v3:

1. **MIX-02 cross-diagnostic outlier signature.** PR 18 D1/D3/D4
   opposite-signed residuals + spike disposition C (kinetics-divergent)
   + PR 20 R_form opposite-signed response. Coherent across four
   independent diagnostics. Mechanism: kinetics divergence × form-face
   coupling = inverse residual response. Routes to hydration sprint series.

2. **Top-surface BC residual is structural, not scalar.** PR 19: none
   of h_conv, α_top, ε_top has authority on the cluster Centerline
   residual within physically reasonable ranges; D3 amplitude ratio
   moves *inversely* under h_conv scaling. Rules out any single-scalar
   top-surface BC knob, not just the three swept. Routes to engine v3
   investigation.

3. **R_form is fully form-face-localized.** PR 20 C3 trivially passed:
   CenterRMS Δ=0.000°F across 25 sweep runs. Strengthens ADR-04's
   "load-bearing physics not legacy calibration" claim. R_form=0.060
   doesn't generalize because of finding #1 (MIX-02 regresses). Steel
   value 0.0862 carries to Sprint 6.

#### Sprint 5 close — gate numbers

Reference set holds S0 5/5 unchanged from sprint-4-complete (no engine
commits). Values from the regenerated `validation/sprint4_baseline.md`
(post-PR-17-defaults, 9f275e4):

| Mix | PeakMax Δ | PeakGrad Δ | FieldRMS | CenterRMS | CornerRMS | S0 |
|---|---|---|---|---|---|---|
| MIX-01 | −0.29°F | −0.29°F | 0.88°F | 0.74°F | 2.22°F | 5/5 |
| MIX-02 | −0.68°F | +1.16°F | 1.23°F | 0.57°F | 1.24°F | 5/5 |
| MIX-03 | −0.30°F | −0.92°F | 0.87°F | 0.79°F | 2.72°F | 5/5 |
| MIX-11 | −0.13°F | +0.07°F | 0.91°F | 0.79°F | 2.12°F | 5/5 |
| MIX-12 | −0.31°F | −0.32°F | 0.87°F | 0.73°F | 2.23°F | 5/5 |

B1 S0 2/5, B2 S0 1/5 — unchanged, no engine commits in Sprint 5.

#### Deliverables

- `validation/diagnostics/pr18_centerline_recon/pr18_summary.md`:
  Centerline RMS recon across 5 Reference mixes. D1–D5 diagnostics.
  Three structural signals (constant offset, amplitude-ratio inversion,
  Mode-A kinetics timing shift).
- `validation/diagnostics/sprint5_mix02_recon/mix02_recon.md`:
  MIX-02 cross-diagnostic field diff (mini-spike). Disposition C:
  kinetics-driven outlier, not form-physics.
- `validation/diagnostics/pr19_top_bc_ablation/pr19_summary.md`:
  Top-surface BC ablation. Negative result — no scalar knob has
  authority on Centerline residual. s2 wait point correctly skipped.
- `validation/diagnostics/pr20_r_form_eval/pr20_r_form_summary.md`:
  R_form Reference-set evaluation. C3 Centerline trivially passed;
  conditional commit deferred (MIX-02 Corner regression at R_form=0.060
  blocked net improvement). R_form 0.0862 confirmed as Sprint 6 steel
  value.
- `docs/engine_v3_release_notes.md` (started, PR 20): Centerline
  residual section drafted. Sprint 6 completes the file.
- Four support commits: 9f275e4 (baseline regeneration), 884471f
  (passdown amendments), 8dc7131 (MIX-02 mini-spike), 1eeb67d
  (MIX-13 row fix).
- This commit (PR 21): §8.2 retrospective, §3 Sprint 6 refinement,
  §8.1 D5 framing fix, `run_all.py` required `--output-md`, §6.10
  follow-on lesson.

#### Floor / Target / Stretch outcome

- **Floor**: met. Reference S0 5/5 maintained. Three structural
  findings committed as durable characterization output. Engine v3
  release notes started.
- **Target**: half-met. Corner improvement *available* (cluster
  CornerRMS ~2.0°F reachable at R_form=0.060 for 4 of 5 Reference
  mixes) but not committed — MIX-02 regresses under R_form=0.060,
  blocking net improvement. Centerline ruled out by PR 19 (structural,
  not scalar). Target required both to close; neither closed.
- **Stretch**: infeasible. Engine v3 ships with documented residuals;
  the residuals were characterized, not closed.

#### What worked (carrying to Sprint 6)

- **Recon-spike-as-PR shape (PR 18).** Characterization-only PRs are
  first-class work, not failed calibration PRs. When diagnostics find
  nothing actionable inside-sprint-scope, landing the recon as a
  committed record is the correct close.
- **Two-phase wait points (PR 19).** The s1 → s2 wait point correctly
  skipped s2 when s1 found no magnitude optimum. Avoids burning scope
  on a Phase 2 investigation when Phase 1 has already falsified the
  hypothesis.
- **Mini-spike sub-PR pattern (MIX-02 recon).** When a recon surfaces
  an unexpected outlier, a small dedicated diagnostic before the next
  PR keeps scope tight without losing the finding. Routes into the
  permanent diagnostic record without inflating the parent PR.
- **Commit-on-main with manual-push review (§6.9).** Five review
  checkpoints across the sprint caught the stale-baseline issue, the
  MIX-13 dropped-row regression, and the h_conv monkey-patch scope
  question before push. Manual review pays.
- **Module-attribute monkey-patching as ablation tool (PR 19, PR 20).**
  `h_forced_convection` and `R_FORM_CONTACT_SI` patched via context
  manager; clean revert, no engine source modifications, uniform
  application across all engine sites. Pattern reuses across future
  ablation work.
- **§6.10's first real application.** PR 17 stale baseline → 9f275e4
  regeneration. The rule caught a downstream regression in 9f275e4
  (MIX-13 row dropped), which is the §6.10 follow-on lesson (Item 5).
  Self-correcting workflow in 2 steps.

#### What didn't work (don't repeat)

- **`run_all.py --output-md` default footgun.** Reproducing the §6.10
  regeneration footgun on every verification run across PRs 18–20 was
  a real ergonomic cost. Fixed in this commit: `--output-md` is now
  required, no default.
- **Loose prompt scoping on ablation breadth (PR 19).** The prompt
  scoped "h_conv top" while prescribing a function-level patch that
  affects both top and side. Ambiguity surfaced as a mid-execution
  question; better-scoped prompt language would have avoided the wait.
  Future ablation prompts should explicitly characterize patch breadth.
- **Durable numbers need verification before drafting.** The `181 →
  122` test count miscount in §8.1's first draft. Even at sprint-close
  time, durable numbers require `pytest` verification before writing.
  Carried from Sprint 4 as a lesson; landed on the wrong side once.

#### Risk register movement

- **R3** (first-48hr hydration shape divergence): unchanged —
  accepted/deferred to hydration series.
- **R5** (R_form form-material parameterization): **scope-narrowed**.
  PR 20 confirmed R_form's authority is form-face-localized AND that
  R_form recalibration on Reference is blocked by MIX-02's kinetics-
  driven contradiction. R5 in Sprint 6 stays scoped to form-material
  parameterization; explicitly does NOT attempt to fix the cluster vs
  MIX-02 split (hydration-series scope). Steel value resolved: 0.0862.
- **R7** (RuntimeWarning): unchanged — Sprint 6.
- **R8** (CW parameter inheritance): **strengthened**. MIX-02 mini-
  spike disposition C is canonical evidence that CW's hydration
  regression consumes inputs not in our `CWMixDesign` (likely cement
  Bogue compounds or admixture chemistry). Inherited by design; not a
  fix item.
- **R9** (Arrhenius rate factor divergence above ~30°C): **strengthened**.
  PR 18 D5 found Mode-A overlap on 5/5 warm-placed Reference mixes
  (placement-temp-shifted time-of-max consistent with kinetics timing).
  PR 19 found D1 DC offset and D5 Mode-A offset cannot be cleanly
  disambiguated until the amplitude residual closes. R9 routing remains
  deferred but now has evidence on both branches.
- **R10** (artifact regeneration): **first application + follow-on
  lesson recorded**. §6.10 rule caught PR 17's stale baseline; 9f275e4
  regeneration fixed it; 9f275e4 dropped MIX-13 row, caught by tests;
  1eeb67d fixed the drop. Follow-on lesson appended to §6.10 in this
  commit.
- **No new risks opened in Sprint 5.** Three structural findings, all
  routing to existing risks (R3/R8/R9) or existing sprint series. The
  thermal sprint series is genuinely converging.

#### Sprint-series-closure framing

Sprint 5 confirms the framework set in §8.1 — the thermal sprint series
ends at engine v3 (`sprint-6-complete`). What Sprint 5 added: the
structural findings explicitly close the hypothesis space *inside* the
thermal series. Top-surface BC authority is ruled out. R_form is
localized. MIX-02's contradiction is routed to the hydration series.
The remaining work for engine v3 is cleanup (R5, R7) and release notes
finalization, not further characterization. Sprint 6 is the last sprint.

---

## Appendix A: Quick reference — Sprint 2 close state

Tag: `sprint-2-complete`
Tests: 92 pass, 0 xfail, 0 xpass
Validation (MIX-01, [48,168]hr window):

| Metric | Result | S0 Gate | S1-aspire | S1 status |
|---|---|---|---|---|
| Peak Max T | −0.3°F | ±1.0°F ✓ | ±0.5°F | ✓ |
| Peak Gradient | −0.3°F | ±2.0°F ✓ | ±1.0°F | ✓ |
| Field RMS | 0.88°F | ≤2.0°F ✓ | ≤1.0°F | ✓ |
| Centerline RMS | 0.74°F | ≤1.0°F ✓ | ≤0.5°F | ✗ |
| Corner RMS | 2.22°F | ≤3.0°F ✓ | ≤3.0°F | ✓ (gate = S0) |

Sprint 3 recalibrates Corner RMS S1-aspire to ≤2.0°F (tightened from
3.0°F, creates pull).
