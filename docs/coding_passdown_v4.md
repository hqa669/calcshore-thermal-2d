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

### Sprint 5 — Dual-peak hydration model (~4–6 weeks)

The deferred first-48hr shape divergence identified in Sprint 2. Current
single-peak α(t_e) Arrhenius-weighted model cannot capture CW's likely
two-peak profile (initial wetting peak + acceleration peak).

Scope candidates:
- Expand hydration model parameterization to support two peaks
- Implement dormant-period onset delay
- Rate-function selection (Schindler vs exponential vs other)
- Re-validate against MIX-01 first, then the 15-mix library

Blocker from Sprint 3 perspective: we need a validated multi-mix
baseline (Sprint 4 deliverable) before investing in a hydration model
upgrade, else we risk overfitting to MIX-01 again.

### Sprint 6 — Production robustness (~3–4 weeks)

- Suppress cosmetic RuntimeWarning at harmonic-mean k-divide
- Parameterize R_form by form material (steel/plywood/plastic liner)
  via new CWConstruction.form_type fields
- Form orientation parsing from CW .dat files if CW stores it
- Error handling and input validation for customer-facing use
- Performance optimization if runtime becomes concern
- API cleanup and customer documentation
- Engine v3 release

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
| R4 | Barber soil parameters (lag, damping) chosen from literature; may not match CW's internal soil model | Sprint 4 PR 17 | Disposition decided: §7.5 found all 8 evaluation-set mixes share soil_lag_hrs=5.0, soil_damping=0.7 in CW exports; combined with Sprint 3 sweep showing no climate variation in library (all TX Austin), parameters are unidentifiable on this dataset. PR 17 sets soil_lag_hrs default → 0.0, deprecates soil_damping (default → 1.0), retains code for future cold-climate data. Updates ADR-08. |
| R5 | R_FORM_CONTACT_SI=0.0862 hardcoded assuming steel form; customer pilots may use plywood or plastic-lined forms | Sprint 6 | Open |
| R6 | Engine tested on half-mat geometry only; behavior on full-mat, slab, or column geometries unvalidated | Sprint 6+ | Open |
| R7 | Cosmetic RuntimeWarning at harmonic-mean k-divide masks potential future numerical issues | Sprint 6 | Mitigated via np.where guard; cleanup deferred |
| R8 | Hydration model uses CW's parameters (τ, β, α_u, Ea, Hu) directly; if those parameters are miscalibrated in CW, our engine inherits the error | Sprint 4–5 | Accepted (inherited by design) |

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

**ADR-08 (Sprint 2 close)**: Ground temperature in form-face LW term
defaults to T_amb (no soil model) as of Sprint 2. Sprint 3 upgrades to
simplified Barber (sinusoidal lag + damping). Full Barber 1957 diffusion
model available as a follow-on if simplified doesn't match CW.

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

## §8 Sprint 4+ retrospectives

(Populated as sprints close.)

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
