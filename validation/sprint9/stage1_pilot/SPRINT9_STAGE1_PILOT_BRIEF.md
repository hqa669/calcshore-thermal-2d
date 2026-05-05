# Sprint 9 Stage 1-pilot — Realistic Mix Validation (Pilot)

**Sprint:** 9 (production calibration validation)
**Stage:** 1-pilot (2-mix go/no-go before full 15-mix validation)
**Type:** Engine validation against CW with realistic mix composition. Wrapper-only — no engine source changes.
**Predecessor:** Sprint 8 closed with two empirical corrections (Hu_residual conditional, c(T_pl) quadratic).
Sprint 9 begins by validating these corrections against realistic production mixes.

---

## 1. Context

Sprint 8 closed with the engine matching CW within the Sprint 7 Structure C gate (R1, R2 ≤ 0.35°F)
across A and F scenarios at α_u ∈ [0.036, 0.80] and T_pl ∈ [40, 110]°F under suppressed-hydration
test conditions (Hu = 1 J/kg). Two corrections apply via wrapper:

1. **Hu_residual conditional:** if user-Hu < 10,000 J/kg → use 12,937 J/kg; else pass through user-Hu.
   (CW floor compensation; only fires when user is in CW's hidden floor regime.)
2. **c(T_pl) quadratic:** `c(T_pl) = -1.3025 + 0.04746·T_pl - 2.081×10⁻⁴·T_pl²` (T_pl in °F).
   Sprint 8 calibrated this under suppressed Hu, where the only active α-dependent paths in the
   engine were k(α) and Cp(α). Therefore c() is interpreted as a shape correction on the α-axis
   used by k(α) and Cp(α), NOT as a scaling of total hydration heat magnitude. Sprint 9 applies it
   accordingly via the wrapper (see §2.2).

These corrections were calibrated under suppressed Hu. Sprint 9 begins with the question: do they
generalize to realistic mixes with full hydration heat release?

This pilot tests 2 mixes (MIX-01 and MIX-07 from a 1–15 production mix set). MIX-01 is the moderate
baseline (Type I/II + 22%FA + 17%slag, τ=29.4 hr, β=0.895, α_u=0.7585, Hu=424,143 J/kg). MIX-07 is
the kinetic-shape extreme (Type I/II + 22%FA + 70%slag, τ=75.1 hr, β=0.516, α_u=0.8935, Hu=463,076
J/kg, E_a=33,523 J/mol). The two mixes bracket the kinetic envelope of the production set.

If both pass, the user will generate the remaining 13 datasets and a follow-up brief will run the
full validation. If either fails, this pilot stops and surfaces the failure mode for diagnosis
before more CW datasets are generated.

---

## 2. Pilot scope

### 2.1 Datasets

User-provided CW datasets at `~/Downloads/pilot_test_mix01,07/`:

- `pilot_test_mix01_73_85/` — MIX-01
- `pilot_test_mix07_73_85/` — MIX-07

Each dataset has the same scenario configuration:

- T_pl = 73°F
- T_soil = 85°F
- T_ambient = 85°F (= T_soil, weather override discipline preserved)
- Geometry: same as Sprint 7/8 baseline (40 ft wide × 80 ft deep half-symmetric submerged)
- Realistic mix kinetic parameters (Hu, α_u, τ, β, E_a) hand-set in input.dat to match the
  user's mix design table. Lines 389–390 of input.dat contain α_u and Hu respectively.

### 2.2 Engine wrapper configuration (B1 — Hu inverse-compensation)

Apply the layered corrections as follows:

```python
# Step 1: read raw values from CW input.dat
alpha_u_raw = parse_input_dat_line_389(...)    # dimensionless
Hu_raw      = parse_input_dat_line_390(...)    # J/kg

# Step 2: Hu_residual conditional (Sprint 8 floor compensation)
HU_FLOOR_THRESHOLD = 10_000.0   # J/kg
HU_RESIDUAL        = 12_937.0   # J/kg
if Hu_raw < HU_FLOOR_THRESHOLD:
    Hu_post_residual = HU_RESIDUAL
else:
    Hu_post_residual = Hu_raw   # passthrough — fires for both pilot mixes

# Step 3: c(T_pl) shape correction
T_pl_F = 73.0
c = -1.3025 + 0.04746 * T_pl_F - 2.081e-4 * T_pl_F**2   # ≈ 1.0532 at T_pl=73

# Step 4: apply c to α_u; inverse-compensate Hu to preserve heat magnitude
mix.alpha_u_effective    = c * alpha_u_raw
mix.Hu_J_kg_effective    = Hu_post_residual / c

# Step 5: pass τ, β, E_a through unmodified (read from input.dat)
# Step 6: k_uc × 0.96 already in engine source (Sprint 7); no wrapper override
```

**Why B1:** Sprint 8 calibrated c() under suppressed Hu, where total heat = Hu × Cc × α_final ≈ 0.
The only place α_u influenced the residual was through k(α) and Cp(α) shape. Applying c to α_u
without compensating Hu would inject ~5.3% extra hydration heat at realistic Hu — a physics change
the calibration data cannot support. Inverse-compensating Hu (Hu_eff = Hu_raw / c) preserves
Q(t) = Hu × Cc × dα/dt at every timestep (the c factor cancels exactly because Q is linear in α_u),
while still routing the c-scaled α through k(α) and Cp(α) as Sprint 8 calibrated.

Other engine settings unchanged from Sprint 7/8 closure:
- `model_soil = False`
- `is_submerged = True`
- `blanket_thickness_m = 0.0`
- 6× grid refinement for concrete domain

Numerical values for the pilot (for verification):

| Mix | α_u_raw | α_u_eff | Hu_raw (J/kg) | Hu_eff (J/kg) |
|---|---|---|---|---|
| MIX-01 | 0.7585 | 0.7989 | 424,143 | 402,719 |
| MIX-07 | 0.8935 | 0.9410 | 463,076 | 439,683 |

c(73) = 1.0532. Hu_residual conditional: passthrough for both.

### 2.3 Residual metric

**Region:** vertical centerline only, from concrete core (di=24, mid-depth) to bottom surface
(di=48), at wi=0. Bottom-half only — surface region (di < 24) explicitly excluded.

**Metric:** `max|R| = max_{di ∈ [24, 48]} |T_engine(di, wi=0, t=168) − T_CW(di, wi=0, t=168)|`

**Gate:** max|R| < 0.5°F per dataset.

The metric is restricted to the bottom half intentionally — to focus on the soil-concrete BC region
(where production thermal cracking analysis is most sensitive) and to avoid contamination from
the top-surface ambient handling. The metric does not include the side profile (R1) or the corner
(R3).

---

## 3. Tasks

### 3.1 Copy and verify CW datasets

Copy from `~/Downloads/pilot_test_mix01,07/` into `validation/sprint9/stage1_pilot/cw_data/`:

- `pilot_test_mix01_73_85/` → `cw_data/mix01/`
- `pilot_test_mix07_73_85/` → `cw_data/mix07/`

Verify file sizes (output.txt should be ~7-8 MB each).

For each input.dat, extract and tabulate:

| Mix | T_pl °F | T_soil °F | T_ambient | α_u (line 389) | Hu (line 390, J/kg) | τ (hr) | β | E_a (J/mol) | Cement (lb/yd³) | w/cm |

Confirm:
- T_pl = 73°F, T_soil = T_ambient = 85°F for both mixes
- MIX-01: α_u=0.7585, Hu=424,143, τ=29.4, β=0.895, E_a=26,458, cement=350, w/cm=0.44
- MIX-07: α_u=0.8935, Hu=463,076, τ=75.1, β=0.516, E_a=33,523, cement=50, w/cm=0.44
- Geometry matches Sprint 7/8 baseline (40 × 80 ft submerged)
- ambient/weather override is in place (lines 519-531 should show flat 85°F climate)

If any inconsistency or unexpected value: flag and pause for user clarification before running.

### 3.2 Extract CW vertical centerline trajectory

For each dataset, extract `T_CW(di, wi=0, t)` for di ∈ [0, 48] at hourly cadence over t ∈ [0, 168 hr]
using the existing Sprint 7 parser (`parse_cw_temp_output`). Apply the 0.02 m depth-origin offset
from Sprint 7 Stage 1 (handled inside the parser; do not double-apply).

Save the full trajectory `T_CW(di, t)` for downstream diagnostic use.

### 3.3 Engine run with all corrections

Run the engine on each dataset with the wrapper config from §2.2.

**Pre-run wrapper log (must appear in stdout for each mix):**
Mix: MIX-01
Raw from input.dat: alpha_u=0.7585, Hu=424143 J/kg, tau=29.4 hr, beta=0.895, Ea=26458 J/mol
Hu_residual: passthrough (Hu_raw=424143 ≥ 10000)
c(T_pl=73) = 1.0532
alpha_u_effective = 1.0532 × 0.7585 = 0.7989
Hu_J_kg_effective = 424143 / 1.0532 = 402719
Engine settings: model_soil=False, is_submerged=True, blanket=0.0, k_uc×0.96

If any computed value differs by more than 0.1% from the table in §2.2, halt and flag.

Capture full T(z, x, t) field at every hour from t=0 to t=168.

For each run, extract `T_engine(di, wi=0, t)` for di ∈ [0, 48] at hourly cadence.

**Sanity checks (perform after each engine run; halt if any fails):**

- `T_engine(di ∈ [24,48], wi=0, t=0)` must equal 73.0°F to within ±0.05°F. Else IC is broken.
- `T_engine(di=48, wi=0, t=168)` must equal 85.0°F to within ±0.5°F. Else bottom Dirichlet broken.
- `T_engine(di=24, wi=0, t=168)`: expect 115–135°F for MIX-01; expect 110–130°F for MIX-07.
  If T_core ≈ 73°F (no rise), hydration heat path didn't fire — wrapper broken. Halt.
- `α(di=24, t=168)`: expect 0.60–0.78 for MIX-01; expect 0.40–0.60 for MIX-07.
  If MIX-07 α(168) > 0.80, kinetics are running too fast — likely τ or β unit error. Halt.

Log all four sanity-check values for each mix to the report.

### 3.4 Compute residual metric

For each dataset:

- `R(di) = T_engine(di, wi=0, t=168) − T_CW(di, wi=0, t=168)` for di ∈ [24, 48]
- `max|R|` over that range
- `di_at_max` (the di index where the maximum occurs)
- Gate verdict: PASS if max|R| < 0.5°F, FAIL otherwise

Tabulate:

| Mix | max\|R\| (°F) | di_at_max | Pass? |
|---|---|---|---|

### 3.5 Decision rule

**Outcome 1 — Both mixes pass (max|R| < 0.5°F).** Run §3.6 spatial diagnostic on BOTH mixes
(diagnostic plots are produced for all outcomes). Report results and recommend proceeding to full
15-mix validation. **Pass margin and shape matter:** distinguish "clean pass" (max|R| < 0.30°F,
residual scattered) from "marginal pass" (max|R| ≥ 0.40°F or residual concentrated at di=48 or
di=24) in §1.4.

**Outcome 2 — Both mixes fail.** Run §3.6 spatial diagnostic on both. Report and stop.

**Outcome 3 — One passes, one fails.** Run §3.6 on both. Report and stop.

In all three cases, do not proceed to additional engine runs, calibration tuning, or correction
modifications. The pilot reports outcome and stops; the user decides next steps.

### 3.6 Spatial diagnostic (run for both mixes regardless of outcome)

For each mix, generate three plots:

**Plot 1 — Centerline profile snapshots.** T_engine(di, wi=0) vs T_CW(di, wi=0) at four timesteps:
t = 0, 24, 84, 168 hr. Two lines per timestep (CW solid, engine dashed). 8 lines total per mix.
**di range: [24, 48]** (bottom-half only). di on x-axis, T on y-axis.

**Plot 2 — Residual depth profile at t=168.** R(di) = T_engine − T_CW at t=168 hr,
for **di ∈ [24, 48]**. Single line per mix. Shows where the residual is concentrated.

**Plot 3 — Residual time evolution at three points.** R(di=36, t), R(di=42, t), R(di=48, t)
for t ∈ [0, 168] at hourly cadence. Three lines per mix. Continuous trace, not discrete snapshots.
Shows whether the residual builds monotonically, peaks and decays, or fluctuates.

**Interpretation guide for the report (under B1 wrapper):**

Under B1, heat source magnitude Q(t) = Hu_eff × Cc × dα_eff/dt = Hu_raw × Cc × dα_raw/dt at every
timestep (the c factor cancels exactly because Q is linear in α_u). Therefore residuals attributable
to heat magnitude should NOT appear under B1; residuals should localize to where k(α) or Cp(α)
shape errors propagate, plus the structural effects below.

- **Residual concentrated near di=48 (bottom surface):** likely the bottom-side BC stencil asymmetry
  documented in Sprint 8 §4.11.8 and Sprint 7 §5 (deferred work). The top-side corner uses a
  quarter-cell + half-cell stencil; the bottom-side uses pure strong Dirichlet write. Magnitude
  scales with α_u and bottom thermal flux. Sprint 8 quantified this at ~0.49°F at α_u=0.80 in
  I-scenario. Under B1, this mechanism is unchanged from prior sprints.

- **Residual concentrated near di=24 (core):** kinetics or k(α)/Cp(α) shape mismatch. Under B1,
  this is the primary signal for c() shape extrapolation failure — heat magnitude is preserved,
  but the c-scaled α produces different k and Cp values than the raw α would. If this fires only
  for MIX-07 and not MIX-01, the failure mechanism is c() not generalizing across kinetic-shape
  regimes (β=0.516 vs β=0.895). If it fires for both, c() is failing at realistic Hu generally.

- **Residual roughly uniform across di ∈ [24, 48]:** bulk mechanism such as ρ × Cp mismatch in
  baseline (not c-induced), or k_uc × 0.96 extrapolation failing at α > 0.04 (Sprint 7 calibrated
  k_uc only at suppressed α; full validation across α ∈ [0.04, 0.80] was deferred).

- **Residual peaks and decays in time:** transient kinetics shape mismatch (engine and CW reach
  the same final α but along different paths). Under B1 with τ, β, E_a passed through unmodified,
  this shouldn't appear unless engine and CW solve the kinetics ODE differently.

- **Residual builds monotonically:** persistent disagreement that doesn't equilibrate by t=168.
  Suggests a structural error (BC, stencil, or thermal property formula) rather than a transient.

**Pre-existing structural conditions to acknowledge in the report:**

- t=24 R values up to ~0.7°F at depth are expected from the boundary-onset convention difference
  (CW applies BC instantaneously at t=0; engine ramps gradually). Sprint 7 §2.2 documented this
  under k×0.96 at suppressed Hu. With realistic Hu, the t=24 transient may differ in magnitude
  but the underlying mechanism is the same. Do not flag t=24 magnitudes as a Sprint 9 finding.

- At t=168, the bottom BC stencil asymmetry contributes ~0.49°F at α_u≈0.80 per Sprint 8 §4.11.8.
  For MIX-01 (α_u_eff≈0.80) this is the relevant scale; for MIX-07 (α_u_eff≈0.94 but α(168)≈0.5)
  the effect may be smaller because hydration is incomplete by t=168.

### 3.7 Reporting

Write `validation/sprint9/stage1_pilot/PILOT_REPORT.md` with sections:

1. **§1.1** Dataset inventory (§3.1 verification table)
2. **§1.2** Engine run setup confirmation: wrapper log for each mix (§3.3 pre-run log), sanity-check
   results (T(t=0), T(di=48,t=168), T_core(di=24,t=168), α(di=24,t=168) for each mix), confirmation
   of which corrections fired and which passed through
3. **§1.3** Residual table (§3.4)
4. **§1.4** Outcome classification (1, 2, or 3 from §3.5), with margin and shape assessment for
   each mix (clean pass / marginal pass / fail; residual concentration region if applicable)
5. **§1.5** Spatial diagnostic plots (§3.6) — 6 plots total (3 per mix)
6. **§1.6** Synthesis paragraph — factual only. Report what residuals are, where they concentrate,
   and which interpretation guide categories from §3.6 they fall under. Do not propose fixes.

Save plots to `validation/sprint9/stage1_pilot/figures/`.

CSV outputs in `data/`:
- `pilot_residual_table.csv` (mix, max|R|, di_at_max, pass/fail)
- `pilot_centerline_profiles.csv` (full T(di, t) at t = 0, 24, 84, 168 for di ∈ [24,48], both mixes,
  both engine and CW)
- `pilot_residual_timetrace.csv` (R(di=36/42/48, t) hourly for both mixes)
- `pilot_dataset_inventory.csv`
- `pilot_wrapper_log.csv` (one row per mix: raw values, effective values, c value, conditional
  outcome, sanity-check results)

Scripts (uncommitted):
- `scripts/run_engine_pilot.py`
- `scripts/analyze_pilot.py`

---

## 4. Constraints

- Same engine state (k_uc × 0.96 from Sprint 7). No engine source changes.
- Hu_residual conditional and c(T_pl) shape correction via wrapper only, applied as B1 (see §2.2).
- τ, β, E_a passed through unmodified from input.dat.
- No new CW runs (the user provides them; the agent does not generate CW datasets).
- No commits to main. All work under `validation/sprint9/stage1_pilot/`.
- No closure recommendation, no calibration tuning, no correction modifications.
- No proceeding to full 15-mix validation in this stage even if pilot passes — that's a separate
  brief with separate dataset preparation.
- Synthesis is factual: report what the residuals are and where they concentrate, do not propose
  fixes.

---

## 5. Sanity checks (consolidated; halt if any fails)

Pre-run (after wrapper config, before engine starts):
- c(73) ≈ 1.0532. If far from this, formula mis-implemented.
- Hu_residual: passthrough for both mixes (Hu_raw >> 10,000). If Hu_eff = 12,937, conditional reversed.
- α_u_eff and Hu_eff match the §2.2 table to within 0.1%. Else wrapper broken.

Post-run (after each engine simulation, before residual analysis):
- T_engine(di ∈ [24,48], t=0) = 73.0 ± 0.05°F. Else IC parsing broken.
- T_engine(di=48, t=168) = 85.0 ± 0.5°F. Else bottom Dirichlet broken.
- T_engine(di=24, t=168): MIX-01 ∈ [115, 135]°F; MIX-07 ∈ [110, 130]°F. Else heat path broken.
- α(di=24, t=168): MIX-01 ∈ [0.60, 0.78]; MIX-07 ∈ [0.40, 0.60]. Else kinetics misconfigured.

Metric region:
- max|R| computed over di ∈ [24, 48] only, NOT full depth. Confirm in code.

---

## 6. Estimated effort

- §3.1 verification: 10 min
- §3.2 CW extraction: 5 min
- §3.3 engine runs: 2 runs × ~30 sec each (realistic Hu, longer CFL-bounded runtime) ≈ 1 min
  compute, ~20 min including wrapper setup and sanity-check logging
- §3.4 residual computation: 10 min
- §3.6 spatial diagnostic (both mixes): 30 min
- §3.7 reporting: 30 min

Total: ~1.5–2 hours.

---

## 7. Exit

Stage 1-pilot closes when:

1. Datasets verified (§3.1)
2. Wrapper config logged with all values (§3.3 pre-run log)
3. Pre-run sanity checks pass for both mixes
4. Engine runs complete with all corrections applied
5. Post-run sanity checks pass for both mixes
6. Residual metric tabulated for both mixes (§3.4)
7. Outcome classified into 1, 2, or 3 (§3.5)
8. Spatial diagnostic produced for both mixes (§3.6)
9. PILOT_REPORT.md written with all six sections
10. CSVs and plots saved
11. Synthesis factual

After exit, the user has a go/no-go signal for the full 15-mix validation.

---

## 8. Notes for the new Claude Code session

This is a fresh session. The agent should be aware of:

- **Sprint 7 calibration record:** `validation/soil_calibration/` has the 9-run baseline, parsers
  (`parse_cw_dat`, `parse_cw_temp_output`), and the documented coordinate convention (di=0 top,
  di=48 bottom, wi=0 centerline, wi=12 form face).
- **Sprint 8 corrections:** Hu_residual = 12,937 J/kg and c(T_pl) quadratic. Documented in
  `validation/sprint8/SPRINT8_CLOSURE_REPORT.md`.
- **Engine source:** `thermal_engine_2d.py` has Sprint 7's `K_UC_CALIBRATION_FACTOR_SPRINT7 = 0.96`.
  Sprint 8's corrections are NOT in the source — applied via wrapper only.
- **Engine α_u flow:** the engine field `mix.alpha_u` (line 1513) feeds three downstream paths:
  (1) hydration kinetics α(t) and dα/dt → Q = Hu × Cc × dα/dt (heat source magnitude),
  (2) k(α) = k_uc × (1.33 − 0.33·α) (thermal conductivity, line 367),
  (3) Cp(α) (specific heat).
  The wrapper sets `mix.alpha_u_effective` and `mix.Hu_J_kg_effective`; the engine reads these.
  Q(t) = Hu_eff × Cc × dα_eff/dt — under B1, c cancels exactly so Q matches the raw user-specified
  physics. k(α) and Cp(α) see the c-scaled α. This is the intended Sprint 9 semantic.
- **CW input.dat schema:** α_u at line 389, Hu at line 390. Climate at lines 519–531.
- **Bottom-side BC stencil asymmetry:** known issue from Sprint 7, deferred. Manifests as residual
  near di=48 with magnitude scaling with α_u and bottom thermal flux. Sprint 8 §4.11.8 quantified
  this at ~0.49°F at α_u=0.80 in I-scenario. Likely candidate if pilot residuals concentrate at di=48.
- **Mix kinetic envelope:** MIX-01 represents the moderate-kinetics regime closer to Sprint 8's
  calibration conditions; MIX-07 is the slow-kinetics, high-slag extreme of the production set.
  τ=75 hr for MIX-07 means α(168) is well below α_u — the simulation does not reach hydration
  completion within the 168 hr window.
