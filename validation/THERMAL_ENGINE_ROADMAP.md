# CalcShore Thermal Engine — Validation Roadmap

**Last updated:** 2026-05-04 (Sprint 9 close)
**Status:** Active — updated at each sprint closure
**Audience:** Internal sprint planning, advisor discussions (Tyler Ley), grant applications (NCHRP, NSF SBIR), customer technical conversations

---

## 1. Current State (post-Sprint 9)

Sprints 7–8–9 produced a thermal solver matching ConcreteWorks (CW) within 0.17–0.76°F at peak
temperature across two production mixes spanning the kinetic envelope (MIX-01: β=0.895,
MIX-07: β=0.516, Hu=424–463 kJ/kg, α_u=0.76–0.89, T_pl=73°F, T_soil=85°F, 80 ft deep
submerged mat). Sprint 8's c(T_pl) shape correction, calibrated under suppressed-Hu test
conditions, transfers cleanly to realistic hydration heat release via B1 inverse-compensation,
confirming that the engine's kinetics shape error and the CW solver's kinetics shape are
decoupled from total heat magnitude. The bottom-side BC stencil asymmetry (di=47 cell,
~1.05°F MIX-01 / ~1.58°F MIX-07) has been isolated as a structural numerical residual driven
by the engine's internal discretization, not a physics mismatch; it is independent of both
calibration knobs (Hu_factor and c_multiplier) within physically meaningful operating ranges
and is deferred to Sprint 10. Calibration is closed for the suppressed-Hu and realistic-Hu
regime at composition-correct Hu_factor, for a single placement scenario (T_pl=73°F,
T_soil=85°F) and rectangular submerged-mat geometry.

---

## 2. Validation Roadmap to CW Replacement

Each phase isolates one physics axis. Estimated timelines assume part-time engineering
effort and depend on CW dataset availability; they are planning estimates, not commitments.

---

### Phase A — Mix × scenario coverage

**Physics axis:** Concrete bulk thermal physics across the production composition and BC
temperature envelope.

**Scope:** 5 production mixes (MIX-01, MIX-04, MIX-06, MIX-07, MIX-10) × 5 placement
scenarios (T_pl/T_soil ∈ {50/50, 73/73, 73/85, 90/60, 60/90}) = 25 CW runs.
Full specification in `validation/sprint9/PHASE_A_SCOPING_BRIEF.md`.

**Engine claim at phase close:** Thermal engine validated against CW within ±2°F at peak
temperature (±gate to be confirmed with Tyler Ley before Phase A launch) for rectangular mass
concrete pours under quiescent BC conditions, across MIX-01 through MIX-10 spanning
β=0.516–1.225. Ready for bulk-temperature and Max ΔT report generation in customer documents.

**Prerequisites:** Sprint 10 stencil fix committed; `compute_hu_factor` specification audited
for all 5 Phase A mixes. (See §7 for ordering discussion.)

**Estimated timeline:** 1 month.

**Deferred at phase close:** Top-surface variation, blanket, curing compound, formwork,
non-rectangular geometry, admixtures, silica fume compositions.

---

### Phase B (Sprint 10) — Bottom-side BC stencil correction

**Physics axis:** Engine numerical accuracy at the cell adjacent to the bottom Dirichlet face.

**Scope:** Engine source work — mirror the top-side quarter-cell + half-cell stencil treatment
to the bottom-side cell (di=47) with appropriate sign conventions. This is not a BC physics
change; both solvers already apply identical strong Dirichlet at di=48. Followed by re-pilot
on MIX-01 and MIX-07 to confirm the di=47 residual drops within gate, and regression check
on Sprint 7's 9-run Structure C suite.

**Engine claim at phase close:** Max|R| < 0.5°F numerical agreement with CW across bulk +
bottom-region in submerged-mat geometry for both pilot mixes. Sprint 8's I-scenario
α_u ∈ {0.60, 0.80} failures (deferred since Sprint 8) expected to close.

**Estimated timeline:** 1–1.5 months including implementation, regression, and re-pilot.

**Deferred at phase close:** Top-surface variation, blanket, curing compound, formwork,
non-rectangular geometry.

---

### Phase C — Top-surface physics

**Physics axis:** Air/concrete interface: convection coefficient, radiation, evaporative
cooling, T_ambient time series (dynamic weather, not constant override).

**Scope:** Sprint 9 explicitly restricted the residual metric to the bottom half of the domain
(di ∈ [24, 48]) to avoid contamination from top-surface handling. Phase C removes that
restriction. Validating the top surface requires characterizing how CW handles convection,
radiation, and evaporative cooling (see `CW_BEHAVIORAL_SPEC.md` §5), and characterizing
whether T_ambient time series inputs change CW's BC application relative to the constant
override used in Sprints 7–9. This is the most physics-rich phase and the one with the most
behavioral uncertainty in CW.

**Engine claim at phase close:** Full-column (di ∈ [0, 48]) validation for quiescent
ambient conditions. T_ambient time series (dynamic weather coupling) validated for standard
placement season profiles.

**Estimated timeline:** 2 months.

**Deferred at phase close:** Blanket, curing compound, formwork, non-rectangular geometry.

---

### Phase D — Insulation blanket physics

**Physics axis:** Thermal resistance and thermal mass of insulation blanket above the
concrete surface.

**Scope:** Validate `blanket_thickness_m` parameter behavior: does CW model the blanket as a
pure resistance layer, as resistance + thermal mass, or as a modified surface BC coefficient?
(See `CW_BEHAVIORAL_SPEC.md` §5.) Sprint 7–9 ran with `blanket_thickness_m=0.0`; Phase D
introduces non-zero blanket thickness as a separate validation axis.

**Engine claim at phase close:** Insulation blanket effect on core temperature and Max ΔT
validated against CW for standard construction blanket thicknesses.

**Estimated timeline:** 1 month.

---

### Phase E — Curing compound physics

**Physics axis:** Surface chemical treatment effects — reduced evaporative cooling, modified
surface emissivity.

**Scope:** Curing compound is physically distinct from a blanket (different layer, different
mechanism). Characterize whether CW models it as a surface emissivity modification, an
evaporation reduction factor, or both. Then validate the engine's implementation against CW's
behavioral response.

**Engine claim at phase close:** Curing compound effect on surface temperature and early-age
hydration validated against CW.

**Estimated timeline:** 1 month.

---

### Phase F — Formwork physics

**Physics axis:** Side-BC thermal mass, formwork removal timing, side-corner stencil.

**Scope:** Sprint 7–9 used `is_submerged=True` (earth on all sides) and symmetric geometry.
Phase F introduces formwork as a distinct side-BC layer, validates formwork thermal mass,
and characterizes the CW behavioral response to formwork removal timing events. Also
validates the engine's side-corner stencil accuracy, which received less scrutiny than the
bottom-corner stencil in Sprints 7–9.

**Engine claim at phase close:** Formed concrete elements (not just earth-formed) validated
against CW. Side-corner numerical accuracy confirmed.

**Estimated timeline:** 1.5 months.

---

### Phase G — Chemical admixtures

**Physics axis:** Accelerators, retarders, water reducers — whether CW has an internal
admixture kinetics model or defers entirely to user-supplied (τ, β, α_u, E_a).

**Scope:** If CW modifies kinetics internally based on admixture concentration inputs
(not captured by user-supplied parameters), then every Sprint 7–9 run that passed through
τ/β/α_u unmodified implicitly bypassed CW's admixture path. Phase G characterizes this.
If CW only uses user-supplied kinetics parameters and has no internal admixture model,
every sprint implicitly validated admixture effects already (they show up in the calibrated
τ/β/α_u values), and Phase G collapses to a documentation exercise.

**Engine claim at phase close:** Engine handles admixture-modified mixes correctly, or
confirmed that admixture effects are fully captured via kinetics parameter pass-through.

**Estimated timeline:** 0.5–1 month (scope contingent on CW's actual internal model).

---

### Phase H — Rectangular-member completeness

**Physics axis:** Integrated — all phases A–G validated together on rectangular mass
concrete pours.

**Scope:** Regression sweep across phases A–G simultaneously — all mix × scenario × BC
combinations in scope. End-to-end confidence that no phase's correction breaks another.

**Engine claim at phase close:** CalcShore thermal engine can replace CW for rectangular
mass concrete members across the production composition envelope, standard placement
scenarios, with insulation blankets, curing compound, and formwork. This is the
commercially meaningful milestone for rectangular-element document generation.

**Estimated timeline:** 1 month.

---

### Phase I — Geometry generalization

**Physics axis:** Non-rectangular element shapes — each geometry type validates separately.

**Scope:** Rectangular validation (Phase H) establishes the methodology and tooling. Each
non-rectangular geometry (cylindrical pier, T-beam, L-wall, custom section) requires its
own CW dataset generation, mesh handling, and residual analysis. Priority order TBD by
user (see §7).

**Engine claim at phase close (per geometry type):** Engine validated against CW for the
specified non-rectangular geometry type. Each validated geometry unlocks a new document
generation vertical (bridge pier calculations, wall calculations, etc.).

**Estimated timeline:** ~1.5 months per geometry type × 3–5 types = 4–7 months.

---

### Total roadmap estimate

**14–18 months** of structured part-time engineering work from current state (post-Sprint 9)
to "engine can replace CW for rectangular + 3–5 non-rectangular geometry types."

Phases A–H (rectangular completeness): ~10–11 months.
Phase I (geometry generalization, 3–5 types): 4–7 months.

Timeline depends on CW dataset availability (user-generated; not parallelizable with sprint
execution until the datasets exist) and engineering bandwidth.

---

## 3. Phase Claim Matrix

| Phase | Internal capability | Customer-facing claim | Deferred items |
|---|---|---|---|
| **Post-Sprint 9** | Bulk core T validated for 2 mixes at T_pl=73°F, T_soil=85°F; bottom-face stencil isolated and quantified | Engine calibrated and characterised for rectangular submerged-mat geometry at nominal conditions | Stencil fix, broader mix/scenario coverage, top surface, blanket, formwork, non-rectangular |
| **Phase A** | 5-mix × 5-scenario bulk T validated, β=0.516–1.225 | "Thermal engine validated against CW within ±2°F at peak temperature for rectangular mass concrete pours across the production composition envelope" | Top-surface variation, blanket, curing compound, formwork, non-rectangular geometry |
| **Phase B** | Bottom-side stencil corrected; max\|R\| < 0.5°F across bulk + bottom-region | "Numerical agreement with CW to ±0.5°F in bulk and bottom-face region for rectangular submerged-mat geometry" | Top surface, blanket, curing compound, formwork, non-rectangular |
| **Phase C** | Full-column (di 0–48) T validated; T_ambient time series validated | "Full-column thermal validation including ambient air interface and dynamic weather coupling" | Blanket, curing compound, formwork, non-rectangular |
| **Phase D** | Blanket thermal resistance + mass validated | "Insulation blanket effect on core temperature and Max ΔT validated against CW" | Curing compound, formwork, non-rectangular |
| **Phase E** | Curing compound surface effects validated | "Curing compound effect on surface temperature and early-age hydration validated against CW" | Formwork, non-rectangular |
| **Phase F** | Formwork thermal mass and removal timing validated; side-corner stencil confirmed | "Formed concrete element thermal behavior validated against CW; formwork removal timing included" | Non-rectangular geometry |
| **Phase G** | Admixture kinetics model characterised and validated (or confirmed as kinetics pass-through) | "Engine handles admixture-modified mixes; admixture effects captured via calibrated kinetics parameters" | Non-rectangular geometry |
| **Phase H** | All phases A–G integrated and regression-verified | "CalcShore thermal engine can replace CW for rectangular mass concrete members across the production composition envelope with standard construction practices" | Non-rectangular geometry |
| **Phase I** (per type) | One non-rectangular geometry type fully validated | "Thermal engine validated against CW for [geometry type] elements" | Other geometry types not yet validated |

---

## 4. Methodology Compounding

Each sprint of CW reverse-engineering produces transferable tooling, conventions, and
diagnostic patterns that accelerate future sprints. The following methodology refinements
were established in Sprints 7–9:

**Parsers and loaders:**
- `parse_cw_dat` — ingests CW's `input.dat`; extracts T_pl, T_soil, α_u, Hu, τ, β, E_a, cement, w/cm, geometry. (Sprint 7 Stage 2)
- `parse_cw_temp_output` — ingests CW's `output.txt` temperature field; handles the 0.02 m depth-origin offset. (Sprint 7 Stage 1)
- `cw_scenario_loader` — combines both parsers with CW-grid alignment for engine comparison runs. (Sprint 7 Stage 2)
- `resample_engine_to_cw` in `stage3_compare.py` — bilinear interpolation of engine output onto CW's (di, wi, t) grid. (Sprint 7 Stage 3)
- `make_neutral_env` in `stage4b_run.py` — deterministic neutral-environment wrapper from a scalar T_pl. (Sprint 7 Stage 4b)
- `align_cw_to_engine` / `ti_near` in sweep scripts — CW cache-once-per-mix pattern for fast sweeps. (Sprint 9)

**Grid and coordinate conventions:**
- di=0 top, di=48 bottom; wi=0 centerline, wi=12 form face. (Sprint 7 Stage 1)
- 0.02 m depth-origin offset in CW depth axis — must be applied during parsing. (Sprint 7 Stage 1)
- 6× grid refinement for concrete domain — confirmed artifact-free at Sprint 7 Stage 5b; carried forward. (Sprint 7)

**Residual metric conventions:**
- Structure C decomposition: R1 side profile (wi slice), R2 bottom profile (di slice), R3 corner RMSE. (Sprint 7 Stage 5d-prep)
- Bottom-half-only metric region (di ∈ [24, 48], wi=0) — isolates BC-region physics from top-surface contamination. (Sprint 9)
- stencil_drop = R(di=47) − R(di=42) — structural decoupling test for the bottom-face stencil dip relative to bulk level. (Sprint 9)
- R_di48 constancy check — range < 0.05°F confirms strong Dirichlet write invariance across interior bulk changes. (Sprint 9)

**Wrapper patterns:**
- B1 inverse-compensation: `α_u_eff = c·α_u`, `Hu_eff = (Hu_raw·Hu_factor) / c` — preserves Q(t) exactly while routing c-scaled α through k(α) and Cp(α). (Sprint 9 formalism, Sprint 8 c() origin)
- Hu_residual conditional passthrough for production mixes (never fires when Hu_raw >> 10,000). (Sprint 8 §4.10)
- Composition-centered sensitivity sweep: ±1% Hu_factor around compute_hu_factor output × ±5% c_multiplier — cleanest signal/noise design for stencil decoupling analysis. (Sprint 9)

**Mix selection strategy:**
- Kinetic envelope bracketing: select mixes at the fast-kinetics extreme (high β, short τ) and slow-kinetics extreme (low β, long τ) as pilot endpoints. MIX-01 (β=0.895) and MIX-07 (β=0.516) bracket the Sprint 9 production set. (Sprint 9 pilot design)

**Future sprint closure reports should include a "Methodology Refinements" section** that adds to this list, so the compounding is documented in one place.

---

## 5. Strategic Context

CalcShore's customer-facing product is document automation — TCP (Thermal Control Plans), mix design reports, cold/hot weather plans, and placement approval documents — supported by a CW-validated thermal engine, not the thermal engine alone. CW replacement builds over 14–18 months; document automation revenue can flow earlier. Phase A validation (5 mixes × 5 scenarios, bulk T and Max ΔT) is sufficient for customer-shippable document accuracy for rectangular mass concrete pours under standard placement conditions — the most common use case. Phases B–H strengthen numerical accuracy and add construction practice coverage (blanket, curing compound, formwork); Phase I unlocks broader element coverage (piers, walls, T-beams) and the corresponding document verticals. The CW Behavioral Specification (`CW_BEHAVIORAL_SPEC.md`) — 14–18 months of observed CW behavioral findings — is a compounding competitive moat: no other party has documented CW's floor artifacts, stencil asymmetries, and kinetics shape mismatch to this level of precision, and that knowledge is embedded in the calibration record.

---

## 6. Update Protocol

This document is updated at each sprint closure. The updating sprint adds:

1. **Actual close date** for the completed phase (replaces the estimate in §2).
2. **Timeline drift** — note actual vs. estimated; update remaining phase estimates if warranted.
3. **Methodology refinements** — append new items to §4.
4. **Claim matrix updates** — update the "Post-Sprint N" row in §3 with the new internal capability and customer-facing claim.
5. **Open questions resolved** — move any resolved items from §7 to the body of the document.

Sprint closure reports should cross-reference this document in their §9 (Next Phase) section.

---

## 7. Open Questions for User Discussion

These items require user decisions outside of any specific sprint:

1. **Phase A gate value.** ±2°F is a placeholder. Confirm the engineering-tolerance threshold
   with Tyler Ley before Phase A launch. This affects how the Phase A pass/fail grid is
   interpreted and what language appears in customer documents.

2. **Phase ordering: A vs. B.** The current brief (PHASE_A_SCOPING_BRIEF.md) lists Phase A
   as requiring Sprint 10's stencil fix (Phase B) as a prerequisite. If Phase A proceeds
   before the stencil fix, the max|R| ≈ 1.05–1.58°F bottom-region residual is present in all
   25 runs; the question is whether the ±2°F gate absorbs it. User to decide whether to run
   Phase A before or after Phase B.

3. **Phase G scope.** Depends on whether CW has an internal admixture kinetics model
   (distinct from user-supplied τ/β/α_u/E_a). Scoping Phase G requires first characterizing
   CW's admixture input fields and whether they modify the kinetics ODE solve or are purely
   for record-keeping.

4. **Phase I geometry priority.** Which non-rectangular geometries first? Candidates based
   on CalcShore's customer pipeline: cylindrical bridge pier (highest concrete volume, most
   common mass concrete classification), L-wall (retaining walls), T-beam. User to decide
   the sequence.

5. **CW Behavioral Specification sharing posture.** See `CW_BEHAVIORAL_SPEC.md` §7.
   Decision needed on whether this document is internal-only, advisor-shared, selectively
   customer-shared, or IP-protected. The document grows in value over 14–18 months and
   is part of the technical moat; sharing posture has IP implications.
