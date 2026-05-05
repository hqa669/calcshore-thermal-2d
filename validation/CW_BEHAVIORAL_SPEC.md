# ConcreteWorks Behavioral Specification

**Last updated:** 2026-05-04 (Sprint 9 close)
**Status:** Active — updated at each sprint closure
**Scope:** CalcShore thermal engine validation program, Sprints 7–9 and forward

---

## 1. Purpose and Scope

This document characterizes ConcreteWorks (CW) behavioral specification as observed during
CalcShore thermal engine validation. It is updated at each sprint closure; new findings are
appended to the appropriate phase section. Findings come only from observed CW outputs
against controlled inputs — no claims are made about CW's internal source code or
implementation details. All cited behaviors are reproducible from the CW datasets and engine
comparison scripts archived in `validation/`.

The purpose is twofold. First, it is an engineering reference: CW behaviors that differ from
naive expectations (hidden floors, coordinate offsets, stencil asymmetries) must be accounted
for in every engine comparison run, and this document is the single place to find them.
Second, it is a compounding competitive asset: over 14–18 months of characterization, this
specification documents CW's behavioral envelope in more detail than any public source, which
is part of CalcShore's technical moat.

---

## 2. CW Behaviors Observed in Sprint 7

*(Sprint 7 closed 2026-05-01. Reference: `docs/sprints/STAGE5e_REPORT.md`)*

---

### 2.1 Depth-axis origin offset

**Observation:** The CW temperature export (`output.txt`) uses a non-zero depth origin.
`cw_depths_m[0]` is not 0.0 m; a 0.02 m offset exists in the depth axis relative to the
concrete domain top face.

**Effect:** If ignored, bilinear resampling of engine output onto the CW grid is biased by
approximately 0.1°F in near-surface rows (di ≈ 0–4).

**Resolution:** `parse_cw_temp_output` accounts for the offset at parse time. Do not
double-apply.

**Cite:** Sprint 7 Stage 1.

---

### 2.2 CW bulk heat-release noise floor (suppressed Hu)

**Observation:** At suppressed hydration (Hu_effective = 1 J/kg), CW produces approximately
0.13°F of spurious temperature rise over 168 hr. This floor is:
- k-independent (present uniformly whether k_uc is swept 0.92–1.04)
- Present in all 9 Sprint 7 calibration runs including the |ΔT|=0 baseline (Run A)
- Spatially uniform (not concentrated at faces or corners)

**Effect:** At suppressed Hu, the engine correctly produces ~0°F rise. The 0.13°F floor makes
the engine appear ~0.13°F colder than CW across all cells.

**Resolution:** Sprint 7 absorbed ~0.13°F floor contribution into the `k_uc × 0.96` calibration
factor. Sprint 8 later characterized this more precisely and introduced the Hu_residual
conditional as the correct mechanism (see §3.1).

**Cite:** Sprint 7 Stage 3 (identified); Sprint 8 §4.10 (quantified).

---

### 2.3 Boundary-onset convention (t=0 → t=24 transient)

**Observation:** CW applies boundary conditions (T_soil, T_ambient) instantaneously at t=0.
The engine applies BCs with a gradual ramp during the first few timesteps.

**Effect:** At t=24 hr, R2 (bottom profile) max|R| exceeds 0.50°F for high-|ΔT| runs even
after `k × 0.96` calibration. By t=84 hr the transient has largely decayed. It is a
structural condition, not a calibration deficiency.

**Resolution:** Treated as a pre-existing acknowledged condition in Sprints 7–9. Not
addressed; early-time agreement at t < 48 hr is not a gate criterion.

**Cite:** Sprint 7 Stage 5d-prep diffusivity quantification; Sprint 7 §4 "t=24 early-time
residuals"; Sprint 9 Stage 1-pilot brief §3.6.

---

### 2.4 Bottom-corner stencil asymmetry (first quantification)

**Observation:** In Sprint 7 runs with |ΔT| > 0, the residual field shows a laterally
distributed bulge near the bottom-side corner (di=43–48, wi=8–12). The top-side corner
uses a quarter-cell energy balance + half-cell BC stencil; the bottom-side corner uses a
pure strong Dirichlet write with no analogous treatment.

**Effect:** At α_u ≈ 0.036 (Sprint 7 suppressed Hu), the asymmetry contributes approximately
0.085°F R3 RMSE; it breaks the F/H sign symmetry (cooling vs. warming at |ΔT|=28°F). Not a
gate failure at Sprint 7's α level.

**Further characterization:** Sprint 8 measured ~0.49°F at α_u=0.80 in I-scenario. Sprint 9
measured 1.05–1.58°F at realistic production-mix α_u (0.80–0.94) in submerged-mat geometry.
Submerged-mat geometry (no formwork, direct earth contact on all sides) exacerbates the
asymmetry relative to Sprint 8's I-scenario configuration.

**Note:** This is a CW vs. engine stencil asymmetry — CW's internal stencil for the
bottom-face cell differs from the engine's. It is NOT a BC physics mismatch; both solvers
apply identical T_soil Dirichlet at the bottom face.

**Cite:** Sprint 7 Stage 5d-prep "Corner BC asymmetry"; Sprint 8 §4.11.8; Sprint 9
SWEEP_V2 §1.3/§1.8.

---

## 3. CW Behaviors Observed in Sprint 8

*(Sprint 8 closed 2026-05-03. Reference: `validation/sprint8/SPRINT8_CLOSURE_REPORT.md`)*

---

### 3.1 Hidden minimum Hu floor

**Observation:** CW has a hidden minimum heat-release floor. When user-supplied Hu falls below
approximately 3,000–10,000 J/kg, CW silently substitutes a floor value of approximately
12,937 J/kg. Above the floor, CW honors the input linearly.

**Characterization:** 14 CW runs at Hu ∈ {1, 10, 100, 1000, 10000, 50000, 100000} J/kg ×
α_u ∈ {0.20, 0.80} in A-scenario. Engine has no floor (clean linear response from Hu→0).
Hu_residual = 12,937 J/kg self-calibrated via probe solve at the c=1 anchor point.

**Effect at production Hu:** For mixes with Hu_raw >> 10,000 J/kg (both Sprint 9 pilot mixes
at 424–463 kJ/kg), the floor never fires. The Hu_residual conditional is a passthrough for
production mixes.

**Effect at suppressed Hu:** The 0.13°F Sprint 7 noise floor (§2.2) is the observable
consequence of CW substituting 12,937 J/kg when the user sets Hu=1 J/kg.

**Resolution:** `Hu_residual` conditional in engine wrapper: if `Hu_raw < 10,000` →
`Hu_J_kg_effective = 12,937`; else passthrough.

**Cite:** Sprint 8 §3.1 and §4.2–§4.9.

---

### 3.2 input.dat parameter schema

**Observation:** CW's `input.dat` encodes key kinetic parameters at specific line numbers:

| Parameter | Line |
|---|---|
| α_u (ultimate hydration degree) | 389 |
| Hu (heat of hydration, J/kg) | 390 |
| Climate temperatures (monthly) | 519–531 |

Weather override: the T_ambient runtime override is at line 440; lines 519–531 are the
underlying monthly baseline and do not represent the run-specific ambient temperature.

**Resolution:** `parse_cw_dat` reads from the correct lines. Hand-setting lines 389–390 in
production mix `input.dat` files bypasses CW's internal interpolation database.

**Cite:** Sprint 8 §4.1; Sprint 9 pilot brief §3.1.

---

### 3.3 High-slag interpolation database underprediction

**Observation:** In CW's 80-mix interpolation database, mixes with ≥50% slag-only replacement
(no hand-set Hu) underpredict heat release by approximately 5–15°F at peak temperature.

**Effect:** If CW derives Hu from its interpolation database rather than from the user-supplied
line 390, high-slag mixes produce lower-than-expected peak temperatures in CW.

**Resolution:** Bypass by hand-setting Hu directly at `input.dat` line 390 to match the
production mix table. Sprint 9's pilot mixes used this approach; CW's database path was not
exercised.

**Open:** Whether the database underprediction reflects an error in the database, a different
mix-chemistry assumption, or a documented CW convention is not yet characterized.

**Cite:** Sprint 8 §4.10.

---

### 3.4 c(T_pl) kinetics shape mismatch

**Observation:** Under suppressed Hu (Hu=12,937 J/kg floor value), CW's α(t) trajectory
diverges from the engine's at temperature extremes. The divergence is:
- Multiplicative on α_u (a scaling, not a shape deformation)
- α_u-independent (cross-α_u slope std/mean = 0.68%)
- Well-modeled as: `α_u_engine_corrected = c(T_pl) × α_u_nominal`

The calibrated quadratic:
```
c(T_pl) = −1.3025 + 0.04746 × T_pl − 2.081×10⁻⁴ × T_pl² (T_pl in °F)
```
Fit quality: R² = 0.99993, max residual 0.0055, anchored at c(70°F) = 1.000.

**Effect:** Engine is "faster" (higher effective α_u) at cold T_pl and "slower" at hot T_pl
relative to CW. The shape is consistent with an Arrhenius kinetics disagreement (engine and
CW likely differ in T_ref or activation energy by a small amount), but the specific parameter
difference has not been directly verified.

**Validity range:** T_pl ∈ [40, 110]°F. Extrapolation outside this range is not validated.

**Resolution:** Apply via B1 inverse-compensation wrapper. See `SPRINT8_CLOSURE_REPORT.md`
§3.2 for implementation detail.

**Cite:** Sprint 8 §3.2 and §4.11.

---

## 4. CW Behaviors Observed in Sprint 9

*(Sprint 9 closed 2026-05-04. Reference: `validation/sprint9/SPRINT9_CLOSURE_REPORT.md`;
`validation/sprint9/stage1_pilot_sweep_v2/SWEEP_V2_REPORT.md`)*

---

### 4.1 No soil element in CW when model_soil=False

**Observation:** With no soil model activated in `input.dat` (CW's `model_soil`-equivalent
flag not set), CW applies T_soil as a direct Dirichlet boundary condition at the bottom face
of the concrete domain. No soil thermal mass, no soil layer, no soil element is present on
the CW side. The concrete domain terminates at di=48.

**Effect:** Both engine (with `model_soil=False`) and CW apply identical T_soil Dirichlet at
di=48. Any residual at or near di=48 is therefore a numerical artifact, not a physics
mismatch.

**Resolution:** Confirmed via user clarification and validated by R_di48 constancy finding
(§4.2). Engine `model_soil=False` setting is the correct comparison mode.

**Cite:** Sprint 9 user clarification (session context); Sprint 9 pilot brief §8.

---

### 4.2 Bottom Dirichlet face residual — small, invariant offset

**Observation:** When both engine and CW apply T_soil=85°F as a Dirichlet BC at di=48, a
residual of +0.206°F ± 0.004°F exists at di=48 across all tested operating points.

**Characterization (Sprint 9 sweep v2):** R_di48 range = 0.00817°F across all 50 grid
points (5 Hu_factor × 5 c_multiplier × 2 mixes). Mean = +0.2061°F. The offset is
independent of mix, Hu_factor, c_multiplier, and bulk physics.

**Effect:** The +0.206°F face offset is a consistent small numerical artifact at the
bottom face itself. It does not propagate strongly into the interior; di=47 and inward are
not dominated by this value.

**Origin:** How each solver enforces the strong Dirichlet write at the face cell. Not
characterized further.

**Cite:** Sprint 9 `SWEEP_V2_REPORT.md` §1.7.

---

### 4.3 Bottom-region stencil asymmetry under realistic Hu — production-scale quantification

**Observation:** At production-mix Hu (424–463 kJ/kg) and realistic α_u (0.80–0.94 effective),
the cell immediately above the bottom Dirichlet face (di=47) shows a residual of:
- MIX-01 (α_u_eff = 0.7988): R_di47 = −1.0507°F, stencil_drop = −1.2023°F
- MIX-07 (α_u_eff = 0.9410): R_di47 = −1.5843°F, stencil_drop = −2.2252°F

(`stencil_drop = R_di47 − R_di42`, measuring the stencil dip relative to bulk level)

This is 2–3× larger than Sprint 8's I-scenario estimate of ~0.49°F at α_u=0.80. The
difference reflects the submerged-mat geometry (bottom heat flux directly into T_soil=85°F
with no soil buffer) vs. Sprint 8's geometry.

**Structural decoupling confirmed:** The stencil dip is independent of c_multiplier within
its physically meaningful ±5% range:
- |∂(stencil_drop)/∂c_mult| = 0.0132°F per 1% (MIX-01), 0.0097°F (MIX-07)
- Threshold: 0.05°F per 1% — both mixes well below threshold

No amount of Hu_factor or c adjustment can close the di=47 gap at production-mix operating
points. The fix must be in the engine's bottom-face stencil.

**Cite:** Sprint 9 `SWEEP_V2_REPORT.md` §1.3, §1.4, §1.8; `SPRINT9_CLOSURE_REPORT.md` §5.3.

---

### 4.4 c(T_pl) correction transfers to realistic Hu under B1 inverse-compensation

**Observation:** The Sprint 8 c(T_pl) shape correction, calibrated under suppressed Hu (Hu
= 12,937 J/kg), transfers to realistic hydration heat release (Hu = 383,000–393,000 J/kg
effective) when applied via the B1 inverse-compensation wrapper.

**Evidence:** Bulk core residuals at composition-correct operating points:
- MIX-01: T_engine(di=24, t=168) = 149.40°F, T_CW = 149.23°F → R = +0.17°F
- MIX-07: T_engine(di=24, t=168) = 148.72°F, T_CW = 147.96°F → R = +0.76°F

MIX-01 bulk is within the 0.5°F gate. MIX-07 bulk exceeds it by a small margin (see §5 open
question on β=0.516 kinetics).

c() sensitivity at composition-correct operating points: 0.010–0.013°F per 1% change in
c_multiplier, both mixes (below the 0.05°F structural-decoupling threshold). The calibrated
c(73°F) = 1.0532 is right-sized; small perturbations do not significantly affect outcomes.

**Cite:** Sprint 9 `PILOT_REPORT.md` §1.6; `SWEEP_V2_REPORT.md` §1.4; `SPRINT9_CLOSURE_REPORT.md`
§5.2.

---

## 5. Open Behavioral Questions for Future Sprints

Items below have not been characterized. Each is a target for a future sprint or phase; the
corresponding entry in `THERMAL_ENGINE_ROADMAP.md` identifies which phase addresses it.

---

**Top-surface physics (Phase C)**
How does CW handle the air/concrete interface? Specifically: what convection coefficient
formula does CW use, does CW apply a radiation term, is there an evaporative cooling model,
and how does CW ingest T_ambient time series (dynamic weather) vs. a constant scalar? The
metric restriction to di ∈ [24, 48] in Sprints 7–9 deliberately avoided this question.

---

**Blanket thermal model (Phase D)**
When `blanket_thickness_m` > 0, does CW model the blanket as a pure thermal resistance
layer (R-value), as resistance + thermal mass, or as a modified surface BC coefficient?
How does the blanket interact with the concrete surface temperature at the interface?

---

**Curing compound physics (Phase E)**
Does CW model curing compound as a surface emissivity modification, an evaporation
reduction factor, or both? Is curing compound effect additive or multiplicative with the
blanket model?

---

**Formwork thermal mass and removal (Phase F)**
How does CW handle formwork as a side-BC layer — as a thermal resistance, as a lumped
thermal mass, or as a temperature-conditioned Dirichlet? How does CW handle the formwork
removal event (time-step BC change) and what happens to the side temperature immediately
after removal?

---

**Chemical admixture kinetics model (Phase G)**
Does CW have an internal admixture kinetics model that modifies the α(t) ODE solve based
on admixture concentration inputs, independent of the user-supplied τ/β/α_u/E_a? Or does
CW treat admixture type and dosage as metadata only, with kinetics fully determined by
the user-supplied parameters? If CW has an internal model, every Sprint 7–9 run that
passed through kinetics parameters directly bypassed it.

---

**Non-rectangular geometry meshing (Phase I)**
For cylindrical, T-beam, L-wall, and custom cross-section elements, does CW use a
structured rectangular mesh with staircase approximation, an unstructured mesh, or a
custom element-specific solver? How are corner and edge cells handled in non-rectangular
meshes, and does the stencil asymmetry observed in the rectangular bottom-corner
(§2.4/§4.3) have analogues in non-rectangular geometries?

---

**MIX-07 bulk offset source (Sprint 9 open item)**
R_di36 = +0.76°F for MIX-07 vs. +0.17°F for MIX-01 at the same scenario. Both mixes
use c(73°F) = 1.0532 via the same wrapper. Possible cause: the Sprint 8 c() calibration
covered β ∈ [0.895 range] primarily; MIX-07's β=0.516 places it in a kinetic-shape regime
not well represented in the 32-point A-scenario sweep. Whether CW handles β=0.516 kinetics
differently than the engine at realistic Hu, and whether this is a c()-form inadequacy
or a c()-range gap, is not yet characterized.

**Cite:** Sprint 9 `SPRINT9_CLOSURE_REPORT.md` §6.3.

---

**compute_hu_factor in SF regime (Sprint 10 target)**
Does `compute_hu_factor` handle silica fume (SF) content correctly? SF's pozzolanic
activity and heat release timing differ from cement, fly ash, and slag. The function was
characterized only for Type I/II + FA + slag compositions in Sprint 9. MIX-13 (SF-bearing)
is deferred from Phase A pending formal specification. See `COMPUTE_HU_FACTOR_NOTE.md` §5(d).

---

## 6. Update Protocol

This document is updated at each sprint closure. New findings are:

1. **Appended to the appropriate sprint section** (§2, §3, §4, or a new §5 for Sprint 10, etc.)
   rather than edited into existing sections. This preserves the observation timeline.
2. **Cited** with the sprint number and report section where the finding was characterized.
3. **Observational only.** Claims about CW's internal implementation are explicitly avoided;
   use "CW produces / CW applies / CW substitutes" rather than "CW's source code does."
4. **Open questions resolved** — move from §5 to the appropriate sprint section with a
   "Resolved in Sprint N" note.

Sprint closure reports should cross-reference this document in their §5 (Key Findings)
when new CW behaviors are characterized.

---

## 7. Confidentiality Posture

Decision deferred to user. Options:

**Internal-only (CalcShore engineering reference).** Restricts access to CalcShore team.
Appropriate during active development; protects specifics of CW's behavioral envelope
from competitors who might use it to calibrate their own solvers against CW.

**Advisor-shared (Tyler Ley, future engineering hires).** Enables technical advisor review
of methodology validity. Advisors with CW expertise could confirm whether observed
behaviors match documented CW features or expose undocumented ones.

**Selectively customer-shared (technical depth signal).** Sharing selected findings
(e.g., the Hu floor characterization, the stencil asymmetry quantification) in technical
conversations demonstrates validation depth. Not necessary to share the full document.

**IP-protected.** The compiled behavioral specification of a commercial solver — built
through 14–18 months of controlled characterization — may have competitive value beyond
CalcShore's internal use. Whether this constitutes a protectable trade secret or technical
know-how is a legal question outside this document.

User decision pending.
