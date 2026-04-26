# CalcShore Thermal Engine v3 — Release Notes

## Overview

Engine v3 is the CalcShore thermal physics engine for mass concrete temperature
prediction. Given a mix design, placement geometry, construction parameters, and
hourly climate data, it solves a 2D explicit finite-difference heat equation over
the half-mat footprint and returns a temperature trajectory at configurable output
resolution. Engine v3 predicts peak temperature, peak temperature differential
(gradient), and RMS agreement with CW reference data across the validated
envelope; it does not predict drying shrinkage, mechanical stress, autogenous
shrinkage, or thermal cracking risk — it produces the temperature inputs those
analyses consume. Engine v3 integrates downstream with the CalcShore TCP document
generator and AI compliance evaluator (see `docs/coding_passdown_v3.md` §1 for
the full product mission).

## Validated envelope

- **Geometry**: 40×60×8 ft half-mat footing (width × length × depth), 2D analysis
  exploiting symmetry at the centerline
- **Climate**: Austin TX summer (CW-exported hourly weather files; library has zero
  climate variation — all 14 mixes are Austin TX)
- **Placement temperature**: 60°F warm placement (cold placement at 45°F is
  out-of-envelope per R9; warm placement at 75°F is out-of-envelope per MIX-14
  findings)
- **Mix range**: 39.1% SCM Reference mixes (single-peak hydration profile; 5 mixes
  — MIX-01, MIX-02, MIX-03, MIX-11, MIX-12)
- **Form material**: steel (R_form = 0.0862 m²·K/W, ADR-04, reinforced PR 20)
- **S0 tolerances**:
  - PeakMax T: ±1.0°F
  - Peak gradient: ±2.0°F
  - Field RMS: ≤2.0°F
  - Centerline RMS: ≤1.0°F
  - Corner RMS: ≤3.0°F
- **Reference set status**: S0 5/5 on all 5 Reference mixes at `sprint-6-complete`.
  Values unchanged from `sprint-4-complete` through Sprint 6; PR 22 was
  bit-identical for steel forms; PR 23 was bit-identical by construction.

## Known residuals

### Top-surface BC amplitude over-oscillation (Centerline RMS)

**Magnitude**: ~0.74°F Centerline RMS systematic across the 5
Reference mixes at sprint-4-complete; passes S0 (≤1.0°F) but
exceeds S1-aspire (≤0.5°F).

**Mechanism (PR 18 / PR 19 finding)**:
- PR 18 D2 + D3 + D4: zero phase lag; engine diurnal amplitude is
  ~1.9× CW's at the centerline mid-depth steady-state window;
  residual is depth-localized at the top surface.
- PR 19 magnitude sweeps on h_conv (0.5×–2.0× ACI 305 default),
  α_top (0.30–0.85), ε_top (0.50–0.95): no scalar BC knob has
  ≥0.24°F authority on cluster-mean CenterRMS within physically
  reasonable ranges. D3 amplitude ratio moves *inversely* to
  expectation under h_conv scaling (higher h → larger ratio,
  not smaller).
- The Centerline RMS variation that *is* observable across the
  PR 19 sweeps (~0.04–0.10°F) comes from DC-offset shift (D1),
  not amplitude correction (D3). The amplitude over-oscillation
  is structural, not a scalar miscalibration.

**Implication**: the residual is not closable by recalibrating
any single top-surface scalar BC parameter. The mechanism lives
in either the BC functional form (how solar G(t), LW T_sky(t),
and convection compose into the surface energy balance), the
coupling structure (e.g., blanket-cell thermal mass, sky
temperature parameterization), or the thermal physics model
itself (effective k(T,α), Cp(T,α), or hydration kinetics
interaction with surface response).

**Engine v3 envelope statement**: engine v3 ships with a
documented ~0.74°F Centerline RMS residual on the validated
Reference geometry (40×60×8ft half-mat footing, Austin summer,
mid-SCM at 60°F placement, steel form). The residual is within
S0 (≤1.0°F) and represents a 32% margin to the S1-aspire
threshold.

**Future investigation candidates** (not Sprint 5/6 scope):
- CW's surface energy balance functional form. If CW uses
  natural-convection-only coefficients or applies wind reduction
  factors not present in our engine, that explains the inverse
  h_conv direction.
- Blanket layer thermal-mass coupling. The engine treats blanket
  as an explicit thermal cell with non-zero properties; CW may
  treat it as a pure thermal resistance.
- Sky temperature parameterization. The engine uses
  Berdahl-Martin T_sky from cloud cover; CW's parameterization
  is undocumented.
- Surface wind coupling vs daily-averaged wind. The library has
  zero wind-speed variation across mixes (`cw_ave_max_wind_m_s`
  = 10.5 m/s identical for all 14 mixes), so the engine's
  wind-coupled h_conv evaluates at a single value across the
  library. CW may average differently.

**R9 implication**: the D1 DC offset (+0.40 to +0.56°F on 4/5
Reference mixes) and the D5 Mode-A overlap finding from PR 18
(+0.77 to +1.21°F in [8,48] hr) remain entangled with the
amplitude residual — they cannot be cleanly separated until the
amplitude story closes. R9's kinetics-vs-inherited-calibration
disambiguation is correspondingly deferred to engine v3
investigation alongside the BC structural work.

### R_form Reference-set evaluation

R_form residual documented — PR 15 candidate (R_form=0.060) did
not generalize to the full Reference set. The cluster mixes
(MIX-01, 11, 12) and MIX-03 show large CornerRMS improvements
at 0.060 vs ADR-04 default 0.0862 (−1.1 to −1.6°F). MIX-02
(kinetics-anomalous) responds oppositely: its CornerRMS and
PeakGrad worsen monotonically as R_form decreases, reaching
S0 failure (PeakGradΔ=+2.679°F > 2.0°F threshold) at 0.060.
The cluster optimum (~0.060–0.070) and MIX-02's optimum
(~0.0862) do not overlap at a globally applicable value.
CenterRMS confirmed fully form-face-decoupled: Δ=0.000°F across
all sweep values on all 5 mixes. ADR-04 remains at 0.0862.
Routing: Sprint 6 R5 (form-material parameterization) evaluated
whether a per-form-type or per-mix-regime R_form seam could
recover the cluster improvement without regressing MIX-02. The
cluster contradiction with MIX-02 was explicitly not addressed
(routes to the hydration series). See
`validation/diagnostics/pr20_r_form_eval/`.

### Corner-side residuals

MIX-02 exhibits a kinetics-divergent signature confirmed across
multiple diagnostics in Sprint 5 (§8.2 finding #1): the engine and
CW disagree on heat generation in ways not closable by thermal-physics
calibration alone. The cross-diagnostic signature — MIX-02 responding
oppositely to R_form decreases while cluster mixes improve, plus
anomalous CornerRMS behavior — is consistent with an upstream
divergence in the hydration kinetics model rather than a thermal-BC
miscalibration.

Engine v3 envelope statement for corner-side accuracy: cluster
CornerRMS of ~2.0°F (range 1.24–2.72°F across Reference mixes at
R_form=0.0862) is the validated envelope statement. This is not a
closed residual — it reflects the accuracy achievable within the
current thermal physics model at the validated R_form value. MIX-02
specifically: PeakGrad S0 holds at R_form=0.0862 (ΔPeakGrad=+1.16°F,
within the ±2.0°F gate), but fails at R_form=0.060. The shipped value
(0.0862) is the safe choice for MIX-02-class mixes.

Routing: MIX-02's kinetics divergence routes to the hydration series
for disambiguation. The thermal sprint series treats the ~2.0°F cluster
CornerRMS as the engine v3 floor — improving it requires the hydration
series, not further thermal recalibration.

## Out-of-envelope conditions

The following conditions are explicitly outside the engine v3 validated
envelope. Each is documented, routed to a successor sprint series, and (where
applicable) emits a runtime warning or error.

- **Cold placement (T_placement < 60°F)**: R9 Mode B Arrhenius rate divergence
  above ~30°C. MIX-15 at 45°F placement: PeakMax Δ=−3.00°F, PeakGrad
  Δ=−8.20°F, FieldRMS 4.06°F, CenterRMS 2.07°F — S0 1/5. Routes to the
  hydration series for kinetics disambiguation.

- **Warm placement (T_placement > 60°F)**: MIX-14 at 75°F placement scores
  S0 0/5 (PeakMax Δ=+1.49°F, PeakGrad Δ=−3.49°F, FieldRMS 3.20°F). The
  mechanism is not closed; routes to the hydration series.

- **High-SCM mixes (≥40% non-Reference SCM, Cluster A: MIX-05/06/07/09/10)**:
  Hydration-shape divergence. Engine systematically over- or under-predicts
  heat generation depending on the SCM fraction and type. Routes to the
  hydration series (dual-peak hydration model).

- **Low-SCM mixes (B1: MIX-04 at 21.7%, MIX-08 at 17.4%)**: S0 2/5 each.
  PeakGrad errors exceed the ±2.0°F gate. Routes to the hydration series.

- **Non-half-mat geometry (full-mat, slab, column, wall)**: R6. Untested.
  Routes to the geometry coverage series.

- **Non-Austin climate**: The 14-mix library has zero climate variation (all
  Austin TX). Behavior on cold-climate placements, high-altitude UV, or marine
  environments is untested. Routes to the multi-climate validation series.

- **Kinetics-divergent mixes (MIX-02-class)**: §8.2 finding #1. Engine and CW
  disagree on heat generation in ways not closable by thermal-physics
  calibration. Routes to the hydration series.

- **Non-steel forms**:
  - *Plywood*: `R_FORM_BY_FORM_TYPE["plywood"] = 0.17 m²·K/W` (ACI 306R-88
    Table 7.3.5/§7.3, bulk panel R — physically distinct from steel's
    0.0862 m²·K/W contact resistance value). `cw_validated=False`. Emits a
    one-time `UserWarning` at runtime via `resolve_r_form()`. Use at
    customer's risk.
  - *Plastic_liner*: not implemented. PR 22 found no peer-reviewed source for
    the interface contact resistance value. Raises `NotImplementedError` at
    `resolve_r_form()` call time.

## Validation summary table

Numbers from `validation/sprint4_baseline.md` (regenerated at
`sprint-4-complete`, confirmed zero diff through Sprint 6). S0 tolerances:
PeakMax ±1.0°F, PeakGrad ±2.0°F, FieldRMS ≤2.0°F, CenterRMS ≤1.0°F,
CornerRMS ≤3.0°F.

| Mix | Group | SCM% | PlaceTemp | PeakMax Δ | PeakGrad Δ | FieldRMS | CenterRMS | CornerRMS | S0 | Routing |
|---|---|---|---|---|---|---|---|---|---|---|
| MIX-01 | Reference | 39.1% | 60°F | −0.29°F | −0.29°F | 0.88°F | 0.74°F | 2.22°F | 5/5 | Validated |
| MIX-02 | Reference | 39.1% | 60°F | −0.68°F | +1.16°F | 1.23°F | 0.57°F | 1.24°F | 5/5 | Validated (kinetics-divergent on R_form sweeps; see Corner-side residuals) |
| MIX-03 | Reference | 39.1% | 60°F | −0.30°F | −0.92°F | 0.87°F | 0.79°F | 2.72°F | 5/5 | Validated |
| MIX-11 | Reference | 39.1% | 60°F | −0.13°F | +0.07°F | 0.91°F | 0.79°F | 2.12°F | 5/5 | Validated |
| MIX-12 | Reference | 39.1% | 60°F | −0.31°F | −0.32°F | 0.87°F | 0.73°F | 2.23°F | 5/5 | Validated |
| MIX-04 | B1 | 21.7% | 60°F | −1.46°F | −2.27°F | 1.37°F | 0.93°F | 3.00°F | 2/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-08 | B1 | 17.4% | 60°F | −1.43°F | −3.42°F | 1.33°F | 0.85°F | 4.03°F | 2/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-05 | Cluster A | 52.2% | 60°F | +0.70°F | +1.02°F | 1.24°F | 1.49°F | 2.05°F | 4/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-06 | Cluster A | 71.3% | 60°F | +2.20°F | +3.04°F | 2.35°F | 2.81°F | 1.84°F | 1/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-07 | Cluster A | 91.3% | 60°F | +4.12°F | +4.97°F | 3.82°F | 4.45°F | 1.83°F | 1/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-09 | Cluster A | 41.7% | 60°F | +3.78°F | +1.33°F | 4.15°F | 4.84°F | 4.57°F | 1/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-10 | Cluster A | 52.2% | 60°F | +0.23°F | +1.43°F | 1.07°F | 1.03°F | 1.40°F | 4/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-15 | B2 | 39.1% | 45°F | −3.00°F | −8.20°F | 4.06°F | 2.07°F | 1.99°F | 1/5 | OUT-OF-ENVELOPE — R9 / hydration series |
| MIX-14 | (warm-OOE) | 39.1% | 75°F | +1.49°F | −3.49°F | 3.20°F | 2.29°F | 4.85°F | 0/5 | OUT-OF-ENVELOPE — hydration series |
| MIX-13 | (sentinel) | — | — | — | — | — | — | — | skipped | No CW output.txt |

## API surface

Customer workflow: construct `CWMixDesign`, `CWGeometry`, `CWConstruction`, and
`CWEnvironment` from TCP parameters, then call `run_one()` with a scenario
directory (or build a scenario directory and use `run_all.py` for batch runs).
Receive the structured comparison dict. Read `s0_overall` and per-gate `s0`
booleans for gate disposition; read `deltas`, `rms`, `engine`, and `cw` for the
raw comparison data. For direct solver access, `solve_hydration_2d()` in
`thermal_engine_2d.py` returns a `HydrationResult` with the full temperature and
hydration field trajectories.

For full API documentation, run `help()` on each class:

```python
>>> from cw_scenario_loader import CWMixDesign, CWGeometry, CWConstruction, CWEnvironment
>>> from compare_to_cw import run_one
>>> help(CWConstruction)  # see form_type enumeration and out-of-envelope caveat
>>> help(run_one)         # see scenario_dir, keyword-only png_path, return shape
```

## Engineering history

Engine v3 is the result of six sprints of thermal-physics development. Sprints
0–2 established the 2D half-mat finite-difference solver and MIX-01 single-mix
validation; Sprint 3 explored Barber soil modeling; Sprint 4 expanded to
multi-mix validation on the 8-mix evaluation set and opened the R_form
characterization thread; Sprint 5 characterized Reference-set residuals and
closed the BC-parameter hypothesis space; Sprint 6 closed cleanup risks
(form-material parameterization R5, RuntimeWarning hygiene R7) and finalized
the customer-facing release documentation. For sprint-by-sprint detail
including hypotheses tested, ablations run, and risks routed, see
`docs/coding_passdown_v3.md` (Sprints 0–2) and `docs/coding_passdown_v4.md`
(Sprints 3–6).
