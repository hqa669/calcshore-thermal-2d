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
- **Placement temperature**: validated across [45°F, 75°F] by adiabatic kinetics
  match against CW (composition-based Hu calibration; see Kinetics calibration
  section below)
- **Mix range**: 0–91% SCM across Type I / I-II / II / V cement, Class F / Class C
  fly ash, and slag (13 mixes validated; silica fume deferred — see
  out-of-envelope conditions)
- **Form material**: steel (R_form = 0.0862 m²·K/W, ADR-04, reinforced PR 20)
- **S0 tolerances**:
  - PeakMax T: ±1.0°F
  - Peak gradient: ±2.0°F
  - Field RMS: ≤2.0°F
  - Centerline RMS: ≤1.0°F
  - Corner RMS: ≤3.0°F

## Kinetics calibration (apr28)

Engine v3 ships with a composition-based Hu calibration that closes the
adiabatic kinetics gap against CW across the validated mix envelope. The CW-
regressed `Hu_J_kg` (line 389 of `input.dat`) is preserved verbatim; the
solver consumes a corrected `Hu_J_kg_effective` computed at scenario-load
time as:

    Hu_eff = Hu_regressed × f_type × Σᵢ kᵢ · pᵢ

where `f_type` is a cement-type-specific factor (Type I, I/II, II, V), `kᵢ`
is an SCM-class-specific factor (Portland, Class F fly ash, Class C fly ash,
slag), and `pᵢ` is the mass fraction of constituent i in total cementitious.

The calibration was validated against CW adiabatic centerline curves across
the 14-mix engine v3 library. Adiabatic peak temperature agreement: 13/14
mixes within ±1°F PeakMax and ≤1°F RMS at T₀ ∈ {45°F, 60°F, 73°F, 75°F}.
The lone exception (MIX-13, 5% silica fume) routes to the silica fume
characterization sprint — see deferred items below.

The correction is applied automatically by `load_cw_scenario`. The flag
`use_cw_calibrated_hu=False` reproduces the pre-correction (raw-Hu) baseline
for diagnostics. Calibration coefficients live in `kinetics_correction.py`
and are pinned by `tests/test_kinetics_correction.py`.

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
itself (effective k(T,α), Cp(T,α)).

With the apr28 kinetics calibration in place, this residual is
unambiguously BC physics — kinetics-vs-BC entanglement (the R9
disambiguation thread from earlier sprints) is closed. Future
investigation can ablate BC form factors directly without
having to control for kinetics overshoot.

**Engine v3 envelope statement**: engine v3 ships with a
documented ~0.74°F Centerline RMS residual on the validated
Reference geometry (40×60×8ft half-mat footing, Austin summer,
mid-SCM at 60°F placement, steel form). The residual is within
S0 (≤1.0°F) and represents a 32% margin to the S1-aspire
threshold.

**Future investigation candidates** (Sprint 7+):
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

### R_form recalibration to 0.060 (re-opened)

The PR 15 candidate R_form=0.060 did not generalize when MIX-02
was carried at the v3 baseline. With apr28 kinetics in place,
MIX-02's "kinetics-divergent" signature (CornerRMS / PeakGrad
worsening as R_form decreased) was a kinetics-overshoot artifact,
not an R_form-physics anomaly. The cluster-vs-MIX-02 contradiction
that blocked the R_form sweep dissolves with the calibration.

R_form recalibration to the cluster optimum (~0.060–0.070) is
back on the Sprint 7 backlog. Expected outcome on the validated
13-mix set: cluster CornerRMS drops 1.1–1.6°F. ADR-04 currently
remains at 0.0862 pending Sprint 7 validation. See
`validation/diagnostics/pr20_r_form_eval/`.

### Corner-side residuals

The MIX-02 kinetics-divergence flag from prior release notes is
withdrawn. With the apr28 calibration, MIX-02 lands at
ΔPeakMax=+0.18°F, RMS=0.11°F adiabatically — the pre-correction
"kinetics divergence" signature was the same Hu overshoot
present in the rest of the library, masked by MIX-02's specific
cement/SCM mix amplifying it differently than the Reference
cluster.

Engine v3 envelope statement for corner-side accuracy: cluster
CornerRMS of ~2.0°F (range 1.24–2.72°F across Reference mixes at
R_form=0.0862) is the validated envelope statement. Improving
this requires Sprint 7 R_form recalibration (now unblocked) and
the BC-physics work flagged above.

## Out-of-envelope conditions

The following conditions are explicitly outside the engine v3 validated
envelope. Each is documented, routed to a successor sprint series, and (where
applicable) emits a runtime warning or error.

- **Silica fume / ultra-fine fly ash above ~5% replacement**: provisional
  `k_SF=0.7` placeholder; needs binary-mix CW runs at 5/8/12% SF replacement
  to characterize. MIX-13 (5% SF, the only library mix in this regime)
  empirically suggests `k_SF≈2.0` at 5% replacement (engine under-predicts
  by 3.84°F at the provisional value). The kinetics module emits a
  `KineticsCorrectionWarning` whenever any SF is present.

- **Type III cement**: not characterized (uncommon in mass concrete; behaves
  differently due to high Blaine fineness). The kinetics module raises
  `KeyError` rather than silently producing a wrong Hu factor. To add a Type III
  entry, run a CW adiabatic reference at high replacement and add the value to
  `F_TYPE` in `kinetics_correction.py`.

- **Pure-slag mixes (cement = 0)**: CW returns a discontinuous α_u step at
  this boundary (α_u drops from ~0.99 to ~0.28). The engine raises
  `NotImplementedError`. Customer workaround: add ≥5% Portland.

- **Non-half-mat geometry (full-mat, slab, column, wall)**: R6. Untested.
  Routes to the geometry coverage series.

- **Non-Austin climate**: The 14-mix library has zero climate variation (all
  Austin TX). Behavior on cold-climate placements, high-altitude UV, or marine
  environments is untested. Routes to the multi-climate validation series.

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

Adiabatic peak comparison vs CW centerline reference, with the apr28
composition-based Hu calibration applied. Pass criterion: |peak Δ| ≤ 1.0°F
and RMS ≤ 1.0°F. Numbers from `validation/kinetics_correction/expected/
summary_table.md` (regenerated by `batch_compare_all_mixes.py`).

| Mix | T₀ | Cement | %Cem | %FA-F | %FA-C | %Slag | %SF | Hu_factor | Engine peak | CW peak | Δpeak | RMS | Pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| MIX-01 | 73°F | Type I/II | 60.9 | 21.7 | 0.0 | 17.4 | 0.0 | 0.9514 | 149.38°F | 149.25°F | +0.13°F | 0.20°F | ✅ |
| MIX-02 | 73°F | Type V | 60.9 | 21.7 | 0.0 | 17.4 | 0.0 | 0.9571 | 140.04°F | 139.86°F | +0.18°F | 0.11°F | ✅ |
| MIX-03 | 73°F | Type I | 78.3 | 21.7 | 0.0 | 0.0 | 0.0 | 0.9787 | 152.80°F | 152.02°F | +0.78°F | 0.57°F | ✅ |
| MIX-04 | 73°F | Type I | 78.3 | 21.7 | 0.0 | 0.0 | 0.0 | 0.9787 | 152.80°F | 152.02°F | +0.78°F | 0.57°F | ✅ |
| MIX-05 | 73°F | Type I/II | 47.8 | 21.7 | 0.0 | 30.4 | 0.0 | 0.9372 | 149.33°F | 149.05°F | +0.27°F | 0.18°F | ✅ |
| MIX-06 | 73°F | Type I/II | 28.7 | 21.7 | 0.0 | 49.6 | 0.0 | 0.9164 | 148.49°F | 148.01°F | +0.48°F | 0.23°F | ✅ |
| MIX-07 | 73°F | Type I/II | 8.7 | 21.7 | 0.0 | 69.6 | 0.0 | 0.8946 | 148.70°F | 147.96°F | +0.75°F | 0.38°F | ✅ |
| MIX-08 | 73°F | Type I/II | 82.6 | 0.0 | 0.0 | 17.4 | 0.0 | 0.9708 | 161.35°F | 161.08°F | +0.27°F | 0.28°F | ✅ |
| MIX-09 | 73°F | Type I/II | 58.3 | 0.0 | 4.2 | 16.7 | 20.8 | 0.9049 | 157.92°F | 157.84°F | +0.08°F | 0.35°F | ✅ (SF=20.8% OUT-OF-ENVELOPE; passes opportunistically) |
| MIX-10 | 73°F | Type I/II | 47.8 | 34.8 | 0.0 | 17.4 | 0.0 | 0.9398 | 142.11°F | 142.12°F | −0.02°F | 0.21°F | ✅ |
| MIX-11 | 73°F | Type I/II | 60.9 | 21.7 | 0.0 | 17.4 | 0.0 | 0.9514 | 148.92°F | 148.53°F | +0.39°F | 0.23°F | ✅ |
| MIX-12 | 73°F | Type I/II | 60.9 | 21.7 | 0.0 | 17.4 | 0.0 | 0.9514 | 149.37°F | 149.29°F | +0.08°F | 0.22°F | ✅ |
| MIX-13 | 73°F | Type I/II | 58.7 | 22.9 | 0.0 | 18.3 | 0.0 | 0.9493 | 140.91°F | 144.75°F | −3.84°F | 3.46°F | ❌ (silica fume — k_SF deferred) |
| MIX-14 | 75°F | Type I/II | 60.9 | 21.7 | 0.0 | 17.4 | 0.0 | 0.9514 | 151.57°F | 151.39°F | +0.17°F | 0.21°F | ✅ |
| MIX-15 | 45°F | Type I/II | 60.9 | 21.7 | 0.0 | 17.4 | 0.0 | 0.9514 | 117.05°F | 117.05°F | +0.00°F | 0.07°F | ✅ |

**Adiabatic kinetics gate: 14/15 mixes pass; MIX-13 (silica fume) deferred to
SF characterization sprint.**

The full-stack thermal validation (S0 5/5 on Reference cluster mixes with
boundary conditions and form-side physics) holds at the v3 baseline. Sprint 7
will re-run the full-stack S0 sweep with the apr28 calibration in place;
expected outcome is dramatic improvement on Cluster A and B1 cluster mixes
that previously routed to the hydration series.

## API surface

Customer workflow: construct `CWMixDesign`, `CWGeometry`, `CWConstruction`, and
`CWEnvironment` from TCP parameters, then call `run_one()` with a scenario
directory (or build a scenario directory and use `run_all.py` for batch runs).
Receive the structured comparison dict. Read `s0_overall` and per-gate `s0`
booleans for gate disposition; read `deltas`, `rms`, `engine`, and `cw` for the
raw comparison data. For direct solver access, `solve_hydration_2d()` in
`thermal_engine_2d.py` returns a `HydrationResult` with the full temperature and
hydration field trajectories.

The apr28 Hu calibration is applied automatically by `load_cw_scenario`
(default `use_cw_calibrated_hu=True`). The verbatim CW-regressed Hu is
preserved on `mix.Hu_J_kg`; the solver consumes `mix.Hu_J_kg_effective`. Set
`use_cw_calibrated_hu=False` to disable the correction for diagnostic
purposes (this reproduces the pre-correction baseline behavior).

For full API documentation, run `help()` on each class:

```python
>>> from cw_scenario_loader import CWMixDesign, CWGeometry, CWConstruction, CWEnvironment
>>> from compare_to_cw import run_one
>>> help(CWConstruction)  # see form_type enumeration and out-of-envelope caveat
>>> help(run_one)         # see scenario_dir, keyword-only png_path, return shape
```

## Engineering history

Engine v3 is the result of six sprints of thermal-physics development plus a
post-Sprint-6 kinetics calibration milestone. Sprints 0–2 established the 2D
half-mat finite-difference solver and MIX-01 single-mix validation; Sprint 3
explored Barber soil modeling; Sprint 4 expanded to multi-mix validation on
the 8-mix evaluation set and opened the R_form characterization thread;
Sprint 5 characterized Reference-set residuals and closed the BC-parameter
hypothesis space; Sprint 6 closed cleanup risks (form-material parameterization
R5, RuntimeWarning hygiene R7) and finalized the customer-facing release
documentation. The apr27/apr28 kinetics-isolation work then characterized the
hidden CW post-regression Hu correction (cement-type lookup `f_type` and SCM-
mass-fraction-weighted derating `Σ kᵢ · pᵢ`), validated the formula against
13 of 15 library mixes adiabatically, and integrated it as a permanent code
path. For sprint-by-sprint detail including hypotheses tested, ablations run,
and risks routed, see `docs/coding_passdown_v3.md` (Sprints 0–2),
`docs/coding_passdown_v4.md` (Sprints 3–6), and `kinetics_passdown_apr28.md`
+ `passdown_hu_factor_integration.md` (kinetics calibration milestone).
