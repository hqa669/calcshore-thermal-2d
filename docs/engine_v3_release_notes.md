# CalcShore Thermal Engine v3 — Release Notes (DRAFT)

This file is in active development through Sprints 5 and 6. Each
section corresponds to a PR's contribution; the final engine v3
release content is consolidated at `sprint-6-complete`.

## Validated envelope

(to be completed in Sprint 6 with the final validation matrix)

## Known residuals

### Top-surface BC amplitude over-oscillation (Reference set, Centerline RMS)

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

**Item 1 outcome (R_form Reference-set evaluation, PR 20)**:
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
Routing: Sprint 6 R5 (form-material parameterization) to
evaluate whether a per-form-type or per-mix-regime R_form
seam can recover the cluster improvement without regressing
MIX-02. See `validation/diagnostics/pr20_r_form_eval/`.
