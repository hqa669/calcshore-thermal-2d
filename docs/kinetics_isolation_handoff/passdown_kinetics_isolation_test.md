# Passdown: thermal_engine_test.py — Kinetics Isolation Validation

## Context for the next chat

You are picking up work on **CalcShore engine v3**, a 2D explicit
finite-difference thermal solver for mass concrete temperature prediction.
The engine is calibrated against ConcreteWorks (CW). Engine v3 ships at
S0 5/5 on the 5 Reference mixes (39.1% SCM, Austin TX summer, 60°F
placement, steel form, 40×60×8 ft half-mat) but has documented residuals:

- ~0.74°F Centerline RMS (BC structural amplitude over-oscillation)
- ~2°F Corner RMS cluster (form-face residual)
- MIX-02 kinetics-divergent on R_form sweeps
- Cluster A high-SCM mixes (52-91% SCM) out-of-envelope
- MIX-15 (45°F) and MIX-14 (75°F) out-of-envelope on placement temp

Per the v3 release notes: "R9 (kinetics-vs-inherited-calibration
disambiguation) is correspondingly deferred to engine v3 investigation
alongside the BC structural work."

**This task closes that disambiguation.** We have decoupled kinetics from
all BC physics by running a CW adiabatic-isolation scenario and extracting
the pure hydration heat curve. We now want to validate that the engine's
kinetics ODE integration matches CW's, independent of any thermal BC.

## Goal

Write `thermal_engine_test.py` that runs the engine in adiabatic mode at a
constant 73°F ambient/initial temperature, extracts the centerline
temperature trajectory, and compares it against a pre-computed CW
reference curve.

If the engine matches CW within tolerance, kinetics is verified clean,
and any residuals in the validated 14-mix envelope are unambiguously in
the BC physics (top BC amplitude, form-face contact resistance, soil
boundary, etc.) — not in kinetics.

## What the test script needs to do

1. Load mix design parameters from a CW input.dat. The relevant fields:
   - `activation_energy_J_mol` (line 385 of input.dat)
   - `tau_hrs` (line 386)
   - `beta` (line 387)
   - `alpha_u` (line 388)
   - `Hu_J_kg` (line 389)
   - `cement_lb_yd3`, `water_lb_yd3`, aggregates, SCMs, air content
   - `thermal_conductivity_BTU_hr_ft_F` (line 404)
   - `aggregate_Cp_BTU_lb_F` (line 405)
   - Use `cw_scenario_loader.load_cw_scenario(input_dat_path,
     weather_dat=None, cw_output_txt=None)` — weather/output are optional
     when only the mix is needed.
2. Build a small grid (any size — adiabatic mode produces uniform field,
   so even a 3×3 grid works). Use the half-mat builder from
   `thermal_engine_2d.py` (`build_grid_2d` or similar).
3. Set initial temperature uniformly to 73°F.
4. Call `solve_hydration_2d(grid, mix, T_initial_C=22.78,  # 73°F
   duration_s=168*3600, boundary_mode="adiabatic", ...)`.
5. Extract any cell trajectory (in adiabatic mode, all cells are
   identical by construction).
6. Compare against `T_center_F_adiabatic` from
   `cw_adiabatic_reference_mix01.csv`:
   - Plot engine vs CW on the same axes
   - Compute peak deviation, RMS error, time-to-100°F, time-to-130°F
   - Report pass/fail against proposed tolerances

## Proposed pass tolerances

- Peak T at t=168 hr: ±1.0°F (CW=149.25°F → engine in [148.25, 150.25])
- RMS over [0, 168] hr: ≤ 1.0°F
- Time to T=100°F: ±2 hr
- Time to T=130°F: ±2 hr

## Files to attach to the new chat

You need 5 files, plus possibly the input.dat from the hydration scenario.

1. **`thermal_engine_2d.py`** — the engine source (~2400 lines). The
   relevant entry points are:
   - `solve_hydration_2d()` — main solver
   - `boundary_mode="adiabatic"` — option that zero-fluxes all four edges
   - `build_grid_2d()` — grid builder
   - Material property functions: `arrhenius_vec()`,
     `hydration_alpha_vec()`, `hydration_rate_vec()`,
     `specific_heat_variable()`, `thermal_conductivity_variable()`
2. **`cw_scenario_loader.py`** — for loading CW input.dat into the
   `CWMixDesign` / `CWGeometry` / `CWConstruction` / `CWEnvironment`
   dataclasses. The relevant entry point is `load_cw_scenario()`.
3. **`cw_adiabatic_reference_mix01.csv`** — CW kinetics ground truth.
   2016 rows, 5-min sampling, 168 hr. Columns:
   `time_hrs, T_center_F_adiabatic, T_max_xs_F, T_ambient_F`.
   **Use the `T_center_F_adiabatic` column** as the comparison reference.
4. **`cw_adiabatic_reference_mix01_README.md`** — full documentation of
   how the CW curve was generated, why it's a clean kinetics reference,
   and the temperature-dependence constraint (the curve is valid for
   T₀=73°F only — cannot be reused at different ambient temps).
5. **`engine_v3_release_notes.md`** — for context on validated envelope,
   known residuals, and out-of-envelope conditions.

**Optional but recommended**: the `input.dat` from the MIX-01 hydration
scenario directory (`HydrationCenter_mix01/input.dat`). This is what was
used to generate the CW reference curve, and reading the kinetics
parameters directly from it ensures the engine consumes the exact same
Ea, τ, β, α_u, Hu values CW used. Without it, you'd have to use a
generic MIX-01 input.dat from the existing 14-mix library — should
produce the same numbers, but worth confirming.

## Key facts the new chat needs to know

### Kinetics math the engine uses

Schindler-Folliard hydration model:

```
α(te) = α_u · exp(-(τ/te)^β)
dα/dte = β · (τ/te)^β · (1/te) · α(te)
te = ∫₀ᵗ arrhenius(T(t')) dt'
arrhenius(T) = exp(-Ea/R · (1/T - 1/T_ref))
T_ref = 296.15 K  (= 23°C, matches CW's Tr=533°R)
Q_hyd = Hu · Cc · dα/dt  [W/m³]
```

Where:
- `Cc` = cementitious content (kg/m³ of mix)
- `Hu` = ultimate heat of hydration per kg cement (J/kg)
- `te` = equivalent age at reference temperature

### Adiabatic mode behavior

In `boundary_mode="adiabatic"`, all four edges have zero flux. With
uniform initial temperature, the field stays spatially uniform for all
time. Every cell follows the local kinetics ODE:

```
ρ Cp dT/dt = Q_hyd
```

So `T_field_C[t, :, :]` is uniform at every t, and any cell trajectory
is the engine's adiabatic temperature curve.

### Why this validation matters

The v3 release notes flag two entangled residuals — the BC structural
amplitude residual (~0.74°F Centerline RMS) and the kinetics-divergent
mixes (MIX-02, Cluster A, cold/warm placement). Per the release notes,
"R9's kinetics-vs-inherited-calibration disambiguation is correspondingly
deferred to engine v3 investigation alongside the BC structural work."

This test cleanly disentangles them:
- If engine matches CW adiabatically → kinetics is correct, all
  residuals live in BC physics. Future work focuses on the
  Berdahl-Martin sky temperature, blanket layer coupling, h_conv
  parameterization, etc.
- If engine diverges from CW adiabatically → kinetics has an issue.
  Different signatures isolate which parameter:
  - Wrong peak amplitude → Hu, Cc, α_u, or Cp(α,T)
  - Wrong inflection timing → τ, β
  - Wrong steepness → Ea (Arrhenius)
  - Wrong asymptote → α_u

### Watch for unit conversions

The engine works in SI (kg, m, s, K, J, W). The CW reference curve and
input.dat use Imperial (lb, ft, hr, °F, BTU). Conversions to verify:

- `T°C = (T°F - 32) × 5/9`
- `T_K = T°C + 273.15`
- 73°F = 22.78°C = 295.93 K
- The engine has `R_IMP_TO_R_SI = 0.1761` for blanket R-values
  (hr·ft²·°F/BTU → m²·K/W) — not relevant for adiabatic mode but worth
  knowing
- `T_REF_K = 296.15` (= 23°C) is hardcoded in the engine to match CW's
  `Tr = 533°R`

### CW's grid convention

CW exports half-mat fields with descending widths:
- `widths_m[0]` = form face (largest width, e.g., ~30 m for 200 ft)
- `widths_m[-1]` = centerline (0 m)

So in the `T_field_F` array of shape `(n_time, nD, nW)`:
- `T_field_F[:, :, 0]` is the form face column
- `T_field_F[:, :, -1]` is the centerline column

In the reference CSV, `T_center_F_adiabatic` comes from
`T_field_F[:, nd//2, 0]`. Wait — that doesn't match what I just said. Let
me clarify: in the CW scenario for the adiabatic isolation, with
`widths_m_unique = sorted(..., reverse=True)`, `widths_m[0]` IS the form
face and `widths_m[-1]` IS the centerline (0 m). So
`T_field_F[:, nd//2, 0]` is actually the form-face mid-depth, not the
centerline.

**HOWEVER**, the c_mid extraction in the original adiabatic_diagnostic.csv
used index 0 for the width axis, but the diagnostic showed c_mid at the
adiabatic core temperature (149.25°F) — meaning the actual extraction was
the centerline. This implies either (a) CW's loader inverts the width
axis after parsing, or (b) the loader's `widths_m_unique` ordering is
different from what's stored in `T_field_F`. Worth checking
`parse_cw_temp_output()` carefully when building the new test.

For practical purposes: **the column labeled `T_center_F_adiabatic` in
the reference CSV is genuinely the centerline mid-depth adiabatic curve**
— that's been verified by the spatial-uniformity check and by the 200×200
ft heatmap visualization showing uniform interior. Trust the CSV column,
not the index notation.

## Validation procedure: step-by-step

1. **Setup**: copy `thermal_engine_2d.py`, `cw_scenario_loader.py`, and
   the reference CSV into a working directory. Place
   `cw_adiabatic_reference_mix01_README.md` alongside for reference.

2. **Load mix parameters** from a MIX-01 input.dat. If using the existing
   14-mix library, point to that scenario. If using the
   adiabatic-isolation scenario's input.dat, point to that.

3. **Build the smallest viable grid**: e.g., 3×3 cells of concrete with
   no soil/blanket extensions. In adiabatic mode the geometry doesn't
   matter — the result is a single uniform trajectory.

4. **Run the engine**:
   ```python
   from thermal_engine_2d import solve_hydration_2d, build_grid_2d
   grid = build_grid_2d(geom, ...)  # any geometry
   result = solve_hydration_2d(
       grid, mix,
       T_initial_C=(73-32)*5/9,  # 22.78°C
       duration_s=168*3600,
       boundary_mode="adiabatic",
       sample_dt_hrs=5/60,  # 5-min sampling to match CW reference
       ...
   )
   T_engine_C = result.T_field_C[:, 0, 0]  # any cell
   T_engine_F = T_engine_C * 9/5 + 32
   t_hrs_engine = result.sample_times_hrs
   ```

5. **Load the CW reference**:
   ```python
   import numpy as np
   ref = np.genfromtxt('cw_adiabatic_reference_mix01.csv',
                       delimiter=',', skip_header=1)
   t_cw = ref[:, 0]
   T_cw_F = ref[:, 1]  # T_center_F_adiabatic column
   ```

6. **Interpolate one onto the other's time grid**:
   ```python
   T_engine_on_cw = np.interp(t_cw, t_hrs_engine, T_engine_F)
   ```

7. **Compute metrics**:
   ```python
   delta = T_engine_on_cw - T_cw_F
   peak_delta = T_engine_on_cw.max() - T_cw_F.max()
   rms = np.sqrt(np.mean(delta**2))
   t_engine_to_100F = np.interp(100, T_engine_on_cw, t_cw)
   t_cw_to_100F = np.interp(100, T_cw_F, t_cw)
   ```

8. **Plot** engine vs CW on the same axes, residual on a secondary
   y-axis or as a separate subplot.

9. **Report pass/fail** against the proposed tolerances.

## Things to anticipate

- **Engine might not have a 5-min sample resolution**. Default is
  probably hourly. May need to add a `sample_dt_hrs` parameter to
  `solve_hydration_2d` if it doesn't exist.
- **Adiabatic mode might not skip the BC code paths cleanly** — verify
  by inspecting `solve_hydration_2d` that when `boundary_mode="adiabatic"`,
  the top, bottom, side BC blocks are actually bypassed. If not, the
  test isn't truly testing kinetics in isolation.
- **The engine reads kinetics parameters from the `CWMixDesign`
  dataclass**. Verify that what the loader extracts from input.dat is
  actually what gets passed to the kinetics functions (no overrides or
  recomputations between).
- **Initial transient**: kinetics ramp up slowly at 73°F due to
  Arrhenius being ~0.78× at that temp vs ref. Both engine and CW should
  show ~5 hr of essentially-flat T before the rise begins.

## After this test passes (or fails)

If the test **passes** (engine matches CW within tolerance):
- Kinetics is verified clean
- Routes the v3 residuals (Centerline RMS 0.74°F, Corner RMS 2°F)
  unambiguously to BC physics work
- Opens up the next sprint: BC structural investigation, blanket layer
  coupling, sky temperature parameterization

If the test **fails** (engine diverges from CW):
- Kinetics has a real issue — track which parameter signature
- Check input.dat parsing: are Ea, τ, β, α_u, Hu values the same in
  the dataclass as in the file?
- Check unit conversions in the engine setup
- Check Arrhenius implementation: T_ref=296.15 K, R=8.314 J/(mol·K)
- Check the integration scheme for `te` (equivalent age)
- This would also explain MIX-02 and Cluster A divergences in the
  validated envelope — they may not be exotic kinetics anomalies but
  baseline kinetics errors that are masked at the validated 60°F
  placement and only show up at extreme conditions

## Out-of-scope for this test

- Comparing different placement temperatures (would require additional
  CW reference curves at T₀ ≠ 73°F)
- Comparing different mix designs (this validates MIX-01 only)
- Re-validating the BC physics (separate work, deferred)
- Changes to engine kinetics parameters (this is a comparison test,
  not a calibration)

## Future expansion

Once MIX-01 passes adiabatic validation, generate additional CW reference
curves at:
- T₀ = 60°F (the validated envelope placement temp — ties this back to
  the existing 14-mix library)
- T₀ = 50°F and 100°F (Arrhenius range coverage)
- For each of MIX-02, MIX-04, MIX-08 (other Reference mixes plus the
  diverged one)

A pass on all of those gives full kinetics coverage of the validated
envelope and the kinetics-divergent mix.
