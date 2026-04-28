# CW Adiabatic Reference Curve — MIX-01 @ 73°F

## Purpose

This reference curve is CW's prediction of the **pure hydration heat
temperature rise** for MIX-01 (39.1% SCM Reference mix) at 73°F initial
temperature, with all boundary conditions neutralized. It is the kinetics
ground truth for validating CalcShore engine v3's hydration ODE
integration.

## What this curve represents

The Schindler-Folliard kinetics ODE integrated forward in time at adiabatic
conditions, starting from T₀ = 73°F:

```
dα/dt = β·(τ/te)^β · (1/te) · α_u · exp(-(τ/te)^β) · arrhenius(T)
dT/dt = (Hu · Cc / (ρ·Cp)) · dα/dt
```

with no heat loss to any boundary. The resulting T(t) curve is a function of
mix-specific parameters (τ, β, α_u, Hu, Ea, Cc, ρ, Cp) and the integrated
Arrhenius temperature acceleration.

## How CW was configured to produce it

To force CW to produce a pure adiabatic curve, all BC heat-loss paths were
neutralized:

### Geometry
- Member: **200 × 200 × 200 ft** rectangular footing
- Analysis: **2D** (third dimension treated as infinite)
- The centerline mid-depth (~30 m from any face) is mathematically
  insulated from any boundary effect over a 168 hr window

### Construction
- **Placement temperature: 73°F** (manual, not calculated from constituents)
- Soil temperature: 73°F (zero gradient at bottom)
- Form type: Steel (matters less than usual since BCs are neutralized)
- Form removal: 168 hr (form on entire run)
- Wet curing blanket on top + sides for full 168 hr
- Blanket R-value: 5.67 (default; doesn't matter due to geometry)
- Footing sides shaded: True (kills side solar)

### Raw material temperatures (Image 1 of CW dialog)
- Cementitious material: 73°F
- Water: 73°F
- Coarse aggregate: 73°F
- Fine aggregate: 73°F

(These don't matter when "Manually enter concrete fresh temperature" is
selected, but were set to 73°F to remove any ambiguity.)

### Environment (manually entered, all 8 days)
- **Ambient air temperature: 73°F constant** (per-day max = min = 73°F)
- Wind speed: 1 mph (minimal forced convection)
- Cloud cover: 100% (suppresses solar AND zeros sky temperature gradient)
- Relative humidity: 50% / 50% (max / min, constant)
- Yearly average temperature: 73.4°F (gives T_gw ≈ 73°F via CW Eq 44:
  T_gw = 0.83·T_aat + 3.7°C)

### Mix design (matches CalcShore engine MIX-01 reference)
- 39.1% SCM, Type I/II cement
- Schindler-Folliard parameters and Hu come from CW's regression
  (read from input.dat lines 385-389 of the resulting scenario)
- Coarse aggregate: Limestone
- Fine aggregate: Siliceous River Sand

## File contents

`cw_adiabatic_reference_mix01.csv` — 2016 rows, 5-min sampling, 168 hr total.

| Column | Units | Description |
|---|---|---|
| `time_hrs` | hours | Time since placement |
| `T_center_F_adiabatic` | °F | **Use this**: centerline mid-depth (clean adiabatic) |
| `T_max_xs_F` | °F | Max temp anywhere in cross-section (BC-contaminated for first ~10 hr) |
| `T_ambient_F` | °F | CW's reported ambient (constant 73°F, sanity check) |

## Which column to use

**Use `T_center_F_adiabatic`** for kinetics validation.

The `T_max_xs_F` column shows the maximum temperature anywhere in the
cross-section. For the first ~10 hr, the maximum occurs at the form face
(width index = corner) where some residual solar heating leaks in despite
100% cloud cover. After ~10 hr, the deep interior overtakes the form face
and `T_max_xs_F` converges to `T_center_F_adiabatic`. The two columns differ
by up to ~5°F in the first hour and converge to within 0.02°F by t=168 hr.

`T_center_F_adiabatic` is the geometric center of the 200×200×200 ft
member — physically too deep for any BC to reach in 168 hr. Its early-time
flatness at 73°F (no rise until ~5 hr) is **physically correct** kinetics
behavior: dα/dt at te ≈ 0.01 hr is essentially zero, and the curve rises
only after Arrhenius-accelerated hydration ramps up.

## Key benchmarks

Use these as reference points when comparing engine vs CW:

| t (hr) | T_center (°F) | Notes |
|---|---|---|
| 0 | 73.00 | Initial |
| 6 | 74.39 | Kinetics ramping |
| 12 | 82.99 | Inflection region |
| 24 | 105.71 | Steep rise |
| 48 | 129.97 | Past most-active hydration |
| 72 | 139.53 | Tail-off begins |
| 168 | 149.25 | Final (still rising slowly) |

**Peak T_center at t=168 hr: 149.25°F** (rise of 76.25°F above placement)

## Important constraint: this curve is for T₀ = 73°F ONLY

The curve is the result of integrating temperature-dependent kinetics. It
**cannot** be reused as a hydration template for runs with different
initial or ambient temperatures.

To validate the engine across a temperature range, you would need to
generate additional adiabatic CW runs at, e.g., T₀ = 50°F, 60°F, 85°F, 100°F.

Each curve at different T₀ has a different shape (timing, slope, peak),
because the integrated Arrhenius effect changes with the temperature
trajectory.

For real-world (non-adiabatic) simulations, the engine must integrate the
kinetics ODE forward in time using the current cell temperature. This
curve is a **validation reference**, not a substitute for that integration.

## Engine-side validation procedure

Procedure for validating engine v3 against this curve:

1. Load the same mix design from the validated MIX-01 input.dat
2. Build a small grid (any size — doesn't matter for adiabatic mode)
3. Set initial temperature to 73°F uniform
4. Call `solve_hydration_2d(grid, mix, T_initial, duration_s=168*3600,
   boundary_mode="adiabatic", ...)`
5. Take any cell from the resulting `T_field_C` (all cells have the same
   trajectory in adiabatic mode by construction)
6. Compare against `T_center_F_adiabatic` from this CSV

Pass criteria (proposed):
- Peak T at t=168 hr: ±1.0°F (CW = 149.25°F → engine in [148.25, 150.25])
- RMS over [0, 168] hr: ≤ 1.0°F
- Time to T = 100°F: ±2 hr
- Time to T = 130°F: ±2 hr

If the engine matches within tolerance, the kinetics ODE integration is
sound and any residuals in real-world runs live in the BC physics
(convection, solar, LW radiation, blanket model).

If the engine diverges, the source is in:
- `arrhenius_vec` (activation energy interpretation)
- `hydration_alpha_vec` / `hydration_rate_vec` (S-F kinetics)
- `specific_heat_variable` (Van Breugel Cp model)
- `thermal_conductivity_variable` (CW Eq 23)
- Imperial-to-SI unit conversions in `solve_hydration_2d` setup
- Or: the input.dat values for Ea, τ, β, α_u, Hu themselves are wrong

## Mix design parameters needed for engine

These come from the input.dat of the CW scenario that produced this curve.
Read these from the .dat file directly (lines from `cw_scenario_loader.py`
CW_DAT_INDEX dict):

- L385: `activation_energy_J_mol`
- L386: `tau_hrs`
- L387: `beta`
- L388: `alpha_u`
- L389: `Hu_J_kg`
- L54: `cement_lb_yd3`
- L55: `water_lb_yd3`
- L56: `coarse_agg_lb_yd3`
- L57: `fine_agg_lb_yd3`
- L58: `air_content_pct`
- L60-63: SCM amounts
- L404: `thermal_conductivity_BTU_hr_ft_F`
- L405: `aggregate_Cp_BTU_lb_F`
- L438: `placement_temp_F` (should be 73)

Or use `cw_scenario_loader.load_cw_scenario(input_dat=..., weather_dat=None,
cw_output_txt=...)` to populate a CWMixDesign dataclass directly.
