# thermal_engine_2d.py — Design Spec

**Status:** Design, not yet implemented. This doc authored 2026-04-22, before coding starts. Update as implementation progresses.

**Purpose:** Define the architecture and decisions for CalcShore's 2D thermal engine, which replaces the current 1D engine (`thermal_engine_v2.py`) as the primary solver. Target: produce output indistinguishable from ConcreteWorks for contractor-facing TCP reports.

---

## When to load this file

Load `thermal_engine_2d_design.md` together with `thermal_engine_2d.py` (once it exists) into any new chat where you want to:

- Extend, modify, or debug the 2D engine
- Add a new member type (column, bent cap, pavement)
- Review architectural decisions before changing them
- Hand off the engine to a new developer or AI agent

This spec documents the *why* behind the code so design intent doesn't get lost.

---

## Strategic context (why 2D, not 1D)

**Business reason:** Contractors trust ConcreteWorks. CalcShore's first users will not adopt a tool that produces different temperature predictions than CW, no matter how "physically correct" we argue the difference is. The engine must produce output that matches CW's 2D model node-for-node to be commercially viable.

**Technical reason:** The 1D engine (`thermal_engine_v2.py`) cannot reproduce:
- CW's `Min Temp in x-section` (corners are 2D phenomena)
- CW's `Max Temp Difference / Gradient` (depends on both max and min)
- Geometry sensitivity (thin mats lose heat laterally, thick mats don't — the 1D fudge factor `Hu_factor` calibrates only for one geometry at a time)

**Decision:** 2D port precedes all other planned physics improvements (solar, longwave, Barber soil, etc.). Rationale: Sprints 1–4 of the original plan all naturally produce more code in 2D (top BC + side BC instead of just top), so doing 1D first would write that code twice.

---

## Authoritative design decisions

These were agreed on 2026-04-22. **Do not change without re-discussion.**

### D1. Half-symmetric geometry (not full mat)

**Decision:** Model only half the mat cross-section, with an adiabatic BC at the centerline (left or right edge).

**Rationale:**
- CW exports half-mat grid (widths 0 → 6.1 m for a 40 ft = 12.2 m full mat)
- 2× faster solve
- Enables direct node-by-node diff against `temp.txt` for validation

**Implication:** Width direction in the solver runs from `corner (x=0)` to `centerline (x=W/2)`. Symmetry BC: `∂T/∂x = 0` at `x=W/2`.

### D2. Grid resolution: 21 × 13 (match CW exactly)

**Decision:** Width nodes = 21, depth nodes = 13, for the concrete portion of a standard 8 ft × 40 ft mat. Match CW's grid node-for-node.

**Rationale:**
- Direct diff of engine output vs. `temp.txt` with no interpolation error
- CW's grid (dx ≈ 1 ft width, 0.67 ft depth) is engineering-adequate for TCP reports
- Once we achieve node-by-node parity, we can optionally refine later

**Implication:** When geometry changes, scale grid proportionally to maintain `dx ≈ 0.305 m` width, `dy ≈ 0.2 m` depth, or use CW's exact export grid when available.

**Note:** CW uses slightly non-uniform spacing in some exports (widths jump 0.30 → 0.31 → 0.30 → 0.30). We use uniform spacing for simplicity; sub-node error from this should be <0.1°F.

### D3. Time step: CFL-bounded internally, sample at 5 min for output

**Decision:** Inner loop uses `dt_s = min(dt_user, CFL_limit)` where CFL limit is computed per time step from max thermal diffusivity across the domain. Output samples at exactly 5-minute intervals to match CW's export cadence.

**Rationale:**
- 5 min matches CW's `temp.txt` exactly (2016 timesteps over 168 hrs = 5 min)
- CFL-bounded inner step keeps us stable
- Explicit scheme, no matrix solve, simple to debug

**Implication:** For MIX-01 grid (dx ≈ 0.305 m, dy ≈ 0.2 m), CFL limit is roughly 90 s, so inner step ≈ 60–90 s. Output sampling every 300 s = every 3–5 inner steps.

### D4. Stencil: 5-point explicit central difference

**Decision:**
```
∂T/∂t = (1/ρCp) × [k × (∂²T/∂x² + ∂²T/∂y²) + Q]

T_new[i,j] = T[i,j] + dt_s × (1/(ρCp)) × [
    k_x_plus × (T[i+1,j] - T[i,j]) / dx² -
    k_x_minus × (T[i,j] - T[i-1,j]) / dx² +
    k_y_plus × (T[i,j+1] - T[i,j]) / dy² -
    k_y_minus × (T[i,j] - T[i,j-1]) / dy² +
    Q[i,j]
]
```

With harmonic mean k at cell interfaces (already used in 1D v2 for blanket/concrete and concrete/soil boundaries).

**Rationale:**
- Same approach as 1D v2, just extended to two dimensions
- No sparse matrix solve required (keep it debuggable)
- CW's time step (~5 min output, likely ~1 min inner step) matches CFL-stable explicit for this grid

**Not chosen:** Crank-Nicolson implicit. Too hard to debug for the marginal benefit.

### D5. Domain composition

**Structure in memory (half-mat, with corner at origin):**
```
   width →                           centerline (symmetry BC)
   ┌─────────────────────────────┐
   │         blanket             │   ← thermal mass + R (covers concrete only)
   ├─────────────────────────────┤   ← blanket-concrete interface
   │                             │
   │          concrete           │   ← hydration heat, variable k and Cp
   │                             │
   ├──────────────────┬──────────┤   ← concrete-soil interface
   │       soil       │   soil   │   ← soil extends 1.5 × mat depth laterally
   │                  │          │      (to give lateral heat loss room)
   │                  │          │
   │       soil       │   soil   │
   └──────────────────┴──────────┘   ← deep ground (fixed T_gw)
                    ↑
                mat edge (soil extends past here)
```

**Key dimensions (for 8 ft × 40 ft mat):**
- Half-mat width = 6.1 m (20 ft) — concrete occupies full width
- Concrete depth = 2.44 m (8 ft)
- Blanket thickness = 0.02 m (CW default) — above concrete only
- Soil depth below concrete = 3.0 m (CW Barber-model equilibration depth)
- Soil lateral extension beyond concrete edge = ~3.7 m (1.5 × mat depth)

**BCs (half-mat):**
- **Top of blanket** (above concrete): convection + radiation to ambient + evaporation (suppressed when plastic cure is in place per `TEST.dat` L478)
- **Top of soil** (beside concrete, where there is no concrete above): convection + radiation to ambient (ground surface BC)
- **Corner side of concrete + blanket** (x=0, between concrete top and bottom): convection + radiation; during form-on period (0 to `form_removal_hrs`), through the steel form as a thermal layer with the side cure blanket behind it
- **Centerline side (x=W/2):** adiabatic `∂T/∂x = 0` (symmetry)
- **Deep ground (bottom of soil domain):** fixed at `T_gw = 0.83 × T_aat + 3.7 °C`
- **Corner side of soil (x=0):** far-field condition — assume `T = T_gw` at x=0 in soil (soil extends to deep ground temperature far from concrete)

### D6. New file, keep v2 alongside for regression

**Decision:** Write `thermal_engine_2d.py` as a new file. Do not modify `thermal_engine_v2.py`. Keep v2 available so we can run the 15-mix validation suite against both engines during the transition.

**Rationale:**
- v2 is production-validated (15/15 mixes pass at ±5°F). Breaking it during 2D development is unacceptable.
- Regression testing: any claim about the 2D engine can be compared to v2's output as a baseline.
- Once 2D is proven, v2 becomes the fallback for 1D-only cases (pavements) or gets deprecated.

---

## Module structure

```
thermal_engine_2d.py
├── Constants (R_GAS, T_REF_K, etc.)
├── Dataclasses (MixDesign2D, Geometry2D, Construction2D, Environment2D)
│   └── Matching cw_scenario_loader.py field names for drop-in compatibility
├── Material property models (imported/shared with v2 where possible)
│   ├── thermal_conductivity_variable(k_uc, alpha_node)     # CW Eq 23
│   ├── specific_heat_variable(...)                          # CW Eq 24-25
│   ├── hydration_rate_vec(te, tau, beta, au)
│   ├── arrhenius_vec(T_K, Ea)
│   └── menzel_evaporation(...)                              # CW Eq 38-39
├── Grid builder
│   ├── build_grid_half_mat(geom, n_concrete_x, n_concrete_y,
│   │                        n_soil_x_ext, n_soil_y, n_blanket)
│   │   → returns x[], y[], material_id[i,j], blanket/concrete/soil masks
│   └── Material ID per cell: 0=blanket, 1=concrete, 2=soil
├── BC layer
│   ├── apply_top_bc(T, t, ...)          # convection + rad + evap on top
│   ├── apply_corner_side_bc(T, t, ...)  # x=0 side: conv + rad, form-aware
│   ├── apply_centerline_bc(T)           # x=W/2 adiabatic (ghost node)
│   └── apply_deep_ground_bc(T)          # y=y_max fixed T_gw
├── Core solver
│   └── solve_thermal_2d(mix, geom, constr, env, ...)
│       → returns {
│           't': (n_time,),
│           'T_field_F': (n_time, ny, nx),     # °F, full domain
│           'T_concrete_F': (n_time, ny_c, nx_c), # concrete subdomain only
│           'T_max_xs_F': (n_time,),            # max over concrete
│           'T_min_xs_F': (n_time,),            # min over concrete
│           'T_diff_xs_F': (n_time,),           # max - min
│           'T_ambient_F': (n_time,),
│           'peak_max_F', 'peak_min_F', 'peak_diff_F', etc.
│         }
├── Diagnostic / plotting helpers
│   └── compare_to_cw(engine_output, cw_validation)
│       → plots and metrics side-by-side vs. cw_scenario_loader output
└── CLI entry point (__main__)
```

---

## Interface contract with cw_scenario_loader.py

The 2D engine must consume `CWScenario` objects directly:

```python
from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import solve_thermal_2d

scn = load_cw_scenario("TEST.dat", "TX__Austin.dat", "temp.txt")
results = solve_thermal_2d(scn.mix, scn.geometry, scn.construction, scn.environment)

# Validation
from thermal_engine_2d import compare_to_cw
compare_to_cw(results, scn.cw_validation)
```

**Requirements on field names:**
- `mix.tau_hrs`, `mix.beta`, `mix.alpha_u`, `mix.activation_energy_J_mol`, `mix.Hu_J_kg` — hydration params
- `mix.thermal_conductivity_BTU_hr_ft_F`, `mix.aggregate_Cp_BTU_lb_F`
- `mix.concrete_density_lb_ft3` (computed property on CWMixDesign)
- `geom.width_ft`, `geom.depth_ft`, `geom.length_ft`
- `constr.placement_temp_F`, `constr.blanket_R_value`, `constr.form_removal_hrs`, etc.
- `env.T_air_F`, `env.RH_pct`, `env.solar_W_m2`, `env.wind_m_s`, `env.cloud_cover` — hourly arrays

This matches what `cw_scenario_loader.py` already produces. No adapter needed.

---

## Output matching against CW — validation strategy

The CW `temp.txt` contains `(nD=13, nW=21)` grid at 5-min intervals. Our 2D engine must produce the matching grid for direct node-by-node comparison.

**Validation metrics (in priority order):**

1. **Peak Max T in x-section** — primary contractor metric. Target ≤ 1°F deviation from CW.
2. **Peak Max Temp Difference (Gradient)** — DOT compliance metric. Target ≤ 2°F deviation.
3. **Node-by-node RMS across full 2D field × all timesteps** — "are we solving the same PDE?" Target ≤ 2°F.
4. **Time-series at centerline core** — catches diurnal phase errors. Target ≤ 1°F RMS.
5. **Time-series at corner (top-surface, x=0)** — catches surface BC errors. Target ≤ 3°F RMS (hardest).

**Reporting template:**
```
MIX-01 Austin 2026-07-15  (168 hr run, 2016 timesteps, 273 nodes)
  Peak Max T:        Engine 129.2°F | CW 129.6°F | Δ = -0.4°F   PASS
  Peak Gradient:     Engine  38.1°F | CW  39.3°F | Δ = -1.2°F   PASS
  Field-wide RMS:    1.4°F  (target ≤ 2.0°F)    PASS
  Centerline core:   0.6°F  (target ≤ 1.0°F)    PASS
  Corner surface:    2.1°F  (target ≤ 3.0°F)    PASS
  OVERALL:                                                     PASS
```

Until these metrics all pass on MIX-01, no other engine changes (solar, longwave, etc.) should be merged.

---

## Performance budget

**Grid size:** 21 width × 13 depth × ~25 blanket/soil extension in y × ~5 soil extension in x ≈ 2000 cells total. 2016 timesteps. Vectorized NumPy.

**Target runtime: < 10 seconds per 7-day run.** Goal is that a full 15-mix validation suite runs in < 3 minutes.

If explicit stencil turns out to be too slow due to CFL-bounded `dt`, we can consider:
- Numba JIT on the inner loop (no dependency changes, minor code restructure)
- ADI (alternating-direction implicit) — doubles the complexity but allows 10× larger `dt`
- Cython for hot paths — last resort

---

## Known design risks

### Risk 1: CW's BC for the ground next to the mat
CW models a chunk of ground beside the mat that's exposed to ambient — this is a 2D phenomenon our 1D engine didn't have to deal with. Open question: does CW use the same top BC (convection + solar + LW) on this exposed ground as on the concrete? Need to read CW V3 manual §3.2 or Riding (2007) Ch. 4 to confirm.

**Mitigation:** For initial 2D port, use the same top BC for all top-surface cells regardless of material (blanket/concrete vs. soil). If validation shows the corner region is off, revisit.

### Risk 2: Form BC during 0–168 hrs
Steel forms have negligible R-value but high emissivity and can conduct significant heat to the side cure blanket. The current CW model likely treats the composite "form + blanket" as a single R-value layer. Need to decode `TEST.dat` L47-53 flag pattern to understand how CW handles this.

**Mitigation:** Start with "side BC = blanket R-value + form R-value (≈0) + convection"; refine based on validation.

### Risk 3: Aspect ratio and convergence
The concrete subdomain is 20 ft wide × 8 ft deep = aspect ratio 2.5. Soil extends further. If we use uniform dx × dy, the soil cells will have poor aspect (very wide compared to deep). Test for convergence issues at the corner region.

**Mitigation:** Allow different dx/dy; choose them to be similar (~0.3 m each) after the grid builder designs the domain.

### Risk 4: Blanket is much thinner than concrete/soil
Blanket thickness is 0.02 m; concrete dy is ~0.2 m. If we give the blanket its own y-layer in the grid, dy_blanket << dy_concrete violates CFL globally.

**Mitigation options:**
- (a) Treat blanket as an R-value surface BC (thermal mass ignored) — simpler, may lose ~1°F accuracy
- (b) Give blanket 1 node, accept small CFL penalty — currently what 1D v2 does
- (c) Nested time stepping (blanket steps faster than concrete/soil) — complex

**Recommend (b) initially**, matching v2's approach. Refine if validation shows blanket is a significant error source.

---

## Known issues / limitations (to fill in during implementation)

(Empty for now. Populate as issues are discovered during coding.)

---

## Milestones

- **M0** — Grid builder + domain composition, unit-tested (no physics yet)
- **M1** — Pure conduction with no hydration, no radiation, constant k/Cp. Validate against analytical solution for a square slab with fixed BCs.
- **M2** — Add hydration heat generation + variable k/Cp. Compare centerline profile against 1D v2 (should match exactly if side BCs are adiabatic).
- **M3** — Add top BC (current v2-style: convection + R-value blanket + Menzel evap). Compare against CW for MIX-01.
- **M4** — Add side BCs (form + side blanket). Validation target: node-by-node RMS ≤ 2°F against CW.
- **M5** — Replace top BC with full radiation balance (Sprint 1 from the original plan). Validation: tighten RMS target.

Sprint 1–4 of the original plan continue from M5.

---

## Changelog

- **2026-04-22** — Initial design spec authored. No code yet. Decisions D1–D6 locked.

### Clarification: vertex-centered interface node ownership (resolved in M0)

At material interfaces between blanket/concrete and concrete/soil, the
interface node is shared between the two materials but `material_id` must
label it as one or the other. Convention adopted in M0:

- **Blanket/concrete interface** (y = t_blanket): node belongs to
  **concrete**. `material_id[iy_concrete_start, ix] = 1`.
- **Concrete/soil interface** (y = t_blanket + D_concrete): node belongs
  to **concrete**. `material_id[iy_concrete_end, ix] = 1`.

As a consequence:
- `iy_concrete_start == iy_blanket_end` (both are the interface row)
- `iy_blanket_end` is best read as "the y-index of the blanket/concrete
  interface" rather than "the last blanket row"
- Blanket has exactly `n_blanket` rows of strictly-blanket cells
  (rows 0 to n_blanket - 1 at concrete columns)
- Concrete has exactly `n_concrete_y` rows of concrete cells including
  both interface nodes

**Implication for M1/M2:** when the stencil applies hydration heat `Q` to
concrete cells, it includes both interface rows. When the harmonic-mean
`k` is computed at the blanket/concrete interface, the interface node's
`k` is the concrete `k`, and the blanket `k` enters only through the
upward face of that cell.
