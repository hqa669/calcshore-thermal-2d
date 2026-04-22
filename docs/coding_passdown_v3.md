# CalcShore Thermal Engine v3 — Coding Passdown (2D, CW-Parity)

**Supersedes:** `coding_passdown_v2.md` (kept only for reference — sprint order has changed).

**Authored:** 2026-04-22

---

## Mission

Build CalcShore's thermal engine to produce output **indistinguishable from ConcreteWorks** for contractor-facing TCP reports. This is a commercial requirement, not an engineering preference: contractors trust CW and will not adopt a tool that gives different numbers.

## Re-sequenced sprint plan

The v2 plan put physics improvements (solar, longwave, etc.) before the 2D port. We've reversed this because:

1. Every physics improvement is easier to validate in a 2D engine that already matches CW structurally
2. Sprint 1's top-surface radiation will be written once in 2D instead of twice (1D then re-adapted)
3. The biggest CW-matching gap is dimensionality, not physics

**New order:**

| Sprint | Goal | Target metric |
|---|---|---|
| **S0 — 2D port** | Reproduce CW's grid and BC structure | Node-by-node RMS ≤ 2°F on MIX-01 |
| **S1 — Solar + longwave** | Full top-surface radiation balance | Peak max T within ±1°F on 3 mixes |
| **S2 — Convection upgrade** | ACI Eq 27 with surface-orientation constants | Diurnal amplitude within ±1°F |
| **S3 — Barber soil** | Physics-correct subgrade initialization | Pavement/column scenarios |
| **S4 — Cleanup** | Remove all remaining calibration fudge | Zero empirical factors |

**Sprint 0 is the biggest lift and the one this passdown focuses on.**

---

## Sprint 0: 2D port (what you are doing now)

### Files involved

| File | Role | Status |
|---|---|---|
| `thermal_engine_2d.py` | The new 2D engine | **TO BUILD** |
| `thermal_engine_2d_design.md` | Design spec (decisions D1–D6) | ✅ Authored |
| `cw_scenario_loader.py` | Loads CW files → dataclasses | ✅ Built & tested |
| `cw_scenario_loader.md` | Loader spec | ✅ Authored |
| `thermal_engine_v2.py` | Current 1D engine | KEEP (regression baseline) |
| `coding_passdown_v2.md` | Old sprint plan | Reference only |

### Milestone sequence (from `thermal_engine_2d_design.md`)

- **M0** — Grid builder + domain composition (no physics)
- **M1** — Pure conduction, constant properties, validate vs. analytical square-slab solution
- **M2** — Add hydration + variable k/Cp, validate centerline vs. 1D v2 (should match with adiabatic sides)
- **M3** — Add top BC (v2-style: convection + blanket R + Menzel evap), validate vs. CW MIX-01
- **M4** — Add side BCs (form + side blanket), validation target ≤ 2°F RMS vs. CW

Each milestone's validation must pass before moving to the next. Do not skip ahead.

### Decision anchors (from the design spec, not to be re-litigated)

- **D1:** Half-symmetric geometry (centerline is right edge, symmetry BC)
- **D2:** Grid = 21 × 13 for 8 ft × 40 ft mat, matching CW exactly
- **D3:** CFL-bounded inner step, 5-min output sampling (matches `temp.txt`)
- **D4:** 5-point explicit central difference, harmonic mean k at interfaces
- **D5:** Domain = blanket + concrete + soil (laterally extended)
- **D6:** New file, keep v2 alive for regression

### Validation harness (must be built alongside the engine)

Not every scenario will be CW-exported. Use these three as the S0 validation set:

1. **MIX-01 Austin Jul 15** (files already on hand: `TEST.dat`, `TX__Austin.dat`, `temp.txt`)
   - Baseline Type I/II with fly ash + slag
   - Peak CW Max T = 129.6°F @ 145.8 hr
   - Peak CW Gradient = 39.3°F @ 146.2 hr
2. **MIX-08 Austin** (high-heat, OPC + slag, no fly ash) — to be exported from CW
   - Tests the heat-generation path under stress
3. **MIX-15 Austin (45°F cold placement)** — to be exported from CW
   - Tests activation energy / early-age behavior under cold conditions

**Acceptance criteria for Sprint 0:**
- All three mixes: peak Max T within ±1°F of CW, peak Gradient within ±2°F, field RMS ≤ 2°F
- Runtime < 10 seconds per mix

### What exists and what needs to be written

**Reusable from `thermal_engine_v2.py`:**
- `arrhenius_vec()`, `hydration_rate_vec()`, `hydration_alpha_vec()` — dimension-agnostic
- `thermal_conductivity_variable()`, `specific_heat_variable()` — work on arrays
- `menzel_evaporation()`, `saturated_vapor_pressure_mmHg()` — scalar, called per surface cell
- `SOIL_PROPERTIES` dict
- Unit conversion helpers

**Must write new:**
- 2D grid builder (`build_grid_half_mat`)
- 2D material ID map (which cell is blanket / concrete / soil)
- 2D stencil application (interior + boundaries)
- BC application functions (`apply_top_bc`, `apply_corner_side_bc`, `apply_centerline_bc`, `apply_deep_ground_bc`)
- Core solver `solve_thermal_2d()`
- `compare_to_cw()` validation harness

**Must adapt:**
- Time-stepping loop structure (mostly carries over, with 2D array instead of 1D)
- Output sampling

### Drop the calibration fudge factors

Current `thermal_engine_v2.py` has:
- `Hu_factor = 1.06` — existed only to compensate for 1D's lack of lateral heat loss. **Eliminate in 2D.** Start with `Hu_factor = 1.0`; if validation is off by a few °F in peak, investigate whether it's truly heat generation or a BC issue before reintroducing a fudge factor.
- `blanket_r_effective()` log-fit — was calibrated to make 1D output match CW's 2D at R=5.67. **Eliminate in 2D.** Use the raw user-input R value directly.

**If 2D validation fails against CW without these fudges**, the next debugging step is NOT to add them back. It's to check:
1. Side BC treatment of the form-on period
2. Soil domain lateral extent (is enough soil modeled to let heat escape like CW does?)
3. Emissivity/absorptivity on the blanket top surface
4. CW's actual time step and output cadence matching ours

### Runtime targets

| Grid | Duration | Target runtime |
|---|---|---|
| 21×13 concrete + soil/blanket extensions (~2000 cells total) | 7 days | < 10 s |
| Full 15-mix validation suite | 7 days each | < 3 min |

If Python/NumPy can't hit this, fall back to Numba JIT (minor code restructure, no API change).

---

## Workflow discipline (critical)

### Claude Code prompts must include verification commands

Previous lesson from the v2 rollout: the cloud sandbox does not update the local checkout. Every prompt should end with:

```
After implementing, run:
  1. python -m pytest tests/
  2. python thermal_engine_2d.py MIX-01  # CLI smoke test
  3. python compare_to_cw.py TEST.dat TX__Austin.dat temp.txt
Verify all pass before claiming completion.
```

### Every code change requires running the validation suite

The 3-mix MIX-01/MIX-08/MIX-15 suite is the gate. Do not merge any change that regresses these.

### After merging, pull locally

```
git pull
```

The sandbox doesn't do this automatically. This was the top source of debugging confusion in v2.

### Milestone completion requires a test

M0, M1, M2, M3, M4 each require a passing test before declaring complete. Do not declare "M1 done" if the analytical validation script hasn't run.

---

## Validation suite context

### The 15-mix library (from v2, still valid)

Pass criteria in v2 was `±5°F`. For v3 the target is **`±1°F` on peak**, realized through exact CW reproduction, not calibration.

Mix library is in `coding_passdown_v2.md` § "Validation Suite (15 Mixes)". Hydration parameters (τ, β, α_u, Ea, Hu) are already embedded there and match what CW writes to `TEST.dat` lines 385–389.

**CW reference values for all 15 mixes** (peak, ΔT) are in v2 passdown. These become the acceptance gate for v3 once the 2D engine is working on MIX-01/08/15.

### New scenario needs

- Export CW `.dat` + `temp.txt` for each mix in the validation library
- Store under `validation/cw_exports/MIX-XX/` with `input.dat`, `weather.dat`, `output.txt`
- One-liner loader per scenario using `cw_scenario_loader.load_cw_scenario()`

---

## Known risks, flagged for implementation

From the design spec, repeated here for emphasis:

1. **CW's ground-beside-mat BC may differ from top-of-mat BC.** Need to confirm from CW V3 manual or just validate empirically.
2. **Form + side blanket composite BC** — decode `TEST.dat` flag pattern when writing side BC.
3. **Aspect ratio** in soil extension may cause stencil issues; check convergence.
4. **Blanket thickness (2 cm) vs. concrete dy (20 cm)** — give blanket 1 node, accept small CFL penalty (matches v2 approach).

---

## Post-Sprint 0: what Sprint 1 looks like

After 2D port matches CW, Sprint 1 replaces the top BC's current `q = h·(T_amb − T_surf) − q_evap` with:
```
q_net = q_solar_absorbed + q_LW_net + q_convection + q_evap
```
Using the hourly `env.solar_W_m2`, `env.T_air_F`, `env.cloud_cover` already loaded by `cw_scenario_loader.py`. This is a *single function replacement* in the new 2D BC module.

Sprint 1 expected gain: once our engine uses actual hourly solar/LW, it will match CW's *shape* but may diverge from CW's *numbers* because CW uses daily-averaged simplifications. This is the first point where we may choose to "surpass CW in correctness, at the cost of exact number-matching." That's a deliberate product decision to defer — for now, match CW's numbers, even if CW's physics is slightly wrong.

---

## Reference documents

Essential to load into any new coding chat for v3:

1. `thermal_engine_2d_design.md` — architecture (this sprint's authority)
2. `cw_scenario_loader.py` + `cw_scenario_loader.md` — data plumbing
3. `thermal_engine_v2.py` — 1D baseline (source of reusable components)
4. `coding_passdown_v3.md` (this file) — workflow and validation discipline
5. `TEST.dat` / `TX__Austin.dat` / `temp.txt` — MIX-01 Austin scenario for M3/M4 validation

Optional but valuable:
- `ConcreteWorks_source_documents.md` — free sources for CW physics equations
- `coding_passdown_v2.md` — legacy sprint plan (reference only; do not follow)

---

## Changelog

- **2026-04-22** — Initial v3 passdown authored. Sprint 0 (2D port) begins. Decisions locked in `thermal_engine_2d_design.md`.
