# PR 16 — H4 Grep Audit (MIX-15 B2)

Phase 2.5 diagnostic for Sprint 4 PR 16.
Audit of `placement_temp_F`, `T_initial`, `T_0`, `T0_C`, `T0_F` usage across
`thermal_engine_2d.py`, `compare_to_cw.py`, `cw_scenario_loader.py`.

---

## What We Know (Phase 2 Preamble)

Before classifying grep matches, two findings from the verification battery sharply
constrain the search space.

**V2: MIX-15 and MIX-01 are byte-identical across all mix-design fields.**
From the verification battery (`pr16_h6_verification.py`), every field compared was
identical to 4+ decimal places: cement_type_I_II (350.0), fly_ash_F (125.0),
fly_ash_C (0.0), ggbfs (100.0), water (253.0), Hu_J_kg (424143.13), tau_hrs (29.401),
beta (0.895), alpha_u (0.7585), Ea (26457.949). Combined with §7.5.2's finding that
`construction.placement_temp_F` is the sole varying construction field, this means:

> **The entire difference between MIX-15 and Reference (MIX-01) is placement_temp_F.**

Any code path where a different value propagates from this single parameter is the
full universe of IC-bug candidates. H6a confirmed that overriding placement_temp_F=60
on MIX-15 inputs produces a bit-identical Reference result — so whatever is wrong, it
manifests only when placement_temp_F < 60°F.

**V3: Blanket-cell temperature is cosmetic, not load-bearing.**
`k_cell[grid.is_blanket] = 0.0` at `thermal_engine_2d.py:1438` (in full_2d mode)
thermally decouples blanket cells from the concrete. Blanket-cell temperature values
are stored in the T array but never enter the heat-transfer stencil. The pure-R
blanket model is fully captured in `_h_top_combined` at the concrete top surface.
Any grep match whose data flow terminates in `grid.is_blanket` is eliminated as a
candidate before analysis.

---

## H4 Grep Results

**Command:**
```
grep -n "placement_temp_F\|T_initial\|T_0\b\|T0_C\|T0_F" \
    thermal_engine_2d.py compare_to_cw.py cw_scenario_loader.py
```

**Total matches: 20**

---

## Classification Table

| # | File | Line | Symbol | Context | Classification | Candidate? |
|---|------|-----:|--------|---------|----------------|:----------:|
| 1 | thermal_engine_2d.py | 1050 | T_initial_C | Function parameter in `solve_conduction_2d` signature | Signature — not a use | No |
| 2 | thermal_engine_2d.py | 1076 | T_initial_C | Docstring for `solve_conduction_2d` | Documentation | No |
| 3 | thermal_engine_2d.py | 1128 | T_initial_C | `T = T_initial_C.astype(np.float64, copy=True)` in `solve_conduction_2d` | Grid init — expected | No |
| 4 | thermal_engine_2d.py | 1249 | T_initial_C | Function parameter in `solve_hydration_2d` signature | Signature — not a use | No |
| 5 | thermal_engine_2d.py | 1288 | T_initial_C | Docstring for `solve_hydration_2d` | Documentation | No |
| 6 | thermal_engine_2d.py | 1415 | T_initial_C | `Cp_cell[grid.is_concrete] = specific_heat_variable(..., T_initial_C[grid.is_concrete], ...)` | Initial Cp seeding at α=0, pre-loop. See §Below. | Negligible |
| 7 | thermal_engine_2d.py | 1540 | T_initial_C | `T = T_initial_C.astype(np.float64, copy=True)` — working array seeded | Grid init — expected | No |
| 8 | thermal_engine_2d.py | 1849 | T_initial_C | `np.where(grid.is_air[_j_top, :], T_initial_C[_j_top, :], <concrete physics>)` inside time loop | Air-cell pin in top concrete row. Applies to soil-extension strip (x < 0), not concrete columns. See §Below. | No |
| 9 | thermal_engine_2d.py | 1854 | T_initial_C | `T_new[0, :] = T_initial_C[0, :]` — "blanket row j=0 is inactive — pin to initial value" | Blanket row pin. V3-eliminated: k=0. | V3-eliminated |
| 10 | thermal_engine_2d.py | 1873 | T_initial_C | `np.where(_top_row_is_air, T_initial_C[0, :], <blanket physics>)` in non-pure-R else branch | Non-full_2d code path. `_use_pure_r_blanket=True` in full_2d mode, so this branch is never reached. | No (dead path) |
| 11 | thermal_engine_2d.py | 2046 | T_initial_C | `T_new[grid.is_air] = T_initial_C[grid.is_air]` inside full_2d branch | Air-cell pin, all-modes. Expected — air has placeholder rho_cp=1. | No |
| 12 | thermal_engine_2d.py | 2085 | T_initial_C | `T_new[grid.is_air] = T_initial_C[grid.is_air]` outside all branches | Air-cell pin, belt-and-suspenders (comment at 2082). Expected. | No |
| 13 | thermal_engine_2d.py | 2088 | T_initial_C | `T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]` | Blanket-cell pin. V3-eliminated: k=0. | V3-eliminated |
| 14 | cw_scenario_loader.py | 88 | placement_temp_F | `'placement_temp_F': 438` in `CW_DAT_INDEX` dict | DAT file column index — data loading | No |
| 15 | cw_scenario_loader.py | 168 | placement_temp_F | `placement_temp_F: float = 60.0` — dataclass field default | Dataclass default | No |
| 16 | cw_scenario_loader.py | 347 | placement_temp_F | `placement_temp_F=_getf(CW_DAT_INDEX['placement_temp_F'])` — parser | Scenario parsing — expected | No |
| 17 | cw_scenario_loader.py | 780 | placement_temp_F | `f"│  Placement temp = {c.placement_temp_F}°F"` | Display/logging | No |
| 18 | compare_to_cw.py | 110 | placement_temp_F / T0_C | `T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0` | T_initial construction — expected entry point | No |
| 19 | compare_to_cw.py | 111 | T_initial | `T_initial = np.full((grid.ny, grid.nx), T0_C)` | T_initial construction — expected | No |
| 20 | compare_to_cw.py | 303 | placement_temp_F | `"placement_temp_F": scn.construction.placement_temp_F` in metadata dict | Logging/metadata | No |

---

## Candidate Analysis

### #6 — Initial Cp seeding at line 1415

**Code:**
```python
Cp_cell[grid.is_concrete] = specific_heat_variable(
    Wc, Wa, Ww, Ca,
    np.zeros(grid.is_concrete.sum()),   # alpha = 0 at t=0
    T_initial_C[grid.is_concrete],      # uses placement temp
    rho_c,
)
```

**Data-flow:** Seeds concrete Cp before the time loop begins, using the initial concrete
temperature at alpha=0. This is the only step where `T_initial_C` enters a concrete
thermal calculation.

**Is it ever updated?** Yes — every timestep in the time loop at line 1606:
```python
Cp_cell[grid.is_concrete] = specific_heat_variable(
    Wc, Wa, Ww, Ca,
    alpha[grid.is_concrete],
    T[grid.is_concrete],   # uses current T, not T_initial_C
    rho_c,
)
```
The first-step seeding is immediately overwritten on the first iteration.

**Cold-placement amplification effect:** For a 15°F IC difference (45°F vs 60°F),
`specific_heat_variable` will return a slightly different Cp value. The error affects
exactly one timestep (CFL-limited dt ≈ 30–60 s out of 604,800 s). Cp variation with
temperature in the Van Breugel model is mild. The contribution to a sustained 3°F
PeakMax deficit is negligible.

**Verdict: Not a plausible root cause.** One-step initialization artifact. Not a candidate.

---

### #8 — T_initial_C in top-row np.where at line 1849

**Code:**
```python
T_new[_j_top, :] = np.where(
    grid.is_air[_j_top, :],
    T_initial_C[_j_top, :],            # for air cells
    _T_surf_c + dt_step * (...),        # concrete physics
)
```
where `_j_top = grid.iy_concrete_start` (the top concrete row, j=1).

**Data-flow:** `grid.is_air[_j_top, :]` is True for columns `< ix_concrete_start`
(the soil-extension strip to the left of x=0). Those cells are material_id=3 (air),
and the concrete columns (`>= ix_concrete_start`) get the physics branch. So
`T_initial_C` enters this line only for non-concrete air cells.

**Cold-placement amplification:** Air cells are set to T_initial_C once at setup
and re-pinned here each step. These cells have k=0 (line 1424) and zero conductance
to the concrete half. No heat transfer crosses into the concrete from them.

**Verdict: Not a candidate.** Air cells, not concrete. Decoupled from concrete physics.

---

### #9, #13 — Blanket pins at lines 1854, 2088

V3-eliminated. k_cell[grid.is_blanket] = 0.0 (line 1438) in full_2d mode. Blanket
cell temperatures are cosmetic. The H6b counter confirmed 16,464 executions of the
patched pin with zero effect on concrete metrics.

---

### #10 — T_initial_C at line 1873 (else branch)

`_use_pure_r_blanket = (boundary_mode == "full_2d") or skip_blanket_node` (line 1436).
In full_2d mode (all production runs), this is always True, and the else branch at
line 1855 is never reached. Dead path for our investigation.

---

## Audit Result

> **Zero non-trivial, non-blanket, non-air candidates identified.**

Every post-initialization use of `T_initial_C` in the time loop routes to either:
- Air cells (k=0, decoupled from concrete)
- Blanket cells (k=0, V3-eliminated)
- A dead code path (non-full_2d branch)

The concrete physics path never re-reads `T_initial_C` after line 1540. The cold IC
propagates through physics only:

1. **Arrhenius / equivalent age**: cold concrete (45°F = 280 K) produces `af ≈ 0.55`
   vs warm concrete (60°F = 289 K) `af ≈ 0.76` (T_REF_K = 296.15 K, Ea = 26,458 J/mol).
   Cold concrete starts with 72% of the warm hydration rate. This is physically correct,
   but the CW model may calibrate or apply this sensitivity differently.

2. **Equivalent age initialization**: `te = np.full(..., 0.01)` — starts at a fixed small
   value, not derived from placement temperature. No differential initialization bug.

3. **No Cp initialization bug**: first-step Cp from T_initial_C is overwritten at line
   1606 on the very first iteration. One-step artifact, negligible.

---

## Physical Interpretation

H6a (warm IC → perfect Reference match) combined with H4 (no code-level IC misuse)
leads to one conclusion:

> The engine's cold-placement discrepancy is a **modeling gap**, not a code bug.

The engine correctly simulates cold concrete given its physics model. The CW reference
for MIX-15 was generated by software that may include cold-placement-specific effects
the engine does not model. Candidates (H5 territory):

- **Different Arrhenius calibration or reference temperature in CW**: if CW calibrates
  the equivalent-age model with a reference temperature or site that implicitly assumes
  ambient conditions near 60–70°F, the fitted τ/β/α_u/Ea may not transpose correctly
  to 45°F concrete, causing the engine to underpredict peak temperature.
- **Condensation latent heat**: at |T₀ − T_amb| ≈ 30°F (MIX-15) vs ≈ 15°F (Reference),
  condensation of humid air on cold concrete releases latent heat not modeled in the
  engine. CW may include this effect, especially when the surface is below dew point.
- **Cold-mix accelerated early hydration from larger thermal gradient**: a 30°F
  differential between concrete and formwork/soil drives a different early heat flux
  pattern, which may interact with CW's calibrated model differently.

None of these is accessible via a single-line code fix in the current engine architecture.
All require Sprint 5/6 scope (new physics term or calibration adjustment).

---

## Blanket-Pin Cosmetic Status (V3 Documentation)

For Sprint 6 session hygiene: the blanket-cell pin at `thermal_engine_2d.py:2088`
(`T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]`) appears to be an IC-propagation
bug at code review but is physically inert. In full_2d mode, `k_cell[grid.is_blanket] = 0.0`
(line 1438) decouples blanket cells entirely — no heat flux crosses the blanket-concrete
interface via conduction, and the BC at the concrete top surface uses `_h_top_combined`
which already incorporates the blanket thermal resistance. Blanket-node temperatures are
stored in the output array but not used in any flux or heat-source calculation.

The air-cell pin at line 2085 is structurally similar but IS load-bearing: air cells have
placeholder `rho_cp=1` and resetting them prevents numerical drift that would corrupt
adjacent flux computations. Do not conflate the two pins when addressing this in Sprint 6.

---

## Decision Routing

- **H4 candidates for code-bug fix: 0**
- Combined with Phase 1 (MARGINAL boundary-physics, primarily hydration-driven PeakMax)
  and Phase 2 (H6a cold-IC causally confirmed, H6b blanket not the mechanism):
- **Proposed decision: G/H routing to Sprint 5** — the cold-IC failure is physics-driven,
  not code-bug-driven. The root cause is in the equivalent-age / hydration-rate model's
  behavior at cold placement temperatures, a domain requiring calibration adjustment or
  new cold-placement physics terms.
- PR 16 closes diagnostic-only. Commit: ablation md, H6 md, H4 audit md, §7.6.4 update.
