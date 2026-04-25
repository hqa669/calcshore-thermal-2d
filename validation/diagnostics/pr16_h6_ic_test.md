# PR 16 — H6 IC Propagation Test (MIX-15 B2)

Phase 2 diagnostic for Sprint 4 PR 16. Two-tier test separating the joint cold-IC
effect (H6a) from the blanket-pin-specifically (H6b).

## Test Design

| Variant | IC | Blanket pin | Reference |
| --- | --- | --- | --- |
| Baseline | MIX-15 cold (45°F) | T_initial (45°F) | MIX-15 CW |
| **H6a** | **Warm (60°F)** | T_initial (60°F) | **MIX-01 CW** |
| **H6b** | MIX-15 cold (45°F) | **T_amb (current)** | MIX-15 CW |

H6a tests joint cold-IC effect — both the concrete grid IC and the blanket pin swap
to 60°F-equivalent. Compare to MIX-01 CW (warm-placed reference, same construction).

H6b isolates the blanket-pin contribution — concrete IC stays cold (45°F), only the
blanket-cell temperature policy changes. Compare to MIX-15 CW (cold-placed; concrete
IC is unchanged so the cold-placed CW is the correct reference).

## Results

### Side-by-Side Metric Table

| Metric | Baseline (MIX-15) | H6a (warm IC) | H6b (blanket→T_amb) |
| --- | ---: | ---: | ---: |
| PeakMax Δ (°F) | -3.01 | -0.29 | -3.01 |
| PeakGrad Δ (°F) | -8.25 | -0.48 | -8.25 |
| FieldRMS (°F) | 4.06 | 0.88 | 4.06 |
| CenterRMS (°F) | 2.07 | 0.74 | 2.07 |
| CornerRMS (°F) | 2.05 | 2.26 | 2.05 |
| S0 pass | 1/5 | 5/5 | 1/5 |

Reference for H6a: MIX-01 CW (warm-placed). Reference for H6b/Baseline: MIX-15 CW (cold-placed).

### FieldRMS Movement

Baseline FieldRMS: 4.06°F
H6a FieldRMS:      0.88°F  (RESOLVED ≤1.5°F)
H6b FieldRMS:      4.06°F  (NOT MOVED)

**H6 verdict: CONCRETE-IC-TRIGGER**

H6a resolves FieldRMS but H6b does not. The concrete grid IC (equivalent-age delay), not the blanket pin, is the primary trigger.

**Proposed Phase 3 decision: F (via H4 audit)** — blanket pin is not the bug;
concrete IC propagation is the trigger. Proceed to H4 grep audit to locate the
specific concrete-IC-dependent term.

## Blanket-Pin Mechanism (§7.6.4 Documentation)

This observation is recorded regardless of H6 outcome as a candidate Sprint 6 cleanup item.

### What the pin does

In `thermal_engine_2d.py`, after every timestep in all boundary modes, the engine
applies two unconditional overrides (`thermal_engine_2d.py:2082-2088`):

```python
# Pin air cells to initial value for ALL modes (air has rho_cp=1 placeholder;
# leaving them at the BC-updated value would corrupt next-step flux calcs).
if grid.is_air.any():
    T_new[grid.is_air] = T_initial_C[grid.is_air]
# When blanket is pure-R (full_2d or skip_blanket_node), pin it to initial value.
if _use_pure_r_blanket:
    T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]
```

The **air-cell pin** is necessary: air cells have a placeholder `rho_cp=1` and are
not physically evolved. Resetting them each step prevents numerical drift from
corrupting flux calculations at adjacent concrete cells.

The **blanket-cell pin** (line 2088) applies the same pattern to the insulating
blanket layer. In a pure-R blanket model, the blanket is a thermal resistance with
no thermal mass, so it has no physical temperature to track. The pin holds blanket
cells at `T_initial_C` — the concrete placement temperature — for the entire 168-hour
simulation.

### Why this is NOT problematic — V3 verification finding

This section originally hypothesized that the cold blanket pin would create a
systematic heat-sink at the top surface for MIX-15. **H6b and V3 falsified this.**

The reason: `k_cell[grid.is_blanket] = 0.0` at `thermal_engine_2d.py:1438` in
full_2d mode sets blanket-cell conductivity to zero. With k=0, no heat flux crosses
the blanket-concrete interface via conduction. The pure-R blanket model is entirely
captured in `_h_top_combined` (the series combination of forced-convection h and the
blanket R-value) applied at the top concrete surface. Blanket-node temperatures are
stored in the output T array but never enter the heat-transfer stencil.

**V3 counter result:** the H6b patch executed 16,464 times (one per timestep).
Blanket T at run end = 23.8°C (tracking ambient) vs. baseline 7.2°C (frozen at IC).
The patched blanket temperature changed completely — yet concrete metrics were
bit-identical to the unpinned baseline. This confirms zero physical coupling.

**The blanket-cell pin at line 2088 is cosmetic.** It stores the IC value in an
inert node for output-array consistency; it does not affect the physics.

The air-cell pin at line 2085 is structurally similar but IS load-bearing (air
cells have placeholder `rho_cp=1`; resetting them prevents stencil corruption at
adjacent cells). Do not conflate the two in a Sprint 6 cleanup.

### Sprint 6 scope note (revised)

The blanket pin is cosmetically misleading — it looks like a cold-sink bug and
would draw attention in a future code review. A Sprint 6 cleanup could either:
(a) document the k=0 decoupling clearly at line 2088, or
(b) replace the pin with a neutral value (e.g., `_T_amb_C`) to prevent confusion,
knowing it has no physical effect.

This is a documentation/clarity issue only, not a correctness fix.

### Sprint 6 scope note

Even if this bug is fixed in PR 16, the air-cell pin at line 2085 applies the same
T_initial_C pattern to air cells. Air cells are non-physical (placeholder rho_cp),
so the pin is correct for them — they should not evolve thermally. The blanket is
different: it represents a physical material with a real thermal interaction with the
ambient environment. A Sprint 6 cleanup could replace the blanket pin with a
quasi-steady boundary equation that accounts for blanket thermal resistance without
freezing blanket temperature at placement IC.

