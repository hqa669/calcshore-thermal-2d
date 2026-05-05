# Sprint 9 Stage 1-pilot — Refined Sensitivity Sweep v2 Report

**Scope:** composition-centered Hu_factor sweep ±1% around each mix's
composition-correct value, × 5 c_multiplier values (±5% around 1.00).
Sensitivity characterization only — no calibration recommendation.

---

## §1.1 Sweep configuration

| Parameter | Value |
|---|---|
| MIX-01 Hu_factor grid | [0.9414, 0.9464, 0.9514, 0.9564, 0.9614] |
| MIX-07 Hu_factor grid | [0.8846, 0.8896, 0.8946, 0.8996, 0.9046] |
| c_multiplier grid     | [0.95, 0.97, 1.0, 1.03, 1.05] |
| c_eff range           | [1.0005, 1.1059] |
| Total runs            | 50 (5 × 5 × 2 mixes) |
| Wall-clock runtime    | 57.5 s (1.15 s/run) |

**Wrapper configuration (B1 inverse-compensation, unchanged from v1):**
```
c_eff  = c_multiplier × 1.0532
Hu_eff = (Hu_raw × Hu_factor_override) / c_eff
α_u_eff = c_eff × α_u_raw
Engine settings: model_soil=False, is_submerged=True, blanket=0.0, k_uc×0.96
τ, β, E_a pass through unmodified.
```

**CW data source:** `validation/sprint9/stage1_pilot/cw_data/{mix01,mix07}/` only.
**No engine source modification.** No prior-sprint CW datasets accessed.

---

## §1.2 MIX-01 surface readout at composition-correct cell

**Reference point:** Hu_factor = 0.9514, c_multiplier = 1.00, c_eff = 1.0532

| Metric | Value |
|---|---|
| max\|R\| (°F)      | 1.0507 |
| di_at_max          | 47 |
| R_di36 (°F)        | +0.1662 |
| R_di42 (°F)        | +0.1516 |
| R_di47 (°F)        | -1.0507 |
| R_di48 (°F)        | +0.2055 |
| stencil_drop (°F)  | -1.2023 |
| T_core engine (°F) | 149.4002 |
| T_core CW (°F)     | 149.2340 |
| Gate (< 0.5°F)     | FAIL |

---

## §1.3 MIX-07 surface readout at composition-correct cell

**Reference point:** Hu_factor = 0.8946, c_multiplier = 1.00, c_eff = 1.0532

| Metric | Value |
|---|---|
| max\|R\| (°F)      | 1.5843 |
| di_at_max          | 47 |
| R_di36 (°F)        | +0.7598 |
| R_di42 (°F)        | +0.6409 |
| R_di47 (°F)        | -1.5843 |
| R_di48 (°F)        | +0.2068 |
| stencil_drop (°F)  | -2.2252 |
| T_core engine (°F) | 148.7157 |
| T_core CW (°F)     | 147.9560 |
| Gate (< 0.5°F)     | FAIL |

---

## §1.4 c_multiplier sensitivity at composition-correct Hu_factor

Linear-fit slope across all 5 c_multiplier points at each mix's correct Hu_factor.
Unit: °F per 1% change in c_multiplier (per 0.01 change).

| Metric | ∂/∂c_mult per 1% (MIX-01) | ∂/∂c_mult per 1% (MIX-07) |
|---|---|---|
| max_R           | -0.0189°F | -0.0130°F |
| R_di36          | +0.0025°F | +0.0004°F |
| R_di47          | +0.0189°F | +0.0130°F |
| stencil_drop    | +0.0132°F | +0.0097°F |

**Diagnostic (stencil_drop sensitivity to c_mult):**
  |∂(stencil_drop)/∂c_mult| = 0.0132°F/1% (MIX-01), 0.0097°F/1% (MIX-07)
  Threshold: < 0.05°F/1%
  → **Both below threshold.** c_multiplier does not substantially affect the
    stencil dip relative to bulk — consistent with structural decoupling from c().

---

## §1.5 Hu_factor sensitivity within ±1% at c_multiplier = 1.00

Linear-fit slope across all 5 Hu_factor points at c_mult = 1.00.
Unit: °F per 0.005 change in Hu_factor (one grid step).

| Metric | ∂/∂Hu per 0.005 (MIX-01) | ∂/∂Hu per 0.005 (MIX-07) |
|---|---|---|
| max_R           | -0.1027°F | -0.0561°F |
| R_di36          | +0.4308°F | +0.5431°F |
| R_di47          | +0.1482°F | +0.1826°F |
| stencil_drop    | -0.2784°F | -0.3530°F |

**Diagnostic (stencil_drop sensitivity to Hu_factor within ±1%):**
  |∂(stencil_drop)/∂Hu_factor| = 0.2784°F/step (MIX-01), 0.3530°F/step (MIX-07)
  Threshold: < 0.05°F/step
  → **One or both mixes above threshold.** stencil_drop has non-trivial Hu dependence within ±1%.

---

## §1.6 Stencil_drop constancy across full restricted grid

mean and std of stencil_drop across all 25 grid points per mix.

| Mix | mean (°F) | std (°F) | range [min, max] | std < 0.2°F? |
|---|---|---|---|---|
| MIX01 | -1.2021 | 0.3968 | [-1.8260, -0.5798] | NO |
| MIX07 | -2.2259 | 0.5005 | [-2.9817, -1.4727] | NO |

**Diagnostic:** std ≥ 0.2°F for at least one mix — stencil_drop shows
non-trivial dependence on grid position within the restricted region.

---

## §1.7 R_di48 constancy check

R_di48 = T_engine − T_CW at the bottom Dirichlet face (di=48, t=168 hr, wi=0).
Both solvers apply T_soil=85°F as a direct Dirichlet at di=48 under model_soil=False.

| Statistic | Value |
|---|---|
| mean R_di48 across 50 runs | +0.2061°F |
| range (max − min)          | 0.00817°F |
| < 0.05°F variation?        | YES |

**Diagnostic:** R_di48 range < 0.05°F across all 50 grid points. BC-face residual
is independent of bulk physics — consistent with strong Dirichlet write at di=48.

---

## §1.8 Synthesis

Sensitivity characterization at physically meaningful operating points complete.
All findings below are factual readouts from the surface data.

### Residual landscape at composition-correct (Hu_factor, c_mult=1.00)

**MIX-01:** max|R| = 1.0507°F at di=47. R_di47 = -1.0507°F, R_di48 = +0.2055°F, R_di36 = +0.1662°F, stencil_drop = -1.2023°F. Gate (< 0.5°F): FAIL.

**MIX-07:** max|R| = 1.5843°F at di=47. R_di47 = -1.5843°F, R_di48 = +0.2068°F, R_di36 = +0.7598°F, stencil_drop = -2.2252°F. Gate (< 0.5°F): FAIL.

### Gate achievability within restricted (Hu_factor ±1%, c_mult ±5%) grid

**MIX-01:** No grid point within the restricted region achieves max|R| < 0.5°F. Minimum max|R| in region = 0.8070°F at (Hu_factor=0.9564, c_mult=1.05).

**MIX-07:** No grid point within the restricted region achieves max|R| < 0.5°F. Minimum max|R| in region = 1.3357°F at (Hu_factor=0.8996, c_mult=1.05).

### Structural-vs-calibration interpretation

**stencil_drop decoupling:** std(stencil_drop) = 0.3968°F (MIX-01), 0.5005°F (MIX-07) across the full composition-meaningful 5×5 region. At least one mix exceeds the 0.2°F threshold: the relative stencil dip shows non-trivial variation over this region.

**R_di48 decoupling:** range = 0.00817°F across all 50 grid points. BC-face residual is constant — strong Dirichlet write at di=48 is insensitive to interior bulk changes, as expected.

Frame: the restricted residual landscape implies that max|R| at the composition-correct operating point is dominated by the stencil mechanism concentrated at di=47, which (if std(stencil_drop) < 0.2°F) is structurally invariant across the calibration knobs available at the wrapper level. The residual cannot be resolved by adjusting Hu_factor or c_multiplier alone within their physically meaningful ranges.

