# Sprint 8 Passdown — k(α) Curve Calibration

**Status as of 2026-05-01:** Sprint 8 in progress. Stage 2 calibration BLOCKED by ambient/weather mismatch finding. Decision pending on path forward.

---

## Sprint 8 goal

Determine whether the Sprint 7 calibration `k_uc × 0.96` (calibrated at α ≈ 0.036) holds across the production-relevant α range (α = 0.04 to α ≈ 0.7), or whether the optimal multiplier varies with α.

Sprint 7 only validated at the very tip of the production hydration curve. Real production runs span α ∈ [0, 0.7] with k changing by ~17% across that range. Sprint 8 covers the rest.

---

## Sprint 7 baseline (closed)

Sprint 7 closed cleanly with the following calibration committed to main:

```
K_UC_CALIBRATION_FACTOR_SPRINT7 = 0.96
```

Applied in `thermal_engine_2d.py` at the k_uc derivation point. Van Breugel formula `k(α) = k_uc × (1.33 − 0.33α)` shape unchanged.

All 9 Sprint 7 validation runs pass Structure C metric (R1, R2 ≤ 0.35°F at α ≈ 0.036). Detailed history archived on `sprint-7-detailed` branch. Single squashed commit `f6e83b3` on main.

Key non-calibration findings preserved in engine:
- `model_soil` flag (default False matches CW)
- `is_submerged` toggle for side-soil cell allocation
- `blanket_thickness_m` parameter with safe handling at zero thickness (CFL guard + dy_minus guard)
- 6× default grid refinement for concrete domain
- CW spurious heat floor identified at Hu=1 (~0.13°F over 168 hr)
- Bottom-side corner BC stencil asymmetry documented and deferred

Full Sprint 7 details: `docs/calibration/k_uc_calibration_sprint7.md` and `STAGE5e_REPORT.md`.

---

## Sprint 8 progress

### Stage 0 — CW formula verification ✅
Read CW V3 manual Section 3.1.4 (page 10). CW Equation 23: `k_c(α) = k_uc × (1.33 − 0.33α)` — same as engine. No formula reverse-engineering needed.

### Stage 1 — CW α-dependence empirical confirmation ✅
Generated 3 CW datasets at varying hydration parameters. Pairwise temperature comparison at t=168 hr showed CW responds to α with magnitude matching Van Breugel predictions (Pair A:B ratio = 2.10 ≈ Δα ratio 2.06). Signal magnitude 22× the constant-k threshold. CW empirically uses Van Breugel.

Datasets at `validation/sprint8/stage1_cw_k_alpha_empirical/`:
- `lowhyd` (α_u=0.10, τ=200, β=0.10) → α reaches 0.036 by t=168
- `highhyd` (α_u=0.70, τ=10, β=0.85) → α reaches 0.638 by t=168
- `highhyd_b010` (α_u=0.70, τ=10, β=0.10) → α reaches 0.329 by t=168

### Stage 2-prep — Engine vs CW α(t) trajectory check ✅ (with caveat)
Engine and CW disagree on α(t) by 5–8% at intermediate times, narrowing to ~1% at t=168. Likely cause: Arrhenius reference temperature mismatch. Empirical curve-fit on Stage 1 data suggests engine uses T_ref=23°C while CW best-fits to T_ref=21°C. **This is a curve-fitting result, not verified from primary documentation.** Engine T_ref hasn't been read directly from `thermal_engine_2d.py` source; CW T_ref not found in V3 manual sections inspected.

**Decision: Don't change T_ref.** Accept residual α(t) discrepancy as part of what the calibration factor compensates for. T_ref alignment deferred.

Artifacts at `validation/sprint8/stage2_prep_alpha_trajectory/` and `validation/sprint8/stage2_prep_alpha_trajectory_v2/`.

### Stage 2 — Multi-α calibration ⚠️ BLOCKED

Generated 12 new CW datasets at engineered α plateau values:
- α_u ∈ {0.20, 0.40, 0.60, 0.80}
- 3 scenarios per α target: A (T_pl=73, T_soil=73), F (T_pl=73, T_soil=45), I (T_pl=100, T_soil=73)
- Same Sprint 7 geometry, mix, BC, Hu=1 suppressed, τ=5 hr, β=0.85

Datasets at `validation/sprint8/stage2_calibration/cw_data/` (uncommitted).

Ran k_uc multiplier sweep (0.92 to 1.04) at each α target.

**Result: 9/21 pass (only Sprint 7 datasets pass).** All 12 new datasets fail by 0.30–3.24°F over the 0.35°F gate.

#### Critical finding — k_uc has near-zero leverage

| α_u | Worst R1/R2 at k×0.96 | Variation across full sweep |
|---|---|---|
| 0.20 | 0.89°F | 0.37°F |
| 0.40 | 1.79°F | 0.18°F |
| 0.60 | 2.69°F | 0.10°F |
| 0.80 | 3.59°F | 0.06°F |

Residuals scale linearly with α_u (~4.4°F per unit α). Optimum k_uc factor pegs at left sweep edge (0.92) for all 4 α targets — sweep-edge artifact, not physical. **No factor in any reasonable range can close the gate.**

#### Root cause hypothesis (under investigation, NOT YET VERIFIED)

The Stage 2 agent proposed: CW used real Austin TX July 2026 weather (diurnal cycling, solar, ambient ~80°F average). Engine uses `make_neutral_env(73°F)` — flat constant ambient with no solar/diurnal.

Proposed mechanism:
- At α ≈ 0.036, high concrete conductivity (k ≈ 1.32 × k_uc) equilibrates to soil BC fast, washing out weather signal — Sprint 7 worked
- At higher α, lower conductivity (k ≈ 1.08 × k_uc at α=0.80), concrete more insulating, weather signature persists
- CW side-face cells read 0.65–2.72°F above engine's Dirichlet-pinned 73°F
- No k_uc factor can close that gap

**This hypothesis has not been verified.** It's the agent's interpretation. Needs verification via residual depth profile inspection — if residual concentrates near top BC (di < 10), ambient is confirmed. If uniformly distributed across depth, additional mechanism at play.

---

## Three options forward (decision pending)

1. **Document narrower validity range** — Sprint 7's k_uc = 0.96 is valid at α ≤ ~0.05 in the high-conductivity regime. Move to other calibration axes (ambient/top-BC, hydration heat) before revisiting high-α k_uc. Sprint 8 closes with documented limitation.

2. **Regenerate the 12 CW datasets with neutral weather** — configure CW with constant T_air = T_pl, no solar, no diurnal. Should isolate k_uc cleanly. Windows-side reruns, same hydration/scenario parameters, only weather setting changes.

3. **Compare at earlier times only** — use t=24 hr instead of t=168 hr. Loses production-relevance horizon. Not recommended.

**Pre-recommendation:** Option 2, but verify the ambient hypothesis first via residual depth profile.

---

## Methodology patterns from Sprint 7 (continue applying)

- **Diagnostic-first discipline** — characterize before fixing, no proposed fixes mid-diagnostic
- **Per-axis scope discipline** — one calibration axis at a time, never co-tune
- **Wrapper-only experiments** — parameter perturbation tests stay uncommitted; engine commits only at calibration close
- **Region-based metric** (R1 side, R2 bottom, R3 corner) — corner reported but not gated due to engine asymmetry
- **Squash-merge to main from sprint branches** — detailed history archived on `sprint-N-detailed`, main stays clean
- **Bidirectional sensitivity sweep before any calibration commit** — confirm calibrated value isn't a sweep-edge artifact

---

## Outstanding uncertainty (worth being explicit about)

The "ambient/weather mismatch" hypothesis is the agent's interpretation of why k_uc has no leverage at high α. It hasn't been:
- Verified by residual depth profile inspection
- Verified by checking what weather setting Sprint 7's CW input.dat used
- Verified by inspecting CW's side-BC implementation

It's plausible but other mechanisms could contribute (e.g., side-BC implementation differences, T_ref propagating through hydration/heat-release coupling at higher α, numerical artifacts at the physics regime where k drops).

The new chat session should treat the ambient hypothesis as a candidate to verify, not a conclusion to act on.

---

## What the new session needs to do

The blocker is: at high α, k_uc cannot be cleanly calibrated against the existing 12 CW datasets. We need to figure out why before deciding the path forward.

**Immediate diagnostic tasks (no new CW runs needed):**

1. **Verify the ambient hypothesis via residual depth profile.** For one or two of the failing datasets (e.g., alpha08_F or alpha08_I), plot residual as a function of depth (di) at t=168. If residual is concentrated near di=0 (top), ambient is confirmed. If relatively uniform, something else is going on.

2. **Check Sprint 7 CW input.dat weather setting.** If Sprint 7's 9 CW datasets used real Austin weather and Sprint 7 still passed, the agent's "high conductivity washes out weather" explanation is supported. If Sprint 7 used neutral weather, that's a mismatch with Sprint 8 dataset configuration.

3. **Inspect CW side-BC handling.** The Stage 2 agent observed "CW side-face cells are 0.65–2.72°F above 73°F due to Austin weather." That implies CW's side-face handling isn't pure Dirichlet — there's some thermal coupling to ambient. Worth understanding the actual CW BC implementation vs the engine's pure Dirichlet.

**Once verification is done:** decide between Options 1, 2, or 3 above based on what the diagnostics show.

---

## Repo state

- Main: clean, last commit `f6e83b3` (Sprint 7 closure)
- Sprint 8 work: all uncommitted in `validation/sprint8/`
- Engine source: untouched since Sprint 7 closure (no Sprint 8 engine changes)

---

## Deferred work (parked for future sprints)

1. k(α) curve calibration at higher α (Sprint 8 — currently blocked)
2. Bottom-side corner BC stencil correction (engine asymmetry from Sprint 7)
3. ρ_concrete and Cp_concrete cross-validation
4. Hydration kinetics calibration (τ, β, α_u, E_a, T_ref)
5. Ambient/top-BC physics calibration (now identified as likely required precursor to Sprint 8 closure)
6. Formwork, blanket physics, curing compound axes
7. T_ref alignment (deferred, accepted as part of Sprint 8 calibration target)
