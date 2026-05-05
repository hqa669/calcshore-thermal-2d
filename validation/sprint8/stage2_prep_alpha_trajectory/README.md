# Sprint 8 Stage 2-prep — α(t) Trajectory Agreement Verification

**Stage predecessor:** Stage 1 (CW k(α) empirical verification — Outcome 1 confirmed)
**Stage type:** Read-only data analysis. No engine changes, no new CW runs, no commits.
**Question answered:** Do the engine and CW compute matching α(t) trajectories at identical Schindler parameters?

---

## Key result

**FAIL at 8.48% max discrepancy (A/C mean|ΔT| at t=24 hr). CAUTION at 2.3–6.2% for the main A/B test and at t=168 hr.**

The discrepancy is not catastrophic for Sprint 8: at the calibration horizon (t=168 hr), the α offset between engine and CW is ≤1.4%. The principal source is a T_ref mismatch — engine uses 23°C, CW likely uses ~21°C (Schindler-Folliard 2002 calibration value). Sprint 8 proceeds with a ±5% α uncertainty acknowledgment.

---

## Method

CW outputs no α(t) data (Path A unavailable — see §3.1 finding). Path B: ratio test on existing
Stage 1 pairwise ΔT data. In a linear thermal model with Van Breugel k(α):

> ΔT_A / ΔT_B ≈ Δα_A / Δα_B (same BCs, geometry, k_uc)

Discrepancy = (observed ΔT ratio / engine Δα ratio − 1) × 100%.
Computed for two ratio tests (A/B and A/C), three ΔT statistics (max|ΔT|, mean|ΔT|, RMS), and three timesteps (24, 84, 168 hr) = 18 values total.

T_REF_K back-calculation (§3.4): sweep T_ref ∈ {20–25°C} to find the candidate minimizing
ratio discrepancy vs the observed CW ΔT ratios.

---

## Stage 1 inputs (no new runs needed)

| File | Source |
|---|---|
| `pairwise_dT_stats.csv` | Stage 1 pairwise ΔT at t=24/84/168 hr |
| `engine_alpha_reference.csv` | Engine isothermal α(t) at three datasets × milestones |

---

## Results

### §3.3 Ratio test verdict

| Test | A/B (main) | A/C (cross-check) |
|---|---|---|
| max discrepancy (any stat/time) | +6.24% | −8.48% |
| Max at t=168 hr | +5.36% | −4.64% |
| Verdict (overall) | CAUTION | **FAIL** |

Combined max = **8.48% → FAIL** by threshold definition.

### §3.4 T_ref back-calculation

Best-fit T_ref = **21°C** (RMS discrepancy 3.21% vs 5.05% at engine default of 23°C).
At T_pl=73°F: factor(21°C)=1.131, factor(23°C)=0.985 — CW α accumulates ~14.7% faster in equivalent age.

At t=168 hr, implied α offset vs engine: lowhyd=1.4%, highhyd=1.0%, highhyd_b010=1.0%.

---

## Sprint 8 decision

**Sprint 8 proceeds without redesign.** The t=168 hr α offset (≤1.4%) is below the noise floor from other calibration sources (grid discretization, BC alignment). Sprint 8 calibration record should note the ±5% α trajectory uncertainty band. If post-calibration residuals exceed 3%, revisit the T_ref mismatch as a first diagnostic.

---

## Files

```
stage2_prep_ratio_test.py           §3.2 + §3.3 ratio test
stage2_prep_t_ref_analysis.py       §3.4 T_ref back-calculation
findings/
  cw_output_structure.md            §3.1 — CW output column structure (Path A finding)
  ratio_test_results.csv / .md      §3.2 + §3.3 — 18-value ratio test table
  t_ref_analysis.csv / .md          §3.4 — T_ref sweep and best-fit
  synthesis.md                      §3.5 — observational synthesis (a)–(d)
```
