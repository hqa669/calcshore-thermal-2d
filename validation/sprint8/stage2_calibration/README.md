# Sprint 8 Stage 2 — Multi-α k_uc Calibration

**Sprint:** 8 (k(α) curve calibration)
**Stage:** 2 — calibrate k_uc multiplier across 4 α targets (α_u ∈ {0.20, 0.40, 0.60, 0.80})

## Purpose

Determine whether the Sprint 7 k_uc × 0.96 calibration holds across the
production-relevant α range, or whether the optimal multiplier varies with α.

- **Outcome A:** constant factor (possibly still 0.96) — Sprint 8 closes.
- **Outcome B:** α-dependent function f(α) — engine needs updated k_uc(α).

This directory contains wrapper scripts only. No engine source modifications.
No T_ref changes. No commits to main during this stage.

## Datasets

- **12 new CW datasets** in `cw_data/` (α_u ∈ {0.20, 0.40, 0.60, 0.80} × 3 BC scenarios)
- **9 Sprint 7 datasets** at `../../soil_calibration/cw_runs/` (α≈0.036 anchor)

## Execution order

```
# §3.1 + §3.2 — copy and validate datasets (must run first)
python stage2_input_consistency.py

# §3.3 — α(t) reference trajectories (fast, no engine runs)
python stage2_alpha_reference.py

# §3.4 — baseline residuals at Sprint 7 factor 0.96 (~12 min)
python stage2_baseline_residuals.py

# §3.5 — k_uc sweep: 4 α × 7 factors × 3 scenarios (~85 min)
python stage2_kuc_sweep.py

# §3.5 refined (conditional — only if sweep hints say "run me")
python stage2_kuc_sweep_refined.py

# §3.6 + §3.7 — calibration decision + curve fit if Outcome B
python stage2_calibration_decision.py

# §3.8 + §3.9 — final 21-point validation + Sprint 7 regression (~21 min)
python stage2_final_validation.py [--factor FLOAT]

# Figures 10, 11, 12
python stage2_figures.py
```

## Override mechanism

`engine_runner.py` monkey-patches `thermal_engine_2d.K_UC_CALIBRATION_FACTOR_SPRINT7`
to the specified absolute factor before each solve and restores the original value
afterward. factor=0.96 exactly replicates the Sprint 7 baked-in calibration.

## Region metric (Structure C, Sprint 7)

| Region | Definition | Gate |
|--------|-----------|------|
| R1 | max\|R\| at di=24 across all widths | ≤ 0.35°F |
| R2 | max\|R\| at wi=0, di=24..48 | ≤ 0.35°F |
| R3 | RMSE at di=43..48, wi=8..12 | reported only |

## Outputs (findings/)

| File | Contents |
|------|----------|
| input_consistency.md | §3.2 dataset validation |
| engine_alpha_reference.csv | §3.3 α(t) trajectories |
| baseline_residuals.csv | §3.4 residuals at factor=0.96 |
| kuc_sweep.csv | §3.5 all sweep results |
| kuc_sweep_optimal.txt | §3.5 optimal factor per α |
| kuc_sweep_refined.csv | §3.5 refined sweep (if run) |
| calibration_decision.md | §3.6 spread + §3.7 curve fit |
| final_validation.csv | §3.8 Sprint 8 21-pt gate |
| sprint7_regression.csv | §3.9 Sprint 7 regression |
| synthesis.md | hand-written synthesis |

## Figures (saved to /mnt/user-data/outputs/ or findings/)

| File | Contents |
|------|----------|
| stage2_calibration_factor_vs_alpha.png | Figure 10 — factor vs α scatter |
| stage2_residuals_per_alpha.png | Figure 11 — R1/R2 before/after calibration |
| stage2_R1_R2_field_at_optimal.png | Figure 12 — residual field heatmaps |
