# Sprint 8 Stage 2-prep v2 — Direct α(t) Comparison via Constant-k Inversion

**Path C: GATE FAIL at §3.1 — inversion degenerate (mid-depth unreachable in 80-ft mat).**
**Path D: COMPLETE — 15-row direct α(t) comparison table delivered.**

## Key results

### Path C finding (stage2_prep_alpha_inversion.py)
Engine T at di=24 (12.2 m depth) = 73.000°F for ALL α_test values — insensitive to k(α). CW shows a 2.16°F uniform bulk offset across di=8–38, produced by a mechanism not in the CalcShore FD engine. Path C is structurally degenerate for this geometry.

### Path D result (stage2_prep_path_d_bulk_inversion.py)
CW's di=24 bulk shift IS the k(α) signal: ΔT_bulk = C(t)·α, linear with CV < 5% across 3 datasets at all 5 timesteps. Inverting directly:

| Metric | Δ at T_ref=23°C | Δ at T_ref=21°C |
|---|---|---|
| max|Δ| | 8.1% | 2.5% |
| mean|Δ| | 2.0% | 1.1% |
| RMS|Δ| | 2.8% | 1.2% |

At t=168hr (Sprint 8 calibration horizon): Δ@21°C ≤ 1.4% — validates Stage 2-prep Path B bound.

## Conclusion

The T_ref=21°C hypothesis (Path B best-fit) is **empirically confirmed** by direct α inversion. Sprint 8 proceeds with ±5% α trajectory uncertainty; the dominant source (T_ref=23°C engine default) contributes 2.0% RMS at t=168hr, reducible to 1.2% RMS with T_ref correction.

## Files

```
stage2_prep_alpha_inversion.py          Path C attempt (degenerate — documented)
stage2_prep_path_d_bulk_inversion.py    Path D — bulk-shift inversion (delivers table)
findings/
  wrapper_validation.md                 §3.1 Path C gate fail
  coarse_sweep.csv                      §3.2 Path C sweep (flat objective)
  bulk_shift_linearity.md               §D.1 linearity gate (PASS)
  t_ref_per_timestep.csv                §D.2 per-timestep best-fit T_ref
  t_ref_validation.md / .csv            §D.4 15-row headline table
  synthesis.md                          full synthesis (Path C + D)
```
