# STAGE5b — Smoke Test: model_soil=True Path at New Default Resolution

## Configuration

- Run: Run B (placement=73°F, soil=60°F, |ΔT|=13°F)
- Flags: `model_soil=True`, `is_submerged=True`

## Check 1 — Native 1× bit-identical to Stage 4b reference

Reference: `stage3_fix2_runs/runB_t168.csv` (Stage 3.5 fix-2 / Stage 4b baseline)  
Engine: `grid_refinement=1`, `model_soil=True`  (same configuration as Stage 4b)

| Metric | Value |
|---|---|
| Engine field shape | (13, 21) |
| Reference field shape | (13, 21) |
| max\|ΔT\| vs reference | **4.443036°F** |
| Tolerance | 0.001°F |
| Wall-clock (1× model_soil=True) | 1.7s |
| Verdict | **FAIL ✗** |

**FAIL**: the 1× engine output differs from the Stage 4b reference by more than 0.001°F.  Investigate: the refactor may have altered the soil-mesh path.

## Check 2 — 6× Default Physically Reasonable vs 1×

Comparison: 6× engine resampled bilinearly onto the 1× node positions.

| Metric | Value |
|---|---|
| 6× field shape | (73, 121) |
| max\|ΔT\| (6× vs 1×, on 1× nodes) | **2.7754°F** |
| mean\|ΔT\| | 0.2896°F |
| Tolerance (sub-degree) | 1.5°F |
| Wall-clock (6× model_soil=True) | 9.4s |
| Wall-clock ratio (6× / 1×) | 5.5× |
| Verdict | **WARN — check residual pattern** |

Max diff 2.7754°F exceeds tolerance.  Inspect residual pattern: if structured (corner concentration, depth gradient) it may indicate a soil-mesh issue at the finer grid.  If diffuse, likely just numerical convergence.

## Overall Verdict: **WARN / FAIL**
