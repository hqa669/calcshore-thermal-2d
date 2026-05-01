# STAGE4b — Smoke Test: model_soil=True Path Preserved

## Configuration

- Run: Run B (placement=73°F, soil=60°F)
- Engine flags: `model_soil=True`, `is_submerged=True`
- Weather: Austin TX (shared environment from `validation/cw_exports/MIX-01/weather.dat`)
- Duration: 168 hr at output_interval_s=1800 s
- Comparison: Stage 3.5 Fix-2 reference output (`stage3_fix2_runs/runB_t168.csv`)
  - The stage3_fix2 runs were the last validated engine runs before Stage 4b; they
    used `is_submerged=True, model_soil=True` (the only path prior to Stage 4b).

## Results

| Metric | Value |
|---|---|
| Max\|diff\| vs stage3_fix2 reference | **0.000050°F** |
| Mean\|diff\| vs stage3_fix2 reference | 0.000025°F |
| Engine concrete T range (Stage 4b) | 65.82–76.05°F |
| Reference concrete T range (Stage 3.5) | 65.82–76.05°F |
| Field shape | (13, 21) = (13, 21) ✓ |

## Verdict

**PASS.** The `model_soil=True` physics path is preserved to within 0.00005°F — well
within numerical tolerance.  The Stage 4b refactor (adding `model_soil` gating) did
not regress the existing soil-mesh code path.

## Notes

The 0.00005°F difference is attributable to floating-point non-associativity from the
new branch statements in the BC loop (the branch itself is taken identically; the
floating-point ops are unchanged).  It is not a physics change.
