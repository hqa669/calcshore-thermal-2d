# Sprint 8 Stage 2-prep v2 — §3.1 Wrapper Validation

## Status: GATE FAIL — Path C inversion is degenerate at specified comparison point

### Engine setup
- Monkey-patch: `te2d.thermal_conductivity_variable` rebound to constant k at α_test
- `K_UC_CALIBRATION_FACTOR_SPRINT7=0.96` active (upstream of patched function)
- `model_soil=False`, baseline grid (no `grid_refinement`), T_amb=100.0°F (constant)
- Comparison point: j_mid=36 (engine CL mid-depth) vs CW di=24 (CL mid-depth)

### Validation results at α_test=0.0361, t=168hr

| Point | Engine T (°F) | CW T (°F) | Gap (°F) | Status |
|---|---|---|---|---|
| CL mid-depth vs lowhyd | 73.0001 | 73.1300 | 0.1299 | TRIVIAL PASS (both ≈ T_placement) |
| CL mid-depth vs highhyd | 73.0001 | 75.2900 | 2.2899 | FAIL (CW bulk shift 2.16°F not reproduced) |
| Top surface (j=0 vs di=0) lowhyd | 77.6356 | 76.4780 | 1.1576 | Mismatch (~1.2°F) |
| Top surface (j=0 vs di=0) highhyd | 77.6356 | 78.4220 | 0.7864 | Mismatch (~0.8°F) |

### Root cause

CW output at t=168hr shows ALL di=8–38 cells (4–18.8 m depth = 15 m interior span) at the SAME temperature: 73.13°F (lowhyd), 75.29°F (highhyd), 74.16°F (highhyd_b010). This uniform bulk interior offset (+2.16°F for highhyd vs lowhyd) cannot be produced by the CalcShore constant-k 2D diffusion engine, because thermal diffusion from either boundary only penetrates ≈2.1 m (1/√(4·α_c·t)) in 168 hours. The engine interior temperature at j_mid remains at T_placement=73.0°F regardless of α_test.

**Path C (constant-k inversion at CL mid-depth) is not feasible for this geometry.** The §3.1 sanity check passes trivially for lowhyd (both near 73°F), but fails for the high-hydration datasets by ~2°F. No α_test sweep will produce a valid inversion at the interior comparison point.
