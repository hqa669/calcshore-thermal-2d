# STAGE2 — Engine vs CW Comparison

## Method

Engine temperature fields (concrete subgrid: 13y × 21x nodes, dx=0.3048 m) are bilinearly resampled onto CW's 49×13 grid using `scipy.interpolate.RegularGridInterpolator`. Engine output is in °C; converted to °F before residual computation.

Engine runs use explicit `T_ground_deep_C = soil_temp_F` (from input.dat), bypassing
the default `compute_T_gw_C(env)` path. This allows the comparison to isolate
model-form differences in how soil is spatially represented, independent of the
known soil-BC plumbing gap. See STAGE2_engine_soil_audit.md §Engine-vs-CW protocol.

Comparison time: t=168 hr.

**Interpolation note**: ~0.1°F noise is expected from grid-resolution mismatch.

## Summary Table

| Run | max\|R\| (°F) | mean\|R\| (°F) | RMS R (°F) | Location of max |
| --- | ------------ | ------------ | ---------- | --------------- |
| A | 1.655 | 0.406 | 0.617 | depth=2.03m, w=0.00m |
| B | 11.436 | 2.764 | 4.511 | depth=0.00m, w=0.00m |
| C | 18.647 | 4.262 | 7.064 | depth=2.03m, w=0.00m |
| D | 14.581 | 3.306 | 5.502 | depth=2.54m, w=0.00m |
| E | 15.292 | 3.677 | 6.021 | depth=0.00m, w=0.00m |
| F | 26.448 | 6.250 | 10.292 | depth=0.00m, w=0.00m |
| G | 28.655 | 6.578 | 10.912 | depth=2.03m, w=0.00m |
| H | 29.498 | 6.749 | 11.246 | depth=4.06m, w=0.00m |
| I | 25.209 | 5.988 | 9.853 | depth=13.21m, w=0.00m |

## Figures

- `plots/runA_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runB_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runC_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runD_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runE_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runF_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runG_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runH_engine_vs_cw.png` — CW | engine | residual at t=168 hr
- `plots/runI_engine_vs_cw.png` — CW | engine | residual at t=168 hr
