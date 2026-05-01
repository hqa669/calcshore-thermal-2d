# STAGE3 Fix 2 — Residual Summary

**Fix applied:** soil_temp_F plumbed + is_submerged=True (concrete sides contact soil)  
**Comparison time:** t=168 hr  
**Method:** bilinear interpolation of engine (21×13) onto CW (49×13) grid.

## Residual Table

| Run | Placement/Soil | max\|R\| (°F) | mean\|R\| (°F) | RMS R (°F) | Location of max |
| --- | -------------- | -------------- | --------------- | ------------ | --------------- |
| A | 73/73 | 5.445 | 0.349 | 0.987 | depth=0.00m, w=4.06m |
| B | 73/60 | 13.673 | 2.418 | 3.881 | depth=0.00m, w=0.51m |
| C | 73/90 | 11.873 | 2.735 | 4.465 | depth=24.38m, w=6.10m |
| D | 60/73 | 9.086 | 2.167 | 3.435 | depth=24.38m, w=6.10m |
| E | 90/73 | 13.217 | 2.878 | 4.655 | depth=0.00m, w=0.51m |
| F | 73/45 | 24.274 | 4.912 | 7.896 | depth=0.00m, w=0.51m |
| G | 73/100 | 18.868 | 4.286 | 7.130 | depth=24.38m, w=6.10m |
| H | 45/73 | 19.566 | 4.409 | 7.284 | depth=24.38m, w=6.10m |
| I | 100/73 | 18.860 | 4.411 | 7.227 | depth=24.38m, w=6.10m |

## Figures

- `plots/STAGE3_runA_engine_vs_cw.png`
- `plots/STAGE3_runB_engine_vs_cw.png`
- `plots/STAGE3_runC_engine_vs_cw.png`
- `plots/STAGE3_runD_engine_vs_cw.png`
- `plots/STAGE3_runE_engine_vs_cw.png`
- `plots/STAGE3_runF_engine_vs_cw.png`
- `plots/STAGE3_runG_engine_vs_cw.png`
- `plots/STAGE3_runH_engine_vs_cw.png`
- `plots/STAGE3_runI_engine_vs_cw.png`
