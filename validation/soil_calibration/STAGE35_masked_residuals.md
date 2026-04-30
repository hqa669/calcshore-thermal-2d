# STAGE 3.5 — Masked-Region Residual Report

**Mask definition:** cells where |residual_RunA| < 0.3°F at t=168 hr.  
**Mask coverage:** 89.8% of grid (572 of 637 cells kept).  
**Engine source:** `stage3_fix2_runs/run{A..I}_t168.csv`  
**Comparison time:** t=168 hr  
**Method:** bilinear interpolation of engine (21×13) onto CW (49×13) grid.

## Residual Table

| Run | Placement/Soil (°F) | Unmasked max\|R\| (°F) | Masked max\|R\| (°F) | Masked mean\|R\| (°F) | Masked max location (width_ft, depth_ft) |
| --- | ------------------- | ----------------------- | --------------------- | ---------------------- | --------------------------------------- |
| A | 73/73 | 5.445 | 0.235 | 0.097 | (3.35, 8.33) |
| B | 73/60 | 13.673 | 9.074 | 2.055 | (20.01, 79.99) |
| C | 73/90 | 11.873 | 11.873 | 2.769 | (20.01, 79.99) |
| D | 60/73 | 9.086 | 9.086 | 2.111 | (20.01, 79.99) |
| E | 90/73 | 13.217 | 11.873 | 2.669 | (20.01, 79.99) |
| F | 73/45 | 24.274 | 19.566 | 4.420 | (20.01, 79.99) |
| G | 73/100 | 18.868 | 18.868 | 4.348 | (20.01, 79.99) |
| H | 45/73 | 19.566 | 19.566 | 4.444 | (20.01, 79.99) |
| I | 100/73 | 18.860 | 18.860 | 4.234 | (20.01, 79.99) |

## Gate Verdict (masked region)

- **Gate A:** Run A masked max\|R\| = 0.235°F ≤ 0.5°F → **PASS**
- **Gate B–I:** worst-case Run H masked max\|R\| = 19.566°F ≤ 2.0°F → **FAIL**
- **Overall:** **FAIL**

Stage 4 calibration needed: **yes**

## Bottom-corner artifact disposition

13 of 13 bottom-row cells (depth ≈ 24.38 m) survive the mask.

### Bottom-row max|R| per run (masked cells only)

| Run | ΔT (soil−pl, °F) | Bottom-row max\|R\| (°F) | Location (width_ft) |
| --- | ----------------- | ------------------------- | ------------------- |
| B | -13 | 9.074 | 20.01 |
| C | +17 | 11.873 | 20.01 |
| D | +13 | 9.086 | 20.01 |
| E | -17 | 11.873 | 20.01 |
| F | -28 | 19.566 | 20.01 |
| G | +27 | 18.868 | 20.01 |
| H | +28 | 19.566 | 20.01 |
| I | -27 | 18.860 | 20.01 |

### ΔT-scaling check

Paired runs share the same |ΔT|; if the bottom-row residual tracks |ΔT| monotonically the artifact is soil-driven.

- |ΔT|=13°F (B/D): avg bottom-row max|R|=9.080°F
- |ΔT|=17°F (C/E): avg bottom-row max|R|=11.873°F
- |ΔT|=27°F (G/I): avg bottom-row max|R|=18.864°F
- |ΔT|=28°F (F/H): avg bottom-row max|R|=19.566°F

**Monotonically scales with |ΔT|: yes**

### Symmetry check (CL w≈6.10m vs outer edge w≈0m at bottom row)

Note: in the CW/engine arrays, column 0 = CL (max width = 6.10 m); column -1 = outer edge (w ≈ 0 m).

- Run B: CL (w=6.10m) |R|=9.074°F, outer edge (w=0.00m) |R|=5.837°F
- Run C: CL (w=6.10m) |R|=11.873°F, outer edge (w=0.00m) |R|=7.640°F
- Run D: CL (w=6.10m) |R|=9.086°F, outer edge (w=0.00m) |R|=5.849°F
- Run E: CL (w=6.10m) |R|=11.873°F, outer edge (w=0.00m) |R|=7.639°F
- Run F: CL (w=6.10m) |R|=19.566°F, outer edge (w=0.00m) |R|=12.594°F
- Run G: CL (w=6.10m) |R|=18.868°F, outer edge (w=0.00m) |R|=12.144°F
- Run H: CL (w=6.10m) |R|=19.566°F, outer edge (w=0.00m) |R|=12.594°F
- Run I: CL (w=6.10m) |R|=18.860°F, outer edge (w=0.00m) |R|=12.136°F

**Note:** Do not attempt a fix in this session. If the artifact is confirmed soil-coupled (large magnitude, scales with |ΔT|), flag as a Stage 4 task.
