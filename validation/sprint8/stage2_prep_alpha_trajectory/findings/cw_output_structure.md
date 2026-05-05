# Sprint 8 Stage 2-prep — §3.1 CW Output Structure Finding

## Result: Path A unavailable. Path B (ratio test) required.

CW's `output.txt` contains exactly **640 columns per data row**:

| Column(s) | Content |
|---|---|
| 1 | Time (hr) |
| 2–638 | T(w,d) temperature field — 637 spatial cells (13 widths × 49 depths), °C |
| 639 | Gradient (max − min T across section, °C) |
| 640 | Ambient temperature (°C) |

CW outputs **no degree of hydration (α), DOH, maturity, or equivalent age** anywhere
in this file. No separate section, no additional file checked. The entire file is a
pure T(x,z,t) tensor plus two scalar summary columns.

**Implication:** CW's internal α(t) trajectory is not directly observable from the
output files. It must be inferred indirectly — §3.2 uses the ratio test approach on
the Stage 1 pairwise ΔT data.

## Raw header (line 2 of output.txt)

```
Time 6.1/0 5.59/0 5.08/0 ... 0/24.38 Gradient Ambient
```

Coordinates are `width_m/depth_m` measured from the upper corner. Widths decrease
6.1→0 (CL is at w=6.1, form face at w=0). Depths increase 0→24.38 (top to bottom).
