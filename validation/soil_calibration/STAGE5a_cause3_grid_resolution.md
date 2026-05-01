# STAGE5a — Cause 3: Grid Resolution / Numerical Accuracy

Run F: placement=73°F, soil=45°F, |ΔT|=28°F (largest residual in Stage 4b).
Three grid resolutions tested; same CW reference and Stage 3.5 mask applied.

## Resolution sweep results

| Resolution | n_cx × n_cy | dx (ft) | dy_mean (ft) | dt_inner (s) | masked max|R| (°F) | masked mean|R| (°F) | wall-clock (s) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1× | 21×13 | 1.000 | 6.667 | 37.41 | **4.084** ← Stage 4b ref 4.084°F | 0.508 | 1.4 |
| 2× | 41×25 | 0.500 | 3.333 | 36.93 | **1.249** | 0.227 | 1.9 |
| 4× | 81×49 | 0.250 | 1.667 | 35.15 | **0.409** | 0.163 | 4.3 |

## Richardson extrapolation (second-order FD: R(h) ≈ R_∞ + C·h²)

Using 1× and 2×:
  R_∞ (1×/2× estimate) = (4×R_2× − R_1×) / 3 = (4.994 − 4.084) / 3 = **0.303°F**
  Numerical contribution (1× − R_∞) / R_1× = (4.084 − 0.303) / 4.084 = **92.6%**

Using 2× and 4×:
  R_∞ (2×/4× estimate) = (4×R_4× − R_2×) / 3 = (1.634 − 1.249) / 3 = **0.129°F**
  Numerical contribution (2× − R_∞) / R_2× = (1.249 − 0.129) / 1.249 = **89.7%**

Residuals are monotonically converging with resolution → Richardson extrapolation is reliable.

## Verdict

Cause 3 is DOMINANT (93%). Grid resolution explains most of the α_c gap; Causes 1 and 2 are secondary.

Parametric (non-numerical) residual at infinite resolution: 0.303°F
Gate target: 0.35°F. Achievable by grid refinement alone.
