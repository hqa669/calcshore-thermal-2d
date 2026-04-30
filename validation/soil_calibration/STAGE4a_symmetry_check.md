# STAGE4a — Side vs Bottom Symmetry Check

## S4a.5 Purpose

Compare α_c extracted via M1 (side-direction penetration, from w=0.00 m Dirichlet) against
M2 (bottom-direction penetration, from depth=24.38 m Dirichlet) using CW data alone.
This tests whether CW itself is geometrically symmetric in its thermal treatment of the
side and bottom boundaries.

- Agreement within 5%  → CW is symmetric; any engine bottom-CL bug is engine-side only.
- Disagreement > 5%    → CW has an asymmetry that must be understood before Stage 4b.

## S4a.5 Data

All eight non-A runs; α_c extracted by `stage4a_extract.py`.

| Run | ΔT (°F) | M1 α_c (m²/hr) | M2 α_c (m²/hr) | M2/M1 ratio |
|---|---|---|---|---|
| B | −13 | 0.00447 | 0.00455 | 1.018 |
| C | +17 | 0.00466 | 0.00473 | 1.015 |
| D | +13 | 0.00475 | 0.00483 | 1.017 |
| E | −17 | 0.00411 | 0.00418 | 1.017 |
| F | −28 | 0.00447 | 0.00452 | 1.011 |
| G | +27 | 0.00488 | 0.00493 | 1.010 |
| H | +28 | 0.00475 | 0.00480 | 1.011 |
| I | −27 | 0.00433 | 0.00439 | 1.014 |

**M2/M1 mean ratio: 1.014 (1.4% higher)**
**M2/M1 std: 0.003**

## S4a.5 Result

M1 and M2 agree to within **1.0–1.8%** across all eight runs — well within the 5%
threshold.

CW is geometrically symmetric: the thermal penetration rate from the side Dirichlet and
from the bottom Dirichlet are consistent to within method precision.  The 1.4% systematic
offset (M2 slightly above M1) is within the quantization noise from the 0.508 m grid
spacing: the bottom direction (24.38 m depth range) has the same 0.508 m pitch as the
side (6.10 m width range), so both are equally quantized.

## S4a.5 Interpretation

The small M2 > M1 bias has a geometric explanation: both directions use a 1°F front
tracked on a 0.508 m grid, but the two directions have different domain lengths.  The
bottom direction (24.38 m, 49 nodes) has more nodes available to track the front, so
the front position estimate benefits from slightly more spatial resolution at early times.
Quantization noise is symmetric and slightly smaller per step in the longer axis.

**This does not represent a material asymmetry.**

## S4a.5 Verdict

**CW is symmetric to within 1.4%.  The engine bottom-CL corner bug is engine-side only.**

CW applies the same Dirichlet mechanism at both the side face (w=0.00 m) and the bottom
face (depth=24.38 m), and the thermal material (concrete) is the same in both directions.
The engine's asymmetry arises because it adds 3 m of soil below the concrete on the
bottom side only — there is no analogous soil extension on the side in the engine's
`is_submerged=False` mode (the Stage 3 fix).

This finding simplifies Stage 4b: the fix need only address the bottom boundary.  No
adjustment to the side boundary treatment is indicated.
