# STAGE1 — Penetration Depth Analysis

## Method

1°F thermal front penetration from each soil-concrete interface vs √t.
Linear fit across full 168hr run → effective slope (m/√hr).
In a half-space with step BC: penetration ≈ C·√(α_eff·t).

**Grid resolution caveat**: CW's output grid has 0.51m node spacing.
The 1°F front position is quantized to 0.51m steps. This makes the step
curves jagged and degrades R² even when the underlying physics is √t-like.
The slope and R² values should be treated as indicative, not precise.

- Threshold: 1.0°F from placement temperature
- Linear fit window: 0 – 168.0 hr (full run)
- Side: walk from edge (0 m) inward at mid-depth (12.19 m)
- Bottom: walk from bottom upward on centerline (w=6.10 m)

## Results

| Run | Placement°F | Soil°F | ΔT | side slope (m/√hr) | side R² | side pen@168hr (m) | bot slope (m/√hr) | bot R² | bot pen@168hr (m) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A | 73 | 73 | 0 | — (baseline) | — | — | — | — | — |
| B | 73 | 60 | -13 | 0.1672 | 0.9138 | 2.03 | 0.1687 | 0.9142 | 2.03 |
| C | 73 | 90 | +17 | 0.1824 | 0.9255 | 2.54 | 0.1837 | 0.9258 | 2.54 |
| D | 60 | 73 | +13 | 0.1725 | 0.9252 | 2.03 | 0.1739 | 0.9255 | 2.03 |
| E | 90 | 73 | -17 | 0.1713 | 0.9260 | 2.03 | 0.1728 | 0.9262 | 2.03 |
| F | 73 | 45 | -28 | 0.1985 | 0.9389 | 2.54 | 0.1997 | 0.9393 | 2.54 |
| G | 73 | 100 | +27 | 0.2060 | 0.9474 | 2.54 | 0.2072 | 0.9476 | 2.54 |
| H | 45 | 73 | +28 | 0.2047 | 0.9461 | 2.54 | 0.2058 | 0.9464 | 2.54 |
| I | 100 | 73 | -27 | 0.1940 | 0.9340 | 2.54 | 0.1954 | 0.9345 | 2.54 |

## Interpretation

- If slope values cluster across runs: single effective diffusivity, model is
  consistent (property is not ΔT-dependent).
- If slopes scatter: model-form non-linearity or material property mismatch.
- R² < 0.95 on early-time linear fit: penetration is not √t-like even at early time
  (e.g., BC-dominated or grid-resolution artifact).
- R² drops at long time (deviation from linearity) is expected when the finite
  soil domain saturates to Dirichlet BC temperature.

## Plots

- `plots/penetration_side.png`
- `plots/penetration_bottom.png`
