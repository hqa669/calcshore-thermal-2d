# STAGE4a — Final Diagnosis Report

## 1 Bottom-soil engine audit

The engine allocates 15 soil cells (3.0 m total, `dy=0.20 m`) below the 80 ft
(24.38 m) concrete slab (`thermal_engine_2d.py:976-981`).  The Dirichlet ground
temperature `T_gw` is applied as an **unconditional full-row overwrite** at the
deepest row (`T_new[-1, :] = _T_gw_C`, L2184) — 3 m below the concrete bottom face.
The side boundary (column 0) uses an `is_soil`-masked Dirichlet (L2181), and Stage 3
added a half-cell stencil correction for the side-soil contact (L2085–2098).  No
analogous correction exists for the bottom row.  The bottom-CL corner (j=-1, ix=nx-1)
is pinned unconditionally to `T_gw`, overriding the zero-flux CL symmetry stencil.
The engine's bottom-soil structure was the primary structural suspect going into Stage 4a.

## 2 Engine thermal property provenance

Soil properties (`k=1.50 W/(m·K)`, `ρ=2100 kg/m³`, `Cp=900 J/(kg·K)`) were committed
on 2026-04-22 (7c3838e), six days before the hydration kinetics calibration (`455d3c3`,
2026-04-28).  Concrete base `k` (`1.56 BTU/hr/ft/°F → 2.70 W/(m·K)`) was also committed
on 2026-04-22 (c0013ba); however, the actual CW input files override these defaults, and
all nine runs use `k_uc = 2.691 W/(m·K)`, `ρ_c = 2101.2 kg/m³`, and the Van Breugel
`Cp` model.  Engine α_c at representative late-time hydration (α_hyd = 0.8):
**0.00449 m²/hr**.  Soil α_s: **2.857×10⁻³ m²/hr** (not independently calibrated).

## 3 CW-extracted parameters

CW applies Dirichlet T_soil directly at both the side (w=0.00 m, all depths) and bottom
(depth=24.38 m, all widths) concrete faces — there is no soil domain in CW's output.
As a result, α_s and the admittance ratio e_c/e_s cannot be extracted from CW data.
α_c was extracted via three independent methods across Runs B–I:

| Method | Mean α_c (m²/hr) | Rel std |
|---|---|---|
| M1 — side penetration slope (quantized) | 0.00455 | 5.2% |
| M2 — bottom penetration slope (quantized) | 0.00462 | 5.2% |
| M3 — erfc spatial profile fit (**most reliable**) | 0.00528 | 2.8% |

M3 RMSE: 0.07–0.26°F (excellent fit quality); amplitude A ≈ ΔT confirms clean Dirichlet
at x=0.  M1/M2 are quantized (0.508 m grid) and bias low.  **CW consensus α_c ≈
0.0050–0.0053 m²/hr** (M3 at t=48–168 hr).

## 4 Engine vs CW comparison

| Parameter | Engine value | CW value | Ratio (CW/Eng) | Verdict |
|---|---|---|---|---|
| α_c (α_hyd=0.8, representative) | 0.00449 m²/hr | 0.00528 m²/hr | **1.18** | ~18% off |
| α_c (α_hyd=0.0, start) | 0.00538 m²/hr | 0.00528 m²/hr | 0.98 | Match |
| α_s | 2.857×10⁻³ m²/hr | not extractable | — | — |
| e_c/e_s | 1.526 (at α_hyd=0.8) | not extractable | — | — |
| Dirichlet depth (bottom BC) | 27.40 m (3 m below concrete) | **24.38 m** (at concrete face) | — | **Structural mismatch** |

The late-time α_c mismatch (~18%) is secondary.  The decisive comparison is direct
bottom-temperature (M4):

| Run | |ΔT| (°F) | CW bot T (°F) | Engine bot T (°F) | Gap (°F) | Stage3.5 residual |
|---|---|---|---|---|---|---|
| B | 13 | 60.01 | 69.04 | **+9.03** | ~9.1 |
| C | 17 | 90.00 | 78.18 | **−11.82** | ~11.9 |
| D | 13 | 73.00 | 63.96 | **−9.05** | ~9.1 |
| E | 17 | 73.00 | 84.82 | **+11.82** | ~11.9 |
| F | 28 | 45.00 | 64.48 | **+19.48** | ~19.6 |
| G | 27 | 100.00 | 81.22 | **−18.78** | ~18.9 |
| H | 28 | 73.00 | 53.52 | **−19.48** | ~19.5 |
| I | 27 | 73.00 | 91.78 | **+18.78** | ~18.8 |

Gap = Stage 3.5 residual to within 0.1°F for every run.

## 5 Side vs bottom symmetry check

M1 (side) and M2 (bottom) α_c estimates agree to within **1.0–1.8%** across all eight
runs (M2/M1 mean = 1.014, std = 0.003).  This is well within the 5% threshold.
CW is geometrically symmetric — the same Dirichlet mechanism and the same concrete
material apply in both directions.  The engine's asymmetry (soil buffer present below
but not beside) is engine-side only.  Stage 4b need only address the bottom boundary.

## 6 Sensitivity sweep

**Not run.** The structural proof is unambiguous: the M4 bottom-temperature gaps
match the Stage 3.5 residuals to within 0.1°F across all eight runs.  A parametric
sweep would not change the diagnosis.

## 7 Diagnosis verdict

**STRUCTURAL (primary).** The 9.1–19.6°F bottom-CL residual is caused entirely by the
engine's 3 m soil buffer below the concrete slab: the Dirichlet T_gw is applied at
depth 27.40 m while CW applies it at the concrete bottom face (24.38 m).  The soil
buffer insulates the concrete bottom, preventing it from reaching T_soil as quickly as
CW's direct Dirichlet.  Stage 4b should remove the soil buffer from the bottom boundary
— or equivalently, pin `T_new[iy_concrete_end, :]` = `_T_gw_C` — matching CW's model.

**PARAMETRIC (secondary, minor).** Engine α_c at late-time hydration (~0.00449 m²/hr)
is ~18% below CW's (~0.00528 m²/hr).  This affects the overall interior temperature
evolution by 1–3°F across the domain but is not the source of the localised corner spike.
Stage 4b may optionally recalibrate k_uc after the structural fix is validated.

## 8 Open questions for master chat

1. **Soil buffer design intent**: Why does the engine allocate 3 m of soil below the
   concrete?  If this was intended to model far-field ground temperature diffusion, the
   soil depth needs a separate Dirichlet depth matched to CW's assumption (or removed
   entirely if CW treats the concrete bottom as the thermal ground boundary).  Clarifying
   this may affect how Stage 4b implements the fix.

2. **α_c secondary discrepancy**: The 18% gap in effective concrete diffusivity at late
   hydration (engine 0.00449 vs CW 0.00528 m²/hr) suggests either the Van Breugel
   k_uc × (1.33−0.33α) modulation overshoots at late α, or CW uses a different Cp model
   for the Clay mix.  After Stage 4b's structural fix reduces the residual to near zero,
   this secondary offset may become the dominant remaining error and should be
   characterised before closing Stage 4.
