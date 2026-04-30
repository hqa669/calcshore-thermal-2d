# STAGE4a — Engine vs CW Comparison

## S4a.4 Parameters: what can and cannot be compared

CW's model geometry (40ft × 80ft) contains **only concrete** — there is no soil domain in
the CW output.  CW applies Dirichlet T_soil directly at the concrete side face (w=0.00 m)
and concrete bottom face (depth=24.38 m).  Consequently:

| Parameter | Extractable from CW? | Reason |
|---|---|---|
| α_c (concrete diffusivity) | **Yes** — three methods | Concrete fills the entire domain |
| α_s (soil diffusivity) | **No** | No soil domain in CW output |
| e_c/e_s (admittance ratio) | **No** | CW uses Dirichlet (not Robin) at all interfaces |

The engine parameter staleness question (predate `455d3c3`) applies only to α_c.

## S4a.4 α_c comparison

### CW-extracted values (Stage4a extract.py, Runs B–I)

| Method | Mean α_c (m²/hr) | Std | Rel std |
|---|---|---|---|
| M1 — side penetration slope | 0.00455 | 0.00024 | 5.2% |
| M2 — bottom penetration slope | 0.00462 | 0.00024 | 5.2% |
| M3 — erfc profile fit (most reliable) | 0.00528 | 0.00015 | 2.8% |

M1/M2 are biased low by the 0.508 m grid quantization of the 1°F front.  M3 (direct
curve_fit to the erfc spatial profile) is the most reliable estimator: RMSE 0.07–0.26°F,
amplitude A ≈ ΔT (confirms clean Dirichlet at x=0), narrow inter-run spread.

The M3 timestamp trend (t=24hr → t=168hr: α_c decreases from ~0.0055 to ~0.0050)
reflects the finite domain: at long times the zero-flux CL BC at w=6.10 m and the top
BC begin to corrupt the pure half-space erfc.  The t=24–48hr M3 values are cleanest.

**CW consensus α_c: ~0.0050–0.0053 m²/hr** (using M3 at t=48–168 hr).

### Engine values (STAGE4a_engine_thermal_params.md)

Concrete properties come from CW input files (identical across all 9 runs):

| Quantity | Engine value | Source |
|---|---|---|
| k_uc (unhydrated concrete k) | 2.691 W/(m·K) | CW input files (BTU→SI converted) |
| ρ_c | 2101.2 kg/m³ | CW mix: `total_solids/27 × (1−air%)` |
| Cp_c(α_hyd) | Van Breugel model | `specific_heat_variable()`, L361 |
| k_c(α_hyd) | k_uc × (1.33 − 0.33α) | `thermal_conductivity_variable()`, L341 |

Engine α_c varies with hydration degree:

| α_hyd | k_c (W/m·K) | Cp_c (J/kg·K) | α_c (m²/hr) |
|---|---|---|---|
| 0.0 (start) | 3.580 | 1139.8 | 0.00538 |
| 0.5 | 3.136 | 1111.4 | 0.00481 |
| 0.8 (representative) | 2.869 | 1094.4 | 0.00449 |
| 1.0 (fully hydrated) | 2.691 | 1083.0 | 0.00426 |

At t=168 hr, most of the concrete has hydrated to α_hyd ≈ 0.7–0.9, so **α_c_engine ≈
0.00449 m²/hr** (at α_hyd=0.8) is the representative value.

### Ratio table

| Metric | Engine | CW (M3 consensus) | Ratio CW/Engine | Within 10%? |
|---|---|---|---|---|
| α_c at α_hyd=0.8 | 0.00449 m²/hr | 0.00528 m²/hr | 1.18 | **No (~18% off)** |
| α_c at α_hyd=0.5 | 0.00481 m²/hr | 0.00528 m²/hr | 1.10 | Borderline |
| α_c at α_hyd=0.0 | 0.00538 m²/hr | 0.00528 m²/hr | 0.98 | Yes |

The engine's effective α_c at late time is ~10–18% below CW's.  This is a secondary
discrepancy — it would affect the **shape** and **rate** of interior temperature evolution,
but would not produce the sharp 9–19°F residual localised at the bottom-CL corner.

## S4a.4 Structural comparison: bottom-temperature gap (M4)

The decisive comparison is direct bottom-temperature:

| Run | ΔT (°F) | CW bot T (°F) | Engine bot T (°F) | Gap (°F) | Stage3.5 residual (°F) |
|---|---|---|---|---|---|
| B | −13 | 60.01 | 69.04 | +9.03 | ~9.1 |
| C | +17 | 90.00 | 78.18 | −11.82 | ~11.9 |
| D | +13 | 73.00 | 63.96 | −9.05 | ~9.1 |
| E | −17 | 73.00 | 84.82 | +11.82 | ~11.9 |
| F | −28 | 45.00 | 64.48 | +19.48 | ~19.6 |
| G | +27 | 100.00 | 81.22 | −18.78 | ~18.9 |
| H | +28 | 73.00 | 53.52 | −19.48 | ~19.5 |
| I | −27 | 73.00 | 91.78 | +18.78 | ~18.8 |

CW bottom temperature = T_soil (Dirichlet applied directly at the concrete bottom face).
Engine bottom temperature ≈ T_placement + fraction(T_soil − T_placement), buffered by 3 m soil.
Gap matches Stage 3.5 masked residual at (w=6.10 m, depth=24.38 m) to within 0.1°F.

## S4a.4 Decision

Applying the brief's verdict rules:

| Ratio band | Criterion | Verdict |
|---|---|---|
| All within ~10% | → Structural only | — |
| One or more > 30% | → Parametric (at least partly) | — |
| **α_c 18% off, structural gap = residual** | → **Both, structural dominant** | ✓ |

**Verdict: STRUCTURAL (primary) + PARAMETRIC (secondary, minor)**

- **Primary (structural)**: The engine places the Dirichlet T_gw 3 m below the concrete
  bottom; CW places it at the concrete bottom face.  The 3 m soil buffer accounts for
  100% of the 9–19°F bottom-CL residual.
- **Secondary (parametric)**: Engine α_c at late-time hydration (≈0.00449 m²/hr) is ~18%
  below CW's (~0.00528 m²/hr).  This would affect interior temperature evolution by
  ~1–3°F across the domain but does not explain the localised bottom-CL spike.

**Stage 4b primary scope**: apply Dirichlet T_gw at the concrete bottom row
(`iy_concrete_end`) rather than at the bottom of the soil buffer (`T_new[-1, :]`).  Or
equivalently: eliminate the soil buffer entirely for the bottom boundary and apply the
Dirichlet directly at the concrete-soil interface, matching CW's boundary formulation.

**Stage 4b secondary scope (optional)**: recalibrate concrete base k to match CW's
effective α_c ≈ 0.00528 m²/hr at representative α_hyd, or revisit the k_uc × (1.33−0.33α)
modulation against CW's actual thermal model.
