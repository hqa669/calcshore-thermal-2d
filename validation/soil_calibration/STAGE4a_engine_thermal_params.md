# STAGE4a — Engine Thermal Property Values and Provenance

## S4a.2 Hardcoded values

### Soil properties

File: `thermal_engine_2d.py:229`

```python
SOIL_PROPERTIES_2D: dict = {'k': 1.50, 'rho': 2100.0, 'Cp': 900.0}
```

| Property | Value | Units |
|---|---|---|
| k_soil | 1.50 | W/(m·K) |
| ρ_soil | 2100.0 | kg/m³ |
| Cp_soil | 900.0 | J/(kg·K) |
| ρ·Cp_soil | 1,890,000 | J/(m³·K) |

Git blame: commit **7c3838e** (2026-04-22). Predates `455d3c3` (2026-04-28, hydration kinetics calibration) by **6 days**. Unverified against CW Clay since the kinetics calibration.

Override at runtime: `solve_hydration_2d(..., soil_props={'k': X, 'rho': Y, 'Cp': Z})`. Merge at `thermal_engine_2d.py:1515`.

### Concrete properties

File: `cw_scenario_loader.py:161–162`

```python
thermal_conductivity_BTU_hr_ft_F: float = 1.56   # → 2.702 W/(m·K)
aggregate_Cp_BTU_lb_F: float = 0.20              # → 837 J/(kg·K)
```

Git blame: commit **c0013ba** (2026-04-22). Also predates `455d3c3` by 6 days.

Note: these are **defaults** only. CW input files specify mix-specific values.
The actual values loaded from `runA_baseline/input.dat` through `runI_100_73/input.dat`
are identical across all 9 runs (same mix):

| Property | Base value | Units |
|---|---|---|
| thermal_conductivity_BTU_hr_ft_F | 1.5551 | BTU/hr/ft/°F |
| k_concrete (converted) | 2.6914 | W/(m·K) |
| aggregate_Cp_BTU_lb_F | 0.2047 | BTU/lb/°F |
| Ca (aggregate Cp, converted) | 856.9 | J/(kg·K) |
| concrete_density_lb_ft3 (computed from mix) | 131.17 | lb/ft³ |
| ρ_concrete (converted) | 2101.2 | kg/m³ |

Concrete `k` varies with hydration degree α via `thermal_conductivity_variable` (L341):
`k_c(α) = k_uc × (1.33 − 0.33 × α)`
— ranges from 1.33 × k_uc (fully unhydrated, α=0) to k_uc (fully hydrated, α=1).

Concrete `Cp` varies with α and T via `specific_heat_variable` (Van Breugel model, L361),
using the full mix composition (cement + aggregate + water weights per m³):

| Hydration degree α | k_c (W/m·K) | Cp_c (J/kg·K) | ρ·Cp_c (MJ/m³·K) |
|---|---|---|---|
| 0.0 (start) | 3.580 | 1139.8 | 2.395 |
| 0.5 | 3.136 | 1111.4 | 2.335 |
| 0.8 (representative) | 2.869 | 1094.4 | 2.299 |
| 1.0 (full) | 2.691 | 1083.0 | 2.276 |

## S4a.2 Implied diffusivities

For soil (constant properties):

```
α_s_engine = k_s / (ρ_s · Cp_s)
           = 1.50 / (2100 × 900)
           = 7.937 × 10⁻⁷ m²/s
           = 2.857 × 10⁻³ m²/hr
```

For concrete (representative at α_hyd = 0.8):

```
α_c_engine = k_c(0.8) / (ρ_c · Cp_c(0.8))
           = 2.869 / (2101.2 × 1094.4)
           = 1.247 × 10⁻⁶ m²/s
           = 4.488 × 10⁻³ m²/hr
```

At α_hyd = 0.0 (start): α_c = 5.382 × 10⁻³ m²/hr (highest, before any hydration)
At α_hyd = 1.0 (full): α_c = 4.261 × 10⁻³ m²/hr (lowest)

## S4a.2 Implied admittance ratio (e_c/e_s)

Thermal effusivity `e = √(k · ρ · Cp)`.

At α_hyd = 0.8:
```
e_c = √(k_c · ρ_c · Cp_c) = √(2.869 × 2101.2 × 1094.4)
    = √(6,603,058) = 2569.6 J/(m²·K·s^0.5)

e_s = √(k_s · ρ_s · Cp_s) = √(1.50 × 2100 × 900)
    = √(2,835,000) = 1683.8 J/(m²·K·s^0.5)

(e_c/e_s)_engine = 2569.6 / 1683.8 = 1.526
```

| Hydration degree α | e_c/e_s |
|---|---|
| 0.0 | 1.700 |
| 0.5 | 1.624 |
| 0.8 | 1.526 |
| 1.0 | 1.457 |

## S4a.2 Summary and staleness assessment

| Parameter | Engine value | Commit date | Prior to 455d3c3? |
|---|---|---|---|
| k_soil | 1.50 W/(m·K) | 2026-04-22 | **Yes (6 days before)** |
| ρ_soil·Cp_soil | 1.89 MJ/(m³·K) | 2026-04-22 | **Yes** |
| k_concrete_base | 2.691 W/(m·K) | 2026-04-22 | **Yes** |
| α_c (at α_hyd=0.8) | 4.49×10⁻³ m²/hr | derived | — |
| α_s | 2.86×10⁻³ m²/hr | derived | — |
| (e_c/e_s) at α_hyd=0.8 | 1.526 | derived | — |

All thermal property values predate the April 28 hydration kinetics calibration.
They have not been verified against CW's effective Clay soil model.

The concrete `k` and `Cp` also predate any verification against CW's concrete thermal
model. The hydration-degree dependence (`thermal_conductivity_variable`) is based on
CW Eq 23/24/25, but the base values were established at project initialization and
may not reflect CW's actual material assumptions for the Clay/Clay validation dataset.
