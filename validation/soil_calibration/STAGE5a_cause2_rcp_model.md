# STAGE5a — Cause 2: ρCp Model Assessment

## Van Breugel specific_heat_variable — formula

```
Cp = (1/ρ) × [Wc·α·(8.4·T + 339) + Wc·(1−α)·Ca + Wa·Ca + Ww·4186]
```

## Run B mix composition (SI units)

| Parameter | Value | Units |
| --- | --- | --- |
| ρ_concrete | 2101.15 | kg/m³ |
| Wc (cement) | 341.13 | kg/m³ |
| Wa (aggregate) | 1720.50 | kg/m³ |
| Ww (water) | 150.10 | kg/m³ |
| Ca (aggregate Cp) | 856.90 | J/(kg·K) |
| Wc+Wa+Ww total | 2211.73 | kg/m³ |

## At α_hyd=0.0359, T_soil_C (cold) (15.6°C / 60.0°F)

c_cef = 8.4×15.6 + 339 = 469.7 J/(kg·K)

| Term | Contribution (J/m³·K) | % of total |
| --- | --- | --- |
| Wc·α·c_cef (hydrated cement) | 5751 | 0.2% |
| Wc·(1−α)·Ca (unhydrated cement) | 281823 | 11.8% |
| Wa·Ca (aggregate) | 1474291 | 61.7% |
| Ww·4186 (water) | 628314 | 26.3% |
| **ρCp total** | **2390179** | 100% |

## At α_hyd=0.0359, T_mid (avg) (19.2°C / 66.5°F)

c_cef = 8.4×19.2 + 339 = 500.0 J/(kg·K)

| Term | Contribution (J/m³·K) | % of total |
| --- | --- | --- |
| Wc·α·c_cef (hydrated cement) | 6123 | 0.3% |
| Wc·(1−α)·Ca (unhydrated cement) | 281823 | 11.8% |
| Wa·Ca (aggregate) | 1474291 | 61.7% |
| Ww·4186 (water) | 628314 | 26.3% |
| **ρCp total** | **2390551** | 100% |

## At α_hyd=0.0359, T_pl (placement) (22.8°C / 73.0°F)

c_cef = 8.4×22.8 + 339 = 530.3 J/(kg·K)

| Term | Contribution (J/m³·K) | % of total |
| --- | --- | --- |
| Wc·α·c_cef (hydrated cement) | 6494 | 0.3% |
| Wc·(1−α)·Ca (unhydrated cement) | 281823 | 11.8% |
| Wa·Ca (aggregate) | 1474291 | 61.7% |
| Ww·4186 (water) | 628314 | 26.3% |
| **ρCp total** | **2390922** | 100% |

## Temperature sensitivity of ρCp

T swing ±10°C around mid-range 19.2°C:
  ρCp(9°C) = 2389522 J/m³·K
  ρCp(19°C)    = 2390551 J/m³·K
  ρCp(29°C) = 2391579 J/m³·K
  Swing: 0.1% — negligible (<5%)

## Comparison to textbook normal-weight concrete

Textbook ρCp range: 2.1–2.6 ×10⁶ J/m³·K
Engine ρCp at operating point (mid-range T): 2390551 J/m³·K = 2.391 ×10⁶ J/m³·K
→ Within textbook range.

## Verdict

If ρCp = 2390551 J/m³·K is accepted as correct, the entire α_c gap must come from k.
If ρCp were as high as 2600000 (upper textbook), α_c would be 0.00491 m²/hr — worse (lower), not better.
Van Breugel ρCp appears within textbook range — Cause 2 is unlikely to be the dominant driver of the α_c gap.
