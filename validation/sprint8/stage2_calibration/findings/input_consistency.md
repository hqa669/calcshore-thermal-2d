# Sprint 8 Stage 2 §3.2 — Input.dat Consistency

## Sprint 7 Reference (runA_baseline)

- width_ft: 40.0
- depth_ft: 80.0
- cement_lb_yd3: 575.0
- water_lb_yd3: 253.0
- coarse_lb_yd3: 1800.0
- fine_lb_yd3: 1100.0
- k_BTU: 1.55509402162759
- alpha_u: 0.1
- tau_hrs: 200.0
- beta: 0.1
- Hu_J_kg: 1.0
- Ea_J_mol: 50000.0
- T_pl_F: 73.0
- T_soil_F: 73.0
- form_type: steel

## New Dataset Verification

| Folder | α_u | T_pl (°F) | T_soil (°F) | τ | β | Hu | width_ft | depth_ft | k_BTU | Status |
|---|---|---|---|---|---|---|---|---|---|---|
| thermal_alpha02_A_73_73 | 0.20 | 73 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha02_F_73_45 | 0.20 | 73 | 45 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha02_I_100_73 | 0.20 | 100 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha04_A_73_73 | 0.40 | 73 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha04_F_73_45 | 0.40 | 73 | 45 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha04_I_100_73 | 0.40 | 100 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha06_A_73_73 | 0.60 | 73 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha06_F_73_45 | 0.60 | 73 | 45 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha06_I_100_73 | 0.60 | 100 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha08_A_73_73 | 0.80 | 73 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha08_F_73_45 | 0.80 | 73 | 45 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |
| thermal_alpha08_I_100_73 | 0.80 | 100 | 73 | 5.0 | 0.85 | 1.0 | 40 | 80 | 1.5551 | OK |

**Verdict:** All 12 consistent ✓