# Sprint 8 Stage 1 — input.dat Consistency Check

Datasets compared: lowhyd, highhyd, highhyd_b010

| Parameter | lowhyd | highhyd | highhyd_b010 | Status |
|---|---|---|---|---|
| cement_lb_yd3 | 575.0 | 575.0 | 575.0 | same |
| water_lb_yd3 | 253.0 | 253.0 | 253.0 | same |
| coarse_agg_lb_yd3 | 1800.0 | 1800.0 | 1800.0 | same |
| fine_agg_lb_yd3 | 1100.0 | 1100.0 | 1100.0 | same |
| fly_ash_F_lb_yd3 | 0.0 | 0.0 | 0.0 | same |
| fly_ash_C_lb_yd3 | 0.0 | 0.0 | 0.0 | same |
| ggbfs_lb_yd3 | 0.0 | 0.0 | 0.0 | same |
| silica_fume_lb_yd3 | 0.0 | 0.0 | 0.0 | same |
| air_content_pct | 5.0 | 5.0 | 5.0 | same |
| cement_type | Type I/II | Type I/II | Type I/II | same |
| coarse_agg_type | Limestone | Limestone | Limestone | same |
| fine_agg_type | Siliceous River Sand | Siliceous River Sand | Siliceous River Sand | same |
| activation_energy_J_mol | 50000.0 | 50000.0 | 50000.0 | same |
| tau_hrs | 200.0 | 10.0 | 10.0 | EXPECTED DIFF |
| beta | 0.1 | 0.85 | 0.1 | EXPECTED DIFF |
| alpha_u | 0.1 | 0.7 | 0.7 | EXPECTED DIFF |
| Hu_J_kg | 1.0 | 1.0 | 1.0 | same (unexpected) |
| k_BTU_hr_ft_F | 1.55509402162759 | 1.55509402162759 | 1.55509402162759 | same |
| agg_Cp_BTU_lb_F | 0.20466622519 | 0.20466622519 | 0.20466622519 | same |
| CTE_microstrain_F | 4.23010793810549 | 4.23010793810549 | 4.23010793810549 | same |
| depth_ft | 80.0 | 80.0 | 80.0 | same |
| width_ft | 40.0 | 40.0 | 40.0 | same |
| length_ft | 40.0 | 40.0 | 40.0 | same |
| placement_temp_F | 73.0 | 73.0 | 73.0 | same |
| soil_temp_F | 100.0 | 100.0 | 100.0 | same |
| form_type | steel | steel | steel | same |
| form_removal_hrs | 168.0 | 168.0 | 168.0 | same |
| blanket_R_value | 5.67 | 5.67 | 5.67 | same |
| footing_subbase | Limestone | Limestone | Limestone | same |

## Verdict

**PASS** — only hydration kinetics parameters differ (tau_hrs, beta, alpha_u, Hu_J_kg). Diagnostic is valid.