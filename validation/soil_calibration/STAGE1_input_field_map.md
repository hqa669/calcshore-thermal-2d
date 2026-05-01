# STAGE1 — input.dat Field Map and Consistency Check

## Field-to-Parameter Mapping

All indices are **0-based** (Python `lines[idx]`). File line = idx + 1.

| Field | idx | 1-indexed line | Description |
| ----- | --- | -------------- | ----------- |
| `length_units` | 1 | 2 | ft |
| `mass_units` | 2 | 3 | lb/yd³ |
| `temp_units` | 4 | 5 | °F |
| `analysis_duration` | 8 | 9 | 7 days |
| `placement_time_str` | 10 | 11 | 5 am |
| `project_location` | 11 | 12 | TX, Austin |
| `placement_date` | 12 | 13 | 2026/7/15 |
| `member_width_ft` | 21 | 22 | 40 |
| `member_depth_ft` | 22 | 23 | 80 |
| `member_length_ft` | 23 | 24 | 40 |
| `cement_lb_yd3` | 54 | 55 | 575 |
| `water_lb_yd3` | 55 | 56 | 253 |
| `coarse_agg_lb_yd3` | 56 | 57 | 1800 |
| `fine_agg_lb_yd3` | 57 | 58 | 1100 |
| `air_content_pct` | 58 | 59 | 5 |
| `silica_fume_lb_yd3` | 59 | 60 | 0 |
| `class_c_fly_ash_lb_yd3` | 60 | 61 | 0 |
| `class_f_fly_ash_lb_yd3` | 61 | 62 | 0 |
| `fly_ash_CaO_pct` | 62 | 63 | 0 |
| `ggbfs_lb_yd3` | 63 | 64 | 0 |
| `cement_type` | 361 | 362 | Type I/II |
| `activation_energy_J_mol` | 385 | 386 | 50000 |
| `tau_hrs` | 386 | 387 | 200 |
| `beta` | 387 | 388 | 0.1 |
| `alpha_u` | 388 | 389 | 0.1 |
| `Hu_J_kg` | 389 | 390 | 1 |
| `coarse_agg_type` | 390 | 391 | Limestone |
| `fine_agg_type` | 393 | 394 | Siliceous River Sand |
| `CTE_microstrain_F` | 403 | 404 | 4.23010793810549 |
| `thermal_conductivity_BTU_hr_ft_F` | 404 | 405 | 1.55509402162759 |
| `aggregate_Cp_BTU_lb_F` | 405 | 406 | 0.20466622519 |
| `maturity_a` | 407 | 408 | -4000 |
| `maturity_b` | 408 | 409 | 2800 |
| `form_color` | 431 | 432 | Red |
| `form_type` | 432 | 433 | Steel |
| `form_removal_hrs` | 433 | 434 | 168 |
| `blanket_R_value` | 434 | 435 | 5.67 |
| `placement_temp_F` | 438 | 439 | 73 |
| `soil_temp_F` | 439 | 440 | 73 |
| `delay_strip_to_cure_hrs` | 445 | 446 | 1 |
| `footing_subbase` | 446 | 447 | Limestone |
| `steel_cover_in` | 447 | 448 | 2 |
| `top_cure_blanket_time_hrs` | 478 | 479 | 2 |
| `Dref_x1e13_m2_s` | 491 | 492 | 73 |
| `chloride_aging_m` | 492 | 493 | 73 |
| `steel_type` | 495 | 496 | 73 |
| `exposure_class` | 508 | 509 | 1 |

### Extra fields NOT in CW_DAT_INDEX but present in input.dat

| idx | 1-indexed line | Content (from runA_baseline) | Notes |
| --- | -------------- | ----------------------------- | ----- |
| 465 | 466 | `Clay` | Side soil type (parsed by CW but not by loader) |
| 466 | 467 | `Clay` | Bottom soil type (parsed by CW but not by loader) |
| 486–494 | 487–495 | soil temp profile (9 values) | Per-depth soil temperatures; mirrors soil_temp_F |

**Known gap**: `footing_subbase` (idx 446) reads `Limestone` — a CW UI default that doesn't
reflect the actual soil used by CW's model. The true soil type (`Clay`) lives at idx 465/466,
which the loader does NOT parse. The engine ignores `footing_subbase` entirely.

---

## Consistency Table

| Label | Folder | Placement °F | Soil °F | ΔT | Geometry | Hydration | Form | Soil strings (idx465/466) | Notes |
| ----- | ------ | ------------ | ------- | -- | -------- | --------- | ---- | ------------------------- | ----- |
| A | runA_baseline | 73 | 73 | 0 | OK | OK | OK | Clay/Clay | — |
| B | runB_73_60 | 73 | 60 | -13 | OK | OK | OK | Clay/Clay | — |
| C | runC_73_90 | 73 | 90 | 17 | OK | OK | OK | Clay/Clay | — |
| D | runD_60_73 | 60 | 73 | 13 | OK | OK | OK | Clay/Clay | — |
| E | runE_90_73 | 90 | 73 | -17 | OK | OK | OK | Clay/Clay | — |
| F | runF_73_45 | 73 | 45 | -28 | OK | OK | OK | Clay/Clay | — |
| G | runG_73_100 | 73 | 100 | 27 | OK | OK | OK | Clay/Clay | — |
| H | runH_45_73 | 45 | 73 | 28 | OK | OK | OK | Clay/Clay | — |
| I | runI_100_73 | 100 | 73 | -27 | OK | OK | OK | Clay/Clay | — |

## Summary

All 9 runs confirm: geometry=40×80×40 ft, Hu_J_kg=1 (suppressed), Ea=50000, τ=200, β=0.1, αu=0.1, form=Steel, soil strings=Clay/Clay. Placement and soil temperatures match labels.

## Mirror Pair Note

| Pair | Run+ | Run- | |ΔT| | Clean mirror? |
| ---- | ---- | ---- | ----- | ------------- |
| B↔D | B (73/60) | D (60/73) | 13 | Yes |
| C↔E | C (73/90) | E (90/73) | 17 | Yes |
| F↔H | F (73/45) | H (45/73) | 28 | Yes |
| G↔I | G (73/100) | I (100/73) | 27 | Yes — note: brief's IMPORTANT warning about Run I having soil=60 is incorrect; actual data shows soil=73 |
