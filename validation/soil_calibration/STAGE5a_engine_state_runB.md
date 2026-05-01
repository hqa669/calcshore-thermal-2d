# STAGE5a — Engine State: Run B at Multiple Time Slices

Run B: placement=73°F, soil=60°F, |ΔT|=13°F.
Grid: model_soil=False, is_submerged=True.
Grid dimensions: 13×21 concrete nodes, dx=0.3048 m = 1.000 ft.

**Sanity check vs Stage 4b:** masked max|R| at t=168 hr = 1.935°F (Stage 4b reference: 1.935°F, tol ±0.01°F) — PASS.

**Note on 'bulk' value:** mean over interior concrete cells only (1-cell erosion from all four concrete-domain edges).

---

## t = 1.0 hr (sample index 1)

Temperature range: 60.0–73.0°F

| Field | min | max | mean | std | bulk mean (interior) |
| --- | --- | --- | --- | --- | --- |
| α_hyd | 0.01677 | 0.01828 | 0.01809 | 0.00049 | **0.01827** |
| k (W/m·K) | 3.56333 | 3.56467 | 3.56349 | 0.00044 | **3.56333** |
| ρCp (J/m³·K) | 2392705.98499 | 2392884.74406 | 2392862.22332 | 58.02048 | **2392883.85058** |
| α_c (m²/hr) | 0.00536 | 0.00536 | 0.00536 | 0.00000 | **0.00536** |

## t = 24.0 hr (sample index 24)

Temperature range: 60.0–73.0°F

| Field | min | max | mean | std | bulk mean (interior) |
| --- | --- | --- | --- | --- | --- |
| α_hyd | 0.02718 | 0.02900 | 0.02874 | 0.00059 | **0.02895** |
| k (W/m·K) | 3.55381 | 3.55542 | 3.55404 | 0.00052 | **3.55385** |
| ρCp (J/m³·K) | 2391329.78959 | 2391690.64847 | 2391627.83172 | 122.41637 | **2391666.88081** |
| α_c (m²/hr) | 0.00535 | 0.00535 | 0.00535 | 0.00000 | **0.00535** |

## t = 84.0 hr (sample index 84)

Temperature range: 60.0–73.0°F

| Field | min | max | mean | std | bulk mean (interior) |
| --- | --- | --- | --- | --- | --- |
| α_hyd | 0.03169 | 0.03355 | 0.03322 | 0.00061 | **0.03342** |
| k (W/m·K) | 3.54977 | 3.55142 | 3.55006 | 0.00055 | **3.54988** |
| ρCp (J/m³·K) | 2390734.90618 | 2391183.84799 | 2391077.86019 | 159.23821 | **2391119.65690** |
| α_c (m²/hr) | 0.00534 | 0.00535 | 0.00534 | 0.00000 | **0.00534** |

## t = 168.0 hr (sample index 168)

Temperature range: 60.0–73.0°F

| Field | min | max | mean | std | bulk mean (interior) |
| --- | --- | --- | --- | --- | --- |
| α_hyd | 0.03422 | 0.03609 | 0.03571 | 0.00062 | **0.03590** |
| k (W/m·K) | 3.54751 | 3.54917 | 3.54785 | 0.00055 | **3.54768** |
| ρCp (J/m³·K) | 2390400.10464 | 2390900.26309 | 2390753.63252 | 180.03540 | **2390794.07671** |
| α_c (m²/hr) | 0.00534 | 0.00535 | 0.00534 | 0.00000 | **0.00534** |

---

**Engine's actual α_hyd at t=168 hr (interior bulk): 0.0359 ± 0.0003**
α_hyd = 0.036 << α_u = 0.100 — hydration strongly suppressed (36% of α_u; kinetics effectively frozen by Hu_J_kg≈1).

**Engine's actual α_c at t=168 hr (interior bulk): 0.00534 m²/hr**
CW M3 reference: ~0.00528 m²/hr. Ratio (engine/CW): 1.012.
