# s1 — h_conv Magnitude Sweep

**Knob**: `h_forced_convection` output scaled by factor `s` (ablation
monkey-patch — context-manager-reverted after each run, not a production pattern).

**Side-channel scope**: patch replaces the module attribute and therefore
affects all 4 call sites in thermal_engine_2d (lines 496, 524, 1354, 1355),
including side form-face combined-h (line 524). Side-h is dominated by
R_FORM_CONTACT_SI in series; cluster CenterRMS authority is top-surface-driven.

Wind trace: v = 10.5 m/s (bit-identical across all Reference mixes).
Default h_eff = 20.30 W/(m²·K)  [ACI 305: `5.6 + 3.5·(0.4·v)`].

## Sweep table

CenterRMS = S0 gate metric, window [48, 168] hr, centerline mid-depth.
Authority threshold: cluster-mean CenterRMS variation ≥ 0.24°F.

| Mix | scale (h_eff) | D1 ΔT̄ (°F) | D2 lag (hr) | D3 ratio | CenterRMS (°F) | D4 TopRMS (°F) | D4 spread (°F) | S0 |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.50× (10.2 W/m²K) | +0.608 | +0.00 | 1.6233 | 0.8083 | 2.5661 | 2.2810 | **5/5** |
| MIX-11 | 0.50× (10.2 W/m²K) | +0.709 | +0.00 | 1.5443 | 0.8843 | 2.7451 | 2.4482 | **5/5** |
| MIX-12 | 0.50× (10.2 W/m²K) | +0.588 | +0.00 | 1.5789 | 0.7939 | 2.5347 | 2.2514 | **5/5** |
| MIX-03 | 0.50× (10.2 W/m²K) | +0.748 | +0.00 | 1.6110 | 0.8974 | 2.5912 | 2.2357 | 4/5 |
| MIX-01 | 0.75× (15.2 W/m²K) | +0.493 | +0.00 | 1.7977 | 0.7599 | 2.1636 | 1.9694 | **5/5** |
| MIX-11 | 0.75× (15.2 W/m²K) | +0.587 | +0.00 | 1.6974 | 0.8224 | 2.3056 | 1.9998 | **5/5** |
| MIX-12 | 0.75× (15.2 W/m²K) | +0.475 | +0.00 | 1.7523 | 0.7481 | 2.1393 | 1.9457 | **5/5** |
| MIX-03 | 0.75× (15.2 W/m²K) | +0.631 | +0.00 | 1.7678 | 0.8306 | 2.2123 | 1.9674 | **5/5** |
| MIX-01 | 1.00× (20.3 W/m²K) | +0.419 | +0.00 | 1.9102 | 0.7381 | 1.9338 | 1.6966 | **5/5** |
| MIX-11 | 1.00× (20.3 W/m²K) | +0.508 | +0.00 | 1.7961 | 0.7915 | 2.0449 | 1.8672 | **5/5** |
| MIX-12 | 1.00× (20.3 W/m²K) | +0.402 | +0.00 | 1.8640 | 0.7281 | 1.9154 | 1.6527 | **5/5** |
| MIX-03 | 1.00× (20.3 W/m²K) | +0.556 | +0.00 | 1.8689 | 0.7950 | 2.0022 | 1.7338 | **5/5** |
| MIX-01 | 1.25× (25.4 W/m²K) | +0.368 | +0.00 | 1.9884 | 0.7277 | 1.7937 | 1.4601 | **5/5** |
| MIX-11 | 1.25× (25.4 W/m²K) | +0.454 | +0.00 | 1.8647 | 0.7747 | 1.8790 | 1.6956 | **5/5** |
| MIX-12 | 1.25× (25.4 W/m²K) | +0.350 | +0.00 | 1.9417 | 0.7190 | 1.7802 | 1.4438 | **5/5** |
| MIX-03 | 1.25× (25.4 W/m²K) | +0.503 | +0.00 | 1.9391 | 0.7743 | 1.8780 | 1.5485 | **5/5** |
| MIX-01 | 1.50× (30.5 W/m²K) | +0.330 | +0.00 | 2.0457 | 0.7227 | 1.7040 | 1.3573 | **5/5** |
| MIX-11 | 1.50× (30.5 W/m²K) | +0.413 | +0.00 | 1.9150 | 0.7651 | 1.7678 | 1.5253 | **5/5** |
| MIX-12 | 1.50× (30.5 W/m²K) | +0.313 | +0.00 | 1.9986 | 0.7149 | 1.6945 | 1.3480 | **5/5** |
| MIX-03 | 1.50× (30.5 W/m²K) | +0.464 | +0.00 | 1.9906 | 0.7614 | 1.8008 | 1.4729 | **5/5** |
| MIX-01 | 2.00× (40.6 W/m²K) | +0.278 | +0.00 | 2.1240 | 0.7196 | 1.6030 | 1.2422 | **5/5** |
| MIX-11 | 2.00× (40.6 W/m²K) | +0.358 | +0.00 | 1.9837 | 0.7558 | 1.6340 | 1.2795 | **5/5** |
| MIX-12 | 2.00× (40.6 W/m²K) | +0.262 | +0.00 | 2.0765 | 0.7130 | 1.5995 | 1.2389 | **5/5** |
| MIX-03 | 2.00× (40.6 W/m²K) | +0.412 | +0.00 | 2.0609 | 0.7471 | 1.7177 | 1.3907 | **5/5** |

## Wait-point verdict

- Cluster-mean CenterRMS variation across sweep: **0.0994°F**
- Authority threshold ≥0.24°F: **NO**
- Cluster optimum: scale **2.00×** → h_eff = 40.6 W/(m²·K)  (cluster-mean CenterRMS = 0.7295°F)
- MIX-03 CenterRMS: baseline=0.7950°F → optimum=0.7471°F  (Δ=-0.0479°F)
- MIX-03 generalizes (Δ ≤ 0.10°F): **YES**

s2 (functional-form check): **SKIP** — no h_conv magnitude authority; form-check moot.
