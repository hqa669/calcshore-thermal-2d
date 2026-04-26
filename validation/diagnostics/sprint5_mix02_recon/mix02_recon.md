# MIX-02 within-Reference field diff

Pre-PR-19 mini-spike. Compares CWMixDesign / CWGeometry / CWConstruction / CWEnvironment fields across the 5 Reference mixes (MIX-01, 02, 03, 11, 12) to characterize MIX-02's outlier behavior in PR 18's recon.

## Disposition: C

**Disposition C — hydration regression divergence.** hydration regression params vary across mixes: ['Hu_J_kg', 'activation_energy_J_mol', 'alpha_u', 'beta', 'tau_hrs']; stress-only fields also vary (['CTE_microstrain_F']) but do not affect the heat equation. Kinetics sub-cluster analysis: cluster {MIX-01, MIX-11, MIX-12} shares Ea=26457.9, tau=29.4010 hr, beta=0.89500, Hu=424143.1 J/kg (alpha_u varies within cluster — see kinetics table). Outlier mixes with divergent (Ea, tau, beta, Hu): {MIX-02, MIX-03}. Note: MIX-03 also has divergent kinetics but was NOT flagged as a diagnostic outlier in PR 18 (D1/D3/D4 within the cluster signature); PR 19 should target cluster {MIX-01, MIX-11, MIX-12} for the BC ablation and report MIX-02 and MIX-03 separately. MIX-02's outlier behavior in PR 18 (D1 negative mean, D3 inverse amplitude, D4 bottom-concentrated profile) is hydration-kinetics-driven, not boundary-physics-driven. MIX-02 routes to the hydration sprint series. R9 implication: differing CW-computed kinetics on mixes with identical parsed composition (MIX-01/02/03 share composition exactly) indicates CW's regression consumes inputs not captured in our parsed CWMixDesign fields — likely cement Bogue compounds or admixture chemistry.

## Construction fields that vary

_No scalar fields vary in `construction` across the 5 Reference mixes._

## Geometry fields that vary

_No scalar fields vary in `geometry` across the 5 Reference mixes._

## Environment scalar fields that vary

_No scalar fields vary in `environment` across the 5 Reference mixes._

## Environment ndarray fields that vary (max abs diff vs MIX-01)

_All hourly weather arrays bit-identical across the 5 Reference mixes._

## Environment list fields that vary

_All daily-summary lists bit-identical across the 5 Reference mixes._

## Mix design — composition (lb/yd³)

| Mix | cement | fly_ash_F | fly_ash_C | ggbfs | silica_fume | water | total | wcm |
|---|---|---|---|---|---|---|---|---|
| MIX-01 | 350.0 | 125.0 | 0.0 | 100.0 | 0.0 | 253.0 | 575.0 | 0.4400 |
| MIX-02 | 350.0 | 125.0 | 0.0 | 100.0 | 0.0 | 253.0 | 575.0 | 0.4400 |
| MIX-03 | 350.0 | 125.0 | 0.0 | 100.0 | 0.0 | 253.0 | 575.0 | 0.4400 |
| MIX-11 | 350.0 | 125.0 | 0.0 | 100.0 | 0.0 | 218.0 | 575.0 | 0.3791 |
| MIX-12 | 350.0 | 125.0 | 0.0 | 100.0 | 0.0 | 259.0 | 575.0 | 0.4504 |

## Mix design — hydration regression parameters

| Mix | Ea_J_mol | tau_hrs | beta | alpha_u | Hu_J_kg |
|---|---|---|---|---|---|
| MIX-01 | 26457.9 | 29.4010 | 0.89500 | 0.75852 | 424143.1 |
| MIX-02 | 25797.2 | 34.1580 | 0.81800 | 0.74452 | 395949.6 |
| MIX-03 | 26036.3 | 33.4520 | 0.99000 | 0.76652 | 432957.0 |
| MIX-11 | 26457.9 | 29.4010 | 0.89500 | 0.72554 | 424143.1 |
| MIX-12 | 26457.9 | 29.4010 | 0.89500 | 0.76342 | 424143.1 |

## Mix design — derived thermal properties

| Mix | k_BTU_hr_ft_F | agg_Cp_BTU_lb_F | CTE_microstrain_F |
|---|---|---|---|
| MIX-01 | 1.5551 | 0.2047 | 4.2476 |
| MIX-02 | 1.5551 | 0.2047 | 4.2476 |
| MIX-03 | 1.5551 | 0.2047 | 4.2476 |
| MIX-11 | 1.5551 | 0.2047 | 4.2070 |
| MIX-12 | 1.5551 | 0.2047 | 4.2544 |
