# Sprint 8 Stage 2-prep v2 — §D.4 Direct α(t) Comparison Table

## Method
CW α(effective) inverted from CW di=24 wi=0 bulk-shift signal via:
  CW α(eff) = ΔT_bulk[d, t] / C(t)
where C(t) is calibrated by fitting the T_ref that makes ΔT_bulk/α constant across
all 3 datasets at each timestep (per-timestep best-fit T_ref from §D.2).

## Per-timestep best-fit T_ref

| t (hr) | Best T_ref (°C) | CV% |
|---|---|---|
| 24 | 20.0 | 3.29% |
| 48 | 22.0 | 0.28% |
| 84 | 19.5 | 2.82% |
| 120 | 23.0 | 4.58% |
| 168 | 19.0 | 0.93% |

## Headline Table: Engine α vs CW α (effective)

| Dataset | t (hr) | Engine α (T_ref=23°C) | Engine α (T_ref=21°C) | CW α (eff) | Δ at 23°C | Δ at 21°C | Δ at 23°C % | Δ at 21°C % |
|---|---|---|---|---|---|---|---|---|
| lowhyd | 24 | 0.0290 | 0.0295 | 0.0297 | +0.0007 | +0.0003 | +2.6% | +0.8% |
| lowhyd | 48 | 0.0315 | 0.0320 | 0.0318 | +0.0003 | -0.0003 | +0.8% | -0.8% |
| lowhyd | 84 | 0.0335 | 0.0341 | 0.0344 | +0.0009 | +0.0004 | +2.7% | +1.1% |
| lowhyd | 120 | 0.0349 | 0.0354 | 0.0349 | +0.0000 | -0.0005 | +0.0% | -1.4% |
| lowhyd | 168 | 0.0361 | 0.0366 | 0.0371 | +0.0010 | +0.0005 | +2.8% | +1.4% |
| highhyd | 24 | 0.4326 | 0.4563 | 0.4676 | +0.0351 | +0.0114 | +8.1% | +2.5% |
| highhyd | 48 | 0.5360 | 0.5520 | 0.5441 | +0.0082 | -0.0079 | +1.5% | -1.4% |
| highhyd | 84 | 0.5930 | 0.6040 | 0.6116 | +0.0186 | +0.0076 | +3.1% | +1.3% |
| highhyd | 120 | 0.6193 | 0.6277 | 0.6193 | +0.0000 | -0.0085 | +0.0% | -1.4% |
| highhyd | 168 | 0.6384 | 0.6450 | 0.6509 | +0.0125 | +0.0060 | +2.0% | +0.9% |
| highhyd_b010 | 24 | 0.2796 | 0.2832 | 0.2850 | +0.0053 | +0.0018 | +1.9% | +0.6% |
| highhyd_b010 | 48 | 0.2974 | 0.3009 | 0.2991 | +0.0018 | -0.0018 | +0.6% | -0.6% |
| highhyd_b010 | 84 | 0.3115 | 0.3150 | 0.3177 | +0.0061 | +0.0026 | +2.0% | +0.8% |
| highhyd_b010 | 120 | 0.3205 | 0.3240 | 0.3205 | +0.0000 | -0.0035 | +0.0% | -1.1% |
| highhyd_b010 | 168 | 0.3289 | 0.3323 | 0.3358 | +0.0069 | +0.0035 | +2.1% | +1.0% |

## Summary

| Metric | Δ at T_ref=23°C | Δ at T_ref=21°C |
|---|---|---|
| max|Δ| (%) | 8.1% | 2.5% |
| mean|Δ| (%) | 2.0% | 1.1% |
| RMS|Δ| (%) | 2.8% | 1.2% |

_Δ is signed: (CW α − Engine α). Positive = CW α > Engine α (CW hydrates faster)._