# PR 18 — Centerline RMS Reconnaissance Summary

**Sprint**: 5  **Source commit**: sprint-4-complete (5596c92)
**Reference mixes**: MIX-01, MIX-02, MIX-03, MIX-11, MIX-12
**Gate window**: t ∈ [48, 168] hr  **Mode-A window**: t ∈ [8, 48] hr

## Per-mix wall times

| Mix | wall_s |
|---|---|
| MIX-01 | 1.8 |
| MIX-02 | 1.7 |
| MIX-03 | 1.7 |
| MIX-11 | 1.8 |
| MIX-12 | 1.7 |

## Consolidated diagnostic table

| Mix | D1 mean ΔT | D1 std ΔT | D2 lag_hr | D3 swing_ratio | D4 mid_rms | D4 spread | D5 mean ΔT | D5 t_at_max | D5 match |
|---|---|---|---|---|---|---|---|---|---|
| MIX-01 | +0.419 | 0.608 | +0.00 | 1.9102 | 0.7867 | 1.6966 | +1.1503 | 41.0 | N |
| MIX-02 | -0.136 | 0.559 | +0.00 | 0.6318 | 0.6022 | 2.2391 | +0.7676 | 39.0 | N |
| MIX-03 | +0.556 | 0.569 | +0.00 | 1.8689 | 0.8576 | 1.7338 | +0.9243 | 46.9 | N |
| MIX-11 | +0.508 | 0.607 | +0.00 | 1.7961 | 0.8391 | 1.8672 | +1.2106 | 42.0 | N |
| MIX-12 | +0.402 | 0.608 | +0.00 | 1.8640 | 0.7767 | 1.6527 | +1.1385 | 41.0 | N |

## Routing analysis

### Per-diagnostic signature summary

- **D1 constant offset**: positive bias on all mixes = False; near-constant (std < 0.5×|mean|) = False
- **D2 phase lag**: any |lag| > 0.5 hr = False; consistent sign = True
- **D3 amplitude**: swing ratio departs >5% from 1.0 on all mixes = True
- **D4 depth profile**: spread > 0.2°F on any mix = True
- **D5 Mode-A**: 0/5 mixes match (need ≥+1.0°F mean + t_at_max ∈ [12,24] hr)

### Routing paragraph

D4 is the dominant finding: residual is strongly depth-localized (spread 1.65–2.24°F) with peak residual at the top concrete surface (d=0: 1.9–2.0°F for MIX-01/03/11/12) and a bottom-surface-concentrated profile for MIX-02 (monotonically rising to 2.60°F at d=12), pointing to surface boundary-condition physics as the primary driver rather than internal thermal properties. D3 reinforces this: engine diurnal amplitude is ~90% larger than CW at [144,168] hr on 4/5 mixes (MIX-02 inverted), confirming a BC amplitude mismatch rather than a bulk thermal property offset. D5 formal Mode-A criteria fail 0/5 (peak time [12,24] hr was calibrated to cold-placed MIX-15; Reference 60°F mixes peak at 39–47 hr), but positive early-hydration bias is present (+0.77 to +1.21°F mean ΔT in [8,48] hr), consistent with R9 Mode-A at warm placement and contributing a DC offset component visible in D1 (4/5 mixes: mean ΔT +0.40 to +0.56°F in gate window). Recommended path: PR 19 fix-conditional — investigate top-surface BC amplitude response (convection coefficient, solar/LW diurnal forcing) and bottom-surface BC (form-face R_form, soil Barber lag/damping) as the primary target; note that the DC component traced to shifted Mode-A (R9) may persist as a residual even after BC correction.

### Key anomaly: MIX-02 diverges from cluster

MIX-02 differs from MIX-01/03/11/12 on all three discriminating diagnostics:
- D1: mean ΔT = −0.14°F (only mix with negative mean; others +0.40–0.56°F)
- D3: swing ratio 0.63 (engine UNDER-oscillates; others ~1.81–1.91)
- D4: bottom-concentrated profile (others are top-concentrated)

MIX-02 has a different SCM blend (different fly_ash_C / ggbfs fractions despite same total SCM%) than the MIX-01 cluster. The inverted amplitude and opposite mean-ΔT suggest MIX-02's hydration heat curve differs in timing or shape in a way that produces different BC-residual interaction. PR 19 should track whether MIX-02 continues to behave as an outlier after BC recalibration, or whether the cluster aligns once BC amplitude is corrected.
