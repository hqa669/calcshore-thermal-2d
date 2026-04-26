# s1 — R_form Reference-Set Sweep

**Knob**: `R_FORM_CONTACT_SI` (module-level constant, `thermal_engine_2d.py:127`).
Monkey-patch via context manager; all 14 in-engine sites perturbed uniformly
(code sites: lines 524, 1913, 1914, 1935, 1944, 1949, 2020, 2192, 2193, 2212, 2221).

**Sweep values**: [0.05, 0.06, 0.07, 0.0862, 0.1]
**Evaluation set**: ['MIX-01', 'MIX-02', 'MIX-03', 'MIX-11', 'MIX-12']
MIX-02 reported alongside cluster mixes; kinetics-anomalous per `mix02_recon.md`.

## Sweep table

Primary gate: **CornerRMS** (PR 15 finding; S0 threshold 3.0°F).
Secondary: **PeakGrad Δ** (also moves with R_form; S0 threshold 2.0°F).
Sanity check: **CenterRMS** (should be form-face-decoupled; S0 threshold 1.0°F).

| Mix | R_form (m²K/W) | PeakGrad Δ (°F) | CenterRMS (°F) | CornerRMS (°F) | S0 |
| --- | ---: | ---: | ---: | ---: | :---: |
| MIX-01 | 0.0500 | +2.324 | 0.7381 | 1.3335 | 4/5 |
| MIX-02 | 0.0500 | +3.346 | 0.5748 | 2.3631 | 4/5 |
| MIX-03 | 0.0500 | +1.881 | 0.7950 | 1.1338 | **5/5** |
| MIX-11 | 0.0500 | +2.627 | 0.7915 | 1.3989 | 4/5 |
| MIX-12 | 0.0500 | +2.300 | 0.7281 | 1.3287 | 4/5 |
| MIX-01 | 0.0600 | +1.534 | 0.7381 | 1.0034 | **5/5** |
| MIX-02 | 0.0600 | +2.679 | 0.5748 | 1.8113 | 4/5 |
| MIX-03 | 0.0600 | +1.032 | 0.7950 | 1.0923 | **5/5** |
| MIX-11 | 0.0600 | +1.854 | 0.7915 | 1.0244 | **5/5** |
| MIX-12 | 0.0600 | +1.507 | 0.7281 | 1.0016 | **5/5** |
| MIX-01 | 0.0700 | +0.799 | 0.7381 | 1.2723 | **5/5** |
| MIX-02 | 0.0700 | +2.064 | 0.5748 | 1.4053 | 4/5 |
| MIX-03 | 0.0700 | +0.244 | 0.7950 | 1.6160 | **5/5** |
| MIX-11 | 0.0700 | +1.136 | 0.7915 | 1.2233 | **5/5** |
| MIX-12 | 0.0700 | +0.771 | 0.7281 | 1.2750 | **5/5** |
| MIX-01 | 0.0862 | -0.288 | 0.7381 | 2.2192 | **5/5** |
| MIX-02 | 0.0862 | +1.160 | 0.5748 | 1.2362 | **5/5** |
| MIX-03 | 0.0862 | -0.925 | 0.7950 | 2.7238 | **5/5** |
| MIX-11 | 0.0862 | +0.073 | 0.7915 | 2.1198 | **5/5** |
| MIX-12 | 0.0862 | -0.318 | 0.7281 | 2.2251 | **5/5** |
| MIX-01 | 0.1000 | -1.129 | 0.7381 | 3.0777 | 4/5 |
| MIX-02 | 0.1000 | +0.603 | 0.5748 | 1.5493 | **5/5** |
| MIX-03 | 0.1000 | -1.830 | 0.7950 | 3.6569 | 4/5 |
| MIX-11 | 0.1000 | -0.750 | 0.7915 | 2.9583 | **5/5** |
| MIX-12 | 0.1000 | -1.161 | 0.7281 | 3.0847 | 4/5 |

## Commit-decision criteria (Item 1)

**C1** — all 5 Reference mixes S0 5/5 at R_form=0.060: **FAIL ✗**
  - MIX-01: S0 5/5 ✓
  - MIX-02: S0 FAIL ✗
  - MIX-03: S0 5/5 ✓
  - MIX-11: S0 5/5 ✓
  - MIX-12: S0 5/5 ✓

**C2** — ≥3/5 mixes improve CornerRMS at 0.060 vs 0.0862: **PASS ✓** (4/5 improve)
  - MIX-01: ✓ improves
  - MIX-02: ✗ no improvement
  - MIX-03: ✓ improves
  - MIX-11: ✓ improves
  - MIX-12: ✓ improves

**C3** — no mix CenterRMS degrades >0.05°F at 0.060 vs 0.0862: **PASS ✓**
  - MIX-01: Δ=-0.0000°F
  - MIX-02: Δ=-0.0000°F
  - MIX-03: Δ=-0.0000°F
  - MIX-11: Δ=-0.0000°F
  - MIX-12: Δ=-0.0000°F

**Outcome**: **DO NOT COMMIT** — C1 failed. Route to engine v3 release notes.

