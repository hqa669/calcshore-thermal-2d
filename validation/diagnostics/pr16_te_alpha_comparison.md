# PR 16 — Phase 2.6: Equivalent-Age and Hydration Trajectory Comparison

Compares engine vs CW for MIX-01 (warm IC, 60°F) and MIX-15 (cold IC, 45°F).
CW-implied te/α derived by integrating engine's Arrhenius model over the CW
centerline temperature trajectory (mid-depth, centerline column).

If the CW and engine diverge in the same place for both mixes: baseline difference.
If only MIX-15 diverges: cold-IC-specific difference.
If neither diverges: cold-IC failure is downstream of equivalent age.

## Mix Parameters (identical for MIX-01 and MIX-15 per V2 finding)

| Parameter | Value |
| --- | --- |
| Ea (J/mol) | 26457.95 |
| tau (hr) | 29.4010 |
| beta | 0.8950 |
| alpha_u | 0.7585 |
| T_REF_K | 296.15 K (23.0°C = 73.4°F) |

## Divergence Summary

| Mix | First |ΔT| > 1°F (hr) | Peak ΔT (°F) | At t (hr) |
| --- | --- | --- | --- |
| MIX-01 (warm IC) | t = 24 | +1.70 | 42 |
| MIX-15 (cold IC) | t = 96 | -3.04 | 168 |

### Finding 3

**Finding 3: Both MIX-01 and MIX-15 show ΔT > 1°F vs CW.**
There is a baseline implementation difference independent of cold IC.
This is a Sprint 5 priority — affects the whole library.

## MIX-01 — Temperature / te / α Comparison (warm IC, 60°F)

CW reference: MIX-01 CW. Engine comparison: engine(MIX-01, 60°F IC).
ΔT = Engine − CW (positive = engine hotter). CW-implied te/α uses engine's Arrhenius model.

| t (hr) | Eng T (°F) | CW T (°F) | ΔT (°F) | Eng te (hr) | CW te (hr) | Δte (hr) | Eng α | CW α | Δα | Eng dα/dt (×10⁻³/hr) | CW dα/dt (×10⁻³/hr) |
|---|---|---|---|---|---|---|---|---|---|---|---|
|    0 |   60.00 |   60.01 |  -0.01 |    0.010 |    0.010 |  +0.000 | 0.0000 | 0.0000 | +0.0000 | +0.0934 | +0.0809 |
|    4 |   60.04 |   60.04 |  -0.00 |    3.042 |    2.979 |  +0.063 | 0.0004 | 0.0003 | +0.0000 | +1.6022 | +1.5413 |
|    8 |   61.47 |   61.41 |  +0.06 |    6.109 |    6.045 |  +0.064 | 0.0128 | 0.0123 | +0.0005 | +5.7968 | +5.6883 |
|   12 |   65.35 |   65.10 |  +0.25 |    9.356 |    9.282 |  +0.074 | 0.0467 | 0.0458 | +0.0009 | +10.1880 | +10.0764 |
|   16 |   70.80 |   70.27 |  +0.53 |   12.939 |   12.838 |  +0.101 | 0.0943 | 0.0929 | +0.0014 | +12.6435 | +12.5097 |
|   20 |   76.94 |   76.12 |  +0.82 |   16.976 |   16.820 |  +0.156 | 0.1479 | 0.1459 | +0.0020 | +13.5205 | +13.3550 |
|   24 |   83.21 |   82.09 |  +1.11 |   21.547 |   21.305 |  +0.242 | 0.2025 | 0.1998 | +0.0027 | +13.1997 | +13.0253 |
|   30 |   92.07 |   90.61 |  +1.47 |   29.501 |   29.064 |  +0.438 | 0.2799 | 0.2762 | +0.0037 | +12.0978 | +11.9525 |
|   36 |   99.77 |   98.11 |  +1.65 |   38.800 |   38.097 |  +0.703 | 0.3477 | 0.3432 | +0.0044 | +10.3813 | +10.2952 |
|   42 |  106.11 |  104.41 |  +1.70 |   49.373 |   48.354 |  +1.019 | 0.4045 | 0.3997 | +0.0048 | +8.6170 | +8.5913 |
|   48 |  111.19 |  109.54 |  +1.64 |   61.080 |   59.718 |  +1.362 | 0.4511 | 0.4463 | +0.0048 | +6.4106 | +6.4422 |
|   60 |  118.31 |  116.92 |  +1.39 |   87.233 |   85.188 |  +2.045 | 0.5199 | 0.5157 | +0.0042 | +4.7910 | +4.8460 |
|   72 |  122.64 |  121.57 |  +1.07 |  116.025 |  113.386 |  +2.639 | 0.5660 | 0.5626 | +0.0034 | +2.8192 | +2.8755 |
|   96 |  126.82 |  126.32 |  +0.50 |  178.076 |  174.649 |  +3.427 | 0.6213 | 0.6192 | +0.0022 | +1.7935 | +1.8365 |
|  120 |  128.20 |  128.10 |  +0.10 |  242.840 |  239.111 |  +3.729 | 0.6521 | 0.6508 | +0.0014 | +1.0423 | +1.0692 |
|  144 |  128.39 |  128.57 |  -0.18 |  308.363 |  304.688 |  +3.675 | 0.6714 | 0.6705 | +0.0009 | +0.6735 | +0.6902 |
|  168 |  128.01 |  128.37 |  -0.36 |  373.745 |  370.373 |  +3.372 | 0.6845 | 0.6839 | +0.0006 | +0.5450 | +0.5579 |

## MIX-15 — Temperature / te / α Comparison (cold IC, 45°F)

CW reference: MIX-15 CW. Engine comparison: engine(MIX-15, 45°F IC).
ΔT = Engine − CW (positive = engine hotter). CW-implied te/α uses engine's Arrhenius model.

| t (hr) | Eng T (°F) | CW T (°F) | ΔT (°F) | Eng te (hr) | CW te (hr) | Δte (hr) | Eng α | CW α | Δα | Eng dα/dt (×10⁻³/hr) | CW dα/dt (×10⁻³/hr) |
|---|---|---|---|---|---|---|---|---|---|---|---|
|    0 |   45.00 |   45.00 |  +0.00 |    0.010 |    0.010 |  +0.000 | 0.0000 | 0.0000 | +0.0000 | +0.0070 | +0.0058 |
|    4 |   45.00 |   45.00 |  +0.01 |    2.195 |    2.149 |  +0.046 | 0.0000 | 0.0000 | +0.0000 | +0.3912 | +0.3714 |
|    8 |   45.36 |   45.34 |  +0.02 |    4.385 |    4.339 |  +0.046 | 0.0031 | 0.0030 | +0.0002 | +2.1242 | +2.0725 |
|   12 |   46.95 |   46.87 |  +0.08 |    6.622 |    6.573 |  +0.048 | 0.0170 | 0.0166 | +0.0004 | +4.8645 | +4.7999 |
|   16 |   49.82 |   49.60 |  +0.21 |    8.973 |    8.918 |  +0.056 | 0.0420 | 0.0414 | +0.0007 | +7.2310 | +7.1630 |
|   20 |   53.59 |   53.22 |  +0.36 |   11.506 |   11.434 |  +0.071 | 0.0749 | 0.0739 | +0.0010 | +8.7924 | +8.7145 |
|   24 |   57.91 |   57.42 |  +0.49 |   14.272 |   14.176 |  +0.096 | 0.1124 | 0.1111 | +0.0013 | +9.7585 | +9.6738 |
|   30 |   64.83 |   64.24 |  +0.59 |   18.955 |   18.805 |  +0.150 | 0.1725 | 0.1706 | +0.0018 | +9.9987 | +9.9266 |
|   36 |   71.73 |   71.17 |  +0.56 |   24.367 |   24.153 |  +0.214 | 0.2324 | 0.2302 | +0.0022 | +9.6908 | +9.6549 |
|   42 |   78.17 |   77.77 |  +0.40 |   30.566 |   30.293 |  +0.273 | 0.2887 | 0.2865 | +0.0022 | +8.9493 | +8.9600 |
|   48 |   83.94 |   83.77 |  +0.17 |   37.562 |   37.250 |  +0.312 | 0.3398 | 0.3377 | +0.0020 | +7.5005 | +7.5628 |
|   60 |   93.17 |   93.58 |  -0.40 |   53.806 |   53.535 |  +0.271 | 0.4238 | 0.4226 | +0.0011 | +6.0911 | +6.1731 |
|   72 |   99.65 |  100.62 |  -0.96 |   72.605 |   72.579 |  +0.026 | 0.4859 | 0.4859 | +0.0001 | +3.9340 | +4.0033 |
|   96 |  107.04 |  108.82 |  -1.79 |  115.502 |  116.581 |  -1.079 | 0.5654 | 0.5668 | -0.0014 | +2.6020 | +2.6458 |
|  120 |  110.38 |  112.71 |  -2.33 |  162.479 |  165.318 |  -2.838 | 0.6108 | 0.6129 | -0.0020 | +1.5387 | +1.5578 |
|  144 |  111.82 |  114.55 |  -2.73 |  211.378 |  216.443 |  -5.065 | 0.6392 | 0.6415 | -0.0023 | +0.9896 | +0.9967 |
|  168 |  112.30 |  115.34 |  -3.04 |  261.061 |  268.695 |  -7.634 | 0.6583 | 0.6607 | -0.0024 | +0.7963 | +0.7997 |

## Interpretation Notes

### te divergence interpretation
Δte > 0 at a given clock time means the engine accumulates equivalent age faster than
CW-implied te (i.e., the engine is hotter than CW at that time, driving more Arrhenius
factor). Δte < 0 means the engine is lagging — colder → slower hydration → lower te.

### Arrhenius asymmetry at cold IC
At 45°F (280.4 K), af = 0.5462.
At 60°F (288.7 K), af = 0.7581.
At 73.4°F (296.15 K, T_REF), af = 1.0000.
Cold concrete starts at 72% of the warm hydration rate. If CW uses a different
T_REF or Ea, this ratio changes and te accumulates at a different pace.

### Limitation of this diagnostic
CW-implied te is computed with the ENGINE's Arrhenius model (same Ea, T_REF_K=296.15 K).
If CW uses the same model, the CW-implied te is directly comparable to the engine's te.
If CW uses a different T_REF or a different Arrhenius form, the comparison shows the
divergence attributable to model difference (not just temperature trajectory).
