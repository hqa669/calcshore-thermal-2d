# PR 19 — Top-Surface BC Ablation Summary

**Sprint**: 5  **Evaluation set**: MIX-01, MIX-11, MIX-12 (cluster)  +  MIX-03 (control)
**Total wall time**: 114s
**Baseline residual** (from PR 18): D3 ratio ~1.90×, D4 top-concentrated,
cluster CenterRMS ≈ 0.79–0.84°F (target: <0.50°F S1-aspire).

## Knob 1 — h_conv magnitude

Default: h_eff = 20.30 W/(m²·K) at v=10.5 m/s  [ACI 305 linear-in-wind]
Sweep: 0.5× – 2.0× (10.2 – 40.6 W/m²K)

- Cluster-mean CenterRMS variation: 0.0994°F
- **Authority (≥0.24°F): NO**
- Cluster optimum: scale=2.00×  h_eff=40.6 W/(m²·K)  cluster-mean CenterRMS=0.7295°F
- MIX-03 at optimum: baseline=0.7950°F → 0.7471°F  Δ=-0.0479°F  generalizes=YES

## Knob 1b — h_conv functional-form check

Status: **SKIPPED** — see s2_h_conv_form.md for reason.

## Knob 2 — α_top (solar absorptivity)

Default: α_top = 0.65  Sweep: 0.30 – 0.85

- Cluster-mean CenterRMS variation: 0.0809°F
- **Authority (≥0.24°F): NO**
- Cluster optimum: α_top=0.45  cluster-mean CenterRMS=0.7267°F
- MIX-03 Δ at optimum: -0.0515°F  generalizes=YES

## Knob 3 — ε_top (LW emissivity)

Default: ε_top = 0.88  Sweep: 0.50 – 0.95

- Cluster-mean CenterRMS variation: 0.0399°F
- **Authority (≥0.24°F): NO**
- Cluster optimum: ε_top=0.95  cluster-mean CenterRMS=0.7482°F
- MIX-03 Δ at optimum: -0.0072°F  generalizes=YES

## Routing for PR 20

**PR 20 shape: route to engine v3 release notes as a documented residual.** None of the three top-surface BC knobs (h_conv, α_top, ε_top) produced ≥0.24°F authority on cluster-mean CenterRMS. The BC amplitude mismatch identified in PR 18 (D3 ~1.90×, D4 top-concentrated) is not attributable to any single top-surface scalar BC parameter within physically reasonable ranges. PR 20 documents this finding and escalates to a deeper investigation (engine v3: revised BC physics or additional parameter interactions).

## R9 implication

No BC knob authority found. PR 18's D1 DC offset (+0.40–0.56°F cluster mean) and D3 amplitude mismatch remain entangled. R9 Mode-A DC offset cannot be cleanly separated from the BC amplitude residual without first closing D3. R9 assessment is unchanged from PR 18; deferred to engine v3.

