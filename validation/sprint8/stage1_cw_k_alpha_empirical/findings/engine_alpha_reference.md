# Sprint 8 Stage 1 — Engine α(t) Reference (§3.8)

Placement temperature: 73°F for all three datasets.

- **isothermal**: CW actual regime (Hu=1 J/kg ≈ no self-heating, T fixed at T_pl).

- **adiabatic_ref**: hypothetical with Hu=350 kJ/kg, Cc=400 kg/m³ (representative OPC, for context only).

| Dataset | Regime | α(24hr) | α(48hr) | α(84hr) | α(120hr) | α(168hr) |
|---|---|---|---|---|---|---|
| lowhyd (α_u=0.10, τ=200, β=0.10) | isothermal | 0.0290 | 0.0315 | 0.0335 | 0.0349 | 0.0361 |
| lowhyd (α_u=0.10, τ=200, β=0.10) | adiabatic_ref | 0.0294 | 0.0320 | 0.0340 | 0.0354 | 0.0366 |
| highhyd (α_u=0.70, τ=10,  β=0.85) | isothermal | 0.4326 | 0.5360 | 0.5930 | 0.6193 | 0.6384 |
| highhyd (α_u=0.70, τ=10,  β=0.85) | adiabatic_ref | 0.6388 | 0.6749 | 0.6863 | 0.6904 | 0.6931 |
| highhyd_b010 (α_u=0.70, τ=10,  β=0.10) | isothermal | 0.2796 | 0.2974 | 0.3115 | 0.3205 | 0.3289 |
| highhyd_b010 (α_u=0.70, τ=10,  β=0.10) | adiabatic_ref | 0.3103 | 0.3297 | 0.3451 | 0.3547 | 0.3636 |