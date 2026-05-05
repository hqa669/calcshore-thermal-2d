# Sprint 8 Stage 1 — CW k(α) Empirical Verification

**Stage predecessor:** Stage 0 (documentation pass confirming CW V3 manual states Van Breugel Eq. 23)
**Stage type:** Read-only data analysis. No engine changes, no commits.
**Question answered:** Does CW empirically use α-dependent k(α)?

---

## Setup

Three CW runs at identical geometry/BC/placement (73°F placement, 100°F soil, 80×40 ft mat), varying only hydration kinetics:

| Dataset | α_u | τ (hr) | β | Hu (J/kg) |
|---|---|---|---|---|
| lowhyd | 0.10 | 200 | 0.10 | 1.0 |
| highhyd | 0.70 | 10 | 0.85 | 1.0 |
| highhyd_b010 | 0.70 | 10 | 0.10 | 1.0 |

With Hu=1 J/kg in all three, hydration heat release is negligible. Any T-field divergence across runs must come from k(α).

---

## Key results

**input.dat consistency:** PASS — only {tau_hrs, beta, alpha_u} differ. All mix design, geometry, BC, and placement parameters are identical.

**Pairwise ΔT at t=168 hr** (mask: di=4–47, all wi):

| Pair | max\|ΔT\| | mean\|ΔT\| | RMS ΔT |
|---|---|---|---|
| A: highhyd − lowhyd | **2.16°F** | 1.44°F | 1.63°F |
| B: highhyd_b010 − lowhyd | **1.03°F** | 0.66°F | 0.76°F |
| C: highhyd − highhyd_b010 | **1.13°F** | 0.77°F | 0.87°F |

**Engine α(t) at t=168 hr (isothermal, what CW ran):**
- lowhyd: α = 0.036 → k_c = k_uc × 1.318
- highhyd: α = 0.638 → k_c = k_uc × 1.119 (−15.1% vs lowhyd)
- highhyd_b010: α = 0.329 → k_c = k_uc × 1.221 (−7.3% vs lowhyd)

---

## Verdict

**Outcome 1: CW empirically uses Van Breugel α-dependent thermal conductivity.**

Pair A max|ΔT| = 2.16°F — 20× the Outcome 2 threshold of <0.1°F. The magnitude ordering (A > B, A > C), the ~2:1 ratio matching Δα ratios, and the monotonically growing ΔT over time are all consistent with Van Breugel k(α) = k_uc · (1.33 − 0.33·α).

See `findings/synthesis.md` for the full quantitative analysis.

---

## Files

```
cw_data/          CW reference datasets (read-only copies)
findings/
  input_consistency.md       §3.2 — parameter comparison table
  pairwise_dT_stats.csv/md   §3.4 — pairwise ΔT statistics
  engine_alpha_reference.csv/md  §3.8 — engine α(t) trajectories
  synthesis.md               §7  — observational synthesis
figures/
  stage1_cw_pairwise_dT_fields.png   3×3 ΔT heatmap grid
  stage1_cw_midprofile_t168.png      mid-depth lateral profile at t=168
  stage1_cw_CL_T_evolution.png       centerline T(t) 0–168 hr
stage1_input_consistency.py
stage1_pairwise_dT.py
stage1_figures.py
stage1_engine_alpha_reference.py
```
