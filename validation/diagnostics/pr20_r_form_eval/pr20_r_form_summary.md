# PR 20 — R_form Reference-Set Evaluation Summary

## Outcome: DO NOT COMMIT R_form=0.060

R_form=0.060 fails the Reference-set generalization test. Commit-decision
criterion C1 is not met; the recalibration is routed to engine v3 release
notes rather than ADR-04.

## Per-criterion pass/fail

| Criterion | Result | Evidence |
| --- | :---: | --- |
| C1: all 5 Reference mixes S0 5/5 at R_form=0.060 | **FAIL** | MIX-02 S0=4/5 (PeakGradΔ=+2.679°F > 2.0°F threshold) |
| C2: ≥3/5 mixes improve CornerRMS at 0.060 vs 0.0862 | PASS | 4/5 improve (MIX-01, 03, 11, 12); MIX-02 does not |
| C3: no mix CenterRMS degrades >0.05°F at 0.060 vs 0.0862 | PASS | Δ=0.000°F on all 5 mixes — CenterRMS is fully form-face-decoupled |

## CornerRMS per mix: R_form=0.060 vs 0.0862

| Mix | CornerRMS at 0.0862 (°F) | CornerRMS at 0.060 (°F) | Δ (°F) | Direction |
| --- | ---: | ---: | ---: | --- |
| MIX-01 | 2.2192 | 1.0034 | −1.216 | improves ✓ |
| MIX-02 | 1.2362 | 1.8113 | +0.575 | **worsens ✗** |
| MIX-03 | 2.7238 | 1.0923 | −1.632 | improves ✓ |
| MIX-11 | 2.1198 | 1.0244 | −1.095 | improves ✓ |
| MIX-12 | 2.2251 | 1.0016 | −1.224 | improves ✓ |

## Generalization finding: kinetics-coupled BC heterogeneity

The cluster mixes (MIX-01, 11, 12) and MIX-03 show a consistent,
large-magnitude improvement at R_form=0.060: CornerRMS drops from the
2.1–2.7°F range to the 1.0–1.1°F range. PeakGrad also moves toward pass
for these mixes (e.g., MIX-01 PeakGradΔ: −0.288°F at 0.0862 → +1.534°F at
0.060; all still below the 2.0°F threshold). The 4/5 cluster+MIX-03 result
is a robust, high-confidence improvement.

MIX-02 behaves oppositely across the entire sweep. As R_form decreases from
0.100 to 0.050, MIX-02's CornerRMS and PeakGrad worsen monotonically:

| R_form | MIX-02 CornerRMS | MIX-02 PeakGradΔ | MIX-02 S0 |
| ---: | ---: | ---: | :---: |
| 0.0500 | 2.3631°F | +3.346°F | 4/5 |
| 0.0600 | 1.8113°F | +2.679°F | 4/5 |
| 0.0700 | 1.4053°F | +2.064°F | 4/5 |
| **0.0862** | **1.2362°F** | **+1.160°F** | **5/5** |
| 0.1000 | 1.5493°F | +0.603°F | 5/5 |

MIX-02's response has a minimum CornerRMS at R_form≈0.0862 (the current
ADR-04 default) and is at S0 5/5 only from 0.0862 upward. The cluster's
CornerRMS minimum is around R_form=0.060–0.070. The two regimes do not
overlap at a common globally optimal value.

MIX-02 is kinetics-anomalous per `mix02_recon.md`: its hydration heat
profile diverges from the cluster's. The physical mechanism behind the
opposite R_form response is likely that MIX-02's phase relationship between
internal heat generation and the ambient diurnal cycle differs — lower form-
face contact resistance changes the direction of thermal gradient at the
corner rather than just its magnitude. This is a fundamentally different
sensitivity regime, not a quantitative deviation.

**Implication for Sprint 6 / R5**: the heterogeneity is likely parameterizable.
If R_form is allowed to vary by form material or mix design type (R5's scope),
a per-form-type optimal can be found without a global compromise. The cluster
optimum (~0.060–0.070) and MIX-02's optimum (~0.0862) would each be
recoverable under a per-type parameterization. This routes cleanly to R5.

## C3 note: CenterRMS form-face decoupling confirmed

CenterRMS is identical (to machine precision) across all 5 R_form sweep
values for every mix. The Δ values at R_form=0.060 vs 0.0862 are all
0.0000°F. This confirms: the centerline mid-depth temperature is fully
decoupled from the form-face contact resistance. The structural centerline
residual documented in PR 18/19 has no R_form path. This is consistent
with R_form being a local BC parameter at the form face (x=0 edge),
separated from the centerline (x=W/2 symmetry) by the full half-width
of the concrete.

## Routing

The R_form=0.060 recalibration cannot be committed globally without a
per-mix or per-form-type override seam. Routes to:

1. **Engine v3 release notes** (`docs/engine_v3_release_notes.md`, Item 2):
   document that the CornerRMS residual at the current R_form default has
   a known recalibration candidate (0.060) that benefits 4/5 Reference mixes
   but cannot be applied globally due to MIX-02's opposite response.

2. **Sprint 6 R5** (form-material parameterization): evaluate whether the
   R_form heterogeneity across mixes is better explained by form material,
   ambient coupling, or kinetics × BC interaction. If the per-form-type seam
   can isolate MIX-02's regime, the R_form=0.060 improvement for the cluster
   can land cleanly.

ADR-04 remains at 0.0862 for now. The load-bearing note in §7.6.3 (0.0862
sits at the high end of MIX-01's 5/5 tolerance window with limited upward
margin) is confirmed by this sweep; the upward margin is because 0.0862 is
*already near MIX-02's optimum*, not because the cluster wants higher
resistance.
