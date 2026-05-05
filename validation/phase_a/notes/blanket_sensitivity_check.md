# Phase A Blanket Sensitivity Check

**Date:** 2026-05-04
**Context:** Phase A closure — verifying that CW input.dat blanket R-value
does not affect the Phase A gate region.

## Question

CW input.dat files for Phase A were generated with blanket R=5.68 deg F.ft^2.hr/Btu
(default). The engine wrapper sets `blanket_thickness_m = 0.0`. Does this
mismatch affect Phase A's validation result?

## Test

Re-ran CW for MIX-01 at scenario S3 (T_pl=73 deg F, T_soil=85 deg F) at two R-values:
- R=5.68 (default used in Phase A)
- R=0.18 (CW minimum)

Compared the two output.txt files cell-by-cell across all 49x13 grid points
and all 2016 timesteps.

## Result

**Lower half (di in [24, 48], all wi, all t): bit-identical.** Maximum
difference between the two runs is 0.000 deg F at every cell in the validated
region. T_max in the lower half is 65.13 deg C (149.234 deg F) in both runs,
occurring at di=24, wi=0, t=167.58 hr in both. delta-T_max in the lower half
is 35.69 deg C (64.24 deg F) in both runs.

**Top region (di < 5):** Large differences — up to 25 deg C (45 deg F) at
di=0, decaying to 0.5 deg C (0.9 deg F) by di=5 and to numerical zero by di=10.

## Implication for Phase A

The blanket R-value is *effectively inert* for Phase A's gate region
(di in [5, 48]). The R=5.68 used in Phase A's CW data does not affect the
+/-0.82 deg F sub-envelope claim or the +/-1.72 deg F worst-case claim. Phase A's
result is independent of blanket configuration in CW.

## Implication for Phase C (top-surface validation)

The blanket-handling difference between the engine (`blanket_thickness_m = 0.0`)
and CW (R-min = 0.18 deg F.ft^2.hr/Btu) DOES affect the top-surface region.
Phase C work must explicitly handle blanket-BC consistency between engine
and CW before validating against di in [0, 4].
