# Sprint 4 — S0 Gate Baseline Spike

| Mix | Location | Climate | SCM% | Family | PeakMax Δ | PeakGrad Δ | FieldRMS | CenterRMS | CornerRMS | S0 pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| MIX-01 | TX, Austin | hot-humid | 39.1% | mid_scm | -0.3°F ✓ | -0.5°F ✓ | 0.88°F ✓ | 0.74°F ✓ | 2.26°F ✓ | 5/5 |
| MIX-02 | TX, Austin | hot-humid | 39.1% | mid_scm | -0.7°F ✓ | ++1.0°F ✓ | 1.23°F ✓ | 0.57°F ✓ | 1.32°F ✓ | 5/5 |
| MIX-03 | TX, Austin | hot-humid | 39.1% | mid_scm | -0.3°F ✓ | -1.1°F ✓ | 0.87°F ✓ | 0.79°F ✓ | 2.76°F ✓ | 5/5 |
| MIX-04 | TX, Austin | hot-humid | 21.7% | mid_scm | -1.5°F ✗ | -2.5°F ✗ | 1.37°F ✓ | 0.93°F ✓ | 3.03°F ✗ | 2/5 |
| MIX-05 | TX, Austin | hot-humid | 52.2% | high_scm | ++0.7°F ✓ | ++0.8°F ✓ | 1.24°F ✓ | 1.49°F ✗ | 2.09°F ✓ | 4/5 |
| MIX-06 | TX, Austin | hot-humid | 71.3% | high_scm | ++2.2°F ✗ | ++2.9°F ✗ | 2.35°F ✗ | 2.81°F ✗ | 1.89°F ✓ | 1/5 |
| MIX-07 | TX, Austin | hot-humid | 91.3% | high_scm | ++4.1°F ✗ | ++4.8°F ✗ | 3.82°F ✗ | 4.45°F ✗ | 1.88°F ✓ | 1/5 |
| MIX-08 | TX, Austin | hot-humid | 17.4% | low_scm | -1.4°F ✗ | -3.6°F ✗ | 1.33°F ✓ | 0.85°F ✓ | 4.05°F ✗ | 2/5 |
| MIX-09 | TX, Austin | hot-humid | 41.7% | high_scm | ++3.8°F ✗ | ++1.1°F ✓ | 4.15°F ✗ | 4.84°F ✗ | 4.58°F ✗ | 1/5 |
| MIX-10 | TX, Austin | hot-humid | 52.2% | high_scm | ++0.2°F ✓ | ++1.2°F ✓ | 1.07°F ✓ | 1.03°F ✗ | 1.46°F ✓ | 4/5 |
| MIX-11 | TX, Austin | hot-humid | 39.1% | mid_scm | -0.1°F ✓ | -0.1°F ✓ | 0.91°F ✓ | 0.79°F ✓ | 2.16°F ✓ | 5/5 |
| MIX-12 | TX, Austin | hot-humid | 39.1% | mid_scm | -0.3°F ✓ | -0.5°F ✓ | 0.87°F ✓ | 0.73°F ✓ | 2.27°F ✓ | 5/5 |
| MIX-13 | — skipped (no_cw_output) — |  |  |  |  |  |  |  |  |  |
| MIX-14 | TX, Austin | hot-humid | 39.1% | mid_scm | ++1.5°F ✗ | -3.7°F ✗ | 3.20°F ✗ | 2.29°F ✗ | 4.86°F ✗ | 0/5 |
| MIX-15 | TX, Austin | hot-humid | 39.1% | mid_scm | -3.0°F ✗ | -8.3°F ✗ | 4.06°F ✗ | 2.07°F ✗ | 2.05°F ✓ | 1/5 |

## Surprises noted

- **All 15 mixes are `TX, Austin`** — climate-bucket column is uniform (`hot-humid`).
  Climate-based clustering is N/A for this dataset.
- **`cw_comparison_MIX-01.png` hardcoded** in `compare_to_cw.py:437` — every run
  overwrites the same PNG in the repo root. PR 13 should fix the output path.
- **MIX-01 constants at `compare_to_cw.py:41-44`** (`CW_PEAK_MAX_F`, etc.) are used
  to display the CW peak *time* label at line 259 for the gradient row. For
  non-MIX-01 mixes the displayed CW peak time is wrong, but the delta/pass values
  are computed from runtime CW data and are trustworthy. Fix belongs in PR 13.
