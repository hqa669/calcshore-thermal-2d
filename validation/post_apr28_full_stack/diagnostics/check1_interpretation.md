# CHECK 1.3 — Interpretation

**Outcome: 1B — the calibration IS being applied on the run_all.py path.**

Empirical evidence from `check1_mix01_output.txt`:

```
Hu_J_kg              = 424143.13
Hu_J_kg_effective    = 403531.0094293787
Hu_factor_calibrated = 0.9514029130434783
ratio (effective/raw) = 0.951403
```

The factor `0.9514` matches the apr28 release-notes Hu_factor for MIX-01 (Type I/II, 60.9% cement, 21.7% FA-F, 17.4% slag) to four decimal places. The dataclass default `Hu_factor_calibrated = 1.0` has been overwritten by `compute_hu_factor()` at `cw_scenario_loader.py:822`, and `Hu_J_kg_effective = 403531` is consistent with `424143 × 0.9514 = 403532` (matches to 1 J/kg, rounding-equivalent).

The post-solve sanity check reproduces the committed run exactly:
- Engine peak: 125.72°F @ 152.0 hr (`run_all_output.md` recorded 125.7°F @ 152.0 hr)
- CW peak: 129.56°F @ 145.8 hr
- Δ = −3.84°F (matches the −3.84°F in `comparison_table.md`)

**Conclusion**: The post-apr28 regression is NOT a wiring bug. `Hu_factor_calibrated = 0.9514` propagates correctly to `Hu_J_kg_effective = 403531`, the solver consumes that calibrated value at `thermal_engine_2d.py:1438`, and the resulting peak (125.72°F) is exactly what was recorded in the committed full-stack run. Alternative A (wiring bug) is ruled out.

**Side observation worth flagging** (not a hard finding): `Hu_calibration_note` is an empty string. The note field exists to surface envelope-violation warnings (e.g., "Silica fume 20.8% exceeds the 15% envelope" for MIX-09). For MIX-01 — well within the validated composition envelope — an empty note is expected, but the loader should arguably emit a positive confirmation note like "in-envelope" so downstream consumers can distinguish "calibration applied with no envelope flags" from "calibration not applied". This is a UX nit, not a correctness issue.

**Proceed to CHECK 2.**
