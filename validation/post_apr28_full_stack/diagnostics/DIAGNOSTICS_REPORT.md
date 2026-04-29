# Post-apr28 full-stack regression — diagnostic checks

**Summary**: All three checks converge on a single explanation. The calibration is correctly wired (CHECK 1). The engine sheds significantly more heat than CW in full-stack mode (CHECK 2). The time-domain residual is DC-offset dominated with weak AC (CHECK 3). Alternative D (BC over-loss) is the only explanation consistent with all three. The pre-apr28 5/14 S0 pass count was error cancellation, not genuine agreement.

---

## CHECK 1 — Calibration wiring

**Outcome: 1B — the apr28 Hu calibration IS applied on the full-stack (run_all.py) path.**

Instrumented MIX-01: `Hu_J_kg = 424143.13`, `Hu_factor_calibrated = 0.9514`, `Hu_J_kg_effective = 403531.01` (ratio 0.951403). The engine peak (125.72°F) matches the committed run_all.py output (Δ = −3.84°F). The call chain `run_all.py → run_one() → load_cw_scenario(default use_cw_calibrated_hu=True) → solve_hydration_2d(Hu = mix.Hu_J_kg_effective)` applies the calibrated value at every step.

**Alternative A (wiring bug) is ruled out.**

---

## CHECK 2 — Adiabatic vs full-stack heat shedding

**Outcome: 2B — the engine sheds significantly MORE heat than CW going from adiabatic to full-stack.**

Across 14 mixes (MIX-13 skipped):

- Mean (engine_drop − cw_drop) = **+3.95°F** (SD 0.95°F)
- Mean as % of adiabatic rise = **+5.33%** (SD 1.42%)
- All 14 mixes positive — the direction is uniform
- Range: +2.19°F (MIX-14, 75°F warm placement) to +6.34°F (MIX-15, 45°F cold placement)

The tight clustering (SD 0.95°F) rules out CW internal mode mismatch as the driver — a composition-dependent or CW-internal discrepancy would produce wider scatter. The cold-regime outlier (MIX-15: CW barely drops 0.07°F at 45°F placement while engine drops 6.27°F) is a secondary modulation of the dominant uniform BC over-loss.

**Alternative B (CW mode mismatch) is ruled out.** The evidence points to a uniform BC-side phenomenon in the engine.

---

## CHECK 3 — Time-domain residual signature (MIX-01 centerline)

**Outcome: 3B — DC offset dominates; AC component is small.**

The MIX-01 full-stack centerline residual (T_engine − T_cw):

| Metric | Value |
|---|---|
| First \|residual\| > 1°F | hr 38.2 |
| \|residual\| max | −3.86°F at hr 167.7 |
| Residual at t = 168 hr | −3.85°F |
| RMS, full trajectory | 2.79°F |
| RMS, [0, 48] hr (early) | 0.69°F |
| RMS, [48, 168] hr (late) | 3.27°F |
| DC offset (late-time mean) | −3.13°F |
| AC component (std after 24-hr detrend) | 0.29°F |

The residual is near-zero in the early phase (early RMS 0.69°F, comparable to the v3-documented structural 0.74°F floor), emerges at hr 38.2 as the concrete begins to cool and BC heat flux grows, then grows monotonically negative through hr 168. The late-phase is dominated by a −3.13°F DC offset; the AC component (0.29°F) is well below the v3-documented top-surface oscillation residual (0.74°F).

**Alternatives C and "kinetics under-supply" are ruled out.** Kinetics under-supply would concentrate the residual in the heating phase (but early RMS is small). Alternative C (kinetics temperature dependence in the cooling tail) would produce a residual concentrated in the final ~20–30 hr after the peak; instead the residual is steady-then-growing from hr 38 onward. Neither matches. A DC-dominant signature with weak AC is the fingerprint of a boundary condition that removes heat at a near-constant rate rather than oscillating with the diurnal cycle.

---

## Combined attribution

**Alternative D (BC over-loss) is the only explanation consistent with all three checks.**

1. The calibration is wired — the engine uses Hu_factor 0.9514, not 1.0. The kinetics under-supply implied by Alternative A does not exist.
2. The engine sheds 3.95°F more heat than CW when transitioning from adiabatic to full-stack boundary conditions, uniformly across all 14 mixes. The BC over-loss is real and substantial.
3. The time-domain residual is DC-dominant (AC = 0.29°F), ruling out the top-surface amplitude oscillation residual that would show strong diurnal AC. The residual emerges post-rise and grows monotonically, consistent with a BC term that accumulates steadily rather than cycling.

The pre-apr28 5/14 S0 pass count on reference mixes (MIX-01, 02, 03, 11, 12) was **error cancellation**: the pre-apr28 engine used Hu_factor ≈ 1.0 (un-calibrated), generating ~3–4°F of extra hydration heat that approximately offset the BC over-loss and produced accidental agreement with CW. The apr28 calibration (Hu_factor ≈ 0.95) correctly removed the over-hydration and thereby unmasked the BC over-loss that was always present. The regression is not a regression in the calibration — it is a correctly calibrated engine exposing a pre-existing BC physics error.

---

## Sprint 7 routing recommendation

**Primary: BC heat-loss recalibration.** Identify which boundary condition term (or combination of terms) causes the engine to shed ~3.95°F more heat than CW in full-stack mode and reduce it to match.

The DC-dominant signature (AC = 0.29°F) provides discriminating information about which BC term is most likely:

| BC term | Expected signature | Compatibility |
|---|---|---|
| Soil-bottom conduction | Steady DC sink to deep thermal mass; minimal diurnal modulation | ✓ (best match) |
| Form-side conduction | Near-DC with form thermal mass moderating diurnal variation | ✓ (good match) |
| LW radiation | Weakly diurnal if dominated by ambient/cloud; moderate DC component | ✓ (reasonable) |
| Top-surface convection | q_conv ∝ (T_conc − T_amb); T_amb has clear diurnal cycling → AC component should be visible | ✗ (AC too small) |

The recommended first ablation experiment for Sprint 7: independently zero out or scale down each BC term and re-run MIX-01, targeting the −3°F DC offset without re-introducing the v3 0.74°F top-surface AC residual. Expected order of magnitude: a 30–50% reduction in soil conductance, OR a 40% reduction in form-side conductance (increase in form R-value), OR a ~20% reduction in LW emissivity, would close the gap. Which is "correct" depends on CW's own parameterization — the next step after identifying the over-loss term is to read CW's BC parameter documentation and align the engine's values to match.

**Runner-up: MIX-15 cold-placement investigation** (+6.34°F drop_diff vs mean +3.95°F). The cold-placement regime appears to activate a BC mode where CW removes almost no heat (−0.07°F CW drop at 45°F) while the engine continues to remove 6.27°F. This is a secondary target — resolving the primary ~4°F uniform offset may or may not close the MIX-15 outlier, and it deserves separate investigation if not.

**Do not begin Sprint 7 work until the diagnostic package is reviewed.** This report is the routing handoff only.
