# CHECK 3.3 — Interpretation

**Outcome: 3B — the residual is dominantly a DC offset that grows slowly through the simulation. AC (diurnal) component is small. The signature is consistent with a uniform BC heat-loss term, not with a top-surface amplitude oscillation issue.**

## Numerical signature (MIX-01 full-stack centerline)

```
First |residual| > 1°F:  hr 38.2
|residual| max:          -3.86°F at hr 167.7  (essentially at end)
Residual at t=168 hr:    -3.85°F
RMS, full trajectory:    2.79°F
RMS, [0, 48] hr (early): 0.69°F
RMS, [48, 168] hr (late):3.27°F
DC offset (late-time mean):                          -3.13°F
AC component (std after 24-hr running-mean detrend):  0.29°F
```

## What the signature tells us

**1. The early phase tracks closely.** RMS over [0, 48] hr is 0.69°F, comparable to the RMS the v3 release notes documented as the structural BC residual (0.74°F). This means the engine and CW reproduce the hydration rise and early-phase boundary response very similarly. The kinetics correction is not breaking the early-phase shape.

**2. The residual emerges past the rise.** The first time |residual| exceeds 1°F is at hr 38.2 — past the hydration peak's onset, when concrete temperature is high enough to drive significant BC heat flux. From hr 38 onward, the residual grows monotonically negative (engine cooler than CW).

**3. The late phase is dominated by DC offset.** The late-time mean residual is -3.13°F. The AC component (after detrending with a 24-hour running mean) has std 0.29°F — a small fraction of the DC offset and well below the v3-documented 0.74°F top-surface amplitude residual. This is diagnostic: **the new post-apr28 residual is not the v3 0.74°F top-surface oscillation issue. A different mechanism, producing a steady DC offset, has been unmasked.**

**4. The residual is still growing at simulation end.** At hr 168, the residual is -3.85°F vs. the late-time mean of -3.13°F — the engine continues to cool relative to CW as the run progresses. The integrated heat loss continues to compound. This is consistent with a BC heat-loss path that the engine over-applies, with the cumulative effect growing throughout the cooling phase.

## Discriminating between candidates

| Candidate | Predicted signature | Match? |
|---|---|---|
| Top-surface BC amplitude over-oscillation (v3 documented 0.74°F) | Strong AC, small DC | ✗ (AC is 0.29°F, small) |
| Kinetics under-supply (uniform Hu shortfall) | Residual concentrated in the heating phase, not late | ✗ (early RMS 0.69°F, late 3.27°F) |
| Kinetics with hidden temperature dependence (Alternative C) | Residual concentrated in the cooling tail when T drops | ✗ (residual steady-then-growing, not cooling-tail-only) |
| Uniform BC over-loss (top conv / form conv / LW / soil) | DC offset emerging post-rise, slow growth, minimal AC | ✓ |

The signature most cleanly matches **Outcome 3B (DC offset dominates with small AC)**, supplemented by the observation that the residual is continuing to grow (not parallel-and-stable). The BC mechanism is unambiguously suggested.

## What this implies for Sprint 7

The DC dominance with weak AC tells us the over-loss is *not* a top-surface convection / radiation amplitude issue (those would carry stronger diurnal AC). Candidates ranked by signature compatibility:

1. **Soil-bottom conduction** — produces steady DC heat flow to a deep thermal sink; minimal diurnal modulation; cumulative loss grows over simulation.
2. **Form-side conduction** — produces near-DC heat flow to the form-soil-air interface; some diurnal but moderated by form thermal mass; cumulative loss grows.
3. **LW radiation** — depends on effective sky temperature; some diurnal but typically small over a 168-hour window if sky temp is dominated by cloud cover and ambient.
4. **Top-surface convection** — produces both DC and AC; the AC component should be larger than 0.29°F if this were the dominant term, since `q_conv ∝ (T_conc − T_amb)` and `T_amb` has clear diurnal cycling.

The cleanest first ablation experiment for Sprint 7: independently zero out / scale down each BC term and re-run MIX-01 to see which closes the -3°F DC offset without re-introducing the 0.74°F AC residual. The expected order of magnitude effect: a 30–50% reduction in soil conduction, OR a 40% reduction in form-side R-value (i.e., higher R = less heat loss), OR a ~20% reduction in convection coefficient, would close the gap. Which of these is "correct" depends on whether CW's own BC physics is parameterized differently than the engine's.

**Proceed to combined attribution and Sprint 7 routing.**
