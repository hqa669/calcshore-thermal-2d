# STAGE5a — Final Report: Diagnosis of 18% α_c Gap

## 1. Framing

Stage 4b closed the structural bottom-soil gap, leaving 0.15–4.1°F masked
residuals that scale linearly with |ΔT|. Stage 4a attributed this to an 18%
α_c underestimate (engine 0.00449 vs CW 0.00528 m²/hr at α_hyd=0.8).
Stage 5a tests three candidate causes before any fix is applied.

**Critical finding up front:** the "18% gap" was a diagnosis artifact.
Stage 4a computed α_c at α_hyd=0.8 ("representative late-time"), but in the
9-run suppressed-hydration scenario, the engine actually operates at
α_hyd≈0.036 — not 0.8. At the correct operating point, engine α_c ≈ CW α_c.
The residuals are almost entirely numerical (grid resolution), not parametric.

---

## 2. Engine Operating Point (S5a.1)

Run B (placement=73°F, soil=60°F, |ΔT|=13°F), model_soil=False, hourly output.
Sanity check: masked max|R| at t=168 hr = 1.935°F ✓ (matches Stage 4b to <0.001°F).

| Time | α_hyd (interior bulk) | k_c (W/m·K) | ρCp (J/m³·K) | α_c (m²/hr) |
| --- | --- | --- | --- | --- |
| t=1 hr   | 0.0183 | 3.564 | 2,392,884 | 0.00536 |
| t=24 hr  | 0.0290 | 3.554 | 2,391,667 | 0.00535 |
| t=84 hr  | 0.0334 | 3.550 | 2,391,120 | 0.00534 |
| t=168 hr | 0.0359 | 3.548 | 2,390,794 | 0.00534 |

**Engine's actual α_hyd at t=168 hr: 0.0359 ± 0.0003** (bulk interior, not 0.8).

Hydration is strongly suppressed (Hu_J_kg≈1, α_u=0.10). The engine reaches
only 36% of α_u over 168 hr. The α_c stays essentially constant throughout
the run because α_hyd barely evolves.

**Engine's actual α_c at t=168 hr: 0.00534 m²/hr** vs CW M3 reference
0.00528 m²/hr — ratio 1.012, a +1.2% overshoot, not an 18% deficit.

---

## 3. Cause 1 Evaluation: k(α) Value (S5a.2)

`k_c(α) = k_uc × (1.33 − 0.33·α)`, k_uc = 2.6914 W/m·K.

| α_hyd | k (W/m·K) | α_c (m²/hr) | α_c / CW_target |
| --- | --- | --- | --- |
| 0.036 (operating) | 3.548 | 0.00534 | **1.012** |
| 0.1  | 3.491 | 0.00527 | 0.998 |
| 0.8  | 2.869 | 0.00447 | 0.846 |

The "18% gap" from Stage 4a (ratio 0.846 at α=0.8) is real but irrelevant to
these runs. At the actual operating point α=0.036, the engine's α_c essentially
matches CW.

**Cause 1 contribution: ~0% of the observed residuals.** The engine's k value
is correct for the operating α_hyd.

---

## 4. Cause 2 Evaluation: ρCp Model (S5a.3)

Van Breugel model at α=0.036, T=19°C (mid-range):

| Term | Contribution (J/m³·K) | % |
| --- | --- | --- |
| Wc·α·c_cef (hydrated cement) | 6,123 | 0.3% |
| Wc·(1−α)·Ca (unhydrated cement) | 281,823 | 11.8% |
| Wa·Ca (aggregate) | 1,474,291 | 61.7% |
| Ww·4186 (water) | 628,314 | 26.3% |
| **ρCp total** | **2,390,551** | 100% |

ρCp = 2.39 × 10⁶ J/m³·K — within textbook normal-weight concrete range
(2.1–2.6 × 10⁶ J/m³·K). Temperature sensitivity of ρCp over the 60–73°F run
range: **0.1% swing** (negligible; cement terms are only 12% of total).

**Cause 2 contribution: ~0%.** Van Breugel ρCp is correct for this mix and
operating point.

---

## 5. Cause 3 Evaluation: Grid Resolution (S5a.4)

Run F (placement=73°F, soil=45°F, |ΔT|=28°F — largest Stage 4b residual).

| Resolution | Grid | dx (ft) | dy_mean (ft) | masked max|R| (°F) | Reduction | Wall-clock |
| --- | --- | --- | --- | --- | --- | --- |
| 1× (native) | 21×13 | 1.000 | 6.667 | **4.084** | — | 1.4 s |
| 2× | 41×25 | 0.500 | 3.333 | **1.249** | 3.3× | 1.9 s |
| 4× | 81×49 | 0.250 | 1.667 | **0.409** | 10.0× | 4.3 s |

Convergence is monotone. Richardson extrapolation (second-order FD):

- **Using 1× and 2×:** R_∞ ≈ 0.303°F; numerical fraction = **92.6%**
- **Using 2× and 4×:** R_∞ ≈ 0.129°F; numerical fraction = 89.7%

The 2×/4× estimate of R_∞ = 0.129°F is more reliable (less polluted by the
coarse 1× grid's large error). Numerical contribution to the gap: **~90–93%**.

Note on cell sizes: the native dx=1.0 ft (horizontal) is the limiting factor
for the horizontal gradients, but the native dy_mean=6.67 ft (vertical) is
very coarse for the 80-ft slab. Both axes must be refined together.

**Cause 3 contribution: ~90–93% of the observed residuals.**

---

## 6. Verdict

**Cause 3 (grid resolution) dominates. Causes 1 and 2 are not responsible.**

The apparent "18% α_c gap" from Stage 4a was a false alarm — it arose from
evaluating α_c at the wrong operating point (α_hyd=0.8 instead of ~0.036).
At the true operating point, engine and CW diffusivities match to within 1.2%.
The entire 2–4°F residual at native grid resolution is numerical discretization
error from the 1.0 ft × 6.67 ft cells.

---

## 7. Recommended Stage 5b Fix

**Grid refinement to n_concrete_x=121, n_concrete_y=73 (~6× refinement).**

Rationale: using the Richardson 2nd-order extrapolation from 4× (R_4×=0.409°F,
R_∞=0.129°F), at 6× (h₀/6):

```
R(6×) ≈ R_∞ + (R_4× − R_∞) × (4/6)² = 0.129 + 0.280 × 0.444 ≈ 0.253°F
```

Expected masked max|R| for Run F at 6×: **~0.25°F** — below the 0.35°F gate.

For all 9 runs (using the Stage 4b |ΔT|-linear scaling, which holds because
the dominant error mechanism is grid-resolution-sensitive gradient steepness):

| |ΔT| | Runs | Stage 4b max|R| | Predicted 6× max|R| | Gate |
| --- | --- | --- | --- | --- | --- |
| 0°F  | A | 0.15°F | ~0.02°F | PASS ✓ |
| 13°F | B, D | 1.9°F | ~0.10°F | PASS ✓ |
| 17°F | C, E | 2.4°F | ~0.13°F | PASS ✓ |
| 27°F | G, I | 3.9°F | ~0.22°F | PASS ✓ |
| 28°F | F, H | 4.1°F | ~0.25°F | PASS ✓ |

All 9 runs predicted to pass the 0.35°F gate at 6× resolution.

As a practical minimum, Stage 5b should also test n_concrete_x=81,
n_concrete_y=49 (4×, giving ~0.41°F for Run F) to confirm the 0.409°F result
across all 9 runs. If all 9 runs scale with |ΔT| as expected, 4× passes runs
B–E but marginally fails F, G, H, I; 6× should clear all.

### Performance impact

The engine's inner time step is NOT governed by CFL (dt≈37s is nearly identical
across all resolutions — likely limited by the hydration rate term). Wall-clock
time scales approximately with cell count: 4× (81×49) runs in ~4s vs 1× (21×13)
in ~1.4s — roughly proportional to cells (factor ~15×), which takes only 4s
total. A full 9-run sweep at 4× takes ~35s. Grid refinement is essentially free
at this scale.

### Side effects to watch

- **k(α) at non-suppressed hydration:** the current k(α) formula (k_uc ×
  (1.33−0.33α)) is correct for suppressed runs but produces k=2.87 W/m·K
  at α=0.8 — ~18% below CW. When Stage 5b runs the full-physics validation
  (Hu active, α_hyd→0.8), a residual parametric gap will reappear. Grid
  refinement will not close that gap. **After Stage 5b validates the
  suppressed-hydration gate, the active-hydration k(α) curve should be
  checked separately.**
- **No engine source modifications** in Stage 5b: only `build_grid_half_mat`
  call-site changes (larger n_concrete_x, n_concrete_y).
- **Test compatibility:** existing tests use `ny == 14` (model_soil=False).
  The 4× grid has ny=50 (1 blanket + 49 concrete). Tests that check exact
  grid shape will need updating; tests that check physics should be unaffected.

---

## 8. Estimated Post-Fix Residual (Gate Assessment)

If Stage 5b implements 6× grid (n_cx=121, n_cy=73):

- **Predicted masked max|R| on all 9 runs: ≤ 0.25°F**
- **Gate (0.35°F): ACHIEVABLE** for all runs, with margin.

If 4× only (n_cx=81, n_cy=49):

- Runs A–E (|ΔT|≤17°F): PASS ✓
- Runs F, G, H, I (|ΔT|=27–28°F): borderline, ~0.38–0.41°F — marginally FAIL.
- **Gate: borderline; recommend 6× to be safe.**

---

## 9. Open Questions for Master Chat

**Q1 (primary):** When full hydration is active (Stage 6, Hu>>1, α_hyd→0.8),
the k(α) curve produces k=2.87 W/m·K at α=0.8 — which gives α_c=0.00447 vs
CW 0.00528 m²/hr (18% gap, Stage 4a's original finding). Grid refinement alone
will not close that gap. Should the k(α) calibration be addressed in Stage 5b
(after the suppressed-hydration gate is closed), or deferred to Stage 6?

**Q2 (minor):** The Richardson R_∞ estimate shifts from 0.303°F (1×/2×) to
0.129°F (2×/4×), indicating the solution hasn't fully converged at 4×. The
vertical dy at 4× is still 1.667 ft — quite coarse for an 80-ft slab. Is the
residual dominated by x-direction gradients (near the submerged side face) or
y-direction gradients (near the submerged bottom face)? Running 4×dx-only
(n_cx=81, n_cy=13) and 4×dy-only (n_cx=21, n_cy=49) would isolate the dominant
axis and may allow a smaller-cost targeted fix.
