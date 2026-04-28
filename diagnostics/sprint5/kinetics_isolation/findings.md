# Kinetics-Isolation Test — Findings

**Sprint 5 | MIX-01 adiabatic vs CW | T₀ = 73°F | 2026-04-27**

## Headline

The kinetics-isolation test **fails by +4.32°F** at t = 168 hr (RMS 3.60°F).
The engine's Schindler-Folliard hydration ODE integration is running hot relative
to the CW adiabatic reference, with all BC physics decoupled.

This cleanly disambiguates the engine v3 residuals. In the validated 14-mix
Reference set (60°F placement, steel form, Austin TX), the engine ships at
PeakMax Δ in the range −0.13 to −0.31°F. Those near-zero residuals were a
**coincidental cancellation**: BC heat losses were compensating for a +4.3°F
kinetics overshoot. Engine v3 is not "BC-residual-limited" — it is
"kinetics-error-and-BC-loss-canceling." Future sprint work routes to
`Hu / Cc / αu / Cp(α,T)` before returning to BC structural investigation.

## Gate dispositions

| Gate | Threshold | Observed | Disposition |
|---|---|---|---|
| \|peak Δ\| at t=168 hr | ≤ 1.0°F | +4.317°F | **FAIL** |
| RMS over [0, 168] hr | ≤ 1.0°F | 3.602°F | **FAIL** |
| \|Δ t→100°F\| | ≤ 2.0 hr | 0.769 hr | **PASS** |
| \|Δ t→130°F\| | ≤ 2.0 hr | 4.680 hr | **FAIL** |

Overall: **FAIL** (3/4 gates failing).

Reference benchmarks from `kinetics_isolation_report.txt`:
- Engine peak at t=168 hr: 153.569°F vs CW 149.252°F
- max|Δ| over full run: 4.333°F
- CW kinetics params consumed: Ea=26457.95 J/mol, τ=29.4010 hr, β=0.8950,
  αu=0.75852, Hu=424143.1 J/kg, Cc=575.0 lb/yd³ (39.1% SCM)

## Residual signature

- **t ∈ [0, ~5] hr**: engine and CW are stationary at T₀ = 73°F. Both curves
  match to <0.01°F. This is physically correct: Schindler-Folliard kinetics at
  te ≈ 0 has dα/dt ≈ 0; Arrhenius at 73°F (22.8°C) is ~0.78× reference rate.
- **t_100°F crossing**: engine reaches 100°F at 20.43 hr vs CW at 21.20 hr
  (Δ = −0.77 hr). The inflection timing is within 1 hour — well within the 2 hr
  gate.
- **t ∈ [10, 168] hr**: residual grows monotonically from ≈0 to +4.33°F and
  effectively saturates by t ≈ 100 hr.
- **t_130°F crossing**: engine reaches 130°F at 43.47 hr vs CW at 48.15 hr
  (Δ = −4.68 hr). The engine reaches a higher absolute temperature earlier
  because it is generating more cumulative heat.

**Signature classification** (per passdown failure-signature glossary):

| Wrong feature | Points to |
|---|---|
| Wrong inflection timing | τ, β |
| Wrong steepness | Ea (Arrhenius) |
| Wrong asymptote | αu |
| **Wrong peak amplitude / cumulative heat** | **Hu, Cc, αu, Cp(α, T)** |

The timing is essentially correct (t_100°F within 0.77 hr), so τ, β, and Ea
are **not** the primary issue. The engine simply generates too much total heat
over the 168 hr window. This points to `Hu`, `Cc`, `αu`, or the Van Breugel
specific-heat model `Cp(α, T)`.

## Routing implications for engine v3

The validated 14-mix Reference set ships at:

| Mix | PeakMax Δ |
|---|---|
| MIX-01 | −0.29°F |
| MIX-02 | −0.68°F |
| MIX-03 | −0.30°F |
| MIX-11 | −0.13°F |
| MIX-12 | −0.31°F |

With a +4.3°F adiabatic kinetics overshoot, BC heat losses at the validated
60°F / Austin TX / steel-form conditions are offsetting approximately +4.5°F
of kinetics error — leaving a residual of only −0.13 to −0.29°F at the
centerline peak. This is a **lucky cancellation**, not genuine kinetics fidelity.

For out-of-envelope conditions where the BC heat-loss balance is different,
this cancellation breaks down unevenly:

- **Cold placement (45°F, MIX-15)**: PeakMax Δ = −3.00°F, PeakGrad Δ = −8.20°F.
  Kinetics run colder at lower Arrhenius acceleration, so the overshoot is
  smaller; BC heat loss dominates and the net result goes deeply negative.
- **Warm placement (75°F, MIX-14)**: PeakMax Δ = +1.49°F. Kinetics overshoot
  is larger (more Arrhenius acceleration); BC losses don't compensate fully;
  net result flips positive.
- **High-SCM (Cluster A)**: MIX-06 (+2.20°F), MIX-07 (+4.12°F), MIX-09
  (+3.78°F) — kinetics divergence grows with SCM fraction, consistent with a
  systematic Hu/Cc/αu error that scales with the cementitious content.

This routing strongly suggests the kinetics fix should come before re-opening
the BC structural investigation, because the BC knobs are currently serving
as a hidden compensating offset.

## Prime suspects in the engine

Read these first when continuing the hydration series:

1. **`thermal_engine_2d.py:361` — `specific_heat_variable`**
   Van Breugel Cp model: `c_cef = 8.4·T_C + 339` (fictitious cement
   specific heat). The `Cp` denominator controls how many degrees each
   joule of hydration heat raises the temperature. If `Cp` is
   underestimated at any phase of hydration, the temperature rise will
   be amplified. Check CW's implementation of CW Eq 24-25 against this
   formula for discrepancies in the α-weighting or temperature coefficient.

2. **`thermal_engine_2d.py:1444-1446` — mix → SI conversion block**
   ```python
   Cc = mix.total_cementitious_lb_yd3 * LB_YD3_TO_KG_M3   # kg_cement/m³
   Ca = mix.aggregate_Cp_BTU_lb_F * BTU_LB_F_TO_J_KG_K    # J/(kg·K)
   ```
   `Cc` drives total heat generation (`Q = Hu·Cc·dα/dt`). If
   `total_cementitious_lb_yd3` includes aggregates or water in CW's
   convention (and the engine interprets it as cement-only), the heat
   source is over-stated. Similarly, `aggregate_Cp_BTU_lb_F` feeds the
   Van Breugel Cp denominator; an underestimate here reduces heat
   capacity and amplifies the temperature response.

3. **`cw_scenario_loader.py:388` — `alpha_u` read from input.dat line 388**
   `mix.alpha_u = 0.75852`. Verify this is the actual value CW used to
   generate the reference curve by cross-checking against the handoff
   `input.dat`. If CW computes an effective αu during its run that differs
   from the stored value, the engine and CW diverge at high equivalent age
   (tail of the hydration curve).

## Cross-section uniformity observation

`kinetics_isolation_xs_snapshots.png` shows 34 snapshot panels (t = 0 to
168 hr at 5 hr intervals), each a 3×3 grid rendered in a single flat colour.
The colour shifts from 73°F (initial) to 153.6°F (final). Every panel is
spatially uniform.

This is the **physically correct outcome** for adiabatic mode with uniform IC,
not a bug. The engine's adiabatic ghost-node copy at all four edges
(`thermal_engine_2d.py:1761-1765`) correctly prevents any spatial gradient
from developing when there is no external BC to drive one. The spatial-spread
check confirms: max(T) − min(T) = 0.000e+00°F at t = 168 hr.

Similarly, `kinetics_isolation_centerline_profiles.png` shows flat horizontal
profiles at every snapshot time — uniform temperature with depth, as expected.

The visualizer now supports both CW temp.txt output and engine-emitted npz
files via `--input-npz`. See the updated `visualize_xs_snapshots.py` docstring.

## What this test is not

- Not a calibration step: no engine parameters were adjusted.
- Not a multi-T₀ sweep: this curve is valid for T₀ = 73°F only. Different
  initial temperatures produce different hydration trajectories due to
  Arrhenius temperature acceleration.
- Not a multi-mix validation: only MIX-01 (39.1% SCM Reference) is tested here.

**Future expansion** (see passdown §Future expansion):
Generate CW adiabatic reference curves at T₀ = 50°F, 60°F, 85°F, 100°F and for
MIX-02, MIX-04, MIX-08. A pass at T₀ = 60°F would confirm kinetics is clean at
the validated-envelope placement temp; passes across the T₀ range would give full
Arrhenius coverage.
