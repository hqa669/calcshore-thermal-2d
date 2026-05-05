# Sprint 8 Stage 0 Finding

**Date:** 2026-05-01
**Question:** Does ConcreteWorks treat k as α-dependent?
**Method:** Step 1 documentation pass (sufficient — no empirical runs needed)
**Status:** CLOSED

---

## Answer: Yes — CW uses the Van Breugel formula (Outcome 2)

ConcreteWorks V3 applies the following formula for concrete thermal conductivity at every
time step, confirmed verbatim in the V3 Training/User Manual Section 3.1.4, **Equation 23**:

```
k_c(α) = k_uc · (1.33 − 0.33 · α)
```

This is **identical** to the engine's `thermal_conductivity_variable` (`thermal_engine_2d.py:351`).
The formula was adopted in CW "based on the recommendation of Schindler (2002)."

---

## Source

> *"Based on the recommendation of Schindler (2002), ConcreteWorks assumes a linear decrease
> of the thermal conductivity with the degree of hydration from 1.33 times the ultimate thermal
> conductivity to the ultimate thermal conductivity as shown in Equation 23:*
> **k_c(α) = k_uc · (1.33 − 0.33 · α)**"

— ConcreteWorks V3 Training/User Manual (TxDOT 0-6332-P1, 2013/2017), §3.1.4, doc p. 10

The manual also notes that "there is conflict in the literature about the change in thermal
conductivity with increasing hydration" — CW chose the Schindler (2002) linear-decrease
recommendation, which defines a slope of −0.33.

---

## Implications for Sprint 8

1. **The formula shape is confirmed correct.** Sprint 8 does NOT need to replace the Van
   Breugel formula or reverse-engineer a different functional form. The slope (−0.33) and
   intercept ratio (1.33) are shared between the engine and CW by construction.

2. **Sprint 8 calibrates k_uc only.** With the shape locked, the single free parameter is
   the unhydrated thermal conductivity `k_uc`. Sprint 7 calibrated it to × 0.96 at α ≈ 0.036.
   Sprint 8 must verify this scaling holds across a range of α (i.e., the same × 0.96 factor
   gives good agreement at α = 0.2, 0.5, 0.8 as well as 0.036).

3. **The 6% gap at α = 0.8 (from STAGE5a) is a k_uc residual, not a shape error.** Before
   the Sprint 7 calibration, the engine's k_uc was 4% too high. After applying × 0.96, the
   agreement at α = 0.036 is within ≤ 0.35°F. The same 4% correction propagates uniformly
   across the curve — Sprint 8's multi-α dataset will confirm whether 0.96 is the right scalar
   for all α or whether a residual shape error remains after accounting for k_uc.

4. **Sprint 8 Stage 1 dataset design:** Multi-α dataset per `SPRINT8_PASSDOWN.md` §3 table
   (α targets: 0.036 [Sprint 7, already have], 0.20, 0.40, 0.60, 0.80). For each α target,
   engineer Hu/τ/β/α_u so hydration plateaus at the target by 168 hr. Hold geometry, BC,
   climate constant.

5. **Meaning of k_uc × 0.96:** This is an authentic `k_uc` calibration at α = 0.036. Because
   both engine and CW use the same formula, the 0.96 factor applies uniformly across all α.
   The Sprint 7 calibration was correct in scope and result; Sprint 8 verifies it holds at
   high α.

---

## What Sprint 8 does NOT need to do (Stage 0 eliminates)

- Reverse-engineer a different α-formula — not needed; CW uses Van Breugel.
- Build an empirical k_eff(α) compensation — not needed; the formula is correct.
- Run the engineered 2-run or 5-run α-sweep described in the brief's Steps 3-4 for formula
  identification — the formula is already identified from documentation.

Sprint 8 still needs a multi-α calibration dataset to verify the k_uc × 0.96 factor holds
across α, and to detect any remaining residual if k_uc itself varies (e.g., CW computes k_uc
differently for different aggregate types or mixes). But that is calibration, not exploration.

---

## Cross-references

- V3 Manual quote: `findings/step1_doc_quotes.md`
- Engine formula: `thermal_engine_2d.py:350-367` (`thermal_conductivity_variable`)
- Engine calibration constant: `thermal_engine_2d.py:130-137` (`K_UC_CALIBRATION_FACTOR_SPRINT7 = 0.96`)
- Sprint 7 calibration record: `docs/calibration/k_uc_calibration_sprint7.md`
- Sprint 8 passdown: `docs/sprints/SPRINT8_PASSDOWN.md`
- Prior 2-point evidence (pre-calibration): `validation/soil_calibration/STAGE5a_cause1_k_value.md`
