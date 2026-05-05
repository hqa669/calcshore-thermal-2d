# Sprint 8 Stage 0 — CW k(α) Diagnostic

**Status: CLOSED — Step 1 definitive**
**Date closed: 2026-05-01**

## Question
Does ConcreteWorks treat concrete thermal conductivity k as a function of degree of hydration α?

## Answer
**Yes — CW uses the Van Breugel formula (Outcome 2):**

```
k_c(α) = k_uc · (1.33 − 0.33 · α)    [CW V3 Manual §3.1.4 Eq. 23]
```

This is identical to the engine's `thermal_conductivity_variable` at `thermal_engine_2d.py:351`.

## How it was determined
Step 1 documentation pass. ConcreteWorks V3 Training/User Manual (TxDOT 0-6332-P1, 2013/2017),
Section 3.1.4 "Concrete Thermal Properties", page 10, states the formula verbatim and cites
Schindler (2002) as the authority. Steps 2-4 (input inspection, empirical runs) were not needed.

## Deliverables
- `findings/step1_doc_quotes.md` — verbatim source quote
- `findings/stage0_finding_1page.md` — 1-page finding + Sprint 8 implications

## Sprint 8 implication
Sprint 8 calibrates k_uc only (the scalar multiplier). The shape (slope −0.33, intercept 1.33)
is theoretically confirmed as matching CW. Multi-α dataset still needed to verify k_uc × 0.96
from Sprint 7 holds across the full α range.
