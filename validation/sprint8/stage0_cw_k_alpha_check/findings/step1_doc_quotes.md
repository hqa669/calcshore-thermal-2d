# Stage 0 Step 1 — Verbatim Documentation Quotes

**Source:** ConcreteWorks V3 Training/User Manual (TxDOT 0-6332-P1, July 2013; Published April 2017)
**File:** `~/Downloads/Deep research/ConcreteWorks V3 Training and User Manual (TxDOT 0-6332 P1P2 2013).pdf`
**Location:** Section 3.1.4 "Concrete Thermal Properties", document page 10 (PDF page 20)

---

## Verbatim quote — Section 3.1.4 Concrete Thermal Properties (doc p. 10)

> **3.1.4. Concrete Thermal Properties**
>
> Because of the constantly changing early age properties of concrete, the concrete thermal
> properties must be updated at every time step. The thermal conductivity is known to be a function
> of "the moisture content, content and type of aggregate, porosity, density and temperature (Van
> Breugel, 1998)." The concrete thermal conductivity increases with increasing moisture content.
> There is conflict in the literature about the change in thermal conductivity with increasing
> hydration. Some suggest that the thermal conductivity increases with the degree of hydration,
> while others report that it decreases up to 30% (Van Breugel, 1998; Schindler, 2002). **Based on
> the recommendation of Schindler (2002), ConcreteWorks assumes a linear decrease of the thermal
> conductivity with the degree of hydration from 1.33 times the ultimate thermal conductivity to the
> ultimate thermal conductivity as shown in Equation 23:**
>
> **k_c(α) = k_uc · (1.33 − 0.33 · α)    [Equation 23]**

---

## Source verification

- The equation number "Eq. 23" in `thermal_engine_2d.py:351` is a direct reference to this
  V3 Training Manual equation.
- The Riding 2007 PhD dissertation (the CW theory source) contains an equivalent section at
  Chapter 9 "ConcreteWorks User Manual", §9.3.1 "Heat Transfer Modeling" (dissertation p. 204),
  which repeats the same thermal modeling framework as the V3 Manual.
- Both documents were authored by the same team (Riding, Schindler, Folliard et al.).

---

## Implication for Stage 0

**Outcome 2 (Van Breugel) confirmed.** CW uses exactly the same formula as the engine.
Steps 2-4 are NOT needed — Step 1 alone closes the k(α) question.

The Sprint 7 k_uc × 0.96 calibration was a `k_uc` correction applied uniformly across the
entire curve — it does not affect the slope (−0.33) or the intercept ratio (1.33). Sprint 8
must calibrate k_uc at multiple α values to verify the slope holds in practice; the formula
shape is now theoretically confirmed as matching CW.
