# STAGE1 — Reciprocity Check

## Method

**Corrected residual (used for linearity verdict):**
`R_corr(x, y, t) = (T_pos − T_place_pos) + (T_neg − T_place_neg)`

This measures whether the DEVIATION of each run from its own placement temperature
is equal and opposite — the true reciprocity condition for a linear system.

**Brief's original formula (shown for reference):**
`R_orig = T_pos + T_neg − 2 × T_baseline`

Note: `R_orig ≈ −ΔT` at all times for any linear diffusion system with mismatched
initial conditions (Run B starts at 73°F, Run D starts at 60°F). This is a
mathematical identity, not non-linearity. The brief's formula conflates the IC
offset with the reciprocity residual; `R_corr` removes that offset.

## Note on Run I

The brief's IMPORTANT warning stated that Run I has soil=60°F and ΔT=−40, making it
a non-mirror of G. **This was incorrect.** Actual `input.dat` inspection confirms:
- Run I (thermal_100_73): placement=100°F, soil=73°F, ΔT=−27°F
- Run G (thermal_73_100): placement=73°F, soil=100°F, ΔT=+27°F

Run G ↔ I is a clean mirror pair at |ΔT|=27°F and is included in the analysis.

## Results

| Pair | Run+ | Run− | |ΔT|°F | max|R_corr|°F | mean|R_corr|°F | mean|R_corr| interior°F | max|R_orig|°F | t of max | Location of max |
| ---- | ---- | ---- | ------ | ------------ | ------------- | ------------------------------ | ------------ | -------- | --------------- |
| BD | B (73/60) | D (60/73) | 13 | 5.1000 | 0.2572 | 0.1947 | 13.7160 | 2.0 hr | depth=0.00m, width=6.10m |
| CE | C (73/90) | E (90/73) | 17 | 5.9400 | 0.3269 | 0.2535 | 17.0640 | 2.0 hr | depth=0.00m, width=6.10m |
| FH | F (73/45) | H (45/73) | 28 | 6.7320 | 0.2366 | 0.1694 | 30.3480 | 2.0 hr | depth=0.00m, width=1.02m |
| GI | G (73/100) | I (100/73) | 27 | 7.2280 | 0.3567 | 0.2746 | 27.0900 | 2.0 hr | depth=0.00m, width=6.10m |

## Interpretation

### Top-surface observation

All pairs show max|R_corr| = 5–7°F concentrated at **depth=0m (top surface), t≈2hr**.
This is from CW's atmospheric top BC (solar + evaporation), which is non-linear in
concrete surface temperature. It does NOT reflect soil-coupling non-linearity.

### Soil-coupling linearity (interior domain)

**CW's soil-concrete coupling is linear within 0.5°F** (mean interior residual). Worst interior mean: 0.2746°F (pair GI). The dataset has 4 effective independent degrees of freedom: A (baseline) plus |ΔT|∈{13, 17, 27, 28}°F.

## Plots

- `plots/reciprocity_BD.png` — residual field at t=168 hr
- `plots/reciprocity_CE.png` — residual field at t=168 hr
- `plots/reciprocity_FH.png` — residual field at t=168 hr
- `plots/reciprocity_GI.png` — residual field at t=168 hr
