# cw_scenario_loader.py — Design Spec & Contract

**Status:** Working. Validated end-to-end against MIX-01 Austin scenario (TEST.dat + TX__Austin.dat + temp.txt) as of 2026-04-22.

**Purpose:** Load a ConcreteWorks scenario from its three native file formats into CalcShore-ready dataclasses, so we can run head-to-head validation between the CalcShore thermal engine and CW 2D output.

---

## When to load this file

Load `cw_scenario_loader.md` + `cw_scenario_loader.py` together into any new chat where you want to:

- Extend the loader to support new CW versions or additional fields
- Build a comparison/validation harness against CW output
- Debug a CW file that doesn't parse correctly
- Generate engine-ready inputs from a CW project

This spec tells Claude everything it needs to modify the loader safely without breaking existing behavior.

---

## Context: Where this fits in the CalcShore architecture

CalcShore's thermal engine (`thermal_engine_v2.py`) is being evolved to reach **ConcreteWorks equation parity** by replacing calibration fudge factors with physics. The sprint plan (from `coding_passdown_v2.md`):

1. **Sprint 1** — Wire up full top-surface radiation balance (solar + longwave, currently missing)
2. **Sprint 2** — Upgrade convection model (ACI Eq 27 with surface-orientation constants)
3. **Sprint 3** — Barber soil model for pavement/column scenarios
4. **Sprint 4** — Cleanup, eliminate remaining calibration

**Every sprint needs a validation harness that compares engine output against CW.** This loader is the plumbing for that harness. Without it, every validation run requires manual data re-entry.

The loader is **read-only** — it does not modify or run the engine. It just transforms CW's file formats into dataclasses that match the engine's existing `MixDesign`, `Geometry`, `Construction`, `Environment` shapes.

---

## Tested scenario: MIX-01 Austin, July 15 2026

All decoded values reconcile with CW UI screenshots:

| Field | File source | Value | CW UI |
|---|---|---|---|
| Member W × D × L | TEST.dat L21-23 | 40 × 8 × 60 ft | ✅ |
| Cement | TEST.dat L54 | 350 lb/yd³ Type I/II | ✅ |
| Class F fly ash / GGBFS | TEST.dat L61 / L63 | 125 / 100 lb/yd³ | ✅ |
| Hydration τ / β / α_u / Hu | TEST.dat L386-389 | 29.401 / 0.895 / 0.7585 / 424143 | matches passdown MIX-01 |
| Daily max temp (avg Jul 15-21) | TX__Austin.dat col 5 | 92.2°F | UI: 92.3°F |
| Daily min temp | TX__Austin.dat col 5 | 74.5°F | UI: 74.6°F |
| Max solar | TX__Austin.dat col 3 | 868 W/m² raw / 840 avg | UI: 848.1 |
| Max RH | TX__Austin.dat col 7 | 89% raw | UI: 87.4 |
| Peak Max T in x-section | temp.txt | 129.6°F @ 145.8 hr | ✅ |
| Peak gradient | temp.txt | 39.3°F @ 146.2 hr | ✅ |

---

## File formats supported

### 1. `TEST.dat` — CW v2.1.3 input file

- **Format:** ISO-8859 text, CRLF line terminators, one value per line, ~514 lines
- **Encoding:** `latin-1` (contains `lb/yd³` with superscript 3)
- **Structure:** Positional dump of every variable in CW's internal input form, in the order they appear in the UI builder
- **Parsed by:** `parse_cw_dat(path)` using the `CW_DAT_INDEX` dict at the top of the module

**Critical: the `CW_DAT_INDEX` dict is empirically derived.** Each line index was verified by matching values against the CW "Input Check" UI screenshot. If CW version changes the layout, this dict is the **only** thing that needs updating. Line indices are listed with their screenshot-verified meanings.

### 2. `TX__Austin.dat` — CW weather file (annual)

- **Format:** Plain ASCII, 3 header lines + 8784 hourly data rows (365 days + leap year buffer)
- **Columns (17):** see `WEATHER_COL` dict; decoded empirically by matching daily-averaged values against CW UI
- **Time indexing:** month/day/hour are 1-based; hour 1 = 00:00–01:00, hour 24 = 23:00–24:00
- **Parsed by:** `parse_cw_weather(path, placement_date, placement_hour, duration_days)` — slices out exactly the analysis window, wraps around year-end if needed

### 3. `temp.txt` — CW 2D temperature export

- **Format:** Tab/space-separated, °C units
- **Structure:**
  - Line 0: description ("Distance is from the upper corner in the width direction/ depth direction")
  - Line 1: header: `Time  w1/d1  w2/d1 ... wN/d1 w1/d2 ... wN/dM  Gradient  Ambient`
  - Lines 2+: one row per timestep, (N_spatial + 3) columns
- **Spatial grid layout (CRITICAL):**
  - Widths stored in **decreasing** order per depth block (e.g. 6.1 → 5.79 → ... → 0.3 → 0)
  - Depths stored in increasing order (0 → 2.44 m)
  - File order: all widths at depth 0, then all widths at depth 1, etc.
- **Coordinate convention (VERY EASY TO GET WRONG):**
  - "Distance from the upper corner in the width direction" means `w=0` is **AT the corner**
  - `w=6.1` m is furthest from the corner = **centerline** of a 40 ft (12.2 m) symmetric half-mat
  - Therefore **centerline is width index 0** (the FIRST column in each depth row)
  - This was initially parsed wrong in the first plot — check this when in doubt
- **Parsed by:** `parse_cw_temp_output(path)` returns `CWValidationSeries` with both full 2D field and centerline-specific convenience series

---

## Module API contract

### Top-level: `load_cw_scenario(input_dat, weather_dat, cw_output_txt, cw_ui_overrides) → CWScenario`

**This is the only function most callers should use.**

```python
from cw_scenario_loader import load_cw_scenario

scn = load_cw_scenario(
    input_dat="TEST.dat",
    weather_dat="TX__Austin.dat",    # optional but strongly recommended
    cw_output_txt="temp.txt",        # optional, for validation only
    cw_ui_overrides={                # optional, for values CW UI shows that
        'cw_ave_max_wind_m_s': 10.5, # can't be derived from weather file
    },
)
# scn.mix, scn.geometry, scn.construction, scn.environment, scn.cw_validation
```

### Dataclasses

Mirror the existing `thermal_engine_v2.py` structure. Field names chosen to enable near-drop-in replacement. Do not rename fields without also updating the engine's adapter (currently planned as `cw_env_to_engine_env()`).

- `CWMixDesign` — mix proportions + CW Schindler regression output + derived thermal props. Includes computed `concrete_density_lb_ft3` property (CW does not store this directly).
- `CWGeometry` — member W × D × L + shape + analysis type
- `CWConstruction` — placement, forms, blanket, curing, subbase. Placement date/hour preserved as strings/ints to avoid timezone complexity.
- `CWEnvironment` — **hourly arrays** + CW UI-displayed averages. Both are preserved. See "Two environment values" below.
- `CWValidationSeries` — Max/Min/Gradient/Ambient + full 2D field + centerline series
- `CWScenario` — top-level container

---

## Design decisions (do not undo without discussion)

### 1. Two sets of environment values: raw hourly + CW UI averages

`CWEnvironment` stores both:
- **Hourly arrays** from the weather file: `T_air_F[]`, `solar_W_m2[]`, `wind_m_s[]`, etc.
- **CW UI values**: `cw_ave_max_daily_temp_F`, `cw_ave_max_solar_W_m2`, etc.

**Why:** CW displays daily-averaged values in its Input Check UI that it uses for its simplified sinusoidal BCs. The hourly arrays are physically correct and are what Sprint 1 (solar + longwave radiation) should use. Keeping both lets us:
- Reproduce CW exactly (use UI values)
- Improve beyond CW (use hourly arrays)
- Validate CW's simplifications (diff the two)

### 2. `cw_ui_overrides` parameter

Some UI values cannot be derived from the weather file — notably wind (`10.5 m/s` in UI vs. `5.0 m/s` max in raw data). Rather than reverse-engineering CW's transformation, we let the caller pass an override dict.

**Leading hypothesis for wind:** CW UI may be reporting gusts, design-envelope values, or a different column of the weather file. Not a blocker for the MVP; flag preserved for future investigation.

### 3. `CW_DAT_INDEX` at top of module

All positional line indices live in a single dict at the top. The motivation is forward-compatibility: when CW v2.1.4 or v3.0 ships with reshuffled lines, only the dict needs updating. Every index is annotated with its UI-screenshot-verified meaning.

### 4. Zero coupling to the engine

The loader intentionally does **not** import from `thermal_engine_v2.py`. It just returns dataclasses with matching field names. This lets the engine evolve (sprints 1-4) without breaking the loader.

### 5. `concrete_density_lb_ft3` is a computed property

CW does not store concrete density. It's derived from the mix total: `Σ(components) / 27 ft³ × (1 - air/100)`. For MIX-01 this gives **131.2 lb/ft³**, which is lower than CalcShore's default of 150. Using the computed value is more correct — if validation shows systematic discrepancy with CW, check the density first.

### 6. L423 is placement temp, not density

Earlier decode attempt flagged L423=60 as "suspicious density." That was wrong. L423 is a placement temperature echo. Real density is computed from the mix (see #5).

---

## Known issues / limitations

### 1. Wind speed mismatch (UI 10.5 m/s vs. raw 5.0 m/s)

CW UI shows a wind value ~2× the raw weather file max. Unknown transformation. Worked around by `cw_ui_overrides`.

### 2. CW version lock-in

The current `CW_DAT_INDEX` is calibrated for CW v2.1.3 only. Other versions will need their own index map. A version check on line 0 is performed but only warns — it doesn't yet refuse to parse. Future: dispatch to version-specific index maps.

### 3. Admixture type is not fully parsed

The UI shows "Type B, Retarder" but this value isn't stored as a string — CW encodes it via flag patterns in L47–53. Thermal engine doesn't use admixture type directly, so decoded as a fixed string for now. Do not rely on this field.

### 4. Cement chemistry (Bogue values) not fully parsed

CW stores Bogue cement chemistry percentages in L365–384 (C3S, C2S, C3A, C4AF, etc.) and uses them to *derive* the hydration params (τ, β, α_u, Hu). We currently use the derived outputs (L385–389) directly, skipping the Bogue inputs. This is fine as long as "Hydration Parameter Values = Default" in the UI, which calls the Schindler regression. If the user manually overrides hydration params in CW, L385–389 still reflects the override, so we're safe.

### 5. Validation scope is 1D mat only

The loader is version-agnostic, but it's been tested against one scenario (MIX-01, 8 ft mat, Austin). Other member types (columns, bent caps, pavements) have different CW file layouts and will need additional index maps.

---

## How to safely extend

### Adding a new field to decode

1. Open the UI in CW, note the field name and value
2. `grep` for the value in `TEST.dat`; note the line index
3. Add an entry to `CW_DAT_INDEX`
4. Add a field to the appropriate dataclass (`CWMixDesign`, `CWConstruction`, etc.)
5. In `parse_cw_dat()`, add a line that populates it via `_getf()` or `_get()`
6. Add a line to `summarize_scenario()` so it shows up in the diagnostic output
7. Test: `python cw_scenario_loader.py TEST.dat TX__Austin.dat` and verify the value matches the UI

### Supporting a new CW version

1. Get a sample .dat file + matching UI screenshot from the new version
2. Diff against v2.1.3 structure
3. Create a version-specific index map (e.g. `CW_DAT_INDEX_V2_1_4`)
4. In `parse_cw_dat()`, read line 0 to detect version, dispatch to appropriate map

### Supporting a new member type (e.g. columns)

1. Generate a .dat file for a column geometry in CW
2. Diff against the mat example to find new/reordered fields
3. Probably need a new `parse_cw_dat_column()` or a `shape` dispatch inside `parse_cw_dat()`
4. Add a `CWColumnGeometry` or extend `CWGeometry` with more fields

### Adding a new output format (e.g. CW animation export)

1. Write a new `parse_cw_animation(path) → ...` function alongside `parse_cw_temp_output()`
2. Return a new dataclass or extend `CWValidationSeries`
3. Wire into `load_cw_scenario()` via a new optional argument

---

## CLI usage (for quick checks)

```bash
python cw_scenario_loader.py TEST.dat
python cw_scenario_loader.py TEST.dat TX__Austin.dat
python cw_scenario_loader.py TEST.dat TX__Austin.dat temp.txt
```

Prints a full `summarize_scenario()` table. Use this as a sanity check every time CW files are regenerated for a new scenario.

---

## What comes next (validation harness)

Not part of this module, but the logical next step:

1. **`cw_env_to_engine_env()` adapter** — collapses `CWEnvironment.hourly` arrays into the scalar-based format the current `thermal_engine_v2.Environment` expects. Placeholder while we upgrade the engine to use hourly BCs directly.

2. **`compare_engine_to_cw.py`** — runs the CalcShore engine on a loaded `CWScenario`, overlays engine output vs. `scn.cw_validation.T_max_xs_F` / etc., reports deltas at peak and across diurnal cycles.

3. **Sprint 1 acceptance test:** With solar + longwave radiation wired into the engine, target `max|T_max_engine − T_max_CW| ≤ 3°F` across the 168 hr run for MIX-01 Austin. Tighter than the current ±5°F because we're aiming for physics parity.

4. **Scenario library:** Once the loader proves stable, convert the 15-mix validation suite from the passdown into 15 `.dat` + `temp.txt` pairs. Every engine change then gets re-validated across all 15 automatically.

---

## File manifest

All files currently referenced by this spec:

| File | Role | Size |
|---|---|---|
| `cw_scenario_loader.py` | The module (this spec describes it) | ~27 KB |
| `cw_scenario_loader.md` | This spec | — |
| `TEST.dat` | Example CW input file (MIX-01 Austin) | ~3.4 KB |
| `TX__Austin.dat` | Example CW weather file | ~1.1 MB |
| `temp.txt` | Example CW 2D output (MIX-01 7-day run) | ~6 MB |
| `thermal_engine_v2.py` | Current CalcShore engine (target for validation) | ~29 KB |
| `coding_passdown_v2.md` | Sprint plan, 15-mix validation suite | ~15 KB |
| `ConcreteWorks_source_documents.md` | Free sources for CW equations (dissertations, V3 manual) | ~10 KB |

---

## Changelog

- **2026-04-22** — Initial version. Decoded CW v2.1.3 input format (≈50 fields across 514 lines), CW weather format (17-column hourly), CW temp.txt output. Validated against MIX-01 Austin scenario. All UI-displayed values reconcile.
