# CHECK 1.1 — Call chain trace

Trace from `run_all.py` to the solver's `Hu` read, identifying every site where `CWMixDesign` is constructed, where the solver is called with a mix, and where Hu fields are mutated.

## Constructor sites (`CWMixDesign` instantiation)

The dataclass `CWMixDesign` is constructed in exactly one place in the production loader path:

- **`cw_scenario_loader.py:414`** — inside `parse_cw_dat()`, called by `load_cw_scenario()`. Constructed with `Hu_J_kg=_getf(CW_DAT_INDEX['Hu_J_kg'])` (the verbatim CW-regressed value), and all other composition fields from the CW input.dat. `Hu_factor_calibrated` and `Hu_J_kg_effective` are NOT passed at construction; they take their dataclass defaults (1.0 and 0.0 respectively).

Other constructor invocations exist only in test fixtures and standalone diagnostic scripts (`plot_engine_manual_mix.py`, `batch_compare_all_mixes.py`); none are reachable from `run_all.py` → `run_one`.

## Mutation sites for Hu fields

Three mutation sites exist for the Hu trio (`Hu_J_kg`, `Hu_J_kg_effective`, `Hu_factor_calibrated`):

1. **`cw_scenario_loader.py:166-170` (`__post_init__`)** — runs at `CWMixDesign` construction time. If `Hu_J_kg_effective == 0.0` (the default), it is set to `Hu_J_kg`. This is the "raw-Hu fallback": ensures any caller that constructs a mix without setting `Hu_J_kg_effective` gets a sane default that the solver can consume.

2. **`cw_scenario_loader.py:819-824` (in `load_cw_scenario`)** — the calibration application. Conditional on `use_cw_calibrated_hu=True` (the default). Calls `compute_hu_factor(mix)` from `kinetics_correction`, sets:
    - `mix.Hu_factor_calibrated = factor`
    - `mix.Hu_J_kg_effective = mix.Hu_J_kg * factor`
    - `mix.Hu_calibration_note = note`

3. No other mutation sites exist in `run_all.py`, `compare_to_cw.py`, `thermal_engine_2d.py`, or `cw_scenario_loader.py`. `Hu_J_kg` (the raw regressed value) is never mutated post-construction.

## Solver call sites

The 2D thermal solver is `solve_hydration_2d` from `thermal_engine_2d`. Production call sites that pass a mix:

- **`compare_to_cw.py:178-186`** — inside `run_one()`. Called with `scn.mix` (the mix returned by `load_cw_scenario`). This is the only path reached from `run_all.py`.

The solver reads Hu at:

- **`thermal_engine_2d.py:1438`** — `Hu = mix.Hu_J_kg_effective` (with a comment "J/kg_cement (apr28 calibrated)"). The solver consumes `Hu_J_kg_effective`, never `Hu_J_kg`.

## End-to-end path summary

```
run_all.py:145         run_one(scenario_dir, png_path=...)
  └── compare_to_cw.py:68    def run_one(...)
        └── compare_to_cw.py:161   load_cw_scenario(input.dat, weather.dat, output.txt)
              # ↑ positional args only — no use_cw_calibrated_hu kwarg passed
              └── cw_scenario_loader.py:769   def load_cw_scenario(..., use_cw_calibrated_hu: bool = True)
                    # default True applies
                    ├── parse_cw_dat → CWMixDesign(Hu_J_kg=…) at line 414
                    │     └── __post_init__: Hu_J_kg_effective = Hu_J_kg  (raw-Hu fallback)
                    └── if use_cw_calibrated_hu:    # line 819, default-true branch
                          factor, note = compute_hu_factor(mix)
                          mix.Hu_J_kg_effective = mix.Hu_J_kg * factor   # CALIBRATION APPLIED
        └── compare_to_cw.py:178   solve_hydration_2d(grid, scn.mix, ...)
              └── thermal_engine_2d.py:1438   Hu = mix.Hu_J_kg_effective   # solver reads calibrated
```

**Conclusion**: there is no path from `run_all.py` to the solver that bypasses the calibration. The default-True branch at line 819 is on the only path. The instrumentation in CHECK 1.2 will confirm this empirically.
