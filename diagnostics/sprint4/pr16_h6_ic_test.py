"""Sprint 4 PR 16 — H6 two-tier IC-propagation test.

H6a: Override MIX-15 placement_temp_F=60 (warm IC). Compare synthetic output to
     MIX-01 CW reference. Tests joint cold-IC effect: concrete grid AND blanket
     pin both swap to 60°F equivalent.

H6b: Run MIX-15 unchanged (placement_temp_F=45) but patch blanket-cell pinning so
     blanket cells track current ambient T_amb instead of T_initial_C. Compare to
     MIX-15 CW reference (concrete IC is unchanged so cold-placed CW is valid).

Verdict logic (on FieldRMS relative to baseline 4.06°F):
  H6a FieldRMS ≤ 1.5°F                 → cold-IC anywhere is the trigger
  H6b FieldRMS ≤ 1.5°F                 → blanket-pin specifically is the bug (Decision E)
  H6a moves metric, H6b doesn't        → concrete-IC propagation, not blanket
  Neither moves FieldRMS significantly  → H6 falsified, escalate to H4 grep audit

Usage (from repo root):
    python diagnostics/sprint4/pr16_h6_ic_test.py
"""
import sys
import importlib.util
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CW_EXPORTS = REPO / "validation" / "cw_exports"
OUT_MD = Path(__file__).parent / "pr16_h6_ic_test.md"

sys.path.insert(0, str(REPO))

from cw_scenario_loader import load_cw_scenario  # noqa: E402
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d  # noqa: E402

# Gate tolerances — must match compare_to_cw.py exactly
TOL_PEAK_MAX_F    = 1.0
TOL_PEAK_GRAD_F   = 2.0
TOL_FIELD_RMS_F   = 2.0
TOL_CENTER_RMS_F  = 1.0
TOL_CORNER_RMS_F  = 3.0
T_START_RMS_HR    = 48.0

# H6 verdict thresholds
FIELD_RMS_BASELINE  = 4.06   # MIX-15 baseline (from Sprint 4 Round 1)
FIELD_RMS_FIXED     = 1.5    # H6 "resolved" threshold per prompt
FIELD_RMS_MOVED     = 2.5    # "significant movement" threshold for partial verdict


# ---------------------------------------------------------------------------
# H6b patch — modified engine where blanket tracks ambient
# ---------------------------------------------------------------------------

def _load_h6b_engine():
    """Return a module equivalent to thermal_engine_2d but with the blanket
    pin changed to track current ambient temperature instead of T_initial_C.

    Original line 2087-2088:
        # When blanket is pure-R (full_2d or skip_blanket_node), pin it to initial value.
        if _use_pure_r_blanket:
            T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]

    Patched version:
        # H6b: blanket tracks current ambient T_amb (diagnostic patch).
        if _use_pure_r_blanket:
            T_new[grid.is_blanket] = _T_amb_C
    """
    engine_path = REPO / "thermal_engine_2d.py"
    source = engine_path.read_text(encoding="utf-8")

    old = (
        "        # When blanket is pure-R (full_2d or skip_blanket_node), pin it to initial value.\n"
        "        if _use_pure_r_blanket:\n"
        "            T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]"
    )

    new = (
        "        # H6b diagnostic: blanket tracks current ambient T_amb instead of T_initial_C.\n"
        "        if _use_pure_r_blanket:\n"
        "            T_new[grid.is_blanket] = _T_amb_C"
    )

    assert old in source, "Could not locate blanket-pin target in thermal_engine_2d.py — check line numbers"

    patched_source = source.replace(old, new, 1)

    # Write to a temp file; import under a unique name so it doesn't shadow the real module.
    tmp = Path(tempfile.mktemp(suffix="_h6b_engine.py"))
    tmp.write_text(patched_source, encoding="utf-8")

    spec = importlib.util.spec_from_file_location("thermal_engine_2d_h6b", str(tmp))
    mod  = importlib.util.module_from_spec(spec)
    sys.modules["thermal_engine_2d_h6b"] = mod  # required for @dataclass to resolve the module
    spec.loader.exec_module(mod)

    tmp.unlink()  # clean up temp file
    return mod


# ---------------------------------------------------------------------------
# Metric computation (accepts separate run scenario + reference scenario)
# ---------------------------------------------------------------------------

def rms_arr(a, b):
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(diff ** 2)))


def compute_metrics(run_scn, ref_scn, grid, result):
    """Compute 5 gate metrics: run_scn provides geometry/grid context;
    ref_scn.cw_validation provides the CW reference data."""
    jslice, islice = grid.concrete_slice()
    T_conc_C = result.T_field_C[:, jslice, islice]
    T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0

    engine_peak_max_F  = float(T_conc_F.max())
    grad_series        = T_conc_F.max(axis=(1, 2)) - T_conc_F.min(axis=(1, 2))
    engine_peak_grad_F = float(grad_series.max())

    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_center_F = result.T_field_C[:, iy_mid, grid.nx - 1] * 9.0 / 5.0 + 32.0

    iy_top    = grid.iy_concrete_start
    ix_corner = grid.ix_concrete_start
    eng_corner_F = result.T_field_C[:, iy_top, ix_corner] * 9.0 / 5.0 + 32.0

    val = ref_scn.cw_validation
    cw_peak_max_F  = float(val.T_max_xs_F.max())
    cw_peak_grad_F = float(val.T_diff_xs_F.max())

    peak_max_delta  = engine_peak_max_F  - cw_peak_max_F
    peak_grad_delta = engine_peak_grad_F - cw_peak_grad_F

    if val.T_field_F is None:
        return {
            "peak_max_delta":  peak_max_delta,
            "peak_grad_delta": peak_grad_delta,
            "field_rms":  None,
            "center_rms": None,
            "corner_rms": None,
            "s0_pass": (
                (1 if abs(peak_max_delta)  < TOL_PEAK_MAX_F  else 0) +
                (1 if abs(peak_grad_delta) < TOL_PEAK_GRAD_F else 0)
            ),
        }

    cw_t_s = val.time_hrs * 3600.0
    n_cw_t, n_cw_d, _ = val.T_field_F.shape

    cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]
    cw_corner_F = val.T_field_F[:, 0, -1]

    eng_center_interp = np.interp(cw_t_s, result.t_s, eng_center_F)
    eng_corner_interp = np.interp(cw_t_s, result.t_s, eng_corner_F)

    rms_mask     = val.time_hrs >= T_START_RMS_HR
    center_rms_F = rms_arr(eng_center_interp[rms_mask], cw_center_F[rms_mask])
    corner_rms_F = rms_arr(eng_corner_interp[rms_mask], cw_corner_F[rms_mask])

    sq_errs = []
    for ti in range(n_cw_t):
        if val.time_hrs[ti] < T_START_RMS_HR:
            continue
        idx = int(np.searchsorted(result.t_s, cw_t_s[ti]))
        idx = min(idx, len(result.t_s) - 1)
        T_eng_col_C = result.T_field_C[idx, jslice, grid.nx - 1]
        T_eng_col_F = T_eng_col_C * 9.0 / 5.0 + 32.0
        n_cmp = min(len(T_eng_col_F), n_cw_d)
        sq_errs.extend((T_eng_col_F[:n_cmp] - val.T_field_F[ti, :n_cmp, 0]).tolist())
    field_rms_F = float(np.sqrt(np.mean(np.array(sq_errs) ** 2)))

    s0 = {
        "peak_max":  abs(peak_max_delta)  < TOL_PEAK_MAX_F,
        "peak_grad": abs(peak_grad_delta) < TOL_PEAK_GRAD_F,
        "field":     field_rms_F  < TOL_FIELD_RMS_F,
        "center":    center_rms_F < TOL_CENTER_RMS_F,
        "corner":    corner_rms_F < TOL_CORNER_RMS_F,
    }
    return {
        "peak_max_delta":  peak_max_delta,
        "peak_grad_delta": peak_grad_delta,
        "field_rms":  field_rms_F,
        "center_rms": center_rms_F,
        "corner_rms": corner_rms_F,
        "s0_pass":    sum(s0.values()),
        "s0":         s0,
    }


def fmt_delta(v):
    return f"{v:+.2f}" if v is not None else "N/A"

def fmt_rms(v):
    return f"{v:.2f}" if v is not None else "N/A"


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------

def run_standard(scn):
    """Run engine with scn as-is (uses scn.construction.placement_temp_F)."""
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )
    return result, grid


def run_h6b_patched(scn, h6b_engine):
    """Run engine with blanket pin patched to track ambient instead of T_initial."""
    grid = h6b_engine.build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    result = h6b_engine.solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )
    return result, grid


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading scenarios ...")
    scn_mix15 = load_cw_scenario(
        str(CW_EXPORTS / "MIX-15" / "input.dat"),
        str(CW_EXPORTS / "MIX-15" / "weather.dat"),
        str(CW_EXPORTS / "MIX-15" / "output.txt"),
    )
    scn_mix01 = load_cw_scenario(
        str(CW_EXPORTS / "MIX-01" / "input.dat"),
        str(CW_EXPORTS / "MIX-01" / "weather.dat"),
        str(CW_EXPORTS / "MIX-01" / "output.txt"),
    )
    print(f"  MIX-15 placement_temp_F = {scn_mix15.construction.placement_temp_F}°F")
    print(f"  MIX-01 placement_temp_F = {scn_mix01.construction.placement_temp_F}°F")
    print()

    # --- H6a: override placement_temp_F = 60.0 on MIX-15 ---
    print("--- H6a: MIX-15 with placement_temp_F overridden to 60°F ---")
    import dataclasses
    scn_mix15_warm = dataclasses.replace(
        scn_mix15,
        construction=dataclasses.replace(scn_mix15.construction, placement_temp_F=60.0),
    )
    print(f"  Running engine (placement_temp_F = {scn_mix15_warm.construction.placement_temp_F}°F) ...")
    result_h6a, grid_h6a = run_standard(scn_mix15_warm)
    assert np.all(np.isfinite(result_h6a.T_field_C)), "H6a: NaN in engine output"
    m_h6a = compute_metrics(scn_mix15_warm, scn_mix01, grid_h6a, result_h6a)
    print(f"  Done. S0={m_h6a['s0_pass']}/5  PeakMaxΔ={m_h6a['peak_max_delta']:+.2f}°F  "
          f"FieldRMS={fmt_rms(m_h6a['field_rms'])}°F  "
          f"CenterRMS={fmt_rms(m_h6a['center_rms'])}°F")
    print()

    # --- H6b: MIX-15 cold IC, blanket pin patched to track ambient ---
    print("--- H6b: MIX-15 unchanged, blanket pin patched to track T_amb ---")
    print("  Loading patched engine module ...")
    h6b_engine = _load_h6b_engine()
    print(f"  Running engine (placement_temp_F = {scn_mix15.construction.placement_temp_F}°F, blanket→T_amb) ...")
    result_h6b, grid_h6b = run_h6b_patched(scn_mix15, h6b_engine)
    assert np.all(np.isfinite(result_h6b.T_field_C)), "H6b: NaN in engine output"
    m_h6b = compute_metrics(scn_mix15, scn_mix15, grid_h6b, result_h6b)
    print(f"  Done. S0={m_h6b['s0_pass']}/5  PeakMaxΔ={m_h6b['peak_max_delta']:+.2f}°F  "
          f"FieldRMS={fmt_rms(m_h6b['field_rms'])}°F  "
          f"CenterRMS={fmt_rms(m_h6b['center_rms'])}°F")
    print()

    # Also get baseline MIX-15 (no changes) for reference
    print("--- Baseline: MIX-15 unmodified ---")
    result_base, grid_base = run_standard(scn_mix15)
    m_base = compute_metrics(scn_mix15, scn_mix15, grid_base, result_base)
    print(f"  Done. S0={m_base['s0_pass']}/5  PeakMaxΔ={m_base['peak_max_delta']:+.2f}°F  "
          f"FieldRMS={fmt_rms(m_base['field_rms'])}°F")
    print()

    # --- Verdict logic ---
    frs_h6a = m_h6a["field_rms"]
    frs_h6b = m_h6b["field_rms"]

    if frs_h6a is not None and frs_h6a <= FIELD_RMS_FIXED:
        if frs_h6b is not None and frs_h6b <= FIELD_RMS_FIXED:
            verdict = "BLANKET-PIN-CONFIRMED"
            detail  = ("Both H6a and H6b resolve FieldRMS. Blanket-pin is the bug; "
                       "concrete-IC effect is secondary.")
        else:
            verdict = "CONCRETE-IC-TRIGGER"
            detail  = ("H6a resolves FieldRMS but H6b does not. The concrete grid IC "
                       "(equivalent-age delay), not the blanket pin, is the primary trigger.")
    elif frs_h6b is not None and frs_h6b <= FIELD_RMS_FIXED:
        verdict = "BLANKET-PIN-CONFIRMED"
        detail  = ("H6b alone resolves FieldRMS. The blanket-pin-at-cold-IC is the bug. "
                   "Concrete IC effect is not the primary driver.")
    elif frs_h6a is not None and frs_h6a <= FIELD_RMS_MOVED:
        verdict = "H6-PARTIAL"
        detail  = ("H6a partially moves FieldRMS but does not fully resolve it. "
                   "IC is a contributing factor; other mechanisms also involved.")
    elif frs_h6b is not None and frs_h6b <= FIELD_RMS_MOVED:
        verdict = "H6-PARTIAL"
        detail  = ("H6b partially moves FieldRMS but does not fully resolve it. "
                   "Blanket-pin is a contributing factor.")
    else:
        verdict = "H6-FALSIFIED"
        detail  = ("Neither H6a nor H6b significantly moves FieldRMS. "
                   "IC propagation is not the primary trigger. Escalate to H4 grep audit.")

    # --- Write markdown ---
    md = [
        "# PR 16 — H6 IC Propagation Test (MIX-15 B2)",
        "",
        "Phase 2 diagnostic for Sprint 4 PR 16. Two-tier test separating the joint cold-IC",
        "effect (H6a) from the blanket-pin-specifically (H6b).",
        "",
        "## Test Design",
        "",
        "| Variant | IC | Blanket pin | Reference |",
        "| --- | --- | --- | --- |",
        "| Baseline | MIX-15 cold (45°F) | T_initial (45°F) | MIX-15 CW |",
        "| **H6a** | **Warm (60°F)** | T_initial (60°F) | **MIX-01 CW** |",
        "| **H6b** | MIX-15 cold (45°F) | **T_amb (current)** | MIX-15 CW |",
        "",
        "H6a tests joint cold-IC effect — both the concrete grid IC and the blanket pin swap",
        "to 60°F-equivalent. Compare to MIX-01 CW (warm-placed reference, same construction).",
        "",
        "H6b isolates the blanket-pin contribution — concrete IC stays cold (45°F), only the",
        "blanket-cell temperature policy changes. Compare to MIX-15 CW (cold-placed; concrete",
        "IC is unchanged so the cold-placed CW is the correct reference).",
        "",
        "## Results",
        "",
        "### Side-by-Side Metric Table",
        "",
        "| Metric | Baseline (MIX-15) | H6a (warm IC) | H6b (blanket→T_amb) |",
        "| --- | ---: | ---: | ---: |",
        f"| PeakMax Δ (°F) | {fmt_delta(m_base['peak_max_delta'])} | {fmt_delta(m_h6a['peak_max_delta'])} | {fmt_delta(m_h6b['peak_max_delta'])} |",
        f"| PeakGrad Δ (°F) | {fmt_delta(m_base['peak_grad_delta'])} | {fmt_delta(m_h6a['peak_grad_delta'])} | {fmt_delta(m_h6b['peak_grad_delta'])} |",
        f"| FieldRMS (°F) | {fmt_rms(m_base['field_rms'])} | {fmt_rms(m_h6a['field_rms'])} | {fmt_rms(m_h6b['field_rms'])} |",
        f"| CenterRMS (°F) | {fmt_rms(m_base['center_rms'])} | {fmt_rms(m_h6a['center_rms'])} | {fmt_rms(m_h6b['center_rms'])} |",
        f"| CornerRMS (°F) | {fmt_rms(m_base['corner_rms'])} | {fmt_rms(m_h6a['corner_rms'])} | {fmt_rms(m_h6b['corner_rms'])} |",
        f"| S0 pass | {m_base['s0_pass']}/5 | {m_h6a['s0_pass']}/5 | {m_h6b['s0_pass']}/5 |",
        "",
        f"Reference for H6a: MIX-01 CW (warm-placed). Reference for H6b/Baseline: MIX-15 CW (cold-placed).",
        "",
        "### FieldRMS Movement",
        "",
        f"Baseline FieldRMS: {fmt_rms(m_base['field_rms'])}°F",
        f"H6a FieldRMS:      {fmt_rms(frs_h6a)}°F  "
        f"({'RESOLVED ≤1.5°F' if frs_h6a and frs_h6a <= FIELD_RMS_FIXED else 'moved but not resolved' if frs_h6a and frs_h6a <= FIELD_RMS_MOVED else 'NOT MOVED' if frs_h6a else 'N/A'})",
        f"H6b FieldRMS:      {fmt_rms(frs_h6b)}°F  "
        f"({'RESOLVED ≤1.5°F' if frs_h6b and frs_h6b <= FIELD_RMS_FIXED else 'moved but not resolved' if frs_h6b and frs_h6b <= FIELD_RMS_MOVED else 'NOT MOVED' if frs_h6b else 'N/A'})",
        "",
        f"**H6 verdict: {verdict}**",
        "",
        f"{detail}",
        "",
    ]

    if verdict == "BLANKET-PIN-CONFIRMED":
        md += [
            "**Proposed Phase 3 decision: E** — single localized term identified.",
            "Fix: change `T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]` to track",
            "current ambient temperature. The fix must preserve the original pin rationale",
            "(air cells have placeholder rho_cp=1 and must remain pinned).",
            "",
        ]
    elif verdict == "CONCRETE-IC-TRIGGER":
        md += [
            "**Proposed Phase 3 decision: F (via H4 audit)** — blanket pin is not the bug;",
            "concrete IC propagation is the trigger. Proceed to H4 grep audit to locate the",
            "specific concrete-IC-dependent term.",
            "",
        ]
    elif verdict == "H6-PARTIAL":
        md += [
            "**Proposed Phase 3 decision: G** — IC is a contributing factor but not the sole",
            "trigger. Run H4 grep audit to find additional mechanism(s). PR 16 closes with",
            "whatever clean partial fix is identified; remainder routes to Sprint 5/6.",
            "",
        ]
    else:
        md += [
            "**Proposed Phase 3 decision: F** — H6 falsified. Proceed to H4 grep audit.",
            "",
        ]

    md += [
        "## Blanket-Pin Mechanism (§7.6.4 Documentation)",
        "",
        "This observation is recorded regardless of H6 outcome as a candidate Sprint 6 cleanup item.",
        "",
        "### What the pin does",
        "",
        "In `thermal_engine_2d.py`, after every timestep in all boundary modes, the engine",
        "applies two unconditional overrides (`thermal_engine_2d.py:2082-2088`):",
        "",
        "```python",
        "# Pin air cells to initial value for ALL modes (air has rho_cp=1 placeholder;",
        "# leaving them at the BC-updated value would corrupt next-step flux calcs).",
        "if grid.is_air.any():",
        "    T_new[grid.is_air] = T_initial_C[grid.is_air]",
        "# When blanket is pure-R (full_2d or skip_blanket_node), pin it to initial value.",
        "if _use_pure_r_blanket:",
        "    T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]",
        "```",
        "",
        "The **air-cell pin** is necessary: air cells have a placeholder `rho_cp=1` and are",
        "not physically evolved. Resetting them each step prevents numerical drift from",
        "corrupting flux calculations at adjacent concrete cells.",
        "",
        "The **blanket-cell pin** (line 2088) applies the same pattern to the insulating",
        "blanket layer. In a pure-R blanket model, the blanket is a thermal resistance with",
        "no thermal mass, so it has no physical temperature to track. The pin holds blanket",
        "cells at `T_initial_C` — the concrete placement temperature — for the entire 168-hour",
        "simulation.",
        "",
        "### Why this is potentially problematic for cold-placed mixes",
        "",
        "For Reference (MIX-01, placement_temp_F=60°F), `T_initial_C` ≈ 15.6°C. The ambient",
        "weather fluctuates but averages near this range during the simulation. The blanket",
        "pin is therefore approximately consistent with ambient conditions.",
        "",
        "For MIX-15 (placement_temp_F=45°F), `T_initial_C` ≈ 7.2°C. As the concrete cures",
        "and ambient warms (e.g., to 25–30°C/77–86°F daytime), the blanket cells remain",
        "frozen at 7.2°C throughout the run. This creates a persistent cold boundary",
        "condition at the top surface that suppresses heat retention in the concrete — a",
        "systematic downward bias on FieldRMS, CenterRMS, and PeakMax for cold-placed pours.",
        "",
        "### Magnitude of effect",
        "",
        "This H6 diagnostic quantifies the blanket-pin contribution. See results table above.",
        "",
        "### Sprint 6 scope note",
        "",
        "Even if this bug is fixed in PR 16, the air-cell pin at line 2085 applies the same",
        "T_initial_C pattern to air cells. Air cells are non-physical (placeholder rho_cp),",
        "so the pin is correct for them — they should not evolve thermally. The blanket is",
        "different: it represents a physical material with a real thermal interaction with the",
        "ambient environment. A Sprint 6 cleanup could replace the blanket pin with a",
        "quasi-steady boundary equation that accounts for blanket thermal resistance without",
        "freezing blanket temperature at placement IC.",
        "",
    ]

    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote: {OUT_MD}")
    print()

    # --- Stdout summary ---
    print("=" * 70)
    print("PHASE 2 H6 SUMMARY")
    print("=" * 70)
    print()
    print(f"Baseline  MIX-15 (cold IC, T_initial blanket):  PeakMaxΔ={fmt_delta(m_base['peak_max_delta'])}  "
          f"FieldRMS={fmt_rms(m_base['field_rms'])}  CenterRMS={fmt_rms(m_base['center_rms'])}")
    print(f"H6a       warm IC 60°F, vs MIX-01 CW:           PeakMaxΔ={fmt_delta(m_h6a['peak_max_delta'])}  "
          f"FieldRMS={fmt_rms(m_h6a['field_rms'])}  CenterRMS={fmt_rms(m_h6a['center_rms'])}")
    print(f"H6b       cold IC, blanket→T_amb, vs MIX-15 CW: PeakMaxΔ={fmt_delta(m_h6b['peak_max_delta'])}  "
          f"FieldRMS={fmt_rms(m_h6b['field_rms'])}  CenterRMS={fmt_rms(m_h6b['center_rms'])}")
    print()
    print(f"FieldRMS movement: baseline {fmt_rms(m_base['field_rms'])}  →  H6a {fmt_rms(frs_h6a)}  /  H6b {fmt_rms(frs_h6b)}")
    print()
    print(f"H6 VERDICT: {verdict}")
    print(f"  {detail}")
    print()
    if verdict == "BLANKET-PIN-CONFIRMED":
        print("Proposed Phase 3: Decision E — fix blanket pin in thermal_engine_2d.py:2088.")
        print("Change: T_new[grid.is_blanket] = T_initial_C[grid.is_blanket]")
        print("    To: T_new[grid.is_blanket] = _T_amb_C")
        print("Validate: run_all --group evaluation_set; pytest; check bit-identical test.")
    elif verdict == "CONCRETE-IC-TRIGGER":
        print("Proposed Phase 3: Decision F — run H4 grep audit to locate IC-dependent term.")
    elif verdict == "H6-PARTIAL":
        print("Proposed Phase 3: Decision G — partial fix plus H4 audit; route remainder Sprint 5/6.")
    else:
        print("Proposed Phase 3: Decision F — H6 falsified; run H4 grep audit.")
    print()
    print("WAIT — review results, confirm Phase 3 decision before proceeding.")


if __name__ == "__main__":
    main()
