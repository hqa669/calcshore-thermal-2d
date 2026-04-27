"""Sprint 4 PR 16 — B2 diagnostic sweeps: F_vert and R_form on MIX-15.

Phase 1 diagnostic — throwaway script, not committed.
Only the output (pr16_b2_ablation.md) is promoted to validation/diagnostics/
on diagnostic-only close.

Mixes: MIX-15 (B2 target) + MIX-01 (Reference control).

Ablation question (binary):
  Does MIX-15 PeakMax Δ vary across the F_vert / R_form sweep,
  or is it invariant the way B1's was?

  ≥0.5°F variation  → boundary-physics-authoritative → Phase 2
  <0.1°F variation  → hydration-routed (B1 precedent) → Decision I
  0.1–0.5°F         → marginal → report, wait for user decision

Sweep 1 — F_vert: overrides construction.vertical_solar_factor, no source edit.
Sweep 2 — R_form: monkey-patches thermal_engine_2d.R_FORM_CONTACT_SI, restored after.

Usage (from repo root):
    python diagnostics/sprint4/pr16_b2_ablation.py
"""
import dataclasses
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CW_EXPORTS = REPO / "validation" / "cw_exports"
OUT_MD = Path(__file__).parent / "pr16_b2_ablation.md"

sys.path.insert(0, str(REPO))

from cw_scenario_loader import load_cw_scenario  # noqa: E402
import thermal_engine_2d  # noqa: E402
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d  # noqa: E402

# Gate tolerances — must match compare_to_cw.py exactly
TOL_PEAK_MAX_F    = 1.0
TOL_PEAK_GRAD_F   = 2.0
TOL_FIELD_RMS_F   = 2.0
TOL_CENTER_RMS_F  = 1.0
TOL_CORNER_RMS_F  = 3.0
T_START_RMS_HR    = 48.0

MIXES = ["MIX-01", "MIX-15"]

F_VERT_SWEEP = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
R_FORM_SWEEP = [0.040, 0.060, 0.0862, 0.100, 0.120, 0.150, 0.200]

# Phase 1 verdict thresholds (on MIX-15 PeakMax Δ variation)
VARIATION_AUTHORITY = 0.5   # ≥ this → boundary-physics-authoritative
VARIATION_MARGINAL  = 0.1   # < this → hydration-routed


def rms_arr(a, b):
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(diff ** 2)))


def compute_metrics(scn, grid, result):
    """Compute the 5 gate metrics, matching compare_to_cw.py exactly."""
    jslice, islice = grid.concrete_slice()
    T_conc_C = result.T_field_C[:, jslice, islice]
    T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0

    engine_peak_max_F = float(T_conc_F.max())
    grad_series = T_conc_F.max(axis=(1, 2)) - T_conc_F.min(axis=(1, 2))
    engine_peak_grad_F = float(grad_series.max())

    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_center_C = result.T_field_C[:, iy_mid, grid.nx - 1]
    eng_center_F = eng_center_C * 9.0 / 5.0 + 32.0

    iy_top    = grid.iy_concrete_start
    ix_corner = grid.ix_concrete_start
    eng_corner_C = result.T_field_C[:, iy_top, ix_corner]
    eng_corner_F = eng_corner_C * 9.0 / 5.0 + 32.0

    val = scn.cw_validation
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
            "s0": {
                "peak_max":  abs(peak_max_delta)  < TOL_PEAK_MAX_F,
                "peak_grad": abs(peak_grad_delta) < TOL_PEAK_GRAD_F,
                "field": None, "center": None, "corner": None,
            },
        }

    cw_t_s = val.time_hrs * 3600.0
    n_cw_t, n_cw_d, _ = val.T_field_F.shape

    cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]
    cw_corner_F = val.T_field_F[:, 0, -1]

    eng_center_interp = np.interp(cw_t_s, result.t_s, eng_center_F)
    eng_corner_interp = np.interp(cw_t_s, result.t_s, eng_corner_F)

    rms_mask = val.time_hrs >= T_START_RMS_HR
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
        "s0_pass": sum(s0.values()),
        "s0": s0,
    }


def run_engine(scn, construction_override=None):
    construction = construction_override if construction_override is not None else scn.construction
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=construction,
        diagnostic_outputs=True,
    )
    return result, grid


def load_scenarios():
    scenarios = {}
    for mix_id in MIXES:
        mix_dir = CW_EXPORTS / mix_id
        scn = load_cw_scenario(
            str(mix_dir / "input.dat"),
            str(mix_dir / "weather.dat"),
            str(mix_dir / "output.txt"),
        )
        scenarios[mix_id] = scn
        print(f"  Loaded {mix_id}")
    return scenarios


# ---------------------------------------------------------------------------
# Sweep 1: F_vert
# ---------------------------------------------------------------------------

def run_fvert_sweep(scenarios):
    print("=== Sweep 1: F_vert (9 values × 2 mixes) ===")
    rows = []
    for mix_id in MIXES:
        scn = scenarios[mix_id]
        for F_vert in F_VERT_SWEEP:
            print(f"  [{mix_id}] F_vert={F_vert:.2f} ... ", end="", flush=True)
            construction = dataclasses.replace(scn.construction, vertical_solar_factor=F_vert)
            result, grid = run_engine(scn, construction_override=construction)
            if not np.all(np.isfinite(result.T_field_C)):
                print("NaN — STOP")
                sys.exit(1)
            m = compute_metrics(scn, grid, result)
            rows.append((mix_id, F_vert, m))
            cr = m["corner_rms"]
            cr_str = f"{cr:.2f}°F" if cr is not None else "N/A"
            print(f"S0={m['s0_pass']}/5  CornerRMS={cr_str}  PeakMaxΔ={m['peak_max_delta']:+.2f}°F")
    return rows


# ---------------------------------------------------------------------------
# Sweep 2: R_form (monkey-patch module constant)
# ---------------------------------------------------------------------------

def run_rform_sweep(scenarios):
    print("=== Sweep 2: R_form (7 values × 2 mixes) ===")
    original = thermal_engine_2d.R_FORM_CONTACT_SI
    rows = []
    try:
        for mix_id in MIXES:
            scn = scenarios[mix_id]
            for R_form in R_FORM_SWEEP:
                print(f"  [{mix_id}] R_form={R_form:.4f} ... ", end="", flush=True)
                thermal_engine_2d.R_FORM_CONTACT_SI = R_form
                result, grid = run_engine(scn)
                if not np.all(np.isfinite(result.T_field_C)):
                    print("NaN — STOP")
                    thermal_engine_2d.R_FORM_CONTACT_SI = original
                    sys.exit(1)
                m = compute_metrics(scn, grid, result)
                rows.append((mix_id, R_form, m))
                cr = m["corner_rms"]
                cr_str = f"{cr:.2f}°F" if cr is not None else "N/A"
                print(f"S0={m['s0_pass']}/5  CornerRMS={cr_str}  PeakMaxΔ={m['peak_max_delta']:+.2f}°F")
    finally:
        thermal_engine_2d.R_FORM_CONTACT_SI = original
    return rows


# ---------------------------------------------------------------------------
# Analysis helpers
# ---------------------------------------------------------------------------

def metric_at(rows, mix_id, val):
    for (mid, v, m) in rows:
        if mid == mix_id and abs(v - val) < 1e-9:
            return m
    return None


def sweep_range(rows, mix_id, metric_key):
    """Return (min, max, range) of a metric across the sweep for one mix."""
    vals = [m[metric_key] for (mid, _, m) in rows if mid == mix_id and m[metric_key] is not None]
    if not vals:
        return None, None, None
    return min(vals), max(vals), max(vals) - min(vals)


def default_metrics(rows, mix_id, default_val):
    """Return metric dict at the default sweep value."""
    return metric_at(rows, mix_id, default_val)


def verdict_str(variation):
    if variation is None:
        return "UNKNOWN (no T_field data)"
    if variation >= VARIATION_AUTHORITY:
        return f"BOUNDARY-PHYSICS-AUTHORITATIVE (variation {variation:.2f}°F ≥ {VARIATION_AUTHORITY:.1f}°F threshold)"
    if variation < VARIATION_MARGINAL:
        return f"HYDRATION-ROUTED (variation {variation:.2f}°F < {VARIATION_MARGINAL:.1f}°F threshold — B1 precedent)"
    return f"MARGINAL (variation {variation:.2f}°F in [{VARIATION_MARGINAL:.1f}, {VARIATION_AUTHORITY:.1f})°F — user decision required)"


# ---------------------------------------------------------------------------
# Markdown builders
# ---------------------------------------------------------------------------

def fmt_delta(v):
    return f"{v:+.2f}" if v is not None else "N/A"

def fmt_rms(v):
    return f"{v:.2f}" if v is not None else "N/A"


def build_table(rows, param_label, fmt_param):
    lines = [
        f"| Mix | {param_label} | PeakMax Δ (°F) | PeakGrad Δ (°F) | FieldRMS (°F) | CenterRMS (°F) | CornerRMS (°F) | S0 pass |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | :---: |",
    ]
    for (mix_id, pval, m) in rows:
        pass_cell = f"**{m['s0_pass']}/5**" if m["s0_pass"] == 5 else f"{m['s0_pass']}/5"
        lines.append(
            f"| {mix_id} | {fmt_param(pval)} "
            f"| {fmt_delta(m['peak_max_delta'])} "
            f"| {fmt_delta(m['peak_grad_delta'])} "
            f"| {fmt_rms(m['field_rms'])} "
            f"| {fmt_rms(m['center_rms'])} "
            f"| {fmt_rms(m['corner_rms'])} "
            f"| {pass_cell} |"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print(f"Loading {len(MIXES)} scenarios ...")
    scenarios = load_scenarios()
    print()

    fvert_rows = run_fvert_sweep(scenarios)
    print()
    rform_rows = run_rform_sweep(scenarios)
    print()

    # --- Analysis: variation ranges for MIX-15 across each sweep ---
    fv_pm_min, fv_pm_max, fv_pm_range = sweep_range(fvert_rows, "MIX-15", "peak_max_delta")
    rf_pm_min, rf_pm_max, rf_pm_range = sweep_range(rform_rows, "MIX-15", "peak_max_delta")

    fv_pg_min, fv_pg_max, fv_pg_range = sweep_range(fvert_rows, "MIX-15", "peak_grad_delta")
    rf_pg_min, rf_pg_max, rf_pg_range = sweep_range(rform_rows, "MIX-15", "peak_grad_delta")

    fv_fr_min, fv_fr_max, fv_fr_range = sweep_range(fvert_rows, "MIX-15", "field_rms")
    rf_fr_min, rf_fr_max, rf_fr_range = sweep_range(rform_rows, "MIX-15", "field_rms")

    fv_cr_min, fv_cr_max, fv_cr_range = sweep_range(fvert_rows, "MIX-15", "corner_rms")
    rf_cr_min, rf_cr_max, rf_cr_range = sweep_range(rform_rows, "MIX-15", "corner_rms")

    # Default metrics (F_vert=0.15, R_form=0.0862)
    dflt_fv  = default_metrics(fvert_rows, "MIX-15", 0.15)
    dflt_rf  = default_metrics(rform_rows, "MIX-15", 0.0862)
    m01_dflt = default_metrics(fvert_rows, "MIX-01", 0.15)

    # Phase 1 verdict: use the larger of the two sweep variations as the authority signal
    max_pm_range = max(
        fv_pm_range if fv_pm_range is not None else 0.0,
        rf_pm_range if rf_pm_range is not None else 0.0,
    )

    if fv_pm_range is None and rf_pm_range is None:
        authority_verdict = "UNKNOWN (no T_field data)"
    elif max_pm_range >= VARIATION_AUTHORITY:
        authority_verdict = "BOUNDARY-PHYSICS-AUTHORITATIVE"
    elif max_pm_range < VARIATION_MARGINAL:
        authority_verdict = "HYDRATION-ROUTED"
    else:
        authority_verdict = "MARGINAL"

    mix01_gate = m01_dflt["s0_pass"] == 5 if m01_dflt else False

    # --- Build markdown ---
    md = [
        "# PR 16 — B2 Diagnostic Ablation (MIX-15)",
        "",
        "Phase 1 diagnostic for Sprint 4 PR 16 (B2 calibration: MIX-15 cold placement).",
        "Mixes: MIX-15 (B2 target), MIX-01 (Reference control).",
        "Gates (S0): PeakMax Δ < ±1.0°F, PeakGrad Δ < ±2.0°F, FieldRMS < 2.0°F, CenterRMS < 1.0°F, CornerRMS < 3.0°F.",
        "",
        "**Ablation question (binary):** does MIX-15 PeakMax Δ vary across the sweep,",
        "or is it invariant the way B1 (MIX-04, MIX-08) was?",
        "",
        "## Phase 1 Verdict",
        "",
        f"| | F_vert sweep | R_form sweep |",
        f"| --- | --- | --- |",
        f"| PeakMax Δ at default values | {fmt_delta(dflt_fv['peak_max_delta'] if dflt_fv else None)}°F | {fmt_delta(dflt_rf['peak_max_delta'] if dflt_rf else None)}°F |",
        f"| PeakMax Δ range (max−min) | {fmt_rms(fv_pm_range)}°F | {fmt_rms(rf_pm_range)}°F |",
        f"| PeakGrad Δ range | {fmt_rms(fv_pg_range)}°F | {fmt_rms(rf_pg_range)}°F |",
        f"| FieldRMS range | {fmt_rms(fv_fr_range)}°F | {fmt_rms(rf_fr_range)}°F |",
        f"| CornerRMS range | {fmt_rms(fv_cr_range)}°F | {fmt_rms(rf_cr_range)}°F |",
        "",
        f"**MIX-01 holds 5/5 at defaults:** {'YES' if mix01_gate else 'NO — gate fail'}",
        "",
        f"**Authority thresholds:** ≥{VARIATION_AUTHORITY:.1f}°F PeakMax variation → boundary-physics-authoritative; "
        f"<{VARIATION_MARGINAL:.1f}°F → hydration-routed (B1 precedent); "
        f"{VARIATION_MARGINAL:.1f}–{VARIATION_AUTHORITY:.1f}°F → marginal.",
        "",
        f"**F_vert verdict:** {verdict_str(fv_pm_range)}",
        f"**R_form verdict:** {verdict_str(rf_pm_range)}",
        f"**Phase 1 verdict (combined):** **{authority_verdict}**",
        "",
    ]

    if authority_verdict == "BOUNDARY-PHYSICS-AUTHORITATIVE":
        md += [
            "**Proposed next step:** Phase 2 — H6 IC-propagation test (§7.6.4).",
            "",
        ]
    elif authority_verdict == "HYDRATION-ROUTED":
        md += [
            "**Proposed next step:** Decision I — diagnostic-only close.",
            "B2 PeakMax is hydration-driven, not boundary-physics. Route alongside B1 to Sprint 5.",
            "Update §7.6.4 to reflect PR 16 routes B2 PeakMax to Sprint 5.",
            "",
        ]
    else:
        md += [
            "**Proposed next step:** WAIT — marginal authority. Report to user for decision before Phase 2.",
            "",
        ]

    # Sweep tables
    md += [
        "## Sweep 1 — F_vert",
        "",
        'Default: `F_VERT_BY_ORIENTATION["unknown"] = 0.15` (Sprint 2 PR 8 calibration).',
        "Override via `construction.vertical_solar_factor` (no source edit required).",
        "",
        build_table(fvert_rows, "F_vert", lambda v: f"{v:.2f}"),
        "",
        "## Sweep 2 — R_form",
        "",
        "Default: `R_FORM_CONTACT_SI = 0.0862` m²·K/W (ADR-04: ACI 347 steel form + wet concrete film).",
        "Override via monkey-patch `thermal_engine_2d.R_FORM_CONTACT_SI` (restored after sweep).",
        "",
        build_table(rform_rows, "R_form (m²·K/W)", lambda v: f"{v:.4f}"),
        "",
    ]

    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote: {OUT_MD}")
    print()

    # --- Stdout Phase 1 summary ---
    print("=" * 70)
    print("PHASE 1 SUMMARY — B2 ABLATION VERDICT")
    print("=" * 70)
    print()

    dflt_pm = dflt_fv["peak_max_delta"] if dflt_fv else None
    print(f"MIX-15 PeakMax Δ at default values:      {fmt_delta(dflt_pm)}°F")
    print(f"MIX-15 PeakMax Δ range — F_vert sweep:   {fmt_rms(fv_pm_range)}°F  "
          f"(min {fmt_delta(fv_pm_min)}, max {fmt_delta(fv_pm_max)})")
    print(f"MIX-15 PeakMax Δ range — R_form sweep:   {fmt_rms(rf_pm_range)}°F  "
          f"(min {fmt_delta(rf_pm_min)}, max {fmt_delta(rf_pm_max)})")
    print()
    print(f"MIX-15 CornerRMS range — F_vert sweep:   {fmt_rms(fv_cr_range)}°F")
    print(f"MIX-15 CornerRMS range — R_form sweep:   {fmt_rms(rf_cr_range)}°F")
    print(f"MIX-15 PeakGrad Δ range — F_vert sweep:  {fmt_rms(fv_pg_range)}°F")
    print(f"MIX-15 PeakGrad Δ range — R_form sweep:  {fmt_rms(rf_pg_range)}°F")
    print(f"MIX-15 FieldRMS range — F_vert sweep:    {fmt_rms(fv_fr_range)}°F")
    print(f"MIX-15 FieldRMS range — R_form sweep:    {fmt_rms(rf_fr_range)}°F")
    print()
    print(f"MIX-01 5/5 at default values:            {'YES — gate passed' if mix01_gate else 'NO — gate FAIL, stop and report'}")
    print()
    print(f"PHASE 1 VERDICT: {authority_verdict}")
    print()

    if authority_verdict == "BOUNDARY-PHYSICS-AUTHORITATIVE":
        print("Proposed Phase 2 path: H6 IC-propagation test first (cheaper).")
        print("  Load MIX-15, override placement_temp_F=60, run engine.")
        print("  Compare synthetic output to MIX-01 CW reference.")
        print("  If H6 confirmed → Decision E candidate.")
        print("  If H6 falsified → H4 grep audit.")
    elif authority_verdict == "HYDRATION-ROUTED":
        print("Proposed path: Decision I — diagnostic-only close.")
        print("  B2 PeakMax is hydration-driven (same finding as B1).")
        print("  Route to Sprint 5 alongside B1.")
        print("  Update §7.6.4 in docs/coding_passdown_v4.md with this finding.")
    else:
        print("Proposed path: WAIT — marginal authority.")
        print("  Variation is in the ambiguous range. User decision required.")
        print("  Options: (a) treat as authoritative and proceed to Phase 2,")
        print("           (b) treat as hydration-routed and close Decision I,")
        print("           (c) run further diagnostics before deciding.")

    print()
    print("WAIT — paste this report, await user confirmation before Phase 2.")


if __name__ == "__main__":
    main()
