"""Sprint 4 PR 15 — B1 diagnostic sweeps: F_vert and R_form.

Phase 1 diagnostic — throwaway script, not committed.
Only the output (pr15_b1_sweeps.md) is promoted to validation/diagnostics/
on Decision D (diagnostic-only close).

Mixes: MIX-04 + MIX-08 (B1 targets) + MIX-01 (Reference control).

Sweep 1 — F_vert: overrides construction.vertical_solar_factor, no source edit.
Sweep 2 — R_form: monkey-patches thermal_engine_2d.R_FORM_CONTACT_SI, restored after.

Usage (from repo root):
    python diagnostics/sprint4/pr15_b1_sweeps.py
"""
import dataclasses
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CW_EXPORTS = REPO / "validation" / "cw_exports"
OUT_MD = Path(__file__).parent / "pr15_b1_sweeps.md"

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

MIXES = ["MIX-01", "MIX-04", "MIX-08"]

F_VERT_SWEEP = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50]
R_FORM_SWEEP = [0.040, 0.060, 0.0862, 0.100, 0.120, 0.150, 0.200]

# Physical plausibility bounds for B.i vs B.ii decision (per prompt)
R_FORM_PLAUSIBLE_LO = 0.05
R_FORM_PLAUSIBLE_HI = 0.12


def rms_arr(a, b):
    diff = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean(diff ** 2)))


def compute_metrics(scn, grid, result):
    """Compute the 5 gate metrics, matching compare_to_cw.py exactly."""
    jslice, islice = grid.concrete_slice()
    T_conc_C = result.T_field_C[:, jslice, islice]
    T_conc_F = T_conc_C * 9.0 / 5.0 + 32.0
    t_hrs_eng = result.t_s / 3600.0

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
    print("=== Sweep 1: F_vert (9 values × 3 mixes) ===")
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
    print("=== Sweep 2: R_form (7 values × 3 mixes) ===")
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

def find_winners(rows):
    """Values where MIX-01 stays 5/5 AND MIX-04 AND MIX-08 both reach 5/5."""
    by_val = defaultdict(dict)
    for (mix_id, param_val, m) in rows:
        by_val[param_val][mix_id] = m
    winners = []
    for val in sorted(by_val):
        mix_map = by_val[val]
        if all(mix_map.get(mid, {}).get("s0_pass", 0) == 5 for mid in MIXES):
            winners.append(val)
    return winners


def b1_optimum(rows):
    """Param value that minimizes mean CornerRMS across MIX-04 + MIX-08."""
    by_val = defaultdict(list)
    for (mix_id, param_val, m) in rows:
        if mix_id in ("MIX-04", "MIX-08") and m["corner_rms"] is not None:
            by_val[param_val].append(m["corner_rms"])
    if not by_val:
        return None
    return min(by_val, key=lambda v: sum(by_val[v]) / len(by_val[v]))


def metric_at(rows, mix_id, val):
    for (mid, v, m) in rows:
        if mid == mix_id and abs(v - val) < 1e-9:
            return m
    return None


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

    # --- Analysis ---
    fvert_winners = find_winners(fvert_rows)
    rform_winners = find_winners(rform_rows)
    fvert_opt     = b1_optimum(fvert_rows)
    rform_opt     = b1_optimum(rform_rows)

    def corner_at(rows, mix_id, val):
        m = metric_at(rows, mix_id, val)
        return m["corner_rms"] if m else None

    # --- Build markdown ---
    md = [
        "# PR 15 — B1 Diagnostic Sweeps",
        "",
        "Phase 1 diagnostic for Sprint 4 PR 15 (B1 calibration: MIX-04, MIX-08).",
        "Mixes: MIX-04 + MIX-08 (B1 targets), MIX-01 (Reference control).",
        "Gates (S0): PeakMax Δ < ±1.0°F, PeakGrad Δ < ±2.0°F, FieldRMS < 2.0°F, CenterRMS < 1.0°F, CornerRMS < 3.0°F.",
        "",
        "## Winner Summary",
        "",
    ]

    # Q1
    md.append("**Q1 — Any F_vert value drives both B1 mixes to S0 PASS without breaking MIX-01?**")
    if fvert_winners:
        md.append(f"YES — F_vert = {[f'{v:.2f}' for v in fvert_winners]} (proposed Decision A).")
    else:
        f04 = corner_at(fvert_rows, "MIX-04", fvert_opt)
        f08 = corner_at(fvert_rows, "MIX-08", fvert_opt)
        opt_str = f"{fvert_opt:.2f}" if fvert_opt is not None else "N/A"
        md.append(
            f"NO — No single F_vert value achieves 5/5 on all three mixes. "
            f"B1 optimum (min mean CornerRMS): F_vert = {opt_str} → "
            f"MIX-04 CornerRMS = {fmt_rms(f04)}°F, MIX-08 CornerRMS = {fmt_rms(f08)}°F."
        )
    md.append("")

    # Q2
    md.append("**Q2 — Any R_form value drives both B1 mixes to S0 PASS without breaking MIX-01?**")
    if rform_winners:
        plaus = all(R_FORM_PLAUSIBLE_LO <= v <= R_FORM_PLAUSIBLE_HI for v in rform_winners)
        path = "Decision B.i (within physical plausibility)" if plaus else "Decision B.ii (outside ACI 347 plausibility)"
        md.append(f"YES — R_form = {[f'{v:.4f}' for v in rform_winners]} (proposed {path}).")
    else:
        r04 = corner_at(rform_rows, "MIX-04", rform_opt)
        r08 = corner_at(rform_rows, "MIX-08", rform_opt)
        opt_str = f"{rform_opt:.4f}" if rform_opt is not None else "N/A"
        md.append(
            f"NO — No single R_form value achieves 5/5 on all three mixes. "
            f"B1 optimum (min mean CornerRMS): R_form = {opt_str} → "
            f"MIX-04 CornerRMS = {fmt_rms(r04)}°F, MIX-08 CornerRMS = {fmt_rms(r08)}°F."
        )
    md.append("")

    # Q3
    md.append("**Q3 — If neither sweep alone fixes both B1 mixes: do their optima agree (direction + magnitude)?**")
    if fvert_winners or rform_winners:
        md.append("At least one sweep found a clean winner — Q3 moot.")
    elif fvert_opt is not None and rform_opt is not None:
        f_dir = "higher" if fvert_opt > 0.15   else ("lower" if fvert_opt < 0.15   else "same as default")
        r_dir = "higher" if rform_opt > 0.0862 else ("lower" if rform_opt < 0.0862 else "same as default")
        agree = (f_dir == r_dir)
        md.append(
            f"F_vert optimum = {fvert_opt:.2f} ({f_dir} than default 0.15). "
            f"R_form optimum = {rform_opt:.4f} ({r_dir} than default 0.0862). "
            f"Direction {'AGREES' if agree else 'DIFFERS'} — both sweeps push the form-face heat balance "
            f"{'in the same direction' if agree else 'in opposite directions'}."
        )
    else:
        md.append("Could not determine optima.")
    md.append("")

    # Q4
    md.append("**Q4 — Is residual CornerRMS at the optimum consistent with H3 catch-all?**")
    if fvert_winners or rform_winners:
        md.append("A clean winner was found — H3 catch-all not needed for this PR.")
    elif fvert_opt is not None:
        f04 = corner_at(fvert_rows, "MIX-04", fvert_opt)
        f08 = corner_at(fvert_rows, "MIX-08", fvert_opt)
        still_fail = all(v is not None and v >= TOL_CORNER_RMS_F for v in [f04, f08])
        if still_fail:
            md.append(
                f"YES — At F_vert = {fvert_opt:.2f} (B1 optimum), MIX-04 CornerRMS = {fmt_rms(f04)}°F and "
                f"MIX-08 CornerRMS = {fmt_rms(f08)}°F remain above the 3.0°F S0 threshold. "
                f"Neither H1 nor H2 fully close the gap — residual consistent with H3 catch-all "
                f"(some boundary physics term beyond F_vert + R_form). Proposed Decision D."
            )
        else:
            md.append(
                f"PARTIAL — At F_vert = {fvert_opt:.2f}, one or both B1 mixes cross the S0 threshold "
                f"(MIX-04={fmt_rms(f04)}°F, MIX-08={fmt_rms(f08)}°F). Check Q3 joint fix potential (Decision C)."
            )
    md.append("")

    # Sweep tables
    md += [
        "## Sweep 1 — F_vert",
        "",
        'Default: `F_VERT_BY_ORIENTATION["unknown"] = 0.15` (Sprint 2 PR 8 calibration). '
        "Override via `construction.vertical_solar_factor` (no source edit required).",
        "",
        build_table(fvert_rows, "F_vert", lambda v: f"{v:.2f}"),
        "",
        "## Sweep 2 — R_form",
        "",
        "Default: `R_FORM_CONTACT_SI = 0.0862` m²·K/W (ADR-04: ACI 347 steel form + wet concrete film). "
        "Override via monkey-patch `thermal_engine_2d.R_FORM_CONTACT_SI` (restored after sweep).",
        "",
        build_table(rform_rows, "R_form (m²·K/W)", lambda v: f"{v:.4f}"),
        "",
    ]

    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote: {OUT_MD}")
    print()

    # --- Stdout four-question summary ---
    print("=" * 70)
    print("PHASE 1 SUMMARY — §7.6.3 FOUR QUESTIONS")
    print("=" * 70)
    print()
    print("Q1 — Any F_vert drives both B1 mixes to S0 PASS, MIX-01 stays 5/5?")
    if fvert_winners:
        print(f"     YES → F_vert = {[f'{v:.2f}' for v in fvert_winners]}")
    else:
        f04 = corner_at(fvert_rows, "MIX-04", fvert_opt)
        f08 = corner_at(fvert_rows, "MIX-08", fvert_opt)
        print(f"     NO  → best F_vert = {fvert_opt:.2f}  "
              f"MIX-04 CornerRMS={fmt_rms(f04)}°F  MIX-08 CornerRMS={fmt_rms(f08)}°F")
    print()
    print("Q2 — Any R_form drives both B1 mixes to S0 PASS, MIX-01 stays 5/5?")
    if rform_winners:
        plaus = all(R_FORM_PLAUSIBLE_LO <= v <= R_FORM_PLAUSIBLE_HI for v in rform_winners)
        path = "B.i (physical)" if plaus else "B.ii (outside ACI 347)"
        print(f"     YES → R_form = {[f'{v:.4f}' for v in rform_winners]}  → proposed {path}")
    else:
        r04 = corner_at(rform_rows, "MIX-04", rform_opt)
        r08 = corner_at(rform_rows, "MIX-08", rform_opt)
        print(f"     NO  → best R_form = {rform_opt:.4f}  "
              f"MIX-04 CornerRMS={fmt_rms(r04)}°F  MIX-08 CornerRMS={fmt_rms(r08)}°F")
    print()

    if not fvert_winners and not rform_winners and fvert_opt is not None and rform_opt is not None:
        f_dir = "higher" if fvert_opt > 0.15   else "lower"
        r_dir = "higher" if rform_opt > 0.0862 else "lower"
        agree = f_dir == r_dir
        print("Q3 — Do optima agree in direction?")
        print(f"     F_vert opt {fvert_opt:.2f} → {f_dir} than default (0.15)")
        print(f"     R_form opt {rform_opt:.4f} → {r_dir} than default (0.0862)")
        print(f"     {'AGREE' if agree else 'DIFFER'}")
        print()
        f04 = corner_at(fvert_rows, "MIX-04", fvert_opt)
        f08 = corner_at(fvert_rows, "MIX-08", fvert_opt)
        print("Q4 — Residual CornerRMS at F_vert optimum vs H3 catch-all?")
        print(f"     MIX-04 CornerRMS = {fmt_rms(f04)}°F  (S0 tol = {TOL_CORNER_RMS_F:.1f}°F)")
        print(f"     MIX-08 CornerRMS = {fmt_rms(f08)}°F")
        still_fail = all(v is not None and v >= TOL_CORNER_RMS_F for v in [f04, f08])
        if still_fail:
            print("     Both remain above S0 threshold → CONSISTENT WITH H3 catch-all")
        else:
            print("     At least one is below S0 threshold → check joint fix (Decision C)")
    print()

    # Decision path proposal
    print("PROPOSED DECISION PATH:")
    if fvert_winners:
        print(f"  A — clean H1 winner (F_vert = {fvert_winners[0]:.2f})")
        print(f"  Action: recalibrate F_VERT_BY_ORIENTATION[\"unknown\"] to {fvert_winners[0]:.2f} in thermal_engine_2d.py")
        print(f"          update ADR-05 in docs/coding_passdown_v4.md")
    elif rform_winners:
        plaus = all(R_FORM_PLAUSIBLE_LO <= v <= R_FORM_PLAUSIBLE_HI for v in rform_winners)
        if plaus:
            print(f"  B.i — clean H2 winner within physical plausibility (R_form = {rform_winners[0]:.4f})")
            print(f"  Action: recalibrate R_FORM_CONTACT_SI to {rform_winners[0]:.4f} in thermal_engine_2d.py")
            print(f"          update ADR-04 in docs/coding_passdown_v4.md")
        else:
            print(f"  B.ii — H2 winner outside ACI 347 physical plausibility (R_form = {rform_winners[0]:.4f})")
            print(f"  Action: diagnostic-only close; route to Sprint 6 (R5 form-material parameterization)")
    else:
        print(f"  D — no single-knob fix; diagnostic-only close")
        print(f"  Action: promote pr15_b1_sweeps.md to validation/diagnostics/")
        print(f"          update §7.6.3 in docs/coding_passdown_v4.md with finding")
        print(f"          route B1 fix to Sprint 5 (dual-peak hydration scope)")

    print()
    print("WAIT — review sweep tables, confirm decision path, then reply with go/no-go on Phase 2.")


if __name__ == "__main__":
    main()
