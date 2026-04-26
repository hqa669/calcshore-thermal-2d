"""
Shared helpers for PR 20 R_form Reference-set evaluation.

Knob: R_FORM_CONTACT_SI (module-level constant, thermal_engine_2d.py:127).
Resolved at 14 in-engine sites (lines 524, 1913, 1914, 1935, 1944, 1949,
2020, 2192, 2193, 2212, 2221; comment sites at 18, 19, 121, 123, 505, 1759,
1901). Module-attr patch perturbs all sites uniformly — correct semantics:
R_form is a single material property across the entire form face.

Ablation-only diagnostic: no engine source modifications.
"""
import contextlib
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _REPO)

import thermal_engine_2d as _te
from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d

# Full Reference set; MIX-02 is kinetics-anomalous (per mix02_recon.md),
# reported separately in tables but evaluated against the same commit criteria.
CLUSTER   = ("MIX-01", "MIX-11", "MIX-12")
REFERENCE = ("MIX-01", "MIX-02", "MIX-03", "MIX-11", "MIX-12")

SCENARIO_ROOT = os.path.join(_REPO, "validation", "cw_exports")

# S0 gate tolerances (mirrors compare_to_cw.py)
TOL_PEAK_MAX_F   = 1.0
TOL_PEAK_GRAD_F  = 2.0
TOL_FIELD_RMS_F  = 2.0
TOL_CENTER_RMS_F = 1.0
TOL_CORNER_RMS_F = 3.0
T_START_RMS_HR   = 48.0

# Item 1 commit-decision constants
R_FORM_DEFAULT        = 0.0862   # ADR-04 anchor value
R_FORM_CANDIDATE      = 0.060    # PR 15 candidate
CENTER_RMS_DEGRADE_MAX_F = 0.05  # criterion 3: CenterRMS must not worsen by more than this


@contextlib.contextmanager
def patched_r_form(value: float):
    """Set R_FORM_CONTACT_SI for one engine run, then restore.

    Ablation tool only. Patches the module attribute so all 14 in-engine
    sites resolve the new value uniformly.
    """
    original = _te.R_FORM_CONTACT_SI
    _te.R_FORM_CONTACT_SI = value
    try:
        yield
    finally:
        _te.R_FORM_CONTACT_SI = original


def load_scenarios():
    """Load all 5 Reference scenarios. Returns dict mix_id -> CWScenario."""
    scenarios = {}
    for mix in REFERENCE:
        mix_dir = os.path.join(SCENARIO_ROOT, mix)
        scn = load_cw_scenario(
            os.path.join(mix_dir, "input.dat"),
            os.path.join(mix_dir, "weather.dat"),
            os.path.join(mix_dir, "output.txt"),
        )
        scenarios[mix] = scn
    return scenarios


def run_engine(scn):
    """Run the 2-D hydration solver. Returns (result, grid)."""
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


def compute_metrics(scn, grid, result):
    """Compute S0 gate metrics and D1-D4 diagnostics.

    Returns dict: peak_max_delta, peak_grad_delta, field_rms, center_rms,
    corner_rms, s0_pass, d1_mean_dt_F, d2_lag_hr, d3_ratio, d4_mid_rms_F,
    d4_spread_F, d4_top_rms_F.

    Identical to PR 19 harness compute_metrics — same D1-D4 framework.
    """
    jslice, islice = grid.concrete_slice()
    T_conc_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0

    engine_peak_max_F  = float(T_conc_F.max())
    grad_series        = T_conc_F.max(axis=(1, 2)) - T_conc_F.min(axis=(1, 2))
    engine_peak_grad_F = float(grad_series.max())

    iy_mid    = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    iy_top    = grid.iy_concrete_start
    ix_corner = grid.ix_concrete_start

    eng_center_F = result.T_field_C[:, iy_mid, grid.nx - 1] * 9.0 / 5.0 + 32.0
    eng_corner_F = result.T_field_C[:, iy_top, ix_corner]   * 9.0 / 5.0 + 32.0

    val            = scn.cw_validation
    cw_peak_max_F  = float(val.T_max_xs_F.max())
    cw_peak_grad_F = float(val.T_diff_xs_F.max())
    peak_max_delta  = engine_peak_max_F  - cw_peak_max_F
    peak_grad_delta = engine_peak_grad_F - cw_peak_grad_F

    if val.T_field_F is None:
        return {
            "peak_max_delta": peak_max_delta, "peak_grad_delta": peak_grad_delta,
            "field_rms": None, "center_rms": None, "corner_rms": None,
            "s0_pass": (int(abs(peak_max_delta) < TOL_PEAK_MAX_F) +
                        int(abs(peak_grad_delta) < TOL_PEAK_GRAD_F)),
            "d1_mean_dt_F": None, "d2_lag_hr": None,
            "d3_ratio": None, "d4_mid_rms_F": None,
            "d4_spread_F": None, "d4_top_rms_F": None,
        }

    cw_t_s = val.time_hrs * 3600.0
    n_cw_t, n_cw_d, _ = val.T_field_F.shape

    cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]
    cw_corner_F = val.T_field_F[:, 0, -1]

    eng_center_at_cw = np.interp(cw_t_s, result.t_s, eng_center_F)
    eng_corner_at_cw = np.interp(cw_t_s, result.t_s, eng_corner_F)

    rms_mask     = val.time_hrs >= T_START_RMS_HR
    diurnal_mask = (val.time_hrs >= 144.0) & (val.time_hrs <= 168.0)

    def _rms(a, b):
        return float(np.sqrt(np.mean((np.asarray(a, dtype=float) -
                                      np.asarray(b, dtype=float)) ** 2)))

    center_rms_F = _rms(eng_center_at_cw[rms_mask], cw_center_F[rms_mask])
    corner_rms_F = _rms(eng_corner_at_cw[rms_mask], cw_corner_F[rms_mask])

    diff_gate  = eng_center_at_cw[rms_mask] - cw_center_F[rms_mask]
    d1_mean_dt = float(np.mean(diff_gate))

    eng_g = eng_center_at_cw[rms_mask]
    cw_g  = cw_center_F[rms_mask]
    n_pts = len(eng_g)
    corr  = np.correlate(eng_g - eng_g.mean(), cw_g - cw_g.mean(), mode="full")
    lags  = np.arange(-(n_pts - 1), n_pts)
    peak_lag  = int(lags[int(np.argmax(corr))])
    stride_hr = float(np.median(np.diff(val.time_hrs)))
    d2_lag_hr = peak_lag * stride_hr

    eng_dm = eng_center_at_cw[diurnal_mask]
    cw_dm  = cw_center_F[diurnal_mask]
    eng_sw = float(eng_dm.max() - eng_dm.min()) if len(eng_dm) > 0 else float("nan")
    cw_sw  = float(cw_dm.max() - cw_dm.min())  if len(cw_dm)  > 0 else float("nan")
    d3_ratio = eng_sw / cw_sw if cw_sw > 0.0 else float("nan")

    n_conc_y = grid.iy_concrete_end - grid.iy_concrete_start + 1
    n_depths  = min(n_conc_y, n_cw_d)
    sq_flat   = []
    sq_per_d  = [[] for _ in range(n_depths)]

    for ti in range(n_cw_t):
        if val.time_hrs[ti] < T_START_RMS_HR:
            continue
        idx   = min(int(np.searchsorted(result.t_s, cw_t_s[ti])), len(result.t_s) - 1)
        col_C = result.T_field_C[idx, jslice, grid.nx - 1]
        col_F = col_C * 9.0 / 5.0 + 32.0
        n_cmp = min(int(col_F.shape[0]), n_depths)
        for d in range(n_cmp):
            err = float(col_F[d]) - float(val.T_field_F[ti, d, 0])
            sq_flat.append(err * err)
            sq_per_d[d].append(err * err)

    field_rms_F = float(np.sqrt(np.mean(sq_flat))) if sq_flat else float("nan")
    depth_rms   = [float(np.sqrt(np.mean(sq))) if sq else float("nan")
                   for sq in sq_per_d]

    mid_idx    = n_depths // 2
    d4_mid_rms = depth_rms[mid_idx]
    d4_top_rms = depth_rms[0]
    d4_spread  = (float(max(depth_rms) - min(depth_rms))
                  if depth_rms and all(np.isfinite(v) for v in depth_rms)
                  else float("nan"))

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
        "field_rms":       field_rms_F,
        "center_rms":      center_rms_F,
        "corner_rms":      corner_rms_F,
        "s0_pass":         sum(s0.values()),
        "d1_mean_dt_F":    d1_mean_dt,
        "d2_lag_hr":       d2_lag_hr,
        "d3_ratio":        d3_ratio,
        "d4_mid_rms_F":    d4_mid_rms,
        "d4_spread_F":     d4_spread,
        "d4_top_rms_F":    d4_top_rms,
    }


def format_sweep_table(rows):
    """Markdown table for the R_form sweep.

    Columns: Mix | R_form | PeakGrad Δ | CenterRMS | CornerRMS | S0
    rows: [(mix_id, r_form_val, metrics_dict), ...]
    """
    header = ("| Mix | R_form (m²K/W) | PeakGrad Δ (°F) "
              "| CenterRMS (°F) | CornerRMS (°F) | S0 |")
    sep    = "| --- | ---: | ---: | ---: | ---: | :---: |"
    lines  = [header, sep]

    def _f(v, fmt):
        return fmt % v if (v is not None and np.isfinite(v)) else "N/A"

    for (mix_id, rval, m) in rows:
        s0_cell = f"**{m['s0_pass']}/5**" if m["s0_pass"] == 5 else f"{m['s0_pass']}/5"
        lines.append(
            f"| {mix_id} | {rval:.4f} "
            f"| {_f(m.get('peak_grad_delta'), '%+.3f')} "
            f"| {_f(m.get('center_rms'),      '%.4f')} "
            f"| {_f(m.get('corner_rms'),      '%.4f')} "
            f"| {s0_cell} |"
        )
    return "\n".join(lines)


def evaluate_commit_criteria(rows):
    """Check the three Item 1 commit-decision criteria.

    C1: all 5 Reference mixes S0 5/5 at R_form=0.060
    C2: ≥3/5 mixes improve CornerRMS at 0.060 vs 0.0862
    C3: no mix CenterRMS degrades >0.05°F at 0.060 vs 0.0862

    Returns dict with per-criterion pass/fail and supporting evidence.
    """
    by_val = {}
    for (mix_id, rval, m) in rows:
        if rval not in by_val:
            by_val[rval] = {}
        by_val[rval][mix_id] = m

    candidate = by_val.get(R_FORM_CANDIDATE, {})
    default_  = by_val.get(R_FORM_DEFAULT,   {})

    # C1: all 5 Reference mixes S0 5/5 at candidate
    c1_per_mix = {mid: candidate.get(mid, {}).get("s0_pass", -1) == 5
                  for mid in REFERENCE}
    c1_pass    = all(c1_per_mix.values())

    # C2: ≥3/5 mixes improve CornerRMS at candidate vs default
    c2_per_mix = {}
    for mid in REFERENCE:
        cand_crms = candidate.get(mid, {}).get("corner_rms")
        def_crms  = default_.get(mid, {}).get("corner_rms")
        if cand_crms is not None and def_crms is not None:
            c2_per_mix[mid] = cand_crms < def_crms
        else:
            c2_per_mix[mid] = None
    c2_pass = sum(1 for v in c2_per_mix.values() if v is True) >= 3

    # C3: no mix CenterRMS degrades >0.05°F at candidate vs default
    c3_per_mix = {}
    stop_c3    = False
    for mid in REFERENCE:
        cand_cen = candidate.get(mid, {}).get("center_rms")
        def_cen  = default_.get(mid, {}).get("center_rms")
        if cand_cen is not None and def_cen is not None:
            delta = cand_cen - def_cen
            failed = delta > CENTER_RMS_DEGRADE_MAX_F
            c3_per_mix[mid] = {"delta": delta, "fail": failed}
            if failed:
                stop_c3 = True
        else:
            c3_per_mix[mid] = {"delta": None, "fail": None}
    c3_pass = not stop_c3

    return {
        "c1_pass":    c1_pass,
        "c1_per_mix": c1_per_mix,
        "c2_pass":    c2_pass,
        "c2_per_mix": c2_per_mix,
        "c3_pass":    c3_pass,
        "c3_per_mix": c3_per_mix,
        "stop_c3":    stop_c3,
        "commit":     c1_pass and c2_pass and c3_pass,
    }
