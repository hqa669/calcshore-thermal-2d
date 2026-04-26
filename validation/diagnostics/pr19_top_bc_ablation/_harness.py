"""
Shared helpers for PR 19 top-surface BC ablation.

Three knobs ablated: h_conv (monkey-patch, lines 496/524/1354/1355),
solar_absorptivity_top (getattr seam, line 1358), emissivity_top
(getattr seam, line 1359).

Ablation-only diagnostic: no engine source modifications.
"""
import contextlib
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _REPO)

import thermal_engine_2d as _te
from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d

# Evaluation sets (per mix02_recon.md disposition C)
CLUSTER = ("MIX-01", "MIX-11", "MIX-12")
CONTROL = ("MIX-03",)
EVAL_MIXES = CLUSTER + CONTROL

SCENARIO_ROOT = os.path.join(_REPO, "validation", "cw_exports")
WIND_M_S = 10.5  # cw_ave_max_wind_m_s — bit-identical across all Reference mixes

# S0 gate tolerances (mirrors compare_to_cw.py)
TOL_PEAK_MAX_F   = 1.0
TOL_PEAK_GRAD_F  = 2.0
TOL_FIELD_RMS_F  = 2.0
TOL_CENTER_RMS_F = 1.0
TOL_CORNER_RMS_F = 3.0
T_START_RMS_HR   = 48.0

AUTHORITY_THRESHOLD_F = 0.24  # ≥0.24°F cluster-mean CenterRMS variation = authoritative


@contextlib.contextmanager
def patched_h_conv(*, scale: float = 1.0, form: str = "A", calibrated_h: float | None = None):
    """Replace thermal_engine_2d.h_forced_convection for one engine run.

    ABLATION TOOL ONLY — not a production pattern. Replaces the module
    attribute so all 4 call sites are affected (lines 496, 524, 1354, 1355),
    including side form-face combined-h (line 524). Side-h authority is
    dominated by R_FORM_CONTACT_SI in series; cluster CenterRMS sensitivity
    is overwhelmingly top-surface-driven.

    Forms:
      A: ACI 305 linear-in-wind × scale  (same shape, different magnitude)
      B: constant = calibrated_h          (no wind coupling; requires calibrated_h)
      C: Duffie-Beckman  h(v) = 5.7 + 3.8·v^0.78  (not scaled by `scale`)
    """
    original = _te.h_forced_convection
    if form == "A":
        _fn = lambda w: original(w) * scale          # noqa: E731
    elif form == "B":
        if calibrated_h is None:
            raise ValueError("patched_h_conv form='B' requires calibrated_h")
        _fn = lambda w: float(calibrated_h)          # noqa: E731
    elif form == "C":
        _fn = lambda w: 5.7 + 3.8 * (w ** 0.78)     # noqa: E731
    else:
        raise ValueError(f"Unknown h_conv form {form!r}")
    _te.h_forced_convection = _fn
    try:
        yield
    finally:
        _te.h_forced_convection = original


def load_scenarios():
    """Load all 4 evaluation scenarios. Returns dict mix_id -> CWScenario."""
    scenarios = {}
    for mix in EVAL_MIXES:
        mix_dir = os.path.join(SCENARIO_ROOT, mix)
        scn = load_cw_scenario(
            os.path.join(mix_dir, "input.dat"),
            os.path.join(mix_dir, "weather.dat"),
            os.path.join(mix_dir, "output.txt"),
        )
        scenarios[mix] = scn
    return scenarios


def run_engine(scn, construction_override=None):
    """Run the 2-D hydration solver. Returns (result, grid)."""
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


def compute_metrics(scn, grid, result):
    """Compute S0 gate metrics and PR 18 D1–D4 diagnostics in one pass.

    Returns a flat dict. Key fields:
      center_rms    — S0 gate / authority target (matches compare_to_cw.py
                      linear-interp formula; same metric as compare_to_cw center_rms)
      d4_mid_rms_F  — D4 mid-depth RMS (searchsorted, matches PR 18 d4 formula)
      d4_top_rms_F  — D4 top-surface RMS (top-concentrated signal indicator)
      d4_spread_F   — D4 max-min depth RMS spread
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

    val = scn.cw_validation
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

    # D1 — constant offset
    diff_gate  = eng_center_at_cw[rms_mask] - cw_center_F[rms_mask]
    d1_mean_dt = float(np.mean(diff_gate))

    # D2 — phase lag (cross-correlation of detrended gate signals)
    eng_g = eng_center_at_cw[rms_mask]
    cw_g  = cw_center_F[rms_mask]
    n     = len(eng_g)
    corr  = np.correlate(eng_g - eng_g.mean(), cw_g - cw_g.mean(), mode="full")
    lags  = np.arange(-(n - 1), n)
    peak_lag   = int(lags[int(np.argmax(corr))])
    stride_hr  = float(np.median(np.diff(val.time_hrs)))
    d2_lag_hr  = peak_lag * stride_hr

    # D3 — amplitude ratio (last full diurnal cycle)
    eng_dm = eng_center_at_cw[diurnal_mask]
    cw_dm  = cw_center_F[diurnal_mask]
    eng_sw = float(eng_dm.max() - eng_dm.min()) if len(eng_dm) > 0 else float("nan")
    cw_sw  = float(cw_dm.max() - cw_dm.min())  if len(cw_dm)  > 0 else float("nan")
    d3_ratio = eng_sw / cw_sw if cw_sw > 0.0 else float("nan")

    # D4 + field_rms — single searchsorted pass over CW time points
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


def analyze_sweep_verdict(rows, *, default_val):
    """Authority + generalization analysis for one knob's magnitude sweep.

    Args:
        rows: [(mix_id, param_val, metrics_dict), ...]
        default_val: the default parameter value (scale=1.0, α=0.65, ε=0.88)

    Returns dict with keys:
        has_authority, variation_F, optimum, cluster_mean_at_optimum,
        mix03_baseline_rms, mix03_optimum_rms, mix03_delta,
        mix03_generalizes, mix03_hurts_substantially
    """
    from collections import defaultdict
    by_val = defaultdict(dict)
    for (mix_id, pval, m) in rows:
        by_val[pval][mix_id] = m

    cluster_mean = {}
    for pval in sorted(by_val):
        crms = [by_val[pval][mid]["center_rms"]
                for mid in CLUSTER
                if mid in by_val[pval]
                and by_val[pval][mid].get("center_rms") is not None]
        if crms:
            cluster_mean[pval] = sum(crms) / len(crms)

    if not cluster_mean:
        return {
            "has_authority": False, "variation_F": 0.0, "optimum": None,
            "cluster_mean_at_optimum": None,
            "mix03_baseline_rms": None, "mix03_optimum_rms": None,
            "mix03_delta": None, "mix03_generalizes": None,
            "mix03_hurts_substantially": None,
        }

    rms_vals      = list(cluster_mean.values())
    variation     = max(rms_vals) - min(rms_vals)
    has_authority = variation >= AUTHORITY_THRESHOLD_F
    optimum       = min(cluster_mean, key=cluster_mean.get)

    # Baseline: param value closest to the published default
    baseline_key = min(by_val.keys(), key=lambda v: abs(v - default_val))

    baseline_mix03 = by_val.get(baseline_key, {}).get("MIX-03", {}).get("center_rms")
    optimum_mix03  = by_val.get(optimum,       {}).get("MIX-03", {}).get("center_rms")

    if baseline_mix03 is not None and optimum_mix03 is not None:
        mix03_delta              = optimum_mix03 - baseline_mix03
        mix03_generalizes        = mix03_delta <= 0.10
        mix03_hurts_substantially = mix03_delta > 0.20
    else:
        mix03_delta = mix03_generalizes = mix03_hurts_substantially = None

    return {
        "has_authority":              has_authority,
        "variation_F":                variation,
        "optimum":                    optimum,
        "cluster_mean_at_optimum":    cluster_mean.get(optimum),
        "mix03_baseline_rms":         baseline_mix03,
        "mix03_optimum_rms":          optimum_mix03,
        "mix03_delta":                mix03_delta,
        "mix03_generalizes":          mix03_generalizes,
        "mix03_hurts_substantially":  mix03_hurts_substantially,
    }


def format_sweep_table(rows, *, param_label, param_fmt):
    """Markdown sweep table matching PR 15's build_table + PR 18 D1–D4 columns.

    rows: [(mix_id, param_val, metrics_dict), ...]
    """
    header = (
        f"| Mix | {param_label} "
        "| D1 ΔT̄ (°F) | D2 lag (hr) | D3 ratio "
        "| CenterRMS (°F) | D4 TopRMS (°F) | D4 spread (°F) | S0 |"
    )
    sep = "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | :---: |"
    lines = [header, sep]

    def _f(v, fmt):
        return fmt % v if (v is not None and np.isfinite(v)) else "N/A"

    for (mix_id, pval, m) in rows:
        s0_cell = f"**{m['s0_pass']}/5**" if m["s0_pass"] == 5 else f"{m['s0_pass']}/5"
        lines.append(
            f"| {mix_id} | {param_fmt(pval)} "
            f"| {_f(m.get('d1_mean_dt_F'), '%+.3f')} "
            f"| {_f(m.get('d2_lag_hr'),    '%+.2f')} "
            f"| {_f(m.get('d3_ratio'),     '%.4f')} "
            f"| {_f(m.get('center_rms'),   '%.4f')} "
            f"| {_f(m.get('d4_top_rms_F'), '%.4f')} "
            f"| {_f(m.get('d4_spread_F'),  '%.4f')} "
            f"| {s0_cell} |"
        )
    return "\n".join(lines)
