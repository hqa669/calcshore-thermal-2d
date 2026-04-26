"""
PR 19 / Knob 1b — h_conv functional-form check.

Conditional on s1's magnitude sweep finding a clean optimum.
At the calibrated magnitude (cluster optimum from s1), compares:
  Form A: ACI 305 linear-in-wind × scale  (already in s1; reference row)
  Form B: constant = calibrated_h         (wind coupling eliminated)
  Form C: Duffie-Beckman h(v)=5.7+3.8·v^0.78  (different functional shape)

2 new forms × 4 mixes = 8 runs.
If s1 finds no clean optimum, this module writes a skip document instead.

Usage:
    python s2_h_conv_form.py  [--calibrated-h VALUE] [--skip [REASON]]
    Or call main(calibrated_h=...) / main(skip=True, skip_reason=...) from driver.
"""
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from _harness import (
    EVAL_MIXES, CLUSTER, WIND_M_S,
    patched_h_conv, load_scenarios, run_engine, compute_metrics,
    format_sweep_table,
)

DEFAULT_H = 5.6 + 3.5 * (0.4 * WIND_M_S)  # 20.3 W/(m²·K) at v=10.5 m/s
DB_H_AT_WIND = 5.7 + 3.8 * (WIND_M_S ** 0.78)  # Duffie-Beckman at 10.5 m/s

_OUT = os.path.join(_HERE, "s2_h_conv_form.md")


def main(*, calibrated_h=None, skip=False, skip_reason=""):
    """Run the functional-form check or write a skip document.

    Args:
        calibrated_h: h_eff (W/m²K) at the cluster optimum from s1 (required if not skip).
        skip: if True, write a skip document and return (None, True).
        skip_reason: explanation for the skip document.

    Returns:
        (rows_or_None, skipped_flag)
    """
    if skip:
        _write_skip_md(skip_reason)
        return None, True

    if calibrated_h is None:
        raise ValueError("s2.main() requires calibrated_h when not skip=True")

    print("=== s2 — h_conv functional-form check ===")
    print(f"  Calibrated h = {calibrated_h:.2f} W/(m²·K)")
    print(f"  Form A (reference, from s1): {calibrated_h:.2f} W/(m²·K)  [ACI 305 × scale]")
    print(f"  Form B (constant):           {calibrated_h:.2f} W/(m²·K)  [no wind coupling]")
    print(f"  Form C (Duffie-Beckman):     {DB_H_AT_WIND:.2f} W/(m²·K)  [5.7+3.8·v^0.78 at v={WIND_M_S}]")

    print("\nLoading scenarios ...", flush=True)
    scenarios = load_scenarios()

    rows = []
    # Form A reference: run at scale=calibrated_h/DEFAULT_H so the patched form
    # evaluates to calibrated_h at the library wind. This gives the Form A baseline
    # at the calibrated magnitude, directly comparable to s1's optimum row.
    form_a_scale = calibrated_h / DEFAULT_H

    for form, label in [("A", f"A-linear×{form_a_scale:.2f}"), ("B", "B-constant"), ("C", "C-DB")]:
        print(f"\n  Form {form} ({label})")
        for mix_id in EVAL_MIXES:
            t0 = time.perf_counter()
            scn = scenarios[mix_id]
            patch_kwargs = dict(form=form)
            if form == "A":
                patch_kwargs["scale"] = form_a_scale
            elif form == "B":
                patch_kwargs["calibrated_h"] = calibrated_h
            with patched_h_conv(**patch_kwargs):
                result, grid = run_engine(scn)
            wall = time.perf_counter() - t0
            if not np.all(np.isfinite(result.T_field_C)):
                print(f"    [{mix_id}] NaN — STOP")
                sys.exit(1)
            m = compute_metrics(scn, grid, result)
            rows.append((mix_id, form, m))
            print(f"    [{mix_id}]  CenterRMS={m['center_rms']:.4f}°F"
                  f"  D3={m['d3_ratio']:.4f}  S0={m['s0_pass']}/5  ({wall:.1f}s)")

    _write_md(rows, calibrated_h=calibrated_h)
    return rows, False


def _write_md(rows, *, calibrated_h):
    form_labels = {
        "A": f"A — ACI 305 linear×scale (h={calibrated_h:.2f} W/m²K at v={WIND_M_S})",
        "B": f"B — constant={calibrated_h:.2f} W/m²K (no wind coupling)",
        "C": f"C — Duffie-Beckman (h≈{DB_H_AT_WIND:.1f} W/m²K at v={WIND_M_S})",
    }
    lines = [
        "# s2 — h_conv Functional-Form Check",
        "",
        "Compares three h_conv functional forms at the calibrated magnitude",
        f"(cluster optimum from s1: h = {calibrated_h:.2f} W/(m²·K)).",
        "",
        "Form A is the ACI 305 linear-in-wind form at the calibrated scale.",
        "Form B eliminates wind coupling (constant h). Tests whether the wind",
        "dependence is the issue vs. the magnitude.",
        "Form C is Duffie-Beckman parallel-flow (h = 5.7 + 3.8·v^0.78).",
        f"At v = {WIND_M_S} m/s, Form C evaluates to ≈{DB_H_AT_WIND:.1f} W/(m²·K).",
        "",
        "The library has zero wind-speed variation across mixes (single trace),",
        "so PR 19 cannot calibrate a wind exponent — form substitution tests",
        "whether form choice changes the residual signature at the library wind.",
        "",
        "## Sweep table",
        "",
    ]
    form_rows = [(mix_id, form, m) for (mix_id, form, m) in rows]
    lines.append(format_sweep_table(
        form_rows,
        param_label="form",
        param_fmt=lambda f: form_labels.get(f, f),
    ))
    lines += [
        "",
        "## Form-sensitivity assessment",
        "",
    ]
    # Summarise cluster-mean CenterRMS per form
    from collections import defaultdict
    by_form = defaultdict(list)
    for (mix_id, form, m) in rows:
        if mix_id in CLUSTER and m.get("center_rms") is not None:
            by_form[form].append(m["center_rms"])

    for form in ("A", "B", "C"):
        vals = by_form.get(form, [])
        if vals:
            mean_cr = sum(vals) / len(vals)
            lines.append(f"- Form {form}: cluster-mean CenterRMS = {mean_cr:.4f}°F")

    lines += [
        "",
        "Interpretation: if cluster-mean CenterRMS is similar across forms A/B/C,",
        "the residual is dominated by magnitude, not functional shape. If Form C",
        "diverges substantially, the wind coupling exponent matters at this library wind.",
    ]

    with open(_OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  Written: {_OUT}")


def _write_skip_md(reason):
    lines = [
        "# s2 — h_conv Functional-Form Check: SKIPPED",
        "",
        "This module was not run.",
        "",
        f"**Reason**: {reason}",
        "",
        "Per PR 19 protocol: functional-form check is conditional on s1",
        "identifying a clean magnitude optimum. Without a calibrated magnitude,",
        "the form-substitution question has no anchor point.",
        "",
        "PR 20 routing is driven by s1/s3/s4 authority findings.",
    ]
    with open(_OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  Written (skip doc): {_OUT}")


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--calibrated-h", type=float, default=None)
    p.add_argument("--skip", action="store_true")
    p.add_argument("--skip-reason", default="no magnitude optimum identified from s1")
    args = p.parse_args()
    main(calibrated_h=args.calibrated_h, skip=args.skip, skip_reason=args.skip_reason)
