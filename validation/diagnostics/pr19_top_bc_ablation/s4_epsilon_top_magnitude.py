"""
PR 19 / Knob 3 — top-surface emissivity (ε_top) magnitude sweep.

Overrides construction.emissivity_top via dataclasses.replace()
(engine getattr seam, thermal_engine_2d.py line 1359).
Default: EMISSIVITY_DEFAULT = 0.88 (blanket outer surface).

4 values × 4 mixes = 16 engine runs.
ε_top affects LW radiation balance (primarily nighttime cooling). Like α_top,
this produces a one-sided amplitude correction. ASHRAE values for blanketed
concrete cluster tightly around 0.88, so sweep range is narrow.
"""
import dataclasses
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from _harness import (
    EVAL_MIXES, CLUSTER, AUTHORITY_THRESHOLD_F,
    load_scenarios, run_engine, compute_metrics,
    analyze_sweep_verdict, format_sweep_table,
)

EPS_SWEEP = [0.50, 0.70, 0.88, 0.95]
DEFAULT_EPS = 0.88

_OUT = os.path.join(_HERE, "s4_epsilon_top_magnitude.md")


def main():
    print("=== s4 — ε_top magnitude sweep ===")
    print(f"  Default ε_top = {DEFAULT_EPS}  (blanketed concrete outer surface)")
    print(f"  Sweep: {EPS_SWEEP}")

    print("\nLoading scenarios ...", flush=True)
    scenarios = load_scenarios()
    print(f"  Loaded: {list(scenarios)}")

    rows = []
    stop_flag = None

    for eps in EPS_SWEEP:
        print(f"\n  ε_top={eps:.2f}")
        for mix_id in EVAL_MIXES:
            t0 = time.perf_counter()
            scn = scenarios[mix_id]
            construction = dataclasses.replace(scn.construction, emissivity_top=eps)
            result, grid = run_engine(scn, construction_override=construction)
            wall = time.perf_counter() - t0
            if not np.all(np.isfinite(result.T_field_C)):
                print(f"    [{mix_id}] NaN — STOP")
                sys.exit(1)
            m = compute_metrics(scn, grid, result)
            rows.append((mix_id, eps, m))
            print(f"    [{mix_id}]  CenterRMS={m['center_rms']:.4f}°F"
                  f"  D3={m['d3_ratio']:.4f}  S0={m['s0_pass']}/5  ({wall:.1f}s)")
            if mix_id in CLUSTER and m["s0_pass"] < 5:
                stop_flag = (mix_id, eps)

        if stop_flag is not None:
            mix_id, e = stop_flag
            print(f"\n  STOP: {mix_id} S0 failure at ε_top={e:.2f}")
            print("  Record and do not commit. Investigate before continuing.")
            break

    verdict = analyze_sweep_verdict(rows, default_val=DEFAULT_EPS)
    _write_md(rows, verdict, stop_flag=stop_flag)
    _print_verdict(verdict)
    return rows, verdict, stop_flag


def _write_md(rows, verdict, *, stop_flag):
    lines = [
        "# s4 — ε_top (Emissivity) Magnitude Sweep",
        "",
        "**Knob**: `construction.emissivity_top` override via `dataclasses.replace()`.",
        "Engine reads this via `getattr(construction, 'emissivity_top',",
        "EMISSIVITY_DEFAULT)` at `thermal_engine_2d.py:1359`.",
        "",
        f"Default: ε_top = {DEFAULT_EPS}  (ASHRAE blanketed concrete outer surface).",
        "ε_top affects LW radiation balance (nighttime cooling primarily).",
        "Expected signature: one-sided amplitude correction (nighttime trough",
        "shifts, daytime peak less affected). ASHRAE values cluster tightly",
        "near 0.88, so a large sweep range is not physically justified.",
        "",
    ]
    if stop_flag:
        mix_id, e = stop_flag
        lines += [
            f"⚠️ **STOP CONDITION**: {mix_id} S0 failure at ε_top={e:.2f}.",
            "Sweep truncated. Do not commit. Investigate.",
            "",
        ]
    lines += [
        "## Sweep table",
        "",
        "CenterRMS = S0 gate metric, window [48, 168] hr.",
        f"Authority threshold: cluster-mean CenterRMS variation ≥ {AUTHORITY_THRESHOLD_F}°F.",
        "",
    ]
    lines.append(format_sweep_table(
        rows,
        param_label="ε_top",
        param_fmt=lambda e: f"{e:.2f}",
    ))
    lines += [
        "",
        "## Authority assessment",
        "",
        f"- Cluster-mean CenterRMS variation: **{verdict['variation_F']:.4f}°F**",
        f"- Authority (≥{AUTHORITY_THRESHOLD_F}°F): **{'YES' if verdict['has_authority'] else 'NO'}**",
    ]
    if verdict["optimum"] is not None:
        lines.append(
            f"- Cluster optimum: ε_top = **{verdict['optimum']:.2f}**"
            f"  (cluster-mean CenterRMS = {verdict['cluster_mean_at_optimum']:.4f}°F)"
        )
    if verdict.get("mix03_delta") is not None:
        sign = "+" if verdict["mix03_delta"] >= 0 else ""
        lines.append(
            f"- MIX-03 Δ at optimum: {sign}{verdict['mix03_delta']:.4f}°F"
            f"  (generalizes: {'YES' if verdict['mix03_generalizes'] else 'NO'})"
        )

    with open(_OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  Written: {_OUT}")


def _print_verdict(verdict):
    print("\n--- s4 verdict ---")
    print(f"  ε_top authority: {'YES' if verdict['has_authority'] else 'NO'}"
          f"  (variation={verdict['variation_F']:.4f}°F)")
    if verdict["optimum"] is not None:
        print(f"  Cluster optimum: ε_top={verdict['optimum']:.2f}"
              f"  cluster-mean CenterRMS={verdict['cluster_mean_at_optimum']:.4f}°F")


if __name__ == "__main__":
    main()
