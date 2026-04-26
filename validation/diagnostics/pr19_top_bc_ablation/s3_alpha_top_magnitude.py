"""
PR 19 / Knob 2 — solar absorptivity (α_top) magnitude sweep.

Overrides construction.solar_absorptivity_top via dataclasses.replace()
(engine getattr seam, thermal_engine_2d.py line 1358).
Default: SOLAR_ABSORPTIVITY_DEFAULT = 0.65 (ASHRAE light-gray concrete).

6 values × 4 mixes = 24 engine runs.
α_top is expected to produce one-sided amplitude error (daytime peak shifts,
nighttime trough less affected) given D3's symmetric over-oscillation.
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

ALPHA_SWEEP = [0.30, 0.45, 0.55, 0.65, 0.75, 0.85]
DEFAULT_ALPHA = 0.65

_OUT = os.path.join(_HERE, "s3_alpha_top_magnitude.md")


def main():
    print("=== s3 — α_top magnitude sweep ===")
    print(f"  Default α_top = {DEFAULT_ALPHA}  (ASHRAE light-gray concrete)")
    print(f"  Sweep: {ALPHA_SWEEP}")

    print("\nLoading scenarios ...", flush=True)
    scenarios = load_scenarios()
    print(f"  Loaded: {list(scenarios)}")

    rows = []
    stop_flag = None

    for alpha in ALPHA_SWEEP:
        print(f"\n  α_top={alpha:.2f}")
        for mix_id in EVAL_MIXES:
            t0 = time.perf_counter()
            scn = scenarios[mix_id]
            construction = dataclasses.replace(scn.construction,
                                               solar_absorptivity_top=alpha)
            result, grid = run_engine(scn, construction_override=construction)
            wall = time.perf_counter() - t0
            if not np.all(np.isfinite(result.T_field_C)):
                print(f"    [{mix_id}] NaN — STOP")
                sys.exit(1)
            m = compute_metrics(scn, grid, result)
            rows.append((mix_id, alpha, m))
            print(f"    [{mix_id}]  CenterRMS={m['center_rms']:.4f}°F"
                  f"  D3={m['d3_ratio']:.4f}  S0={m['s0_pass']}/5  ({wall:.1f}s)")
            if mix_id in CLUSTER and m["s0_pass"] < 5:
                stop_flag = (mix_id, alpha)

        if stop_flag is not None:
            mix_id, a = stop_flag
            print(f"\n  STOP: {mix_id} S0 failure at α_top={a:.2f}")
            print("  Record and do not commit. Investigate before continuing.")
            break

    verdict = analyze_sweep_verdict(rows, default_val=DEFAULT_ALPHA)
    _write_md(rows, verdict, stop_flag=stop_flag)
    _print_verdict(verdict)
    return rows, verdict, stop_flag


def _write_md(rows, verdict, *, stop_flag):
    lines = [
        "# s3 — α_top (Solar Absorptivity) Magnitude Sweep",
        "",
        "**Knob**: `construction.solar_absorptivity_top` override via `dataclasses.replace()`.",
        "Engine reads this via `getattr(construction, 'solar_absorptivity_top',",
        "SOLAR_ABSORPTIVITY_DEFAULT)` at `thermal_engine_2d.py:1358`.",
        "",
        f"Default: α_top = {DEFAULT_ALPHA}  (ASHRAE light-gray concrete/steel).",
        "Expected signature: one-sided amplitude error (daytime peak shifts,",
        "nighttime trough less affected). D3's symmetric over-oscillation makes",
        "α_top a lower-probability dominant knob than h_conv.",
        "",
    ]
    if stop_flag:
        mix_id, a = stop_flag
        lines += [
            f"⚠️ **STOP CONDITION**: {mix_id} S0 failure at α_top={a:.2f}.",
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
        param_label="α_top",
        param_fmt=lambda a: f"{a:.2f}",
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
            f"- Cluster optimum: α_top = **{verdict['optimum']:.2f}**"
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
    print("\n--- s3 verdict ---")
    print(f"  α_top authority: {'YES' if verdict['has_authority'] else 'NO'}"
          f"  (variation={verdict['variation_F']:.4f}°F)")
    if verdict["optimum"] is not None:
        print(f"  Cluster optimum: α_top={verdict['optimum']:.2f}"
              f"  cluster-mean CenterRMS={verdict['cluster_mean_at_optimum']:.4f}°F")
    if verdict.get("mix03_delta") is not None:
        sign = "+" if verdict["mix03_delta"] >= 0 else ""
        print(f"  MIX-03: Δ={sign}{verdict['mix03_delta']:.4f}°F"
              f"  generalizes={'YES' if verdict['mix03_generalizes'] else 'NO'}")


if __name__ == "__main__":
    main()
