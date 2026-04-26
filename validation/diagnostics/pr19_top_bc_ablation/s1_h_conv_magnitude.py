"""
PR 19 / Knob 1 — h_conv magnitude sweep.

Scales h_forced_convection output by factors [0.5, 0.75, 1.0, 1.25, 1.5, 2.0].
At the library wind trace (v=10.5 m/s), default h_eff = 20.3 W/(m²·K).
6 scale values × 4 mixes = 24 engine runs.

Wait-point outputs (two questions per §8.1 two-phase protocol):
  1. Does h_conv have ≥0.24°F authority on cluster-mean CenterRMS?
  2. Does the cluster optimum generalize to MIX-03?
"""
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from _harness import (
    EVAL_MIXES, CLUSTER, CONTROL, WIND_M_S,
    patched_h_conv, load_scenarios, run_engine, compute_metrics,
    analyze_sweep_verdict, format_sweep_table, AUTHORITY_THRESHOLD_F,
)

SCALE_FACTORS = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0]
DEFAULT_H = 5.6 + 3.5 * (0.4 * WIND_M_S)  # 20.3 W/(m²·K) at v=10.5 m/s

_OUT = os.path.join(_HERE, "s1_h_conv_magnitude.md")


def main():
    print("=== s1 — h_conv magnitude sweep ===")
    print(f"  Default h_eff = {DEFAULT_H:.2f} W/(m²·K) at v={WIND_M_S} m/s")
    print(f"  Scales: {SCALE_FACTORS}")

    print("\nLoading scenarios ...", flush=True)
    scenarios = load_scenarios()
    print(f"  Loaded: {list(scenarios)}")

    rows = []
    stop_flag = None  # set to (mix_id, scale) on cluster S0 failure

    for scale in SCALE_FACTORS:
        h_eff = DEFAULT_H * scale
        print(f"\n  scale={scale:.2f}×  h_eff={h_eff:.1f} W/(m²·K)")
        scale_rows = []
        for mix_id in EVAL_MIXES:
            t0 = time.perf_counter()
            scn = scenarios[mix_id]
            with patched_h_conv(scale=scale):
                result, grid = run_engine(scn)
            wall = time.perf_counter() - t0
            if not np.all(np.isfinite(result.T_field_C)):
                print(f"    [{mix_id}] NaN in T_field_C — STOP")
                sys.exit(1)
            m = compute_metrics(scn, grid, result)
            scale_rows.append((mix_id, scale, m))
            cr  = m["center_rms"]
            s0  = m["s0_pass"]
            print(f"    [{mix_id}]  CenterRMS={cr:.4f}°F  D3={m['d3_ratio']:.4f}  S0={s0}/5  ({wall:.1f}s)")
            if mix_id in CLUSTER and s0 < 5:
                stop_flag = (mix_id, scale)

        rows.extend(scale_rows)
        if stop_flag is not None:
            mix_id, sc = stop_flag
            print(f"\n  STOP: {mix_id} S0 failure at scale={sc:.2f}×")
            print("  Record this sweep point; do not commit. Investigate before continuing.")
            break

    # "Too-clean" check (PR 16 H6a lesson)
    too_clean = _check_too_clean(rows)
    if too_clean:
        print("\n  WARNING: result appears 'too clean' (<0.05°F CenterRMS on all 4 mixes).")
        print("  Per PR 16 H6a: run tighter cross-check before committing.")

    verdict = analyze_sweep_verdict(rows, default_val=1.0)
    _write_md(rows, verdict, stop_flag=stop_flag, too_clean=too_clean)
    _print_verdict(verdict)
    return rows, verdict, stop_flag


def _check_too_clean(rows):
    """True if any single scale closes CenterRMS to <0.05°F on all 4 mixes."""
    from collections import defaultdict
    by_val = defaultdict(dict)
    for (mix_id, pval, m) in rows:
        by_val[pval][mix_id] = m
    for pval, mix_map in by_val.items():
        if all(
            mix_map.get(mid, {}).get("center_rms", float("inf")) < 0.05
            for mid in EVAL_MIXES
        ):
            return True
    return False


def _write_md(rows, verdict, *, stop_flag, too_clean):
    lines = [
        "# s1 — h_conv Magnitude Sweep",
        "",
        "**Knob**: `h_forced_convection` output scaled by factor `s` (ablation",
        "monkey-patch — context-manager-reverted after each run, not a production pattern).",
        "",
        "**Side-channel scope**: patch replaces the module attribute and therefore",
        "affects all 4 call sites in thermal_engine_2d (lines 496, 524, 1354, 1355),",
        "including side form-face combined-h (line 524). Side-h is dominated by",
        "R_FORM_CONTACT_SI in series; cluster CenterRMS authority is top-surface-driven.",
        "",
        f"Wind trace: v = {WIND_M_S} m/s (bit-identical across all Reference mixes).",
        f"Default h_eff = {DEFAULT_H:.2f} W/(m²·K)  [ACI 305: `5.6 + 3.5·(0.4·v)`].",
        "",
    ]
    if stop_flag:
        mix_id, sc = stop_flag
        lines += [
            f"⚠️ **STOP CONDITION triggered**: {mix_id} S0 failure at scale={sc:.2f}×.",
            "Sweep truncated at this point per stop-condition protocol.",
            "Do not commit. Investigate before continuing.",
            "",
        ]
    if too_clean:
        lines += [
            "⚠️ **Too-clean warning**: a sweep point produces CenterRMS <0.05°F on",
            "all 4 mixes. Per PR 16 H6a: run tighter cross-check (±10% around optimum)",
            "before committing.",
            "",
        ]
    lines += [
        "## Sweep table",
        "",
        "CenterRMS = S0 gate metric, window [48, 168] hr, centerline mid-depth.",
        "Authority threshold: cluster-mean CenterRMS variation ≥ "
        f"{AUTHORITY_THRESHOLD_F}°F.",
        "",
    ]
    lines.append(format_sweep_table(
        rows,
        param_label="scale (h_eff)",
        param_fmt=lambda s: f"{s:.2f}× ({DEFAULT_H * s:.1f} W/m²K)",
    ))
    lines += [
        "",
        "## Wait-point verdict",
        "",
        f"- Cluster-mean CenterRMS variation across sweep: **{verdict['variation_F']:.4f}°F**",
        f"- Authority threshold ≥{AUTHORITY_THRESHOLD_F}°F: "
        f"**{'YES' if verdict['has_authority'] else 'NO'}**",
    ]
    if verdict["optimum"] is not None:
        lines.append(
            f"- Cluster optimum: scale **{verdict['optimum']:.2f}×**"
            f" → h_eff = {DEFAULT_H * verdict['optimum']:.1f} W/(m²·K)"
            f"  (cluster-mean CenterRMS = {verdict['cluster_mean_at_optimum']:.4f}°F)"
        )
    if verdict["mix03_delta"] is not None:
        sign = "+" if verdict["mix03_delta"] >= 0 else ""
        lines.append(
            f"- MIX-03 CenterRMS: baseline={verdict['mix03_baseline_rms']:.4f}°F"
            f" → optimum={verdict['mix03_optimum_rms']:.4f}°F"
            f"  (Δ={sign}{verdict['mix03_delta']:.4f}°F)"
        )
        lines.append(
            f"- MIX-03 generalizes (Δ ≤ 0.10°F): "
            f"**{'YES' if verdict['mix03_generalizes'] else 'NO'}**"
        )
    if verdict.get("mix03_hurts_substantially"):
        lines += [
            "",
            "⚠️ **OPEN QUESTION**: MIX-03 regresses >0.2°F at the cluster optimum.",
            "Surface to user before writing the routing paragraph (per PR 19 brief).",
        ]
    lines += [""]
    if not verdict["has_authority"]:
        s2_decision = "**SKIP** — no h_conv magnitude authority; form-check moot."
    elif verdict.get("mix03_hurts_substantially"):
        s2_decision = "**SKIP** — open-question MIX-03 regression; resolve before form-check."
    else:
        s2_decision = (
            f"**RUN** at calibrated h = {DEFAULT_H * verdict['optimum']:.1f} W/(m²·K)"
            f"  (scale={verdict['optimum']:.2f}×)"
        )
    lines.append(f"s2 (functional-form check): {s2_decision}")

    with open(_OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n  Written: {_OUT}")


def _print_verdict(verdict):
    print("\n--- Wait-point verdict ---")
    print(f"  Authority: {'YES' if verdict['has_authority'] else 'NO'}"
          f"  (variation={verdict['variation_F']:.4f}°F, threshold={AUTHORITY_THRESHOLD_F}°F)")
    if verdict["optimum"] is not None:
        print(f"  Cluster optimum: scale={verdict['optimum']:.2f}×"
              f"  h_eff={DEFAULT_H * verdict['optimum']:.1f} W/(m²·K)")
    if verdict.get("mix03_delta") is not None:
        sign = "+" if verdict["mix03_delta"] >= 0 else ""
        print(f"  MIX-03: Δ={sign}{verdict['mix03_delta']:.4f}°F"
              f"  generalizes={'YES' if verdict['mix03_generalizes'] else 'NO'}")
    if verdict.get("mix03_hurts_substantially"):
        print("  *** OPEN QUESTION: MIX-03 >+0.2°F regression — surface to user ***")


if __name__ == "__main__":
    main()
