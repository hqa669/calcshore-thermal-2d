"""
PR 20 / Item 1 — R_form Reference-set sweep.

Evaluates R_FORM_CONTACT_SI in [0.050, 0.060, 0.070, 0.0862, 0.100] against
the full Reference set (5 mixes × 5 values = 25 runs).

Primary decision gate: CornerRMS (PR 15 found R_form has clean authority;
S0 threshold 3.0°F). Secondary: PeakGrad. Sanity check: CenterRMS.

Commit-decision criteria (Item 1):
  C1: all 5 Reference mixes S0 5/5 at R_form=0.060
  C2: ≥3/5 mixes improve CornerRMS at 0.060 vs 0.0862
  C3: no mix CenterRMS degrades >0.05°F at 0.060 vs 0.0862

C3 violation → STOP condition (unexpected; R_form should be form-face-localized).
Surface to user before committing if C3 fails.
"""
import os
import sys
import time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

from _harness import (
    REFERENCE, CLUSTER, R_FORM_DEFAULT, R_FORM_CANDIDATE,
    patched_r_form, load_scenarios, run_engine, compute_metrics,
    format_sweep_table, evaluate_commit_criteria,
    CENTER_RMS_DEGRADE_MAX_F,
)

R_FORM_VALUES = [0.050, 0.060, 0.070, 0.0862, 0.100]
_OUT_MD = os.path.join(_HERE, "s1_r_form_sweep.md")


def main():
    print("=== s1 — R_form Reference-set sweep ===")
    print(f"  Values: {R_FORM_VALUES}")
    print(f"  Mixes:  {list(REFERENCE)}")
    print(f"  Total:  {len(R_FORM_VALUES) * len(REFERENCE)} runs\n")

    print("Loading scenarios ...", flush=True)
    scenarios = load_scenarios()
    print(f"  Loaded: {list(scenarios)}\n")

    rows = []
    for rval in R_FORM_VALUES:
        print(f"  R_form = {rval:.4f} m²K/W", flush=True)
        for mix_id in REFERENCE:
            t0 = time.perf_counter()
            scn = scenarios[mix_id]
            with patched_r_form(rval):
                result, grid = run_engine(scn)
            wall = time.perf_counter() - t0

            if not np.all(np.isfinite(result.T_field_C)):
                print(f"    [{mix_id}] NaN in T_field_C — STOP")
                sys.exit(1)

            m = compute_metrics(scn, grid, result)
            rows.append((mix_id, rval, m))

            crms   = m.get("corner_rms",   float("nan"))
            cenrms = m.get("center_rms",   float("nan"))
            pgrad  = m.get("peak_grad_delta", float("nan"))
            s0     = m["s0_pass"]
            print(
                f"    [{mix_id}]  CornerRMS={crms:.4f}°F  "
                f"CenterRMS={cenrms:.4f}°F  PeakGradΔ={pgrad:+.3f}°F  "
                f"S0={s0}/5  ({wall:.1f}s)"
            )
        print()

    decision = evaluate_commit_criteria(rows)
    _write_md(rows, decision)
    _print_decision(decision)

    if decision["stop_c3"]:
        print("\n*** STOP — C3 violation: unexpected CenterRMS degradation >0.05°F at R_form=0.060.")
        print("    R_form should be form-face-localized; this signals something unexpected.")
        print("    Do NOT commit. Surface to user before deciding.")
        sys.exit(1)

    return rows, decision


def _write_md(rows, decision):
    d = decision
    lines = [
        "# s1 — R_form Reference-Set Sweep",
        "",
        "**Knob**: `R_FORM_CONTACT_SI` (module-level constant, `thermal_engine_2d.py:127`).",
        "Monkey-patch via context manager; all 14 in-engine sites perturbed uniformly",
        "(code sites: lines 524, 1913, 1914, 1935, 1944, 1949, 2020, 2192, 2193, 2212, 2221).",
        "",
        f"**Sweep values**: {R_FORM_VALUES}",
        f"**Evaluation set**: {list(REFERENCE)}",
        "MIX-02 reported alongside cluster mixes; kinetics-anomalous per `mix02_recon.md`.",
        "",
        "## Sweep table",
        "",
        "Primary gate: **CornerRMS** (PR 15 finding; S0 threshold 3.0°F).",
        "Secondary: **PeakGrad Δ** (also moves with R_form; S0 threshold 2.0°F).",
        "Sanity check: **CenterRMS** (should be form-face-decoupled; S0 threshold 1.0°F).",
        "",
        format_sweep_table(rows),
        "",
        "## Commit-decision criteria (Item 1)",
        "",
    ]

    # C1
    c1_tag = "PASS ✓" if d["c1_pass"] else "FAIL ✗"
    lines.append(f"**C1** — all 5 Reference mixes S0 5/5 at R_form=0.060: **{c1_tag}**")
    for mid in REFERENCE:
        v = d["c1_per_mix"].get(mid)
        tag = "5/5 ✓" if v else "FAIL ✗"
        lines.append(f"  - {mid}: S0 {tag}")
    lines.append("")

    # C2
    c2_tag = "PASS ✓" if d["c2_pass"] else "FAIL ✗"
    n_improved = sum(1 for v in d["c2_per_mix"].values() if v is True)
    lines.append(
        f"**C2** — ≥3/5 mixes improve CornerRMS at 0.060 vs 0.0862: "
        f"**{c2_tag}** ({n_improved}/5 improve)"
    )
    for mid in REFERENCE:
        v = d["c2_per_mix"].get(mid)
        tag = "✓ improves" if v is True else ("✗ no improvement" if v is False else "N/A")
        lines.append(f"  - {mid}: {tag}")
    lines.append("")

    # C3
    c3_tag = "PASS ✓" if d["c3_pass"] else "FAIL ✗ — STOP"
    lines.append(
        f"**C3** — no mix CenterRMS degrades >0.05°F at 0.060 vs 0.0862: **{c3_tag}**"
    )
    for mid in REFERENCE:
        info = d["c3_per_mix"].get(mid, {})
        delta = info.get("delta")
        if delta is not None:
            stop_note = " ← STOP" if info.get("fail") else ""
            lines.append(f"  - {mid}: Δ={delta:+.4f}°F{stop_note}")
        else:
            lines.append(f"  - {mid}: N/A")
    lines.append("")

    # Outcome
    if d["commit"]:
        outcome = "**COMMIT** R_form=0.060 — all three criteria pass."
    elif d["stop_c3"]:
        outcome = ("**STOP — DO NOT COMMIT** — C3 violation (unexpected CenterRMS degradation). "
                   "Surface to user.")
    else:
        fails = []
        if not d["c1_pass"]:
            fails.append("C1")
        if not d["c2_pass"]:
            fails.append("C2")
        outcome = f"**DO NOT COMMIT** — {', '.join(fails)} failed. Route to engine v3 release notes."
    lines += [f"**Outcome**: {outcome}", ""]

    with open(_OUT_MD, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"  Written: {_OUT_MD}")


def _print_decision(d):
    print("--- Commit-decision summary ---")
    print(f"  C1 (all 5 S0 5/5 at 0.060):                  {'PASS' if d['c1_pass'] else 'FAIL'}")
    print(f"  C2 (≥3/5 improve CornerRMS vs 0.0862):        {'PASS' if d['c2_pass'] else 'FAIL'}")
    print(f"  C3 (no CenterRMS degrade >0.05°F vs 0.0862):  {'PASS' if d['c3_pass'] else 'FAIL (STOP)'}")
    if d["commit"]:
        print("  Outcome: COMMIT R_form=0.060")
    else:
        print("  Outcome: DO NOT COMMIT")


if __name__ == "__main__":
    main()
