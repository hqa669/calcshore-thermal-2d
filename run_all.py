#!/usr/bin/env python3
"""
Multi-mix thermal validation driver.

Usage:
    python run_all.py [--root PATH] [--mixes MIX-01,MIX-02,...]
                      [--group reference|b1|b2|cluster_a|all|evaluation_set]
                      [--png-dir PATH] [--output-md PATH] [--quiet]

Defaults:
    --root         validation/cw_exports
    --group        evaluation_set
    --output-md    PATH  (required; see --help)
    PNGs are off by default; pass --png-dir to enable.
"""

import sys
import os
import argparse
from pathlib import Path

from compare_to_cw import run_one, print_gate_table

# §7.5.3 group definitions — static; do not add to dataclasses
GROUPS = {
    "reference":      ["MIX-01", "MIX-02", "MIX-03", "MIX-11", "MIX-12"],
    "b1":             ["MIX-04", "MIX-08"],
    "b2":             ["MIX-15"],
    "cluster_a":      ["MIX-05", "MIX-06", "MIX-07", "MIX-09", "MIX-10"],
    "evaluation_set": ["MIX-01", "MIX-02", "MIX-03", "MIX-04", "MIX-08",
                       "MIX-11", "MIX-12", "MIX-15"],
    "all":            [f"MIX-{n:02d}" for n in range(1, 16) if n != 13],
}

# Reverse lookup: mix name → group label (for annotation, first matching group wins)
_MIX_TO_GROUP = {}
for _grp_name in ("reference", "b1", "b2", "cluster_a"):
    for _mix in GROUPS[_grp_name]:
        _MIX_TO_GROUP.setdefault(_mix, _grp_name)


def _group_label(mix: str) -> str:
    return _MIX_TO_GROUP.get(mix, "ungrouped")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default="validation/cw_exports",
                   help="Scenario library root directory (default: validation/cw_exports)")
    p.add_argument("--mixes", default=None,
                   help="Explicit comma-separated mix list, e.g. MIX-01,MIX-02 (overrides --group)")
    p.add_argument("--group", default="evaluation_set",
                   choices=list(GROUPS.keys()),
                   help="Named mix group (default: evaluation_set)")
    p.add_argument("--png-dir", default=None,
                   help="Directory for PNG outputs; if omitted, no PNGs written")
    p.add_argument("--output-md", default=None,
                   help="Path for markdown gate table (REQUIRED). Pass "
                        "validation/sprint4_baseline.md only when intentionally "
                        "regenerating the committed baseline per §6.10 protocol; "
                        "otherwise pass /tmp/<name>.md or another non-committed path.")
    p.add_argument("--quiet", action="store_true",
                   help="Suppress per-mix stdout from print_gate_table")
    args = p.parse_args()
    if args.output_md is None:
        sys.exit("error: --output-md is required. For verification runs, "
                 "pass --output-md /tmp/run_all_output.md. To regenerate "
                 "the committed baseline (§6.10), pass --output-md "
                 "validation/sprint4_baseline.md.")
    return args


def _build_md_table(results: list) -> str:
    header = (
        "| Mix | Group | SCM% | PlaceTemp | PeakMax Δ | PeakGrad Δ"
        " | FieldRMS | CenterRMS | CornerRMS | S0 pass |\n"
        "|---|---|---|---|---|---|---|---|---|---|\n"
    )
    rows = []
    for r in results:
        mix = os.path.basename(r["scenario_dir"])
        grp = r.get("group", "—")
        if r.get("skipped"):
            rows.append(
                f"| {mix} | {grp} | — | — | — | — | — | — | — | skipped ({r.get('skip_reason', '')}) |"
            )
            continue
        meta = r.get("meta", {})
        scm = f"{meta.get('scm_pct', 0.0):.1f}%" if meta else "—"
        pt  = f"{meta.get('placement_temp_F', 0.0):.0f}°F" if meta else "—"
        d   = r["deltas"]
        rm  = r["rms"]
        pm  = f"{d['peak_max_F']:+.2f}°F"
        pg  = f"{d['peak_grad_F']:+.2f}°F"
        frms = f"{rm['field_F']:.2f}°F"  if rm["field_F"]  is not None else "N/A"
        crms = f"{rm['center_F']:.2f}°F" if rm["center_F"] is not None else "N/A"
        corr = f"{rm['corner_F']:.2f}°F" if rm["corner_F"] is not None else "N/A"
        s0   = f"{r['s0_pass_count']}/{r['s0_pass_total']}"
        rows.append(f"| {mix} | {grp} | {scm} | {pt} | {pm} | {pg} | {frms} | {crms} | {corr} | {s0} |")
    return header + "\n".join(rows) + "\n"


def _count_s0_pass(results: list, group: str) -> tuple[int, int]:
    """Return (passing, total) for non-skipped mixes in the given group."""
    passing = total = 0
    for r in results:
        if r.get("group") != group or r.get("skipped"):
            continue
        total += 1
        if r.get("s0_overall"):
            passing += 1
    return passing, total


def _gate_failure_counts(results: list) -> dict:
    counts = {"peak_max": 0, "peak_grad": 0, "field": 0, "center": 0, "corner": 0}
    for r in results:
        if r.get("skipped"):
            continue
        s0 = r["s0"]
        for gate, val in s0.items():
            if val is False:
                counts[gate] += 1
    return counts


def main():
    args = parse_args()

    if args.mixes:
        selected = [m.strip() for m in args.mixes.split(",") if m.strip()]
    else:
        selected = GROUPS[args.group]

    if args.png_dir:
        os.makedirs(args.png_dir, exist_ok=True)

    results = []
    for mix in selected:
        scenario_dir = os.path.join(args.root, mix)
        png_path = (
            os.path.join(args.png_dir, f"cw_comparison_{mix}.png")
            if args.png_dir else None
        )
        result = run_one(scenario_dir, png_path=png_path)
        result["group"] = _group_label(mix)
        if not args.quiet:
            print_gate_table(result)
        results.append(result)

    # MIX-13 has no CW output (awaits re-export, §7.5.4) and is excluded from
    # GROUPS["all"] to keep the group definition clean. When running --group all
    # the committed baseline must still include MIX-13 as a skipped sentinel row
    # (PR 14 contract; tests test_baseline_has_15_rows, test_mix13_skipped_row_present).
    # run_one()'s sentinel handler produces the correct skipped dict; we just need
    # to include it here and insert it in mix-number order.
    if args.group == "all" and not args.mixes:
        mix13_result = run_one(os.path.join(args.root, "MIX-13"))
        mix13_result["group"] = "ungrouped"
        insert_pos = next(
            (i for i, r in enumerate(results)
             if os.path.basename(r["scenario_dir"]) >= "MIX-14"),
            len(results),
        )
        results.insert(insert_pos, mix13_result)

    # ------------------------------------------------------------------ #
    # Write markdown output
    # ------------------------------------------------------------------ #
    md_path = args.output_md
    os.makedirs(os.path.dirname(md_path) if os.path.dirname(md_path) else ".", exist_ok=True)
    Path(md_path).write_text(_build_md_table(results))

    # ------------------------------------------------------------------ #
    # Stdout summary
    # ------------------------------------------------------------------ #
    non_skipped = [r for r in results if not r.get("skipped")]
    n_pass = sum(1 for r in non_skipped if r.get("s0_overall"))
    n_total = len(non_skipped)
    print(f"\nSprint 4 baseline: {n_pass}/{n_total} mixes pass all S0 gates")

    for grp_label, grp_mixes in [("Reference", "reference"), ("B1", "b1"), ("B2", "b2")]:
        grp_results = [r for r in results if r.get("group") == grp_mixes]
        if not grp_results:
            continue
        gp, gt = _count_s0_pass(results, grp_mixes)
        print(f"  {grp_label}: {gp}/{gt}")

    fc = _gate_failure_counts(results)
    if any(fc.values()):
        print("\nPer-gate failure counts (executed set, non-skipped mixes):")
        print(f"  Peak Max:   {fc['peak_max']} fail{'s' if fc['peak_max'] != 1 else ''}")
        print(f"  Peak Grad:  {fc['peak_grad']} fail{'s' if fc['peak_grad'] != 1 else ''}")
        print(f"  Field RMS:  {fc['field']} fail{'s' if fc['field'] != 1 else ''}")
        print(f"  Center RMS: {fc['center']} fail{'s' if fc['center'] != 1 else ''}")
        print(f"  Corner RMS: {fc['corner']} fail{'s' if fc['corner'] != 1 else ''}")

    print(f"\nMarkdown written to {md_path}")
    sys.exit(0 if n_pass == n_total else 1)


if __name__ == "__main__":
    main()
