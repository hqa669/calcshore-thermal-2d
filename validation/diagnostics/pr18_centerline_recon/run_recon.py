"""
PR 18 reconnaissance driver.

Runs the engine ONCE per Reference mix, dispatches to all 5 diagnostics,
writes individual dN_*.md files and a consolidated pr18_summary.md.

Stop condition (§6.6): if any Reference mix breaks S0 5/5, raises
RuntimeError before any output is written.
"""
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _HERE)
sys.path.insert(0, _REPO)

from _harness import load_mix_cache, REFERENCE_MIXES
import d1_constant_offset as d1
import d2_phase_lag as d2
import d3_amplitude_scaling as d3
import d4_depth_profile as d4
import d5_mode_a_overlap as d5

from compare_to_cw import run_one


def _s0_sanity(mix_name: str) -> None:
    scenario_dir = os.path.join(_REPO, "validation", "cw_exports", mix_name)
    r = run_one(scenario_dir)
    if r.get("skipped"):
        raise RuntimeError(f"STOP: {mix_name} skipped — {r.get('skip_reason')}")
    if not r.get("s0_overall"):
        raise RuntimeError(
            f"STOP: {mix_name} broke S0 ({r['s0_pass_count']}/{r['s0_pass_total']} pass). "
            "Reference set must hold S0 5/5. Do not commit."
        )


def _write(path: str, content: str) -> None:
    with open(path, "w") as f:
        f.write(content)
    print(f"  Written: {path}")


def _build_summary(mix_results: dict, wall_times: dict) -> str:
    # Consolidated header
    lines = [
        "# PR 18 — Centerline RMS Reconnaissance Summary",
        "",
        "**Sprint**: 5  **Source commit**: sprint-4-complete (5596c92)",
        f"**Reference mixes**: {', '.join(REFERENCE_MIXES)}",
        "**Gate window**: t ∈ [48, 168] hr  **Mode-A window**: t ∈ [8, 48] hr",
        "",
        "## Per-mix wall times",
        "",
        "| Mix | wall_s |",
        "|---|---|",
    ]
    for mix in REFERENCE_MIXES:
        lines.append(f"| {mix} | {wall_times[mix]:.1f} |")
    lines.append("")

    # Consolidated signature table
    lines += [
        "## Consolidated diagnostic table",
        "",
        "| Mix | D1 mean ΔT | D1 std ΔT | D2 lag_hr | D3 swing_ratio"
        " | D4 mid_rms | D4 spread | D5 mean ΔT | D5 t_at_max | D5 match |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for mix in REFERENCE_MIXES:
        r1 = mix_results[mix]["d1"]
        r2 = mix_results[mix]["d2"]
        r3 = mix_results[mix]["d3"]
        r4 = mix_results[mix]["d4"]
        r5 = mix_results[mix]["d5"]
        match_str = "**Y**" if r5["mode_a_match"] else "N"
        lines.append(
            f"| {mix}"
            f" | {r1['mean_dt_F']:+.3f}"
            f" | {r1['std_dt_F']:.3f}"
            f" | {r2['lag_hr']:+.2f}"
            f" | {r3['ratio']:.4f}"
            f" | {r4['mid_rms_F']:.4f}"
            f" | {r4['spread_F']:.4f}"
            f" | {r5['mean_dt_F']:+.4f}"
            f" | {r5['t_at_max_hr']:.1f}"
            f" | {match_str} |"
        )
    lines.append("")

    # Routing analysis
    lines += ["## Routing analysis", ""]

    # D1: near-constant offset? std < 0.5 * |mean|
    d1_const = all(
        abs(mix_results[m]["d1"]["std_dt_F"]) < 0.5 * abs(mix_results[m]["d1"]["mean_dt_F"])
        for m in REFERENCE_MIXES
    )
    d1_positive = all(mix_results[m]["d1"]["mean_dt_F"] > 0.0 for m in REFERENCE_MIXES)

    # D2: non-trivial phase lag? |lag| > 0.5 hr on any mix
    d2_any_lag = any(abs(mix_results[m]["d2"]["lag_hr"]) > 0.5 for m in REFERENCE_MIXES)
    d2_consistent = (
        len(set(mix_results[m]["d2"]["lag_samples"] > 0 for m in REFERENCE_MIXES)) == 1
    )

    # D3: amplitude ratio departing from 1.0 by > 5% on all mixes
    d3_scaled = all(
        abs(mix_results[m]["d3"]["ratio"] - 1.0) > 0.05 for m in REFERENCE_MIXES
    )

    # D4: depth-localized? spread > 0.2°F on any mix
    d4_localized = any(mix_results[m]["d4"]["spread_F"] > 0.2 for m in REFERENCE_MIXES)

    # D5: Mode-A overlap counts
    mode_a_matches = [m for m in REFERENCE_MIXES if mix_results[m]["d5"]["mode_a_match"]]
    mode_a_n = len(mode_a_matches)

    lines.append("### Per-diagnostic signature summary")
    lines.append("")
    lines.append(f"- **D1 constant offset**: positive bias on all mixes = {d1_positive};"
                 f" near-constant (std < 0.5×|mean|) = {d1_const}")
    lines.append(f"- **D2 phase lag**: any |lag| > 0.5 hr = {d2_any_lag};"
                 f" consistent sign = {d2_consistent}")
    lines.append(f"- **D3 amplitude**: swing ratio departs >5% from 1.0 on all mixes = {d3_scaled}")
    lines.append(f"- **D4 depth profile**: spread > 0.2°F on any mix = {d4_localized}")
    lines.append(f"- **D5 Mode-A**: {mode_a_n}/5 mixes match"
                 f" (need ≥+1.0°F mean + t_at_max ∈ [12,24] hr)")
    if mode_a_n > 0:
        lines.append(f"  - Matching: {', '.join(mode_a_matches)}")
    if mode_a_n < 5 and mode_a_n > 0:
        non_match = [m for m in REFERENCE_MIXES if not mix_results[m]["d5"]["mode_a_match"]]
        lines.append(f"  - Not matching: {', '.join(non_match)}")
    lines.append("")

    # Routing paragraph
    lines.append("### Routing paragraph")
    lines.append("")

    if mode_a_n == 5:
        routing = (
            "All 5 Reference mixes show Mode-A overlap (D5): the centerline residual "
            "carries a positive early-hydration bias of ≥+1.0°F with peak ΔT in the "
            "12–24 hr window, matching the Mode-A signature traced in PR 16 Phase 2.6 "
            "(§7.6.4). "
        )
        d1_means = ", ".join(
            f"{mix_results[m]['d1']['mean_dt_F']:+.2f}°F" for m in REFERENCE_MIXES
        )
        routing += (
            "D1 shows a consistent positive mean offset "
            f"({d1_means}) "
            "that persists into the gate window [48, 168] hr, suggesting Mode-A's early "
            "heat-generation drift carries forward as a DC bias. "
        )
        routing += (
            "PR 18 strengthens R9 routing: the Centerline RMS residual is dominated by "
            "inherited kinetics drift (R9 Mode A), not a distinct thermal-physics term "
            "with an independent knob. "
            "Recommended path: route Centerline RMS to engine v3 release notes as a "
            "known limitation; skip PR 19 fix-conditional unless a targeted early-hydration "
            "recalibration is scoped separately."
        )
    elif 2 <= mode_a_n <= 4:
        routing = (
            f"D5 Mode-A overlap is partial: {mode_a_n}/5 Reference mixes match the "
            "+1.0°F mean / 12–24 hr peak-time signature. "
            f"Matching: {', '.join(mode_a_matches)}. "
        )
        routing += (
            "This is a meaningful finding but does not cleanly route to either "
            "'PR 19 fix-conditional' or 'route to engine v3 release notes.' "
            "D1 shows consistent positive mean offset across all 5 mixes, suggesting a "
            "common DC bias even where Mode-A formal criteria are not met. "
            "USER DECISION REQUIRED: surface per-mix D5 breakdown to user before "
            "finalizing the routing paragraph in the commit message (per PR 18 prompt "
            "open question)."
        )
    elif mode_a_n == 1:
        routing = (
            f"D5 Mode-A overlap is weak: only {mode_a_n}/5 Reference mixes formally match. "
            "The residual does not cleanly route to R9 Mode-A. "
        )
        routing += _fallback_routing(d1_const, d2_any_lag, d3_scaled, d4_localized)
    else:
        routing = (
            "D5 Mode-A overlap is absent: no Reference mix meets the ≥+1.0°F / "
            "12–24 hr peak-time criteria. "
        )
        routing += _fallback_routing(d1_const, d2_any_lag, d3_scaled, d4_localized)

    lines.append(routing)
    lines.append("")

    if 2 <= mode_a_n <= 4:
        lines.append("---")
        lines.append("**Action required**: per-mix D5 breakdown above is ambiguous.")
        lines.append("Consult user before writing final routing paragraph in commit message.")
        lines.append("")

    return "\n".join(lines)


def _fallback_routing(d1_const, d2_any_lag, d3_scaled, d4_localized):
    candidates = []
    if d1_const:
        candidates.append("constant DC offset (D1)")
    if d2_any_lag:
        candidates.append("phase lag (D2)")
    if d3_scaled:
        candidates.append("amplitude scaling (D3)")
    if d4_localized:
        candidates.append("depth-localized residual (D4)")
    if candidates:
        return (
            f"Thermal-physics candidates with clean signatures: "
            f"{', '.join(candidates)}. "
            "Route to PR 19 fix-conditional knob hunt targeting the identified term(s)."
        )
    return (
        "No clean signature found across D1-D5. Characterization is ambiguous. "
        "Route Centerline RMS to engine v3 release notes as a known unresolved residual."
    )


def main():
    t0_total = time.perf_counter()
    mix_results = {}
    wall_times = {}

    print("PR 18 reconnaissance — S0 sanity check + engine runs")
    print("=" * 60)

    for mix in REFERENCE_MIXES:
        print(f"\n{mix}: S0 sanity check...")
        _s0_sanity(mix)
        print(f"{mix}: loading engine cache...")
        t0 = time.perf_counter()
        cache = load_mix_cache(mix)
        wall_times[mix] = time.perf_counter() - t0

        mix_results[mix] = {
            "d1": d1.compute_for_mix(cache),
            "d2": d2.compute_for_mix(cache),
            "d3": d3.compute_for_mix(cache),
            "d4": d4.compute_for_mix(cache),
            "d5": d5.compute_for_mix(cache),
        }
        r5 = mix_results[mix]["d5"]
        print(f"  mode_a window: {r5['window_min_hr']:.1f} hr → {r5['window_max_hr']:.1f} hr"
              f"  (expected [8.0, 48.0])")
        print(f"  wall: {wall_times[mix]:.1f}s")

    print("\nWriting diagnostic tables...")

    diags = [
        ("d1", d1.format_table, "d1_constant_offset.md"),
        ("d2", d2.format_table, "d2_phase_lag.md"),
        ("d3", d3.format_table, "d3_amplitude_scaling.md"),
        ("d4", d4.format_table, "d4_depth_profile.md"),
        ("d5", d5.format_table, "d5_mode_a_overlap.md"),
    ]
    for key, fmt_fn, fname in diags:
        rows = [mix_results[m][key] for m in REFERENCE_MIXES]
        _write(os.path.join(_HERE, fname), fmt_fn(rows))

    summary = _build_summary(mix_results, wall_times)
    _write(os.path.join(_HERE, "pr18_summary.md"), summary)

    print(f"\nTotal wall time: {time.perf_counter() - t0_total:.1f}s")
    print("\nDone. Review pr18_summary.md routing paragraph before committing.")


if __name__ == "__main__":
    main()
