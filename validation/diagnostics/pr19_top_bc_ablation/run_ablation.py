"""
PR 19 — Top-surface BC ablation driver.

Runs s1 → s3 → s4 unconditionally; invokes s2 conditionally on s1's
wait-point verdict. Writes consolidated pr19_summary.md.

Usage (from repo root or from this directory):
    python validation/diagnostics/pr19_top_bc_ablation/run_ablation.py
"""
import os
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)

import s1_h_conv_magnitude as s1
import s2_h_conv_form as s2
import s3_alpha_top_magnitude as s3
import s4_epsilon_top_magnitude as s4

from s1_h_conv_magnitude import DEFAULT_H as _DEFAULT_H

_OUT = os.path.join(_HERE, "pr19_summary.md")


def main():
    t_start = time.perf_counter()
    print("=" * 60)
    print("PR 19 — Top-surface BC ablation")
    print("Knobs: h_conv × α_top × ε_top  |  Cluster: MIX-01/11/12  |  Control: MIX-03")
    print("=" * 60)

    # ---- s1: h_conv magnitude ----
    print("\n" + "─" * 40)
    s1_rows, s1_verdict, s1_stop = s1.main()
    print()

    # ---- s2: h_conv form (conditional) ----
    print("─" * 40)
    s2_skip_reason = None
    if s1_stop:
        s2_skip_reason = f"s1 stop condition triggered ({s1_stop[0]} S0 failure at scale={s1_stop[1]:.2f}×)"
    elif not s1_verdict["has_authority"]:
        s2_skip_reason = "no h_conv magnitude authority from s1 (variation < 0.24°F)"
    elif s1_verdict.get("mix03_hurts_substantially"):
        s2_skip_reason = (
            "open-question MIX-03 regression (>+0.2°F at cluster optimum) — "
            "resolve before form-check"
        )

    if s2_skip_reason:
        s2_rows, s2_skipped = s2.main(skip=True, skip_reason=s2_skip_reason)
    else:
        calibrated_h = _DEFAULT_H * s1_verdict["optimum"]
        s2_rows, s2_skipped = s2.main(calibrated_h=calibrated_h)
    print()

    # ---- s3: α_top magnitude ----
    print("─" * 40)
    s3_rows, s3_verdict, s3_stop = s3.main()
    print()

    # ---- s4: ε_top magnitude ----
    print("─" * 40)
    s4_rows, s4_verdict, s4_stop = s4.main()
    print()

    # ---- Consolidated summary ----
    t_total = time.perf_counter() - t_start
    _write_summary(
        s1_rows=s1_rows, s1_verdict=s1_verdict, s1_stop=s1_stop,
        s2_rows=s2_rows, s2_skipped=s2_skipped,
        s3_rows=s3_rows, s3_verdict=s3_verdict, s3_stop=s3_stop,
        s4_rows=s4_rows, s4_verdict=s4_verdict, s4_stop=s4_stop,
        t_total_s=t_total,
    )
    print(f"\nTotal wall time: {t_total:.0f}s")
    print(f"Written: {_OUT}")

    # Check any stop condition
    any_stop = s1_stop or s3_stop or s4_stop
    if any_stop:
        print("\n*** STOP CONDITION TRIGGERED — do not commit. "
              "Investigate before continuing. ***")
    else:
        print("\nNo stop conditions. Review pr19_summary.md routing paragraph "
              "before committing.")


def _write_summary(*, s1_rows, s1_verdict, s1_stop,
                   s2_rows, s2_skipped,
                   s3_rows, s3_verdict, s3_stop,
                   s4_rows, s4_verdict, s4_stop,
                   t_total_s):

    from s1_h_conv_magnitude import DEFAULT_H
    from s3_alpha_top_magnitude import DEFAULT_ALPHA
    from s4_epsilon_top_magnitude import DEFAULT_EPS

    any_stop = s1_stop or s3_stop or s4_stop

    lines = [
        "# PR 19 — Top-Surface BC Ablation Summary",
        "",
        "**Sprint**: 5  **Evaluation set**: MIX-01, MIX-11, MIX-12 (cluster)  +  MIX-03 (control)",
        f"**Total wall time**: {t_total_s:.0f}s",
        "**Baseline residual** (from PR 18): D3 ratio ~1.90×, D4 top-concentrated,",
        "cluster CenterRMS ≈ 0.79–0.84°F (target: <0.50°F S1-aspire).",
        "",
    ]

    if any_stop:
        triggers = []
        if s1_stop:
            triggers.append(f"s1: {s1_stop[0]} S0 failure at scale={s1_stop[1]:.2f}×")
        if s3_stop:
            triggers.append(f"s3: {s3_stop[0]} S0 failure at α_top={s3_stop[1]:.2f}")
        if s4_stop:
            triggers.append(f"s4: {s4_stop[0]} S0 failure at ε_top={s4_stop[1]:.2f}")
        lines += [
            "⚠️ **STOP CONDITION(S) TRIGGERED — do not commit:**",
        ]
        for t in triggers:
            lines.append(f"  - {t}")
        lines += [
            "Investigate these sweep points before proceeding to PR 20.",
            "",
        ]

    # ── Knob 1: h_conv ──
    lines += [
        "## Knob 1 — h_conv magnitude",
        "",
        f"Default: h_eff = {DEFAULT_H:.2f} W/(m²·K) at v=10.5 m/s  [ACI 305 linear-in-wind]",
        f"Sweep: 0.5× – 2.0× ({DEFAULT_H*0.5:.1f} – {DEFAULT_H*2.0:.1f} W/m²K)",
        "",
        f"- Cluster-mean CenterRMS variation: {s1_verdict['variation_F']:.4f}°F",
        f"- **Authority (≥0.24°F): {'YES' if s1_verdict['has_authority'] else 'NO'}**",
    ]
    if s1_verdict["optimum"] is not None:
        lines.append(
            f"- Cluster optimum: scale={s1_verdict['optimum']:.2f}×"
            f"  h_eff={DEFAULT_H*s1_verdict['optimum']:.1f} W/(m²·K)"
            f"  cluster-mean CenterRMS={s1_verdict['cluster_mean_at_optimum']:.4f}°F"
        )
    if s1_verdict.get("mix03_delta") is not None:
        sign = "+" if s1_verdict["mix03_delta"] >= 0 else ""
        lines.append(
            f"- MIX-03 at optimum: baseline={s1_verdict['mix03_baseline_rms']:.4f}°F"
            f" → {s1_verdict['mix03_optimum_rms']:.4f}°F"
            f"  Δ={sign}{s1_verdict['mix03_delta']:.4f}°F"
            f"  generalizes={'YES' if s1_verdict['mix03_generalizes'] else 'NO'}"
        )
    if s1_verdict.get("mix03_hurts_substantially"):
        lines.append(
            "- ⚠️ OPEN QUESTION: MIX-03 regression >+0.2°F at cluster optimum "
            "— surface to user."
        )
    lines.append("")

    # ── Knob 1b: h_conv form ──
    lines += ["## Knob 1b — h_conv functional-form check", ""]
    if s2_skipped:
        lines.append("Status: **SKIPPED** — see s2_h_conv_form.md for reason.")
    else:
        lines.append("Status: **RUN** — see s2_h_conv_form.md for per-form CenterRMS table.")
        if s2_rows:
            from collections import defaultdict
            by_form = defaultdict(list)
            for (mix_id, form, m) in s2_rows:
                from _harness import CLUSTER
                if mix_id in CLUSTER and m.get("center_rms") is not None:
                    by_form[form].append(m["center_rms"])
            for form_key in ("A", "B", "C"):
                vals = by_form.get(form_key, [])
                if vals:
                    lines.append(
                        f"  Form {form_key}: cluster-mean CenterRMS = {sum(vals)/len(vals):.4f}°F"
                    )
    lines.append("")

    # ── Knob 2: α_top ──
    lines += [
        "## Knob 2 — α_top (solar absorptivity)",
        "",
        f"Default: α_top = {DEFAULT_ALPHA}  Sweep: 0.30 – 0.85",
        "",
        f"- Cluster-mean CenterRMS variation: {s3_verdict['variation_F']:.4f}°F",
        f"- **Authority (≥0.24°F): {'YES' if s3_verdict['has_authority'] else 'NO'}**",
    ]
    if s3_verdict["optimum"] is not None:
        lines.append(
            f"- Cluster optimum: α_top={s3_verdict['optimum']:.2f}"
            f"  cluster-mean CenterRMS={s3_verdict['cluster_mean_at_optimum']:.4f}°F"
        )
    if s3_verdict.get("mix03_delta") is not None:
        sign = "+" if s3_verdict["mix03_delta"] >= 0 else ""
        lines.append(
            f"- MIX-03 Δ at optimum: {sign}{s3_verdict['mix03_delta']:.4f}°F"
            f"  generalizes={'YES' if s3_verdict['mix03_generalizes'] else 'NO'}"
        )
    lines.append("")

    # ── Knob 3: ε_top ──
    lines += [
        "## Knob 3 — ε_top (LW emissivity)",
        "",
        f"Default: ε_top = {DEFAULT_EPS}  Sweep: 0.50 – 0.95",
        "",
        f"- Cluster-mean CenterRMS variation: {s4_verdict['variation_F']:.4f}°F",
        f"- **Authority (≥0.24°F): {'YES' if s4_verdict['has_authority'] else 'NO'}**",
    ]
    if s4_verdict["optimum"] is not None:
        lines.append(
            f"- Cluster optimum: ε_top={s4_verdict['optimum']:.2f}"
            f"  cluster-mean CenterRMS={s4_verdict['cluster_mean_at_optimum']:.4f}°F"
        )
    if s4_verdict.get("mix03_delta") is not None:
        sign = "+" if s4_verdict["mix03_delta"] >= 0 else ""
        lines.append(
            f"- MIX-03 Δ at optimum: {sign}{s4_verdict['mix03_delta']:.4f}°F"
            f"  generalizes={'YES' if s4_verdict['mix03_generalizes'] else 'NO'}"
        )
    lines.append("")

    # ── Routing paragraph ──
    lines += ["## Routing for PR 20", ""]
    routing = _build_routing(
        s1_verdict=s1_verdict, s2_skipped=s2_skipped, s2_rows=s2_rows,
        s3_verdict=s3_verdict, s4_verdict=s4_verdict,
        any_stop=any_stop,
    )
    lines.append(routing)
    lines.append("")

    # ── R9 implication ──
    lines += ["## R9 implication", ""]
    r9 = _build_r9_implication(s1_verdict, s3_verdict, s4_verdict)
    lines.append(r9)
    lines.append("")

    with open(_OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def _build_routing(*, s1_verdict, s2_skipped, s2_rows,
                   s3_verdict, s4_verdict, any_stop):
    if any_stop:
        return (
            "**Routing DEFERRED — stop condition(s) triggered.** Investigate the "
            "S0 failures before determining PR 20 shape. Do not commit until resolved."
        )

    # Find which knobs have authority + generalization
    knobs_with_authority = []
    knobs_authority_no_gen = []

    for label, v in [("h_conv", s1_verdict), ("α_top", s3_verdict), ("ε_top", s4_verdict)]:
        if v["has_authority"]:
            if v.get("mix03_generalizes"):
                knobs_with_authority.append((label, v))
            elif v.get("mix03_generalizes") is False:
                knobs_authority_no_gen.append((label, v))
            # else: mix03 data unavailable

    open_q = s1_verdict.get("mix03_hurts_substantially")

    if open_q:
        return (
            "**OPEN QUESTION — surface to user before routing.** h_conv has cluster "
            "authority but MIX-03 regresses >+0.2°F at the cluster optimum. "
            "The BC mechanism appears to interact with kinetics variation. "
            "User must decide whether MIX-03's regression is acceptable cost for "
            "cluster S1-aspire closure, or whether the BC fix should be deferred "
            "until kinetics calibration generalizes. "
            "If deferral: PR 20 routes h_conv recalibration to engine v3 release notes "
            "as a documented residual."
        )

    if knobs_with_authority:
        names = " + ".join(label for label, _ in knobs_with_authority)
        details = "; ".join(
            f"{label}: optimum={v['optimum']:.2f}, "
            f"cluster CenterRMS {v['cluster_mean_at_optimum']:.4f}°F"
            for label, v in knobs_with_authority
        )
        return (
            f"**PR 20 shape: fix-conditional commit.** Knob(s) with cluster authority "
            f"+ MIX-03 generalization: {names}. "
            f"Details — {details}. "
            f"PR 20 recalibrates the leading knob and runs the full S0 gate to confirm. "
            f"If the recalibrated value passes gate + achieves cluster S1-aspire on ≥2/3 "
            f"mixes, commit to engine defaults."
        )

    if knobs_authority_no_gen:
        names = ", ".join(label for label, _ in knobs_authority_no_gen)
        return (
            f"**PR 20 shape: route to engine v3 release notes.** Knob(s) with cluster "
            f"authority but WITHOUT MIX-03 generalization: {names}. "
            f"The BC mechanism interacts with kinetics variation; a simple scalar "
            f"recalibration cannot close the residual without regressing kinetics-"
            f"divergent mixes. PR 20 documents this as a known residual attributable "
            f"to kinetics-BC interaction and flags it for engine v3."
        )

    # No knob has authority
    return (
        "**PR 20 shape: route to engine v3 release notes as a documented residual.** "
        "None of the three top-surface BC knobs (h_conv, α_top, ε_top) produced "
        f"≥{0.24}°F authority on cluster-mean CenterRMS. The BC amplitude mismatch "
        "identified in PR 18 (D3 ~1.90×, D4 top-concentrated) is not attributable to "
        "any single top-surface scalar BC parameter within physically reasonable ranges. "
        "PR 20 documents this finding and escalates to a deeper investigation "
        "(engine v3: revised BC physics or additional parameter interactions)."
    )


def _build_r9_implication(s1_verdict, s3_verdict, s4_verdict):
    any_authority = (
        s1_verdict["has_authority"] or
        s3_verdict["has_authority"] or
        s4_verdict["has_authority"]
    )
    if any_authority:
        return (
            "If a BC knob closes the amplitude residual (D3 ratio → 1.0, "
            "D4 spread → 0), the Mode-A DC offset (D1 positive mean ΔT +0.40–0.56°F, "
            "D5 peak timing 39–47 hr) will remain as a separate residual component. "
            "This would strengthen the R9 evidence: the warm-placement Mode-A shift "
            "(PR 18 D5 partial-match finding) is a kinetics-driven DC offset independent "
            "of the BC amplitude mechanism. PR 20 should report D1 mean ΔT at the "
            "recalibrated BC setting as the post-fix R9 residual estimate."
        )
    return (
        "No BC knob authority found. PR 18's D1 DC offset (+0.40–0.56°F cluster mean) "
        "and D3 amplitude mismatch remain entangled. R9 Mode-A DC offset cannot be "
        "cleanly separated from the BC amplitude residual without first closing D3. "
        "R9 assessment is unchanged from PR 18; deferred to engine v3."
    )


if __name__ == "__main__":
    main()
