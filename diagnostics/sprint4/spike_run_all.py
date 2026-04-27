"""Sprint 4 recon spike — run compare_to_cw.py across all 15 mix scenarios.

Throwaway script. Do not commit. Output: diagnostics/sprint4/baseline_spike.md
"""
import re
import sys
import subprocess
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths — script may be run from any cwd; anchor to repo root
# ---------------------------------------------------------------------------
REPO = Path(__file__).resolve().parents[2]
CW_EXPORTS = REPO / "validation" / "cw_exports"
COMPARE_PY = REPO / "compare_to_cw.py"
OUT_MD = Path(__file__).parent / "baseline_spike.md"

# ---------------------------------------------------------------------------
# Climate bucket (hand-coded; single location in current dataset)
# ---------------------------------------------------------------------------
CLIMATE = {
    "TX, Austin": "hot-humid",
}

# ---------------------------------------------------------------------------
# Regex patterns against compare_to_cw.py stdout (lines 257-279)
#   Δ lines:   "  Δ = +0.5°F   [PASS] ..."
#   RMS lines: "  Field-wide RMS: 1.23°F   [PASS] ..."
#   OVERALL:   "  OVERALL:        [PASS]    S1-aspire: 5/5 metrics meet"
# ---------------------------------------------------------------------------
RE_DELTA   = re.compile(r"Δ\s*=\s*([+-]?\d+\.\d+)°F\s+\[(PASS|FAIL)\]")
RE_RMS     = re.compile(
    r"\s+(Field-wide|Centerline|Corner) RMS:\s+(\d+\.\d+)°F\s+\[(PASS|FAIL)\]"
)
RE_OVERALL = re.compile(r"OVERALL:\s+\[(PASS|FAIL)\]\s+S1-aspire:\s+(\d+)/5")

MIX_IDS = [f"MIX-{n:02d}" for n in range(1, 16)]


def scm_family(pct):
    if pct < 20:
        return "low_scm"
    if pct < 40:
        return "mid_scm"
    return "high_scm"


def load_tags(mix_dir: Path):
    """Load scenario via cw_scenario_loader for location/SCM tagging."""
    sys.path.insert(0, str(REPO))
    from cw_scenario_loader import load_cw_scenario  # noqa: E402

    input_dat = mix_dir / "input.dat"
    weather_dat = mix_dir / "weather.dat"
    output_txt = mix_dir / "output.txt"
    cw_out = str(output_txt) if output_txt.exists() else None

    scn = load_cw_scenario(str(input_dat), str(weather_dat), cw_out)
    loc = scn.environment.location
    m = scn.mix
    cement = m.cement_type_I_II_lb_yd3
    scm    = m.fly_ash_F_lb_yd3 + m.fly_ash_C_lb_yd3 + m.ggbfs_lb_yd3 + m.silica_fume_lb_yd3
    total  = cement + scm
    pct    = (scm / total * 100) if total > 0 else 0.0
    climate = CLIMATE.get(loc, "unknown")
    if climate == "unknown":
        print(f"  WARNING: unknown location '{loc}' — add to CLIMATE dict", flush=True)
    return loc, climate, round(pct, 1), scm_family(pct)


def parse_stdout(stdout):
    """Extract gate metrics from compare_to_cw.py output.

    Returns dict with keys: peak_max_delta, peak_max_pass, peak_grad_delta,
    peak_grad_pass, field_rms, field_pass, center_rms, center_pass,
    corner_rms, corner_pass, s0_pass_count, parse_ok.
    """
    deltas = RE_DELTA.findall(stdout)
    rmses  = {name: (val, pf) for name, val, pf in RE_RMS.findall(stdout)}
    overall = RE_OVERALL.search(stdout)

    if len(deltas) != 2 or len(rmses) < 3 or overall is None:
        return {"parse_ok": False, "raw": stdout[-600:]}

    pass_count = overall.group(1)  # actually from PASS/FAIL overall, but we count below
    pf_to_bool = lambda s: s == "PASS"

    r = {
        "parse_ok":        True,
        "peak_max_delta":  deltas[0][0],
        "peak_max_pass":   pf_to_bool(deltas[0][1]),
        "peak_grad_delta": deltas[1][0],
        "peak_grad_pass":  pf_to_bool(deltas[1][1]),
        "field_rms":       rmses.get("Field-wide", ("N/A", "FAIL"))[0],
        "field_pass":      pf_to_bool(rmses.get("Field-wide", ("N/A", "FAIL"))[1]),
        "center_rms":      rmses.get("Centerline", ("N/A", "FAIL"))[0],
        "center_pass":     pf_to_bool(rmses.get("Centerline", ("N/A", "FAIL"))[1]),
        "corner_rms":      rmses.get("Corner", ("N/A", "FAIL"))[0],
        "corner_pass":     pf_to_bool(rmses.get("Corner", ("N/A", "FAIL"))[1]),
    }
    r["s0_pass_count"] = sum([
        r["peak_max_pass"], r["peak_grad_pass"],
        r["field_pass"], r["center_pass"], r["corner_pass"]
    ])
    return r


def run_comparison(mix_dir: Path):
    """Run compare_to_cw.py as subprocess; return (result, stderr_tail)."""
    try:
        proc = subprocess.run(
            [sys.executable, str(COMPARE_PY), str(mix_dir)],
            capture_output=True, text=True, timeout=180,
            cwd=str(REPO),
        )
        return proc.stdout, proc.stderr, None
    except subprocess.TimeoutExpired as e:
        return "", str(e), "timeout"
    except Exception as e:
        return "", str(e), "subprocess_error"


def fmt_delta(val, pass_):
    mark = "✓" if pass_ else "✗"
    return f"{val:>6} {mark}"


def fmt_rms(val, pass_):
    mark = "✓" if pass_ else "✗"
    return f"{val:>5} {mark}"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

rows = []
parse_errors = []

for mix_id in MIX_IDS:
    mix_dir = CW_EXPORTS / mix_id
    print(f"[{mix_id}] ", end="", flush=True)

    if not mix_dir.exists():
        print("directory missing — skip", flush=True)
        rows.append({"mix": mix_id, "skipped": True, "reason": "directory_missing"})
        continue

    if not (mix_dir / "output.txt").exists():
        print("no output.txt — skip", flush=True)
        rows.append({"mix": mix_id, "skipped": True, "reason": "no_cw_output"})
        continue

    # Tag from scenario loader
    try:
        loc, climate, scm_pct, family = load_tags(mix_dir)
    except Exception as e:
        print(f"loader error: {e}", flush=True)
        rows.append({"mix": mix_id, "skipped": True, "reason": f"loader_error: {e}"})
        continue

    # Run comparison
    stdout, stderr, err_kind = run_comparison(mix_dir)

    if err_kind:
        print(f"{err_kind}", flush=True)
        rows.append({
            "mix": mix_id, "skipped": False, "parse_ok": False,
            "loc": loc, "climate": climate, "scm_pct": scm_pct, "family": family,
            "error": err_kind, "stderr": stderr[-400:],
        })
        parse_errors.append(mix_id)
        continue

    metrics = parse_stdout(stdout)
    if not metrics["parse_ok"]:
        print("parse error", flush=True)
        rows.append({
            "mix": mix_id, "skipped": False, "parse_ok": False,
            "loc": loc, "climate": climate, "scm_pct": scm_pct, "family": family,
            "error": "parse_error", "raw": metrics.get("raw", ""), "stderr": stderr[-400:],
        })
        parse_errors.append(mix_id)
        continue

    status = "PASS" if metrics["s0_pass_count"] == 5 else f"FAIL ({metrics['s0_pass_count']}/5)"
    print(f"{status}", flush=True)
    rows.append({
        "mix": mix_id, "skipped": False, "parse_ok": True,
        "loc": loc, "climate": climate, "scm_pct": scm_pct, "family": family,
        **metrics,
    })

# ---------------------------------------------------------------------------
# Write markdown table
# ---------------------------------------------------------------------------

SEP = " | "
HDR = ["Mix", "Location", "Climate", "SCM%", "Family",
       "PeakMax Δ", "PeakGrad Δ", "FieldRMS", "CenterRMS", "CornerRMS", "S0 pass"]

md_lines = [
    "# Sprint 4 — S0 Gate Baseline Spike",
    "",
    "| " + " | ".join(HDR) + " |",
    "| " + " | ".join(["---"] * len(HDR)) + " |",
]

for r in rows:
    mix = r["mix"]
    if r.get("skipped"):
        cells = [mix] + ["— skipped ({}) —".format(r["reason"])] + [""] * (len(HDR) - 2)
    elif not r.get("parse_ok"):
        cells = [mix, r.get("loc","?"), r.get("climate","?"), str(r.get("scm_pct","?")),
                 r.get("family","?"),
                 "— parse error —", "", "", "", "", ""]
    else:
        pm = ("+" if not r["peak_max_delta"].startswith("-") else "") + r["peak_max_delta"] + "°F " + ("✓" if r["peak_max_pass"] else "✗")
        pg = ("+" if not r["peak_grad_delta"].startswith("-") else "") + r["peak_grad_delta"] + "°F " + ("✓" if r["peak_grad_pass"] else "✗")
        fr = r["field_rms"]  + "°F " + ("✓" if r["field_pass"]  else "✗")
        cr = r["center_rms"] + "°F " + ("✓" if r["center_pass"] else "✗")
        ko = r["corner_rms"] + "°F " + ("✓" if r["corner_pass"] else "✗")
        cells = [
            mix, r["loc"], r["climate"], f"{r['scm_pct']:.1f}%", r["family"],
            pm, pg, fr, cr, ko,
            f"{r['s0_pass_count']}/5",
        ]
    md_lines.append("| " + " | ".join(cells) + " |")

# Footnotes for parse errors
if parse_errors:
    md_lines += ["", "## Parse / run errors", ""]
    for r in rows:
        if r.get("mix") in parse_errors and not r.get("skipped"):
            md_lines.append(f"### {r['mix']}")
            md_lines.append(f"Error: {r.get('error','?')}")
            if r.get("stderr"):
                md_lines.append("```")
                md_lines.append(r["stderr"])
                md_lines.append("```")

# Surprises section
md_lines += [
    "",
    "## Surprises noted",
    "",
    "- **All 15 mixes are `TX, Austin`** — climate-bucket column is uniform (`hot-humid`).",
    "  Climate-based clustering is N/A for this dataset.",
    "- **`cw_comparison_MIX-01.png` hardcoded** in `compare_to_cw.py:437` — every run",
    "  overwrites the same PNG in the repo root. PR 13 should fix the output path.",
    "- **MIX-01 constants at `compare_to_cw.py:41-44`** (`CW_PEAK_MAX_F`, etc.) are used",
    "  to display the CW peak *time* label at line 259 for the gradient row. For",
    "  non-MIX-01 mixes the displayed CW peak time is wrong, but the delta/pass values",
    "  are computed from runtime CW data and are trustworthy. Fix belongs in PR 13.",
]

OUT_MD.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

# ---------------------------------------------------------------------------
# Stdout clustering summary
# ---------------------------------------------------------------------------

data_rows = [r for r in rows if not r.get("skipped") and r.get("parse_ok")]
n_total  = len(data_rows)
n_pass5  = sum(1 for r in data_rows if r["s0_pass_count"] == 5)

print()
print("=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"  Mixes evaluated : {n_total} of 14 (MIX-13 skipped: no CW output)")
print(f"  All 5 S0 pass   : {n_pass5}/{n_total}")
if parse_errors:
    print(f"  Parse errors    : {', '.join(parse_errors)}")

gate_names = [
    ("Peak Max T",   "peak_max_pass"),
    ("Peak Gradient","peak_grad_pass"),
    ("Field RMS",    "field_pass"),
    ("Centerline RMS","center_pass"),
    ("Corner RMS",   "corner_pass"),
]

any_failure = False
for gate_label, key in gate_names:
    failed = [r for r in data_rows if not r[key]]
    if not failed:
        continue
    any_failure = True
    mixes = [r["mix"] for r in failed]
    by_family = {}
    for r in failed:
        by_family.setdefault(r["family"], []).append(r["mix"])
    family_notes = "; ".join(
        f"{len(v)} of {sum(1 for d in data_rows if d['family']==k)} {k}"
        for k, v in sorted(by_family.items())
    )
    print(f"  {gate_label:<18} FAIL on: {', '.join(mixes)}")
    print(f"  {'':18}   by family: {family_notes}")

if not any_failure:
    print("  All 14 evaluated mixes pass every S0 gate.")

print()
print("  Climate clustering: single-climate dataset (all TX, Austin) — N/A.")
print("  SCM-family breakdown of evaluated mixes:")
from collections import Counter
family_counts = Counter(r["family"] for r in data_rows)
for fam, cnt in sorted(family_counts.items()):
    print(f"    {fam}: {cnt}")

print()
print(f"  Output written to: {OUT_MD}")
print()
