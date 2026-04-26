"""
D5 — Mode-A overlap.
Does the Reference centerline residual share the Mode-A signature from
PR 16 Phase 2.6: engine warmer than CW by +1.0-1.7°F during active
hydration t = 8-48 hr, with peak ΔT in the t = 12-24 hr range?

Window mask is t ∈ [8, 48] hr — BELOW the gate window. The script
prints the actual mask range to stdout for verification.
"""
import os
import sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _HERE)
sys.path.insert(0, _REPO)

from _harness import load_mix_cache, REFERENCE_MIXES

# Mode-A match thresholds (per §7.6.4 Phase 2.6)
MODE_A_MEAN_MIN_F = 1.0     # mean ΔT must be ≥ +1.0°F
MODE_A_PEAK_LO_HR = 12.0    # time of max ΔT must be in [12, 24] hr
MODE_A_PEAK_HI_HR = 24.0


def compute_for_mix(cache) -> dict:
    mask = cache.mode_a_mask
    t_min = float(cache.cw_time_hrs[mask].min())
    t_max = float(cache.cw_time_hrs[mask].max())

    diff = (cache.eng_center_F_at_cw_t[mask]
            - cache.cw_center_F[mask])
    mean_dt   = float(np.mean(diff))
    max_dt    = float(np.max(diff))
    t_at_max  = float(cache.cw_time_hrs[mask][int(np.argmax(diff))])

    match = (mean_dt >= MODE_A_MEAN_MIN_F
             and MODE_A_PEAK_LO_HR <= t_at_max <= MODE_A_PEAK_HI_HR)

    return {
        "mix": cache.mix_name,
        "window_min_hr": t_min,
        "window_max_hr": t_max,
        "mean_dt_F": mean_dt,
        "max_dt_F": max_dt,
        "t_at_max_hr": t_at_max,
        "mode_a_match": match,
    }


def format_table(rows: list) -> str:
    lines = [
        "# D5 — Mode-A Overlap",
        "",
        "Window: t ∈ [8, 48] hr, centerline mid-depth (below gate window).",
        f"Match criteria: mean ΔT ≥ +{MODE_A_MEAN_MIN_F:.1f}°F AND"
        f" peak ΔT time ∈ [{MODE_A_PEAK_LO_HR:.0f}, {MODE_A_PEAK_HI_HR:.0f}] hr.",
        "Source: PR 16 Phase 2.6 Mode-A signature (§7.6.4).",
        "",
        "| Mix | window_min_hr | window_max_hr | mean ΔT (°F) | max ΔT (°F) | t_at_max (hr) | Mode-A match |",
        "|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        match_str = "**Y**" if r["mode_a_match"] else "N"
        lines.append(
            f"| {r['mix']} | {r['window_min_hr']:.1f} | {r['window_max_hr']:.1f}"
            f" | {r['mean_dt_F']:+.4f} | {r['max_dt_F']:+.4f}"
            f" | {r['t_at_max_hr']:.1f} | {match_str} |"
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    rows = []
    for mix in REFERENCE_MIXES:
        print(f"  {mix}: loading...", end=" ", flush=True)
        cache = load_mix_cache(mix)
        r = compute_for_mix(cache)
        rows.append(r)
        print(f"mode_a window: {r['window_min_hr']:.1f} hr → {r['window_max_hr']:.1f} hr"
              f"  mean ΔT={r['mean_dt_F']:+.4f}°F  max={r['max_dt_F']:+.4f}°F"
              f"  t_at_max={r['t_at_max_hr']:.1f} hr  match={r['mode_a_match']}")

    md = format_table(rows)
    out = os.path.join(_HERE, "d5_mode_a_overlap.md")
    with open(out, "w") as f:
        f.write(md)
    print(f"\nWritten: {out}")
