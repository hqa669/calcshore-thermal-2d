"""
D1 — Constant offset.
Is the centerline residual a near-constant ΔT throughout [48, 168] hr?
Computed at centerline mid-depth; window matches the S0 RMS gate.
"""
import os
import sys
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _HERE)
sys.path.insert(0, _REPO)

from _harness import load_mix_cache, REFERENCE_MIXES


def compute_for_mix(cache) -> dict:
    diff = (cache.eng_center_F_at_cw_t[cache.gate_mask]
            - cache.cw_center_F[cache.gate_mask])
    mean_dt = float(np.mean(diff))
    std_dt  = float(np.std(diff, ddof=1))
    rms_check = float(np.sqrt(mean_dt**2 + std_dt**2))
    return {
        "mix": cache.mix_name,
        "mean_dt_F": mean_dt,
        "std_dt_F": std_dt,
        "rms_check_F": rms_check,
    }


def format_table(rows: list) -> str:
    lines = [
        "# D1 — Constant Offset",
        "",
        "Window: t ∈ [48, 168] hr, centerline mid-depth.",
        "RMS_check = sqrt(mean²+std²); compare to run_all.py Centerline RMS.",
        "Near-constant offset: mean ≈ ±RMS_check and std is small.",
        "",
        "| Mix | mean ΔT (°F) | std ΔT (°F) | RMS_check (°F) |",
        "|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['mix']} | {r['mean_dt_F']:+.4f} | {r['std_dt_F']:.4f}"
            f" | {r['rms_check_F']:.4f} |"
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    rows = []
    for mix in REFERENCE_MIXES:
        print(f"  {mix}: loading...", end=" ", flush=True)
        cache = load_mix_cache(mix)
        r = compute_for_mix(cache)
        rows.append(r)
        print(f"mean ΔT={r['mean_dt_F']:+.4f}°F  std={r['std_dt_F']:.4f}°F"
              f"  RMS_check={r['rms_check_F']:.4f}°F")

    md = format_table(rows)
    out = os.path.join(_HERE, "d1_constant_offset.md")
    with open(out, "w") as f:
        f.write(md)
    print(f"\nWritten: {out}")
