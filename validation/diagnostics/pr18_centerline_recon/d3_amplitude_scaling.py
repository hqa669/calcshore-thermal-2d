"""
D3 — Amplitude scaling.
Is the engine's diurnal swing larger or smaller than CW's?
Evaluated over the last full diurnal cycle in the gate window: t ∈ [144, 168] hr.
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
    dm = cache.diurnal_mask
    eng_sw = float(cache.eng_center_F_at_cw_t[dm].max()
                   - cache.eng_center_F_at_cw_t[dm].min())
    cw_sw  = float(cache.cw_center_F[dm].max() - cache.cw_center_F[dm].min())
    ratio  = eng_sw / cw_sw if cw_sw > 0.0 else float("nan")
    n_pts  = int(dm.sum())
    return {
        "mix": cache.mix_name,
        "eng_swing_F": eng_sw,
        "cw_swing_F": cw_sw,
        "ratio": ratio,
        "n_pts": n_pts,
    }


def format_table(rows: list) -> str:
    lines = [
        "# D3 — Amplitude Scaling",
        "",
        "Window: t ∈ [144, 168] hr (last full diurnal cycle), centerline mid-depth.",
        "ratio = eng_swing / cw_swing; 1.0 = matched amplitude.",
        "",
        "| Mix | eng_swing (°F) | cw_swing (°F) | ratio | n_pts |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['mix']} | {r['eng_swing_F']:.4f} | {r['cw_swing_F']:.4f}"
            f" | {r['ratio']:.4f} | {r['n_pts']} |"
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    rows = []
    for mix in REFERENCE_MIXES:
        print(f"  {mix}: loading...", end=" ", flush=True)
        cache = load_mix_cache(mix)
        r = compute_for_mix(cache)
        rows.append(r)
        print(f"eng_swing={r['eng_swing_F']:.4f}°F  cw_swing={r['cw_swing_F']:.4f}°F"
              f"  ratio={r['ratio']:.4f}")

    md = format_table(rows)
    out = os.path.join(_HERE, "d3_amplitude_scaling.md")
    with open(out, "w") as f:
        f.write(md)
    print(f"\nWritten: {out}")
