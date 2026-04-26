"""
D4 — Depth profile.
Is the residual uniform with depth or concentrated at specific depths?
Per-depth RMS over [48, 168] hr at the centerline column.
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
    gate = cache.gate_mask
    n_eng = cache.eng_centerline_col_F_at_cw_t.shape[1]
    n_cw  = cache.cw_centerline_col_F.shape[1]
    n_depths = min(n_eng, n_cw)

    depth_rms = []
    for d in range(n_depths):
        eng_col = cache.eng_centerline_col_F_at_cw_t[gate, d]
        cw_col  = cache.cw_centerline_col_F[gate, d]
        rms_d   = float(np.sqrt(np.mean((eng_col - cw_col) ** 2)))
        depth_rms.append(rms_d)

    mid_idx  = n_depths // 2
    mid_rms  = depth_rms[mid_idx]
    spread   = float(max(depth_rms) - min(depth_rms))

    return {
        "mix": cache.mix_name,
        "depth_rms": depth_rms,
        "n_depths": n_depths,
        "mid_rms_F": mid_rms,
        "spread_F": spread,
    }


def format_table(rows: list) -> str:
    lines = [
        "# D4 — Depth Profile",
        "",
        "Window: t ∈ [48, 168] hr, all depths at the centerline column.",
        "spread = max(depth_rms) - min(depth_rms); near 0 → uniform, large → localized.",
        "",
    ]
    for r in rows:
        lines.append(f"## {r['mix']}  (spread={r['spread_F']:.4f}°F)")
        lines.append("")
        lines.append("| depth_idx | label | RMS (°F) |")
        lines.append("|---|---|---|")
        n = r["n_depths"]
        for d, rms_d in enumerate(r["depth_rms"]):
            if d == 0:
                label = "top"
            elif d == n - 1:
                label = "bot"
            elif d == n // 2:
                label = "mid"
            else:
                label = f"{d}"
            lines.append(f"| {d} | {label} | {rms_d:.4f} |")
        lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    rows = []
    for mix in REFERENCE_MIXES:
        print(f"  {mix}: loading...", end=" ", flush=True)
        cache = load_mix_cache(mix)
        r = compute_for_mix(cache)
        rows.append(r)
        print(f"n_depths={r['n_depths']}  mid_rms={r['mid_rms_F']:.4f}°F"
              f"  spread={r['spread_F']:.4f}°F")

    md = format_table(rows)
    out = os.path.join(_HERE, "d4_depth_profile.md")
    with open(out, "w") as f:
        f.write(md)
    print(f"\nWritten: {out}")
