"""
D2 — Phase lag.
Is the residual a time-shift of the same diurnal amplitude?
Cross-correlates detrended engine vs CW centerline over [48, 168] hr.
Sign convention: positive lag_hr = engine leads CW.
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
    eng_g = cache.eng_center_F_at_cw_t[cache.gate_mask]
    cw_g  = cache.cw_center_F[cache.gate_mask]

    eng_d = eng_g - eng_g.mean()
    cw_d  = cw_g  - cw_g.mean()

    # xcorr[k+N-1] = sum_t eng_d[t] * cw_d[t-k]; peak k = lag where eng leads cw
    n = len(eng_d)
    corr = np.correlate(eng_d, cw_d, mode="full")
    lags = np.arange(-(n - 1), n)

    peak_idx = int(np.argmax(corr))
    lag_samples = int(lags[peak_idx])
    stride_hr = float(np.median(np.diff(cache.cw_time_hrs)))
    lag_hr = lag_samples * stride_hr

    corr_peak = float(corr[peak_idx])
    corr_zero = float(corr[n - 1])  # lag=0 is at index n-1
    ratio = corr_zero / corr_peak if corr_peak != 0.0 else float("nan")

    return {
        "mix": cache.mix_name,
        "lag_hr": lag_hr,
        "lag_samples": lag_samples,
        "stride_hr": stride_hr,
        "corr_peak": corr_peak,
        "corr_zero": corr_zero,
        "corr_zero_over_peak": ratio,
    }


def format_table(rows: list) -> str:
    lines = [
        "# D2 — Phase Lag",
        "",
        "Window: t ∈ [48, 168] hr, centerline mid-depth, detrended.",
        "Positive lag_hr = engine leads CW.",
        "corr_zero/peak near 1.0 → lag is negligible; near 0 → significant lag.",
        "",
        "| Mix | lag_hr | lag_samples | stride_hr | corr_zero/peak |",
        "|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['mix']} | {r['lag_hr']:+.2f} | {r['lag_samples']:+d}"
            f" | {r['stride_hr']:.4f} | {r['corr_zero_over_peak']:.4f} |"
        )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    rows = []
    for mix in REFERENCE_MIXES:
        print(f"  {mix}: loading...", end=" ", flush=True)
        cache = load_mix_cache(mix)
        r = compute_for_mix(cache)
        rows.append(r)
        print(f"lag={r['lag_hr']:+.2f} hr  corr_zero/peak={r['corr_zero_over_peak']:.4f}")

    md = format_table(rows)
    out = os.path.join(_HERE, "d2_phase_lag.md")
    with open(out, "w") as f:
        f.write(md)
    print(f"\nWritten: {out}")
