#!/usr/bin/env python3
"""§4.11.6 §4 — Build merged 8×4 tables (v2) from §4.11 + extended-range sweep.

Combines:
  - T_pc ∈ {60, 70, 80}: from data/sweep_results.npz (unchanged from §4.11)
  - T_pc ∈ {40, 50, 90, 100, 110}: from data/sweep_results_extended.npz

Also writes:
  - data/extended_sweep_comparison.csv — predicted vs actual c per (T_pc, α_u)
  - data/table_c_optimal_v2.csv
  - data/table_max_dT_at_optimal_v2.csv
  - data/table_max_dT_at_c1_v2.csv
  - data/table_endpoint_dT_at_optimal_v2.csv
"""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"

T_PCS_INTERIOR = [60, 70, 80]
T_PCS_EDGE     = [40, 50, 90, 100, 110]
T_PCS_ALL      = sorted(T_PCS_INTERIOR + T_PCS_EDGE)
ALPHAS         = [0.20, 0.40, 0.60, 0.80]
ALPHA_LABELS   = [f"α_u={au:.2f}" for au in ALPHAS]


def lookup(au_arr, T_arr, vals, au, T_pc):
    n = len(au_arr)
    idx = next((i for i in range(n)
                if abs(au_arr[i] - au) < 0.01 and abs(T_arr[i] - T_pc) < 0.5), None)
    return float(vals[idx]) if idx is not None else float("nan")


def merged_value(key_v1, key_ext, au, T_pc, sweep0, sweep1):
    """Pick from §4.11 sweep_results for interior T_pc, from extended for edge T_pc."""
    if T_pc in T_PCS_INTERIOR:
        return lookup(sweep0["alpha_u"], sweep0["T_pc_F"], sweep0[key_v1], au, T_pc)
    return lookup(sweep1["alpha_u"], sweep1["T_pc_F"], sweep1[key_ext], au, T_pc)


def build_table(label_to_key_v1, label_to_key_ext, sweep0, sweep1):
    rows = []
    for T_pc in T_PCS_ALL:
        row = {"T_pc_F": T_pc}
        for au, label in zip(ALPHAS, ALPHA_LABELS):
            row[label] = merged_value(label_to_key_v1, label_to_key_ext,
                                      au, T_pc, sweep0, sweep1)
        rows.append(row)
    return rows


def write_csv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {path}")


def main():
    sweep0 = np.load(DATA / "sweep_results.npz",          allow_pickle=True)
    sweep1 = np.load(DATA / "sweep_results_extended.npz", allow_pickle=True)

    fields = ["T_pc_F"] + ALPHA_LABELS

    # 4 merged tables: c_optimal, max_dT_at_optimal, max_dT_at_c1, endpoint_dT_at_optimal
    table_specs = [
        ("table_c_optimal_v2.csv",
         "c_optimal", "c_optimal_extended"),
        ("table_max_dT_at_optimal_v2.csv",
         "max_dT_at_optimal_F", "max_dT_at_optimal_F"),
        ("table_max_dT_at_c1_v2.csv",
         "max_dT_at_c1_F", "max_dT_at_c1_F"),
        ("table_endpoint_dT_at_optimal_v2.csv",
         "endpoint_dT_at_optimal_F", "endpoint_dT_at_optimal_F"),
    ]
    for fname, k_v1, k_ext in table_specs:
        rows = build_table(k_v1, k_ext, sweep0, sweep1)
        write_csv(DATA / fname, rows, fields)

    # Comparison: predicted (per-α_u linear extrapolation from §4.11.3 70°F & 80°F)
    # to actual c_optimal_extended at the 5 edge T_pc.
    print("\n--- Linear extrapolation comparison (extended T_pc) ---")
    cmp_rows = []
    print(f"{'T_pc':>5} {'α_u':>5} {'c_pred':>7} {'c_actual':>9} {'Δ':>7} {'%Δ':>7} {'edge':>5}")
    for au in ALPHAS:
        c70 = lookup(sweep0["alpha_u"], sweep0["T_pc_F"], sweep0["c_optimal"], au, 70.0)
        c80 = lookup(sweep0["alpha_u"], sweep0["T_pc_F"], sweep0["c_optimal"], au, 80.0)
        slope = c80 - c70  # per 10 °F
        for T_pc in T_PCS_EDGE:
            c_pred   = c70 + slope * (T_pc - 70.0) / 10.0
            c_actual = lookup(sweep1["alpha_u"], sweep1["T_pc_F"],
                              sweep1["c_optimal_extended"], au, T_pc)
            edge     = bool(lookup(sweep1["alpha_u"], sweep1["T_pc_F"],
                                   sweep1["edge_flag"], au, T_pc))
            delta    = c_actual - c_pred
            pct      = 100.0 * delta / c_pred if c_pred != 0 else float("nan")
            print(f"{T_pc:>5} {au:>5.2f} {c_pred:>7.3f} {c_actual:>9.3f} "
                  f"{delta:>+7.3f} {pct:>+6.1f}% {('YES' if edge else '   '):>5}")
            cmp_rows.append({
                "T_pc_F":       T_pc,
                "alpha_u":      au,
                "c_predicted":  c_pred,
                "c_actual":     c_actual,
                "delta":        delta,
                "delta_pct":    pct,
                "edge_flag":    edge,
            })

    write_csv(DATA / "extended_sweep_comparison.csv", cmp_rows,
              ["T_pc_F", "alpha_u", "c_predicted", "c_actual",
               "delta", "delta_pct", "edge_flag"])

    # Print merged c_optimal table
    print("\n--- Merged c_optimal[T_pc, α_u] (v2) ---")
    print(f"  {'T_pc':>5}" + "".join(f"  {'α='+str(au):>8}" for au in ALPHAS))
    rows = build_table("c_optimal", "c_optimal_extended", sweep0, sweep1)
    for row in rows:
        cells = [f"  {row[l]:>8.4f}" for l in ALPHA_LABELS]
        src = " (interior)" if int(row["T_pc_F"]) in T_PCS_INTERIOR else " (extended)"
        print(f"  {row['T_pc_F']:>5.0f}" + "".join(cells) + src)

    # Trend classification heuristic
    print("\n--- Trend classification heuristic ---")
    pcts = [r["delta_pct"] for r in cmp_rows]
    max_abs_pct = max(abs(p) for p in pcts)
    n_within_5 = sum(1 for p in pcts if abs(p) <= 5.0)
    n_total    = len(pcts)
    print(f"  max |%Δ| across 20 points: {max_abs_pct:.1f}%")
    print(f"  points within ±5%: {n_within_5}/{n_total}")

    # Sublinear sign check: at T_pc < 70, c_pred < 1 so actual−predicted > 0 means actual closer to 1.
    # At T_pc > 70, c_pred > 1 so actual−predicted < 0 means actual closer to 1.
    # Define toward_one = +1 if Δ moves toward 1, −1 if away.
    sublinear_signal = 0
    superlinear_signal = 0
    for r in cmp_rows:
        if r["T_pc_F"] < 70:
            if r["delta"] > 0: sublinear_signal += 1
            elif r["delta"] < 0: superlinear_signal += 1
        else:
            if r["delta"] < 0: sublinear_signal += 1
            elif r["delta"] > 0: superlinear_signal += 1
    print(f"  sublinear signals (Δ moves toward 1):  {sublinear_signal}/{n_total}")
    print(f"  superlinear signals (Δ moves away):    {superlinear_signal}/{n_total}")

    if n_within_5 == n_total:
        verdict = "LINEAR"
    elif sublinear_signal >= 0.7 * n_total:
        verdict = "SUBLINEAR (asymptoting toward 1.0 at extremes)"
    elif superlinear_signal >= 0.7 * n_total:
        verdict = "SUPERLINEAR (diverging from linear at extremes)"
    else:
        verdict = "MIXED (no dominant pattern; possibly non-monotonic — inspect plots)"
    print(f"  → Verdict: {verdict}")

    print("\nDone.")


if __name__ == "__main__":
    main()
