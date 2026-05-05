#!/usr/bin/env python3
"""§3.6 — Build 4 result tables from sweep_results.npz.

Tables (8 rows = T_pc, 4 cols = α_u):
  table_c_optimal.csv          — c_optimal[T_pc, α_u]
  table_max_dT_at_optimal.csv  — max|T_eng − T_CW| at c_optimal (°F)
  table_max_dT_at_c1.csv       — max|T_eng − T_CW| at c=1.0 (baseline gap, °F)
  table_endpoint_dT_at_optimal.csv — (T_eng − T_CW) at t=168 hr (°F)
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"

T_PCS   = [40, 50, 60, 70, 80, 90, 100, 110]
ALPHAS  = [0.20, 0.40, 0.60, 0.80]
ALPHA_LABELS = ["α_u=0.20", "α_u=0.40", "α_u=0.60", "α_u=0.80"]


def build_table(keys: list[str], alpha_u_arr, T_pc_F_arr, values_arr):
    """Return (rows_2d, T_pcs, alphas) where rows_2d[i] = [val(Tpc_i, au_j) for au_j in ALPHAS]."""
    n = len(alpha_u_arr)
    table = []
    for T_pc in T_PCS:
        row = {"T_pc_F": T_pc}
        for au, label in zip(ALPHAS, ALPHA_LABELS):
            idx = next((i for i in range(n)
                        if abs(alpha_u_arr[i] - au) < 0.01 and abs(T_pc_F_arr[i] - T_pc) < 0.5),
                       None)
            row[label] = float(values_arr[idx]) if idx is not None else float("nan")
        table.append(row)
    return table


def write_csv(path, table, fieldnames):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(table)
    print(f"Wrote {path}")


def main():
    npz = np.load(DATA / "sweep_results.npz", allow_pickle=True)
    alpha_u_arr  = npz["alpha_u"]
    T_pc_F_arr   = npz["T_pc_F"]
    c_optimal    = npz["c_optimal"]
    max_dT_opt   = npz["max_dT_at_optimal_F"]
    max_dT_c1    = npz["max_dT_at_c1_F"]
    endpoint_dT  = npz["endpoint_dT_at_optimal_F"]

    fields = ["T_pc_F"] + ALPHA_LABELS

    print(f"\n§3.6 Building tables from {len(alpha_u_arr)} sweep results")
    print(f"  Hu_eff = {float(npz['Hu_eff_J_kg']):.1f} J/kg")

    for fname, arr in [
        ("table_c_optimal.csv",             c_optimal),
        ("table_max_dT_at_optimal.csv",     max_dT_opt),
        ("table_max_dT_at_c1.csv",          max_dT_c1),
        ("table_endpoint_dT_at_optimal.csv", endpoint_dT),
    ]:
        table = build_table(fields, alpha_u_arr, T_pc_F_arr, arr)
        write_csv(DATA / fname, table, fields)

    # Print c_optimal table for quick review
    print("\nc_optimal[T_pc, α_u]:")
    t = build_table(fields, alpha_u_arr, T_pc_F_arr, c_optimal)
    print(f"  {'T_pc':>5}" + "".join(f"  {'α='+str(au):>8}" for au in ALPHAS))
    for row in t:
        print(f"  {row['T_pc_F']:>5.0f}" +
              "".join(f"  {row[l]:>8.4f}" for l in ALPHA_LABELS))

    print("\nDone.")


if __name__ == "__main__":
    main()
