#!/usr/bin/env python3
"""§4.11.7 Task 2 — Dense local refinement at converged points.

For each (T_pc, α_u) in the converged set (T_pc ∈ {50, 60, 70, 80, 90, 100, 110}, 4 α_u = 28 pts):
  - Look up c_v2 from the merged v2 table.
  - Scan c ∈ [c_v2 - 0.10, c_v2 + 0.10] at step 0.01 (21 points).
  - Save c_optimal_t2 = argmin L2; record L2 at v2 value, L2_min, L2_max for diagnostic.
  - Edge check: if c_optimal_t2 hits scan window edge, expand to ±0.20 and re-run.

Reads:
  data/cw_trajectories.npz                — CW reference
  data/sweep_results.npz                  — Hu_eff and §4.11 c_optimal (T_pc=60,70,80)
  data/sweep_results_extended.npz         — §4.11.6 c_optimal (T_pc=50,90,100,110)

Writes:
  data/sweep_results_t2_refined.npz / .csv
  data/resolution_floor_diagnostic.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from sweep_alpha_factor import run_engine_Tcore, compute_loss

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"

T_PCS_INTERIOR = [60.0, 70.0, 80.0]              # from §4.11 sweep_results.npz
T_PCS_EXTENDED = [50.0, 90.0, 100.0, 110.0]      # from §4.11.6 sweep_results_extended.npz
T_PCS_ALL      = sorted(T_PCS_INTERIOR + T_PCS_EXTENDED)
ALPHAS         = [0.20, 0.40, 0.60, 0.80]
SCAN_HALFWIDTH = 0.10
SCAN_STEP      = 0.01


def lookup(au_arr, T_arr, vals, au, T_pc):
    n = len(au_arr)
    idx = next(i for i in range(n)
               if abs(au_arr[i] - au) < 0.01 and abs(T_arr[i] - T_pc) < 0.5)
    return float(vals[idx])


def get_c_v2(au, T_pc, sweep0, sweep1):
    if T_pc in T_PCS_INTERIOR:
        return lookup(sweep0["alpha_u"], sweep0["T_pc_F"], sweep0["c_optimal"], au, T_pc)
    return lookup(sweep1["alpha_u"], sweep1["T_pc_F"], sweep1["c_optimal_extended"], au, T_pc)


def find_folder(folder_names, alpha_u_arr, T_pc_F_arr, au, T_pc):
    n = len(folder_names)
    idx = next(i for i in range(n)
               if abs(alpha_u_arr[i] - au) < 0.01 and abs(T_pc_F_arr[i] - T_pc) < 0.5)
    return CW_DATA / str(folder_names[idx]), idx


def scan_window(folder, T_pc, au, c_lo, c_hi, step, t_cw_hr, T_cw_F, Hu_eff):
    """Return list of (c, L2) pairs for the scan window, with c at SCAN_STEP precision."""
    n = int(round((c_hi - c_lo) / step)) + 1
    cs = [round(c_lo + i * step, 4) for i in range(n)]
    losses = []
    for c in cs:
        t_e, T_e = run_engine_Tcore(folder, T_pc, au * c, Hu_eff)
        L2 = compute_loss(t_e, T_e, t_cw_hr, T_cw_F)
        losses.append((c, L2, t_e, T_e))
    return losses


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    traj = np.load(DATA / "cw_trajectories.npz", allow_pickle=True)
    t_cw_hr      = traj["t_hrs"]
    T_cw_F_all   = traj["T_core_CW_F"]
    alpha_u_arr  = traj["alpha_u"]
    T_pc_F_arr   = traj["T_pc_F"]
    folder_names = traj["folder_names"]

    sweep0 = np.load(DATA / "sweep_results.npz",          allow_pickle=True)
    sweep1 = np.load(DATA / "sweep_results_extended.npz", allow_pickle=True)
    HU_EFF = float(sweep0["Hu_eff_J_kg"])

    print(f"\n§4.11.7 Task 2 — Dense local refinement (28 points)")
    print(f"  Hu_eff: {HU_EFF:.2f} J/kg")
    print(f"  Scan window: ±{SCAN_HALFWIDTH:.2f} around c_v2, step {SCAN_STEP:.2f} (21 pts/scan)")
    print(f"  T_pc set: {T_PCS_ALL}")
    print(f"  α_u set: {ALPHAS}")
    print("=" * 100)

    rows = []
    diag_rows = []

    for au in ALPHAS:
        for T_pc in T_PCS_ALL:
            folder, src_idx = find_folder(folder_names, alpha_u_arr, T_pc_F_arr, au, T_pc)
            T_cw_F = T_cw_F_all[src_idx]
            c_v2 = get_c_v2(au, T_pc, sweep0, sweep1)

            c_lo = round(c_v2 - SCAN_HALFWIDTH, 4)
            c_hi = round(c_v2 + SCAN_HALFWIDTH, 4)

            print(f"  α_u={au:.2f}  T_pc={T_pc:.0f}°F  c_v2={c_v2:.4f}  "
                  f"scan[{c_lo:.2f},{c_hi:.2f}] ...", end="", flush=True)

            losses = scan_window(folder, T_pc, au, c_lo, c_hi, SCAN_STEP,
                                 t_cw_hr, T_cw_F, HU_EFF)

            cs, L2s = [r[0] for r in losses], [r[1] for r in losses]
            i_min = int(np.argmin(L2s))
            c_opt_t2 = cs[i_min]
            L2_min   = L2s[i_min]
            L2_max   = max(L2s)
            # closest c to c_v2 in scan
            i_v2 = int(np.argmin([abs(c - c_v2) for c in cs]))
            L2_at_v2 = L2s[i_v2]

            edge_window = (i_min == 0 or i_min == len(cs) - 1)

            # If at edge of window, expand to ±0.20 once.
            if edge_window:
                print(f"  EDGE — expanding to ±0.20 ...", end="", flush=True)
                c_lo_e = round(c_v2 - 2 * SCAN_HALFWIDTH, 4)
                c_hi_e = round(c_v2 + 2 * SCAN_HALFWIDTH, 4)
                losses = scan_window(folder, T_pc, au, c_lo_e, c_hi_e, SCAN_STEP,
                                     t_cw_hr, T_cw_F, HU_EFF)
                cs, L2s = [r[0] for r in losses], [r[1] for r in losses]
                i_min = int(np.argmin(L2s))
                c_opt_t2 = cs[i_min]
                L2_min   = L2s[i_min]
                L2_max   = max(L2s)
                i_v2 = int(np.argmin([abs(c - c_v2) for c in cs]))
                L2_at_v2 = L2s[i_v2]
                edge_window = (i_min == 0 or i_min == len(cs) - 1)

            # Metrics at c_opt_t2 (use cached trajectory at i_min)
            t_o, T_o = losses[i_min][2], losses[i_min][3]
            T_o_i = np.interp(t_cw_hr, t_o, T_o)
            max_dT = float(np.max(np.abs(T_o_i - T_cw_F)))
            end_dT = float(T_o_i[-1] - T_cw_F[-1])

            delta_v2 = c_opt_t2 - c_v2
            print(f"  c_opt_t2={c_opt_t2:.3f}  Δ={delta_v2:+.3f}  L2={L2_min:.5f}")

            rows.append({
                "T_pc_F":             T_pc,
                "alpha_u":            au,
                "c_v2":               c_v2,
                "c_optimal_t2":       c_opt_t2,
                "delta_from_v2":      delta_v2,
                "L2_residual_t2":     L2_min,
                "max_dT_at_optimal_t2": max_dT,
                "endpoint_dT_at_optimal_t2": end_dT,
                "edge_window":        edge_window,
                "n_scan":             len(cs),
            })
            diag_rows.append({
                "T_pc_F":         T_pc,
                "alpha_u":        au,
                "c_v2":           c_v2,
                "c_optimal_t2":   c_opt_t2,
                "L2_min":         L2_min,
                "L2_at_v2":       L2_at_v2,
                "L2_max":         L2_max,
                "L2_relative_improvement_pct":
                    (100.0 * (L2_at_v2 - L2_min) / L2_at_v2) if L2_at_v2 > 0 else 0.0,
            })

    # Save
    out_csv = DATA / "sweep_results_t2_refined.csv"
    fields = ["T_pc_F", "alpha_u", "c_v2", "c_optimal_t2", "delta_from_v2",
              "L2_residual_t2", "max_dT_at_optimal_t2",
              "endpoint_dT_at_optimal_t2", "edge_window", "n_scan"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_csv}")

    diag_csv = DATA / "resolution_floor_diagnostic.csv"
    with open(diag_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(diag_rows[0].keys()))
        w.writeheader()
        w.writerows(diag_rows)
    print(f"Wrote {diag_csv}")

    np.savez(
        DATA / "sweep_results_t2_refined.npz",
        T_pc_F=np.array([r["T_pc_F"] for r in rows]),
        alpha_u=np.array([r["alpha_u"] for r in rows]),
        c_v2=np.array([r["c_v2"] for r in rows]),
        c_optimal_t2=np.array([r["c_optimal_t2"] for r in rows]),
        delta_from_v2=np.array([r["delta_from_v2"] for r in rows]),
        L2_residual_t2=np.array([r["L2_residual_t2"] for r in rows]),
        max_dT_at_optimal_t2=np.array([r["max_dT_at_optimal_t2"] for r in rows]),
        endpoint_dT_at_optimal_t2=np.array([r["endpoint_dT_at_optimal_t2"] for r in rows]),
        edge_window=np.array([r["edge_window"] for r in rows]),
    )

    # Summary
    print("\nTask 2 refined c_optimal (Δ from v2):")
    print(f"{'α_u':>6}  " + "  ".join(f"{int(T):>5}°F" for T in T_PCS_ALL))
    for au in ALPHAS:
        cells = []
        for T in T_PCS_ALL:
            r = next(r for r in rows if abs(r["alpha_u"] - au) < 0.01
                     and abs(r["T_pc_F"] - T) < 0.5)
            cells.append(f"{r['c_optimal_t2']:>5.3f}({r['delta_from_v2']:+.3f})")
        print(f"  {au:.2f}  " + "  ".join(cells))

    # Resolution floor diagnostic summary
    pcts = [d["L2_relative_improvement_pct"] for d in diag_rows]
    print(f"\nResolution floor diagnostic:")
    print(f"  L2 improvement (median): {np.median(pcts):.2f}%")
    print(f"  L2 improvement (max):    {np.max(pcts):.2f}%")
    print(f"  Points with >5% improvement: {sum(1 for p in pcts if p > 5.0)}/{len(pcts)}")

    print("\nDone.")


if __name__ == "__main__":
    main()
