#!/usr/bin/env python3
"""§4.11.7 Task 1 — T_pc=40°F range extension to c ∈ [0.20, 0.60].

For each α_u ∈ {0.20, 0.40, 0.60, 0.80}:
  1. 9-point coarse scan: c ∈ {0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60}
  2. Golden-section refinement within bracket of coarse min, tol=0.005
  3. If coarse_best is at edge ({0.20} or {0.60}): flag and stop for that α_u

Reuses run_engine_Tcore, compute_loss, golden_section from sweep_alpha_factor.py.
Hu_eff fixed at 12936.96 J/kg (from §4.11 calibration).

Writes:
  data/sweep_results_t1_t40_extended.csv
  data/sweep_results_t1_t40_extended.npz
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from sweep_alpha_factor import run_engine_Tcore, compute_loss, golden_section

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"

T_PC = 40.0
ALPHAS = [0.20, 0.40, 0.60, 0.80]
C_COARSE_T1 = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60]
GSS_TOL = 0.005
EDGE_LO = C_COARSE_T1[0]
EDGE_HI = C_COARSE_T1[-1]


def find_folder(folder_names, alpha_u_arr, T_pc_F_arr, au):
    n = len(folder_names)
    idx = next(i for i in range(n)
               if abs(alpha_u_arr[i] - au) < 0.01 and abs(T_pc_F_arr[i] - T_PC) < 0.5)
    return CW_DATA / str(folder_names[idx]), idx


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    traj = np.load(DATA / "cw_trajectories.npz", allow_pickle=True)
    t_cw_hr      = traj["t_hrs"]
    T_cw_F_all   = traj["T_core_CW_F"]
    alpha_u_arr  = traj["alpha_u"]
    T_pc_F_arr   = traj["T_pc_F"]
    folder_names = traj["folder_names"]

    sweep0 = np.load(DATA / "sweep_results.npz", allow_pickle=True)
    HU_EFF = float(sweep0["Hu_eff_J_kg"])

    print(f"\n§4.11.7 Task 1 — T_pc={T_PC:.0f}°F range extension")
    print(f"  Hu_eff: {HU_EFF:.2f} J/kg")
    print(f"  Coarse grid (n={len(C_COARSE_T1)}): {C_COARSE_T1}")
    print(f"  GSS tol: {GSS_TOL}")
    print("=" * 100)

    rows = []
    for au in ALPHAS:
        folder, src_idx = find_folder(folder_names, alpha_u_arr, T_pc_F_arr, au)
        T_cw_F = T_cw_F_all[src_idx]

        print(f"\n  α_u={au:.2f}  T_pc={T_PC:.0f}°F  coarse[9] ...", end="", flush=True)
        coarse = {}
        for c in C_COARSE_T1:
            t_e, T_e = run_engine_Tcore(folder, T_PC, au * c, HU_EFF)
            coarse[c] = compute_loss(t_e, T_e, t_cw_hr, T_cw_F)
        c_best = min(coarse, key=lambda k: coarse[k])
        print(f"  coarse_best={c_best:.2f} (loss={coarse[c_best]:.5f})")

        edge = False
        if c_best == EDGE_LO or c_best == EDGE_HI:
            edge = True
            c_opt = c_best
            n_gss = 0
            print(f"    EDGE-CLAMPED at {c_best:.2f}  (no GSS, flagged)")
        else:
            ix = C_COARSE_T1.index(c_best)
            a, b = C_COARSE_T1[ix - 1], C_COARSE_T1[ix + 1]

            def loss_fn(c):
                t_e, T_e = run_engine_Tcore(folder, T_PC, au * c, HU_EFF)
                return compute_loss(t_e, T_e, t_cw_hr, T_cw_F)

            c_opt, n_gss = golden_section(loss_fn, a, b, tol=GSS_TOL)
            print(f"    refine[{a:.2f},{b:.2f}] → c_opt={c_opt:.4f}  ({n_gss} GSS evals)")

        # Metrics at c_opt
        t_o, T_o = run_engine_Tcore(folder, T_PC, au * c_opt, HU_EFF)
        T_o_i = np.interp(t_cw_hr, t_o, T_o)
        max_dT = float(np.max(np.abs(T_o_i - T_cw_F)))
        end_dT = float(T_o_i[-1] - T_cw_F[-1])
        L2 = compute_loss(t_o, T_o, t_cw_hr, T_cw_F)

        rows.append({
            "T_pc_F":               T_PC,
            "alpha_u":              au,
            "c_optimal_t1":         c_opt,
            "edge_flag":            edge,
            "L2_residual_t1":       L2,
            "max_dT_at_optimal_t1": max_dT,
            "endpoint_dT_at_optimal_t1": end_dT,
            "n_gss":                n_gss,
        })

    # Save
    out_csv = DATA / "sweep_results_t1_t40_extended.csv"
    fields = ["T_pc_F", "alpha_u", "c_optimal_t1", "edge_flag",
              "L2_residual_t1", "max_dT_at_optimal_t1",
              "endpoint_dT_at_optimal_t1", "n_gss"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_csv}")

    np.savez(
        DATA / "sweep_results_t1_t40_extended.npz",
        T_pc_F=np.array([r["T_pc_F"] for r in rows]),
        alpha_u=np.array([r["alpha_u"] for r in rows]),
        c_optimal_t1=np.array([r["c_optimal_t1"] for r in rows]),
        edge_flag=np.array([r["edge_flag"] for r in rows]),
        L2_residual_t1=np.array([r["L2_residual_t1"] for r in rows]),
        max_dT_at_optimal_t1=np.array([r["max_dT_at_optimal_t1"] for r in rows]),
        endpoint_dT_at_optimal_t1=np.array([r["endpoint_dT_at_optimal_t1"] for r in rows]),
    )

    print(f"\nTask 1 summary @ T_pc={T_PC:.0f}°F:")
    print(f"  {'α_u':>6}  {'c_opt':>8}  {'edge':>5}  {'L2':>10}  {'max_dT':>8}")
    for r in rows:
        flag = "YES" if r["edge_flag"] else ""
        print(f"  {r['alpha_u']:>6.2f}  {r['c_optimal_t1']:>8.4f}  {flag:>5}  "
              f"{r['L2_residual_t1']:>10.5f}  {r['max_dT_at_optimal_t1']:>8.4f}")
    print("\nDone.")


if __name__ == "__main__":
    main()
