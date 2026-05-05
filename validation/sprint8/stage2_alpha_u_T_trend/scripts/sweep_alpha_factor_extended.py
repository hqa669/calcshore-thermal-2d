#!/usr/bin/env python3
"""§4.11.6 §1 — Extended-range α_u factor re-sweep at edge T_pc.

Re-runs the 5 edge T_pc values × 4 α_u values = 20 (T_pc, α_u) points with an
extended search range c ∈ [0.40, 1.80] (vs original [0.80, 1.20] in §4.11).
Reuses Hu_eff from data/sweep_results.npz to keep the 70°F anchor identical.

Per-point procedure:
  1. 15-point coarse scan over c ∈ {0.40, 0.50, ..., 1.80}
  2. If coarse minimum is at the new edge ({0.40} or {1.80}): flag and stop for that point
  3. Otherwise: golden-section refinement within the bracketing pair (tol=0.005)

Reads:
  data/cw_trajectories.npz       — CW reference trajectories
  data/sweep_results.npz         — Hu_eff_J_kg (from §4.11 calibration)
  cw_data/thermal_temp_alpha*/   — input.dat per dataset

Writes:
  data/sweep_results_extended.npz
  data/sweep_results_extended.csv
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from sweep_alpha_factor import (
    run_engine_Tcore,
    compute_loss,
    golden_section,
)

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"

T_PCS_EDGE = [40.0, 50.0, 90.0, 100.0, 110.0]
ALPHAS     = [0.20, 0.40, 0.60, 0.80]
C_COARSE_EXT = [0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00,
                1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80]
GSS_TOL = 0.005
EDGE_LO = C_COARSE_EXT[0]
EDGE_HI = C_COARSE_EXT[-1]


def find_folder(folder_names, alpha_u_arr, T_pc_F_arr, au, T_pc):
    n = len(folder_names)
    idx = next(i for i in range(n)
               if abs(alpha_u_arr[i] - au) < 0.01 and abs(T_pc_F_arr[i] - T_pc) < 0.5)
    return CW_DATA / str(folder_names[idx]), idx


def sweep_one_extended(folder, alpha_nom, T_pc_F, t_cw_hr, T_cw_F, Hu_eff):
    tag = f"α_u={alpha_nom:.2f}  T_pc={T_pc_F:.0f}°F"
    print(f"    {tag}  coarse[15] ...", end="", flush=True)

    coarse = {}
    for c in C_COARSE_EXT:
        t_e, T_e = run_engine_Tcore(folder, T_pc_F, alpha_nom * c, Hu_eff)
        coarse[c] = (compute_loss(t_e, T_e, t_cw_hr, T_cw_F), t_e, T_e)

    c_best = min(coarse, key=lambda k: coarse[k][0])
    print(f" coarse_best={c_best:.2f}", end="", flush=True)

    edge_flag = False
    if c_best == EDGE_LO or c_best == EDGE_HI:
        edge_flag = True
        c_opt = c_best
        n_gss = 0
        print(f"  EDGE-CLAMPED at {c_best:.2f}  (no GSS, flagged)")
    else:
        idx_b = C_COARSE_EXT.index(c_best)
        a = C_COARSE_EXT[idx_b - 1]
        b = C_COARSE_EXT[idx_b + 1]

        def loss_fn(c):
            t_e, T_e = run_engine_Tcore(folder, T_pc_F, alpha_nom * c, Hu_eff)
            return compute_loss(t_e, T_e, t_cw_hr, T_cw_F)

        c_opt, n_gss = golden_section(loss_fn, a, b, tol=GSS_TOL)
        print(f"  refine[{a:.2f},{b:.2f}] → c_opt={c_opt:.4f}  ({n_gss} GSS evals)")

    # Metrics at c_opt
    t_opt, T_opt = run_engine_Tcore(folder, T_pc_F, alpha_nom * c_opt, Hu_eff)
    loss_opt = compute_loss(t_opt, T_opt, t_cw_hr, T_cw_F)
    T_opt_i = np.interp(t_cw_hr, t_opt, T_opt)
    max_dT_opt = float(np.max(np.abs(T_opt_i - T_cw_F)))
    end_dT_opt = float(T_opt_i[-1] - T_cw_F[-1])

    # Metrics at c=1.0 (reuse from coarse)
    t_c1, T_c1 = coarse[1.00][1], coarse[1.00][2]
    T_c1_i = np.interp(t_cw_hr, t_c1, T_c1)
    max_dT_c1 = float(np.max(np.abs(T_c1_i - T_cw_F)))

    return {
        "folder":                    folder.name,
        "alpha_nom":                 alpha_nom,
        "T_pc_F":                    T_pc_F,
        "c_optimal_extended":        c_opt,
        "edge_flag":                 edge_flag,
        "loss_at_optimal_F":         loss_opt,
        "max_dT_at_optimal_F":       max_dT_opt,
        "endpoint_dT_at_optimal_F":  end_dT_opt,
        "max_dT_at_c1_F":            max_dT_c1,
        "n_coarse":                  len(C_COARSE_EXT),
        "n_gss":                     n_gss,
        "t_opt_hr":                  t_opt,
        "T_opt_F":                   T_opt,
        "t_c1_hr":                   t_c1,
        "T_c1_F":                    T_c1,
    }


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
    print(f"\n§4.11.6 Extended-range re-sweep")
    print(f"  Hu_eff (from §4.11 calibration): {HU_EFF:.2f} J/kg")
    print(f"  Coarse grid (n={len(C_COARSE_EXT)}): {C_COARSE_EXT}")
    print(f"  GSS tol: {GSS_TOL}")
    print(f"  Edge T_pc to re-sweep: {T_PCS_EDGE}")
    print(f"  α_u set: {ALPHAS}")
    print("=" * 100)

    results = []
    for au in ALPHAS:
        for T_pc in T_PCS_EDGE:
            folder, src_idx = find_folder(folder_names, alpha_u_arr, T_pc_F_arr, au, T_pc)
            r = sweep_one_extended(folder, au, T_pc, t_cw_hr,
                                   T_cw_F_all[src_idx], HU_EFF)
            r["src_idx"] = src_idx
            results.append(r)

    # Save NPZ
    n = len(results)
    T_opt_all = np.zeros((n, len(t_cw_hr)))
    T_c1_all  = np.zeros((n, len(t_cw_hr)))
    for i, r in enumerate(results):
        T_opt_all[i] = np.interp(t_cw_hr, r["t_opt_hr"], r["T_opt_F"])
        T_c1_all[i]  = np.interp(t_cw_hr, r["t_c1_hr"],  r["T_c1_F"])

    out_npz = DATA / "sweep_results_extended.npz"
    np.savez(
        out_npz,
        alpha_u=np.array([r["alpha_nom"] for r in results]),
        T_pc_F=np.array([r["T_pc_F"] for r in results]),
        folder_names=np.array([r["folder"] for r in results]),
        c_optimal_extended=np.array([r["c_optimal_extended"] for r in results]),
        edge_flag=np.array([r["edge_flag"] for r in results]),
        loss_at_optimal_F=np.array([r["loss_at_optimal_F"] for r in results]),
        max_dT_at_optimal_F=np.array([r["max_dT_at_optimal_F"] for r in results]),
        endpoint_dT_at_optimal_F=np.array([r["endpoint_dT_at_optimal_F"] for r in results]),
        max_dT_at_c1_F=np.array([r["max_dT_at_c1_F"] for r in results]),
        Hu_eff_J_kg=np.float64(HU_EFF),
        t_cw_hr=t_cw_hr,
        T_core_eng_at_optimal_F=T_opt_all,
        T_core_eng_at_c1_F=T_c1_all,
        c_coarse_grid=np.array(C_COARSE_EXT),
    )
    print(f"\nWrote {out_npz}")

    # CSV
    out_csv = DATA / "sweep_results_extended.csv"
    fields = ["folder", "alpha_nom", "T_pc_F", "c_optimal_extended", "edge_flag",
              "loss_at_optimal_F", "max_dT_at_optimal_F",
              "endpoint_dT_at_optimal_F", "max_dT_at_c1_F",
              "n_coarse", "n_gss"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow({k: r[k] for k in fields})
    print(f"Wrote {out_csv}")

    # Summary table
    print("\nExtended c_optimal table (edge T_pc only):")
    print(f"{'α_u':>6}  " + "  ".join(f"{int(T):>5}°F" for T in T_PCS_EDGE))
    for au in ALPHAS:
        row = [next(r for r in results
                    if abs(r["alpha_nom"] - au) < 0.01 and abs(r["T_pc_F"] - T) < 0.5)
               for T in T_PCS_EDGE]
        cells = []
        for r in row:
            mark = "*" if r["edge_flag"] else " "
            cells.append(f"{r['c_optimal_extended']:>6.3f}{mark}")
        print(f"{au:>6.2f}  " + "  ".join(cells))
    print("(* = hit search-range edge {0.40, 1.80} — true optimum outside extended bracket.)")
    print("\nDone.")


if __name__ == "__main__":
    main()
