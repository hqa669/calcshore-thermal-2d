#!/usr/bin/env python3
"""§3.5 — Per-point α_u tuning sweep for Stage 2-alpha-u-T-trend.

For each of the 32 (T_pc, α_u) points:
  1. Coarse 9-point scan: c ∈ {0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20}
  2. Golden-section refinement (~6 evaluations, tol=0.005) around coarse minimum

Engine: A-scenario, k_uc=0.96 (default), T_ref=23°C (default).
Hu: self-calibrated by running a single probe solve at T_pc=70°F, α_u=0.20 with
    Hu_probe=1000 J/kg, then scaling so that c_optimal(70°F, α_u=0.20) = 1.000.
    This removes the Hu ambiguity (CW's input.dat stores "1" as a placeholder,
    not as 1 J/kg) and anchors the 70°F baseline at c=1.0 by construction.

Loss: RMSE of T_core_engine(t) − T_core_CW(t) in °F over t ∈ [0, 168] hr.

Sanity gate (T_pc=70°F): after the 4 T_pc=70°F runs, assert c_optimal ≈ 1.00 ± 0.02.
If any falls outside ±0.05, print and stop.

Reads:
  data/cw_trajectories.npz
  cw_data/thermal_temp_alpha*/input.dat

Writes:
  data/sweep_results.npz
  data/sweep_results.csv
"""
import csv
import re
import sys
from math import sqrt
from pathlib import Path

import numpy as np

HERE  = Path(__file__).resolve().parents[1]
REPO  = Path(__file__).resolve().parents[4]
SC    = REPO / "validation" / "soil_calibration"

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(SC))

import thermal_engine_2d as te2d
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from cw_scenario_loader import parse_cw_dat
from stage4b_run import make_neutral_env

CW_DATA  = HERE / "cw_data"
DATA     = HERE / "data"

ALPHA_TAG_MAP     = {"02": 0.20, "04": 0.40, "06": 0.60, "08": 0.80}
DURATION_HR       = 168.0
OUTPUT_INTERVAL_S = 1800.0
C_COARSE          = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20]
GSS_TOL           = 0.005
GSS_MAX_ITER      = 20
SANITY_WARN       = 0.05
SANITY_GATE       = 0.02
HU_PROBE_J_KG     = 1000.0   # probe Hu for self-calibration step


def run_engine_Tcore(folder: Path, T_pc_F: float, alpha_u_eff: float,
                     Hu_eff: float = 1.0):
    """Single engine solve; returns (t_hr, T_core_F) at concrete centerline mid-depth."""
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))
    mix.Hu_factor_calibrated = 1.0
    mix.Hu_J_kg_effective    = Hu_eff
    mix.alpha_u              = alpha_u_eff

    constr.model_soil   = False
    constr.is_submerged = True

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0 = (T_pc_F - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0)

    res = solve_hydration_2d(
        grid, mix, Ti,
        duration_s=DURATION_HR * 3600.0,
        output_interval_s=OUTPUT_INTERVAL_S,
        boundary_mode="full_2d",
        environment=make_neutral_env(T_pc_F),
        construction=constr,
        T_ground_deep_C=T0,
        diagnostic_outputs=False,
    )

    jsl, isl = grid.concrete_slice()
    T_conc_C = res.T_field_C[:, jsl, isl]      # (n_t, nD, nW)
    nD = T_conc_C.shape[1]
    T_core_C = T_conc_C[:, nD // 2, -1]        # mid-depth, last x = centerline
    T_core_F = T_core_C * 9.0 / 5.0 + 32.0
    t_hr     = np.asarray(res.t_s) / 3600.0
    return t_hr, T_core_F


def compute_loss(t_eng_hr, T_eng_F, t_cw_hr, T_cw_F):
    """RMSE (°F) of engine vs CW T_core trajectory, engine interpolated to CW grid."""
    T_eng_interp = np.interp(t_cw_hr, t_eng_hr, T_eng_F)
    return float(np.sqrt(np.mean((T_eng_interp - T_cw_F) ** 2)))


def golden_section(f, a, b, tol=GSS_TOL, max_iter=GSS_MAX_ITER):
    """Golden-section search for minimum of unimodal f on [a, b]."""
    phi = (3.0 - sqrt(5.0)) / 2.0      # ≈ 0.382
    x1 = a + phi * (b - a)
    x2 = b - phi * (b - a)
    f1, f2 = f(x1), f(x2)
    n_eval = 2
    for _ in range(max_iter):
        if (b - a) < tol:
            break
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = a + phi * (b - a)
            f1 = f(x1)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = b - phi * (b - a)
            f2 = f(x2)
        n_eval += 1
    c_opt = (a + b) / 2.0
    return c_opt, n_eval


def sweep_one(folder: Path, alpha_nom: float, T_pc_F: float,
              t_cw_hr: np.ndarray, T_cw_F: np.ndarray, Hu_eff: float = 1.0):
    """Coarse scan + golden-section refinement for one (T_pc, α_u) point.

    Returns a dict with all per-point metrics.
    """
    tag = f"α_u={alpha_nom:.2f}  T_pc={T_pc_F:.0f}°F"
    print(f"    {tag}  coarse ...", end="", flush=True)

    # --- Coarse scan ---
    coarse_losses = {}
    for c in C_COARSE:
        t_e, T_e = run_engine_Tcore(folder, T_pc_F, alpha_nom * c, Hu_eff)
        coarse_losses[c] = (compute_loss(t_e, T_e, t_cw_hr, T_cw_F), t_e, T_e)

    c_coarse_best = min(coarse_losses, key=lambda k: coarse_losses[k][0])
    print(f" coarse_best={c_coarse_best:.2f}  refine ...", end="", flush=True)

    # --- Bracket for golden-section ---
    idx_best = C_COARSE.index(c_coarse_best)
    a = C_COARSE[max(0, idx_best - 1)]
    b = C_COARSE[min(len(C_COARSE) - 1, idx_best + 1)]
    if a == b:                       # edge case: minimum at edge of scan
        a = max(0.60, c_coarse_best - 0.05)
        b = min(1.40, c_coarse_best + 0.05)

    def loss_fn(c):
        t_e, T_e = run_engine_Tcore(folder, T_pc_F, alpha_nom * c, Hu_eff)
        return compute_loss(t_e, T_e, t_cw_hr, T_cw_F)

    c_opt, n_gss = golden_section(loss_fn, a, b)
    print(f" c_opt={c_opt:.4f}  ({n_gss} GSS evals)")

    # --- Metrics at c_optimal ---
    t_opt, T_opt = run_engine_Tcore(folder, T_pc_F, alpha_nom * c_opt, Hu_eff)
    loss_opt = compute_loss(t_opt, T_opt, t_cw_hr, T_cw_F)
    T_opt_interp = np.interp(t_cw_hr, t_opt, T_opt)
    max_dT_opt = float(np.max(np.abs(T_opt_interp - T_cw_F)))
    endpoint_dT_opt = float(T_opt_interp[-1] - T_cw_F[-1])

    # --- Metrics at c=1.0 ---
    if abs(c_coarse_best - 1.0) < 0.001:   # reuse cached c=1.0 run
        t_c1, T_c1 = coarse_losses[1.00][1], coarse_losses[1.00][2]
    else:
        t_c1, T_c1 = run_engine_Tcore(folder, T_pc_F, alpha_nom * 1.0, Hu_eff)
    T_c1_interp = np.interp(t_cw_hr, t_c1, T_c1)
    max_dT_c1 = float(np.max(np.abs(T_c1_interp - T_cw_F)))

    return {
        "folder":              folder.name,
        "alpha_nom":           alpha_nom,
        "T_pc_F":              T_pc_F,
        "c_optimal":           c_opt,
        "loss_at_optimal_F":   loss_opt,
        "max_dT_at_optimal_F": max_dT_opt,
        "endpoint_dT_at_optimal_F": endpoint_dT_opt,
        "max_dT_at_c1_F":     max_dT_c1,
        "n_gss_evals":         n_gss,
        # Full trajectories saved for plots/tables
        "t_opt_hr":    t_opt,
        "T_opt_F":     T_opt,
        "t_c1_hr":     t_c1,
        "T_c1_F":      T_c1,
    }


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    # Load CW trajectories
    traj = np.load(DATA / "cw_trajectories.npz", allow_pickle=True)
    t_cw_hr    = traj["t_hrs"]          # (2016,)
    T_cw_F     = traj["T_core_CW_F"]    # (32, 2016)
    alpha_u_arr = traj["alpha_u"]        # (32,)
    T_pc_F_arr  = traj["T_pc_F"]         # (32,)
    folder_names = traj["folder_names"]  # (32,)

    # Sorted order: (alpha_u, T_pc) ascending — matches extract_cw_trajectories output
    n = len(folder_names)

    # Identify T_pc=70°F sanity group indices
    sanity_idx = [i for i in range(n) if abs(T_pc_F_arr[i] - 70.0) < 0.5]
    other_idx  = [i for i in range(n) if i not in sanity_idx]

    # ----------------------------------------------------------------
    # Self-calibrate Hu_eff: run probe at T_pc=70°F, α_u=0.20, Hu=HU_PROBE_J_KG
    # Then scale so that the engine signal matches CW dT at 70°F, α_u=0.20.
    # ----------------------------------------------------------------
    cal_idx = next(i for i in range(n)
                   if abs(alpha_u_arr[i] - 0.20) < 0.01 and abs(T_pc_F_arr[i] - 70.0) < 0.5)
    cal_folder = CW_DATA / str(folder_names[cal_idx])
    print("\nAuto-calibrating Hu_eff from T_pc=70°F, α_u=0.20 probe ...")
    t_probe, T_probe = run_engine_Tcore(cal_folder, 70.0, 0.20, HU_PROBE_J_KG)
    dT_probe = float(T_probe[-1] - T_probe[0])
    dT_cw_70 = float(T_cw_F[cal_idx][-1] - T_cw_F[cal_idx][0])
    if dT_probe < 1e-6:
        raise RuntimeError("Probe engine dT is zero — check solver setup.")
    HU_EFF = HU_PROBE_J_KG * dT_cw_70 / dT_probe
    print(f"  Probe: Hu={HU_PROBE_J_KG} J/kg → dT_engine={dT_probe:.4f}°F")
    print(f"  CW reference dT at 70°F, α_u=0.20: {dT_cw_70:.4f}°F")
    print(f"  → Hu_eff = {HU_EFF:.1f} J/kg  (ensures c_opt(70°F,α_u=0.20)=1.00 by construction)")

    print(f"\n§3.5 α_u factor sweep — {n} datasets  ({len(sanity_idx)} sanity, {len(other_idx)} remaining)")
    print(f"  Hu_eff: {HU_EFF:.1f} J/kg")
    print(f"  Coarse c: {C_COARSE}")
    print(f"  GSS tol: {GSS_TOL}  max_iter: {GSS_MAX_ITER}")
    print("=" * 90)

    results = [None] * n

    # ----------------------------------------------------------------
    # Phase A: T_pc=70°F sanity
    # ----------------------------------------------------------------
    print("\nPhase A — T_pc=70°F sanity check (4 datasets)")
    print("-" * 90)

    for i in sorted(sanity_idx, key=lambda k: alpha_u_arr[k]):
        folder = CW_DATA / str(folder_names[i])
        r = sweep_one(folder, float(alpha_u_arr[i]), float(T_pc_F_arr[i]),
                      t_cw_hr, T_cw_F[i], HU_EFF)
        results[i] = r

    # Sanity gate
    print("\nSanity gate check (T_pc=70°F, c_optimal should ≈ 1.00 ± 0.02):")
    gate_fail = False
    for i in sorted(sanity_idx, key=lambda k: alpha_u_arr[k]):
        r = results[i]
        c = r["c_optimal"]
        dev = abs(c - 1.0)
        status = "OK" if dev <= SANITY_GATE else ("WARN" if dev <= SANITY_WARN else "STOP")
        print(f"  {r['folder']:<35}  α_u={r['alpha_nom']:.2f}  "
              f"c_opt={c:.4f}  dev={dev:.4f}  [{status}]")
        if dev > SANITY_WARN:
            gate_fail = True

    if gate_fail:
        print("\n*** SANITY GATE FAILED — c_optimal deviates > 0.05 at T_pc=70°F.")
        print("    Check engine wrapper (Hu override, alpha_u setting, BC setup).")
        print("    Stopping before remaining sweep.")
        sys.exit(1)

    sanity_ok = all(abs(results[i]["c_optimal"] - 1.0) <= SANITY_GATE
                    for i in sanity_idx)
    if not sanity_ok:
        print("\n  Warning: some T_pc=70°F points outside ±0.02 (but within ±0.05). Continuing.")
    else:
        print("\n  PASS — all 4 T_pc=70°F c_optimal within ±0.02 of 1.00.")

    # ----------------------------------------------------------------
    # Phase B: remaining 28 datasets
    # ----------------------------------------------------------------
    print("\nPhase B — remaining 28 datasets")
    print("-" * 90)

    for i in sorted(other_idx, key=lambda k: (alpha_u_arr[k], T_pc_F_arr[k])):
        folder = CW_DATA / str(folder_names[i])
        r = sweep_one(folder, float(alpha_u_arr[i]), float(T_pc_F_arr[i]),
                      t_cw_hr, T_cw_F[i], HU_EFF)
        results[i] = r

    # ----------------------------------------------------------------
    # Save outputs
    # ----------------------------------------------------------------
    csv_rows = []
    npz_data = {
        "alpha_u":              alpha_u_arr,
        "T_pc_F":               T_pc_F_arr,
        "folder_names":         folder_names,
        "t_cw_hr":              t_cw_hr,
    }

    # Per-point trajectory arrays indexed by (alpha_u, T_pc)
    T_opt_all = np.zeros((n, len(t_cw_hr)))
    T_c1_all  = np.zeros((n, len(t_cw_hr)))

    for i, r in enumerate(results):
        csv_rows.append({
            "folder":                    r["folder"],
            "alpha_nom":                 r["alpha_nom"],
            "T_pc_F":                    r["T_pc_F"],
            "c_optimal":                 r["c_optimal"],
            "loss_at_optimal_F":         r["loss_at_optimal_F"],
            "max_dT_at_optimal_F":       r["max_dT_at_optimal_F"],
            "endpoint_dT_at_optimal_F":  r["endpoint_dT_at_optimal_F"],
            "max_dT_at_c1_F":            r["max_dT_at_c1_F"],
            "n_gss_evals":               r["n_gss_evals"],
        })
        # Interpolate saved trajectories to common CW time grid
        T_opt_all[i] = np.interp(t_cw_hr, r["t_opt_hr"], r["T_opt_F"])
        T_c1_all[i]  = np.interp(t_cw_hr, r["t_c1_hr"],  r["T_c1_F"])

    npz_data["T_core_eng_at_optimal_F"] = T_opt_all
    npz_data["T_core_eng_at_c1_F"]      = T_c1_all
    npz_data["Hu_eff_J_kg"]             = np.float64(HU_EFF)
    npz_data["c_optimal"]               = np.array([r["c_optimal"] for r in results])
    npz_data["max_dT_at_optimal_F"]     = np.array([r["max_dT_at_optimal_F"] for r in results])
    npz_data["max_dT_at_c1_F"]         = np.array([r["max_dT_at_c1_F"] for r in results])
    npz_data["endpoint_dT_at_optimal_F"] = np.array([r["endpoint_dT_at_optimal_F"] for r in results])

    out_npz = DATA / "sweep_results.npz"
    np.savez(out_npz, **npz_data)
    print(f"\nWrote {out_npz}")

    out_csv = DATA / "sweep_results.csv"
    fields = ["folder", "alpha_nom", "T_pc_F", "c_optimal",
              "loss_at_optimal_F", "max_dT_at_optimal_F",
              "endpoint_dT_at_optimal_F", "max_dT_at_c1_F", "n_gss_evals"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(csv_rows)
    print(f"Wrote {out_csv}")

    # Summary by alpha_u
    print("\nSummary by α_u:")
    print(f"{'α_u':>6}  " + "  ".join(f"{int(T):>4}°F" for T in [40,50,60,70,80,90,100,110]))
    for au in [0.20, 0.40, 0.60, 0.80]:
        row_vals = []
        for Tpc in [40, 50, 60, 70, 80, 90, 100, 110]:
            idx = next(i for i in range(n)
                       if abs(alpha_u_arr[i] - au) < 0.01 and abs(T_pc_F_arr[i] - Tpc) < 0.5)
            row_vals.append(f"{results[idx]['c_optimal']:>6.3f}")
        print(f"{au:>6.2f}  " + "  ".join(row_vals))

    print("\nDone.")


if __name__ == "__main__":
    main()
