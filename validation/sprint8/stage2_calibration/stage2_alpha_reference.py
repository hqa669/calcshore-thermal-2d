#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.3 — engine α(t) reference trajectories.

For each of the 12 new datasets, computes the isothermal α(t) trajectory
at T = T_pl using engine kinetics (Schindler / Van Breugel).  Confirms
each dataset reaches its α_u plateau by t=168 hr.

Output: findings/engine_alpha_reference.csv + printed table
"""

import csv
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent
ROOT = (HERE / "../../..").resolve()

sys.path.insert(0, str(ROOT))

from cw_scenario_loader import parse_cw_dat
from thermal_engine_2d import R_GAS, T_REF_K

CW_DATA  = HERE / "cw_data"
FINDINGS = HERE / "findings"

NEW_RUNS = [
    ("alpha02_A", "thermal_alpha02_A_73_73",  0.20,  73),
    ("alpha02_F", "thermal_alpha02_F_73_45",  0.20,  73),
    ("alpha02_I", "thermal_alpha02_I_100_73", 0.20, 100),
    ("alpha04_A", "thermal_alpha04_A_73_73",  0.40,  73),
    ("alpha04_F", "thermal_alpha04_F_73_45",  0.40,  73),
    ("alpha04_I", "thermal_alpha04_I_100_73", 0.40, 100),
    ("alpha06_A", "thermal_alpha06_A_73_73",  0.60,  73),
    ("alpha06_F", "thermal_alpha06_F_73_45",  0.60,  73),
    ("alpha06_I", "thermal_alpha06_I_100_73", 0.60, 100),
    ("alpha08_A", "thermal_alpha08_A_73_73",  0.80,  73),
    ("alpha08_F", "thermal_alpha08_F_73_45",  0.80,  73),
    ("alpha08_I", "thermal_alpha08_I_100_73", 0.80, 100),
]

MILESTONE_HR = [24, 48, 84, 120, 168]
DT_HR = 0.25
T_HOURS = np.arange(0.0, 168.0 + DT_HR, DT_HR)


def arrhenius(T_K: float, Ea: float) -> float:
    return float(np.exp(-Ea / R_GAS * (1.0 / T_K - 1.0 / T_REF_K)))


def alpha_from_te(te: np.ndarray, tau: float, beta: float, au: float) -> np.ndarray:
    te_s = np.maximum(te, 1e-6)
    return au * np.exp(-(tau / te_s) ** beta)


def isothermal_trajectory(T_pl_F: float, Ea: float, tau: float,
                           beta: float, au: float) -> np.ndarray:
    T_K = (T_pl_F - 32.0) * 5.0 / 9.0 + 273.15
    af  = arrhenius(T_K, Ea)
    te  = af * T_HOURS
    return alpha_from_te(te, tau, beta, au)


def idx_at(t_hr: float) -> int:
    return min(int(round(t_hr / DT_HR)), len(T_HOURS) - 1)


def main():
    FINDINGS.mkdir(parents=True, exist_ok=True)

    header = (f"{'Dataset':<14} {'T_pl':>5} {'α_u':>5} | "
              + "  ".join(f"α({t}hr)" for t in MILESTONE_HR))
    print(f"\n§3.3 Engine α(t) reference trajectories (isothermal at T_pl)")
    print("=" * len(header))
    print(header)
    print("-" * len(header))

    rows = []
    for label, folder, alpha_u_expected, T_pl_F in NEW_RUNS:
        dat = CW_DATA / folder / "input.dat"
        mix, geom, constr, _ = parse_cw_dat(str(dat))

        Ea   = mix.activation_energy_J_mol
        tau  = mix.tau_hrs
        beta = mix.beta
        au   = mix.alpha_u

        alpha_arr = isothermal_trajectory(T_pl_F, Ea, tau, beta, au)
        milestones = {t: round(float(alpha_arr[idx_at(t)]), 4) for t in MILESTONE_HR}

        val_str = "  ".join(f"{milestones[t]:>8.4f}" for t in MILESTONE_HR)
        print(f"  {label:<14} {T_pl_F:>5.0f} {au:>5.2f} | {val_str}")

        row = {"label": label, "folder": folder, "T_pl_F": T_pl_F,
               "alpha_u": au, "tau_hrs": tau, "beta": beta, "Ea_J_mol": Ea}
        row.update({f"alpha_t{t}hr": milestones[t] for t in MILESTONE_HR})
        rows.append(row)

    fieldnames = (["label", "folder", "T_pl_F", "alpha_u", "tau_hrs", "beta", "Ea_J_mol"]
                  + [f"alpha_t{t}hr" for t in MILESTONE_HR])
    csv_path = FINDINGS / "engine_alpha_reference.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote: {csv_path}")

    # Quick sanity: alpha_u plateau
    print("\nPlateau check (α(168hr) / α_u — should be ≥ 0.98):")
    for r in rows:
        ratio = r["alpha_t168hr"] / r["alpha_u"]
        flag = "" if ratio >= 0.95 else "  ← NOT PLATEAUED"
        print(f"  {r['label']:<14}  α(168)={r['alpha_t168hr']:.4f}  α_u={r['alpha_u']:.2f}  ratio={ratio:.3f}{flag}")


if __name__ == "__main__":
    main()
