"""Sprint 8 Stage 1 §3.8 — engine-side α(t) reference trajectories.

Computes α(t) under two thermal regimes for each of the 3 Stage 1 datasets:
  1. Isothermal at T_pl = 73°F   — what CW actually experienced (Hu=1 J/kg ≈ no heating)
  2. Adiabatic with typical Hu   — hypothetical reference (not what CW ran)

Reports α at t = 24, 48, 84, 120, 168 hr per dataset.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import numpy as np
import csv

from cw_scenario_loader import parse_cw_dat
from thermal_engine_2d import R_GAS, T_REF_K

BASE = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE, 'cw_data')
FINDINGS_DIR = os.path.join(BASE, 'findings')

DATASETS = [
    ('lowhyd',       'thermal_a_cal_73_100_lowhyd',          'α_u=0.10, τ=200, β=0.10'),
    ('highhyd',      'thermal_a_cal_73_100_highhyd',         'α_u=0.70, τ=10,  β=0.85'),
    ('highhyd_b010', 'thermal_a_cal_73_100_highhyd_beta010', 'α_u=0.70, τ=10,  β=0.10'),
]

MILESTONE_HR = [24, 48, 84, 120, 168]

DT_HR = 0.25
T_HOURS = np.arange(0, 168 + DT_HR, DT_HR)
RHO_CP_CONCRETE = 2350 * 900     # J/(m³·K), concrete density × Cp
TYPICAL_Hu_J_kg  = 350_000       # representative OPC, for adiabatic reference only
TYPICAL_Cc_kg_m3 = 400.0         # representative cementitious content kg/m³


def arrhenius(T_K: np.ndarray, Ea: float) -> np.ndarray:
    return np.exp(-Ea / R_GAS * (1.0 / T_K - 1.0 / T_REF_K))


def alpha_from_te(te: np.ndarray, tau: float, beta: float, au: float) -> np.ndarray:
    te_s = np.maximum(te, 1e-6)
    return au * np.exp(-(tau / te_s) ** beta)


def isothermal_trajectory(T_pl_F: float, Ea: float, tau: float, beta: float, au: float) -> np.ndarray:
    T_K = (T_pl_F - 32) * 5.0 / 9.0 + 273.15
    af  = float(arrhenius(np.array([T_K]), Ea)[0])
    te  = af * T_HOURS
    return alpha_from_te(te, tau, beta, au)


def adiabatic_trajectory(T_pl_F: float, Ea: float, tau: float, beta: float, au: float,
                          Hu_J_kg: float, Cc_kg_m3: float) -> np.ndarray:
    T_K  = (T_pl_F - 32) * 5.0 / 9.0 + 273.15
    te   = 0.0
    alph = 0.0
    alpha_arr = np.zeros_like(T_HOURS)
    for i in range(1, len(T_HOURS)):
        af      = float(arrhenius(np.array([T_K]), Ea)[0])
        te     += af * DT_HR
        alph    = float(alpha_from_te(np.array([te]), tau, beta, au)[0])
        alpha_arr[i] = alph
        dT = Hu_J_kg * Cc_kg_m3 / RHO_CP_CONCRETE * (alph - alpha_arr[i - 1])
        T_K += dT
    return alpha_arr


def idx_at(t_hr: float) -> int:
    return min(int(round(t_hr / DT_HR)), len(T_HOURS) - 1)


def main():
    os.makedirs(FINDINGS_DIR, exist_ok=True)

    rows = []
    print(f'\n{"Dataset":<18} {"Regime":<14}  ' +
          '  '.join(f'α(t={t}hr)' for t in MILESTONE_HR))
    print('-' * (18 + 14 + len(MILESTONE_HR) * 12))

    for label, folder, desc in DATASETS:
        dat_path = os.path.join(DATA_DIR, folder, 'input.dat')
        mix, geom, constr, _ = parse_cw_dat(dat_path)

        T_pl = constr.placement_temp_F
        Ea   = mix.activation_energy_J_mol
        tau  = mix.tau_hrs
        beta = mix.beta
        au   = mix.alpha_u

        # Isothermal: CW actual regime (Hu=1 J/kg → no self-heating)
        alpha_iso = isothermal_trajectory(T_pl, Ea, tau, beta, au)
        # Adiabatic with typical OPC Hu (hypothetical reference)
        alpha_adi = adiabatic_trajectory(T_pl, Ea, tau, beta, au, TYPICAL_Hu_J_kg, TYPICAL_Cc_kg_m3)

        for regime, alpha_arr in [('isothermal', alpha_iso), ('adiabatic_ref', alpha_adi)]:
            milestones = {t: round(float(alpha_arr[idx_at(t)]), 4) for t in MILESTONE_HR}
            row = {
                'dataset': label, 'desc': desc,
                'T_pl_F': T_pl, 'Ea_J_mol': Ea,
                'tau_hrs': tau, 'beta': beta, 'alpha_u': au,
                'regime': regime,
            }
            row.update({f'alpha_t{t}hr': milestones[t] for t in MILESTONE_HR})
            rows.append(row)

            val_str = '  '.join(f'{milestones[t]:>9.4f}' for t in MILESTONE_HR)
            print(f'{label:<18} {regime:<14}  {val_str}')

    # CSV
    fieldnames = ['dataset', 'desc', 'T_pl_F', 'Ea_J_mol', 'tau_hrs', 'beta', 'alpha_u', 'regime'] + \
                 [f'alpha_t{t}hr' for t in MILESTONE_HR]
    csv_path = os.path.join(FINDINGS_DIR, 'engine_alpha_reference.csv')
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f'\nWrote: {csv_path}')

    # Markdown
    md_lines = [
        '# Sprint 8 Stage 1 — Engine α(t) Reference (§3.8)\n',
        f'Placement temperature: 73°F for all three datasets.\n',
        f'- **isothermal**: CW actual regime (Hu=1 J/kg ≈ no self-heating, T fixed at T_pl).\n',
        f'- **adiabatic_ref**: hypothetical with Hu={TYPICAL_Hu_J_kg/1000:.0f} kJ/kg, '
        f'Cc={TYPICAL_Cc_kg_m3:.0f} kg/m³ (representative OPC, for context only).\n',
        '| Dataset | Regime | α(24hr) | α(48hr) | α(84hr) | α(120hr) | α(168hr) |',
        '|---|---|---|---|---|---|---|',
    ]
    for r in rows:
        vals = ' | '.join(f"{r[f'alpha_t{t}hr']:.4f}" for t in MILESTONE_HR)
        md_lines.append(f"| {r['dataset']} ({r['desc']}) | {r['regime']} | {vals} |")

    md_path = os.path.join(FINDINGS_DIR, 'engine_alpha_reference.md')
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    print(f'Wrote: {md_path}')
    print('\n' + '\n'.join(md_lines))


if __name__ == '__main__':
    main()
