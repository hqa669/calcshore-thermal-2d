"""Sprint 8 Stage 2-prep §3.4 — T_REF_K back-calculation from observed ratio discrepancy.

Engine uses T_REF_K=296.15 K (23°C). CW's T_ref is undocumented.
At T_pl=73°F (22.78°C), different T_ref values produce different Arrhenius factors
and thus different α(t) trajectories despite identical Schindler parameters.

Because the Schindler function is nonlinear, T_ref differences do NOT cancel in the
ratio test — they produce pair-specific Δα distortions that explain the ratio
discrepancies observed in §3.2/§3.3.

Method:
  For each T_ref candidate, recompute α(t) for all 3 datasets, derive Δα ratios,
  compare to observed CW ΔT ratios (Stage 1 data), report discrepancy.
  Best-fit T_ref = minimizes RMS discrepancy across all ratio tests, stats, timesteps.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import csv
import numpy as np

from thermal_engine_2d import R_GAS

BASE    = os.path.dirname(__file__)
STAGE1  = os.path.join(BASE, '../stage1_cw_k_alpha_empirical/findings')
FINDINGS = os.path.join(BASE, 'findings')

T_DIAG_HR  = [24, 84, 168]
DT_HR      = 0.25
T_HOURS    = np.arange(0, 168.0 + DT_HR, DT_HR)
STATS      = ['max_abs_dT_F', 'mean_abs_dT_F', 'rms_dT_F']
STAT_LABELS = {'max_abs_dT_F': 'max|ΔT|', 'mean_abs_dT_F': 'mean|ΔT|', 'rms_dT_F': 'RMS ΔT'}

# Schindler parameters for each dataset (from engine_alpha_reference.csv)
DATASETS = {
    'lowhyd':       {'tau': 200.0, 'beta': 0.10, 'au': 0.10, 'Ea': 50000.0},
    'highhyd':      {'tau': 10.0,  'beta': 0.85, 'au': 0.70, 'Ea': 50000.0},
    'highhyd_b010': {'tau': 10.0,  'beta': 0.10, 'au': 0.70, 'Ea': 50000.0},
}
T_PL_F   = 73.0
T_PL_K   = (T_PL_F - 32) * 5/9 + 273.15  # 295.928 K

# T_ref candidates (K)
T_REF_CANDIDATES = [
    (293.15, '20.00°C (ACI 308 / Schindler-Folliard 2002 calibration)'),
    (294.15, '21.00°C'),
    (294.26, '21.11°C (70°F — Schindler-Folliard paper placement temp)'),
    (295.15, '22.00°C'),
    (296.15, '23.00°C (engine default T_REF_K)'),
    (297.15, '24.00°C'),
    (298.15, '25.00°C'),
]


def arrhenius_tref(T_K: float, Ea: float, T_ref_K: float) -> float:
    return float(np.exp(-Ea / R_GAS * (1.0 / T_K - 1.0 / T_ref_K)))


def alpha_from_te(te_arr: np.ndarray, tau: float, beta: float, au: float) -> np.ndarray:
    te_s = np.maximum(te_arr, 1e-9)
    return au * np.exp(-(tau / te_s) ** beta)


def isothermal_alpha(T_pl_K: float, Ea: float, tau: float, beta: float, au: float,
                     T_ref_K: float) -> np.ndarray:
    af = arrhenius_tref(T_pl_K, Ea, T_ref_K)
    te = af * T_HOURS
    return alpha_from_te(te, tau, beta, au)


def alpha_at_t(alpha_arr: np.ndarray, t_hr: int) -> float:
    idx = int(round(t_hr / DT_HR))
    return float(alpha_arr[idx])


def load_dT_stats():
    rows = {}
    with open(os.path.join(STAGE1, 'pairwise_dT_stats.csv')) as f:
        for r in csv.DictReader(f):
            key = (r['pair'], int(float(r['t_hr'])))
            rows[key] = {s: float(r[s]) for s in STATS}
    return rows


def compute_ratio_discrepancies(dT_stats, T_ref_K):
    """For a given T_ref, compute ratio discrepancies for all tests/stats/timesteps."""
    alphas = {}
    for ds, p in DATASETS.items():
        arr = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_K)
        alphas[ds] = {t: alpha_at_t(arr, t) for t in T_DIAG_HR}

    ratio_tests = [
        ('A/B', 'A', 'B', 'highhyd', 'lowhyd', 'highhyd_b010', 'lowhyd'),
        ('A/C', 'A', 'C', 'highhyd', 'lowhyd', 'highhyd',      'highhyd_b010'),
    ]

    discrep = []
    detail = []
    for test_name, pA, pB, da_num, db_num, da_den, db_den in ratio_tests:
        for t in T_DIAG_HR:
            d_alpha_num = alphas[da_num][t] - alphas[db_num][t]
            d_alpha_den = alphas[da_den][t] - alphas[db_den][t]
            exp_ratio = d_alpha_num / d_alpha_den if d_alpha_den != 0 else float('nan')
            for stat in STATS:
                val_num = dT_stats[(pA, t)][stat]
                val_den = dT_stats[(pB, t)][stat]
                obs_ratio = val_num / val_den if val_den != 0 else float('nan')
                d_pct = (obs_ratio / exp_ratio - 1.0) * 100.0
                discrep.append(d_pct)
                detail.append({
                    'test': test_name, 'stat': STAT_LABELS[stat], 't_hr': t,
                    'obs_ratio': obs_ratio, 'exp_ratio': exp_ratio, 'discrepancy_pct': d_pct,
                })
    return discrep, detail


def main():
    os.makedirs(FINDINGS, exist_ok=True)
    dT = load_dT_stats()

    # Engine reference α values (T_ref=23°C)
    engine_tref = 296.15
    engine_alphas = {}
    for ds, p in DATASETS.items():
        arr = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], engine_tref)
        engine_alphas[ds] = arr

    print(f'\n=== Stage 2-prep §3.4 T_REF_K Back-calculation ===')
    print(f'\nT_pl = {T_PL_F}°F = {T_PL_K:.3f} K\n')

    # Arrhenius factor table
    print(f'{"T_ref":>8}  {"T_ref (K)":>10}  {"Arrhenius factor":>18}  Description')
    print('-' * 78)
    for T_ref_K, label in T_REF_CANDIDATES:
        af = arrhenius_tref(T_PL_K, 50000.0, T_ref_K)
        marker = ' ← engine' if abs(T_ref_K - 296.15) < 0.01 else ''
        print(f'{label.split("(")[0].strip():>8}  {T_ref_K:>10.2f}  {af:>18.6f}  {label}{marker}')

    # Alpha values at diagnostic timesteps for each T_ref
    print(f'\n--- α(t) at T_pl={T_PL_F}°F for each T_ref candidate ---')
    for ds, p in DATASETS.items():
        print(f'\nDataset: {ds} (α_u={p["au"]}, τ={p["tau"]}, β={p["beta"]})')
        print(f'  {"T_ref":>22}  {"α(24hr)":>10}  {"α(84hr)":>10}  {"α(168hr)":>10}')
        for T_ref_K, label in T_REF_CANDIDATES:
            arr = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_K)
            marker = ' ←' if abs(T_ref_K - 296.15) < 0.01 else ''
            print(f'  {label[:22]:>22}  {alpha_at_t(arr,24):>10.4f}  '
                  f'{alpha_at_t(arr,84):>10.4f}  {alpha_at_t(arr,168):>10.4f}{marker}')

    # Ratio discrepancy vs observed CW ΔT ratios
    print(f'\n--- Ratio discrepancy vs observed CW ΔT ratios ---')
    print(f'{"T_ref":>22}  {"max|disc|":>10}  {"RMS disc":>10}  {"mean disc":>10}')
    print('-' * 60)

    results = []
    for T_ref_K, label in T_REF_CANDIDATES:
        discrep, detail = compute_ratio_discrepancies(dT, T_ref_K)
        max_d  = max(abs(d) for d in discrep)
        rms_d  = float(np.sqrt(np.mean([d**2 for d in discrep])))
        mean_d = float(np.mean(discrep))
        marker = ' ← engine' if abs(T_ref_K - 296.15) < 0.01 else ''
        print(f'{label[:22]:>22}  {max_d:>9.2f}%  {rms_d:>9.2f}%  {mean_d:>+9.2f}%{marker}')
        results.append({'T_ref_K': T_ref_K, 'label': label, 'max_d': max_d, 'rms_d': rms_d,
                        'mean_d': mean_d, 'discrep': discrep, 'detail': detail})

    # Best-fit T_ref
    best = min(results, key=lambda r: r['rms_d'])
    print(f'\nBest-fit T_ref (min RMS discrepancy): {best["label"]}')
    print(f'  Max |discrepancy|: {best["max_d"]:.2f}%')
    print(f'  RMS discrepancy:   {best["rms_d"]:.2f}%')

    # Implied max |Δα|/α between engine (23°C) and best-fit T_ref
    print(f'\n--- Implied |Δα|/α: engine (23°C) vs best-fit T_ref ({best["label"].split("(")[0].strip()}) ---')
    print(f'Dataset            t=24hr  t=84hr  t=168hr  max(0-168hr)')
    for ds, p in DATASETS.items():
        arr_engine = engine_alphas[ds]
        arr_bestfit = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], best['T_ref_K'])
        # Avoid div by zero at t=0
        valid = arr_engine > 1e-6
        frac = np.where(valid, np.abs(arr_bestfit - arr_engine) / arr_engine, 0.0)
        at24  = abs(alpha_at_t(arr_bestfit, 24)  - alpha_at_t(arr_engine, 24))  / max(alpha_at_t(arr_engine, 24), 1e-9)
        at84  = abs(alpha_at_t(arr_bestfit, 84)  - alpha_at_t(arr_engine, 84))  / max(alpha_at_t(arr_engine, 84), 1e-9)
        at168 = abs(alpha_at_t(arr_bestfit, 168) - alpha_at_t(arr_engine, 168)) / max(alpha_at_t(arr_engine, 168), 1e-9)
        peak  = float(np.max(frac[int(1/DT_HR):]))  # skip t<1 hr
        print(f'{ds:<18} {at24*100:>5.1f}%  {at84*100:>5.1f}%  {at168*100:>6.1f}%  {peak*100:>9.1f}%')

    # --- Write Markdown ---
    md_lines = ['# Sprint 8 Stage 2-prep — §3.4 T_REF_K Back-calculation\n',
                '## Method\n',
                f'Engine uses T_REF_K = 296.15 K (23°C). CW\'s reference temperature is undocumented.',
                f'At T_pl = {T_PL_F}°F ({T_PL_K:.2f} K), different T_ref values shift the Arrhenius',
                'factor, stretching/compressing the equivalent-age axis nonlinearly across datasets.',
                'Best-fit T_ref minimizes RMS discrepancy across all ratio tests, statistics, and timesteps.\n',
                '## Arrhenius Factor at T_pl = 73°F\n',
                '| T_ref | T_ref (K) | Arrhenius factor |',
                '|---|---|---|']
    for T_ref_K, label in T_REF_CANDIDATES:
        af = arrhenius_tref(T_PL_K, 50000.0, T_ref_K)
        marker = ' ← engine' if abs(T_ref_K - 296.15) < 0.01 else ''
        md_lines.append(f'| {label} | {T_ref_K:.2f} | {af:.6f}{marker} |')

    md_lines.append('\n## Ratio Discrepancy Summary\n')
    md_lines.append('| T_ref | max\\|disc\\| | RMS disc | mean disc |')
    md_lines.append('|---|---|---|---|')
    for r in results:
        marker = ' ← engine' if abs(r['T_ref_K'] - 296.15) < 0.01 else ''
        md_lines.append(f'| {r["label"]} | {r["max_d"]:.2f}% | {r["rms_d"]:.2f}% | {r["mean_d"]:+.2f}%{marker} |')

    md_lines.append(f'\n## Best-fit T_ref\n')
    md_lines.append(f'**{best["label"]}** — minimizes RMS discrepancy across all 18 ratio test values.')
    md_lines.append(f'Max |discrepancy| at best-fit: **{best["max_d"]:.2f}%**\n')

    md_lines.append('## Implied |Δα|/α vs Engine (23°C)\n')
    md_lines.append('| Dataset | t=24hr | t=84hr | t=168hr | max (0–168hr) |')
    md_lines.append('|---|---|---|---|---|')
    for ds, p in DATASETS.items():
        arr_engine = engine_alphas[ds]
        arr_bestfit = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], best['T_ref_K'])
        valid = arr_engine > 1e-6
        frac = np.where(valid, np.abs(arr_bestfit - arr_engine) / arr_engine, 0.0)
        at24  = abs(alpha_at_t(arr_bestfit,24)  - alpha_at_t(arr_engine,24))  / max(alpha_at_t(arr_engine,24),1e-9)
        at84  = abs(alpha_at_t(arr_bestfit,84)  - alpha_at_t(arr_engine,84))  / max(alpha_at_t(arr_engine,84),1e-9)
        at168 = abs(alpha_at_t(arr_bestfit,168) - alpha_at_t(arr_engine,168)) / max(alpha_at_t(arr_engine,168),1e-9)
        peak  = float(np.max(frac[int(1/DT_HR):]))
        md_lines.append(f'| {ds} | {at24*100:.1f}% | {at84*100:.1f}% | {at168*100:.1f}% | {peak*100:.1f}% |')

    md_path = os.path.join(FINDINGS, 't_ref_analysis.md')
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    print(f'\nWrote: {md_path}')

    # Write CSV for best-fit detail
    csv_path = os.path.join(FINDINGS, 't_ref_analysis.csv')
    fieldnames = ['T_ref_K', 'label', 'max_d_pct', 'rms_d_pct', 'mean_d_pct']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow({'T_ref_K': r['T_ref_K'], 'label': r['label'],
                        'max_d_pct': round(r['max_d'],4), 'rms_d_pct': round(r['rms_d'],4),
                        'mean_d_pct': round(r['mean_d'],4)})
    print(f'Wrote: {csv_path}')


if __name__ == '__main__':
    main()
