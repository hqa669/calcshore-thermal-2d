"""Sprint 8 Stage 2-prep v2 — Path D: CW bulk-shift inversion for direct α(t) comparison.

Path C (constant-k engine inversion at CL mid-depth) was degenerate: engine T at
di=24 (12.2 m depth) is insensitive to k(α) in an 80-ft mat over 168 hours. Path D
instead uses the CW mid-depth bulk-shift signal directly.

Key observation from Path C (stage2_prep_alpha_inversion.py §3.A): CW di=24 wi=0
shows a uniform temperature offset from T_placement (73°F) that is linear in α across
all 3 datasets at each timestep:

  ΔT_bulk[dataset, t] = C(t) · α[dataset, t]     (CV ≈ 1–5% across datasets)

C(t) is a per-timestep proportionality coefficient calibrated by fitting a single T_ref
such that Schindler-isothermal α values make ΔT_bulk/α constant across 3 datasets.
CW α(effective) = ΔT_bulk[d, t] / C(t) = α_iso[d, T_ref(t), t].

Delivers the §3.6 table:
  Dataset | t (hr) | Engine α (23°C) | Engine α (21°C) | CW α (eff) | Δ@23°C | Δ@21°C

No new engine runs required. Pure analysis over cached CW data.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import csv
import pickle
import numpy as np

from thermal_engine_2d import R_GAS

BASE    = os.path.dirname(__file__)
STAGE1  = os.path.join(BASE, '../stage1_cw_k_alpha_empirical/findings')
FINDINGS = os.path.join(BASE, 'findings')

T_DIAG_HR = [24, 48, 84, 120, 168]
DT_HR     = 0.25
T_HOURS   = np.arange(0, 168.0 + DT_HR, DT_HR)

T_PL_F   = 73.0
T_PL_K   = (T_PL_F - 32) * 5 / 9 + 273.15   # 295.928 K
T_PLACEMENT = 73.0                             # anchor for ΔT_bulk

DATASETS = {
    'lowhyd':       {'tau': 200.0, 'beta': 0.10, 'au': 0.10, 'Ea': 50000.0},
    'highhyd':      {'tau': 10.0,  'beta': 0.85, 'au': 0.70, 'Ea': 50000.0},
    'highhyd_b010': {'tau': 10.0,  'beta': 0.10, 'au': 0.70, 'Ea': 50000.0},
}

# T_ref sweep for per-timestep best-fit (°C → K)
T_REF_SWEEP_C = np.arange(19.0, 23.25, 0.25)   # 17 candidates: 19.0, 19.25, ..., 23.0

LINEARITY_CV_GATE = 5.0    # % — linearity sanity gate


# ---------------------------------------------------------------------------
# Helpers: Schindler-isothermal α at parameterised T_ref
# ---------------------------------------------------------------------------

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


def alpha_at_t(alpha_arr: np.ndarray, t_hr: float) -> float:
    idx = int(round(t_hr / DT_HR))
    return float(alpha_arr[min(idx, len(alpha_arr) - 1)])


# ---------------------------------------------------------------------------
# CW bulk shift extraction
# ---------------------------------------------------------------------------

def load_cw_bulk_shifts(cw_fields: dict) -> dict:
    """Return ΔT_bulk[dataset][t_hr] = T_CW(di=24, wi=0, t) − 73.0°F."""
    bulk = {}
    for ds in DATASETS:
        v = cw_fields[ds]
        bulk[ds] = {}
        for t in T_DIAG_HR:
            idx = int(np.argmin(np.abs(v.time_hrs - t)))
            T_cw = float(v.T_field_F[idx, 24, 0])
            bulk[ds][t] = T_cw - T_PLACEMENT
    return bulk


# ---------------------------------------------------------------------------
# §D.1  Linearity sanity check
# ---------------------------------------------------------------------------

def linearity_check(bulk_shifts: dict) -> dict:
    """Compute CV of ΔT_bulk/α(21°C) across datasets, per timestep."""
    T_ref_21_K = 294.15
    results = {}
    for t in T_DIAG_HR:
        ratios = []
        for ds, p in DATASETS.items():
            arr = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_21_K)
            a = alpha_at_t(arr, t)
            dT = bulk_shifts[ds][t]
            r = dT / a if a > 1e-9 else float('nan')
            ratios.append(r)
        ratios = [r for r in ratios if not np.isnan(r)]
        mean_r = float(np.mean(ratios))
        std_r  = float(np.std(ratios, ddof=0))
        cv_pct = (std_r / abs(mean_r) * 100.0) if abs(mean_r) > 1e-9 else float('nan')
        results[t] = {'ratios': ratios, 'mean': mean_r, 'std': std_r, 'cv_pct': cv_pct}
    return results


# ---------------------------------------------------------------------------
# §D.2  Per-timestep best-fit T_ref
# ---------------------------------------------------------------------------

def fit_tref_per_timestep(bulk_shifts: dict) -> list:
    """For each t, sweep T_ref and find T_ref minimizing CV of ΔT_bulk/α."""
    rows = []
    for t in T_DIAG_HR:
        best_tref = None
        best_cv   = float('inf')
        best_C    = None
        for T_ref_C in T_REF_SWEEP_C:
            T_ref_K = T_ref_C + 273.15
            ratios = []
            for ds, p in DATASETS.items():
                arr = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_K)
                a = alpha_at_t(arr, t)
                dT = bulk_shifts[ds][t]
                r = dT / a if a > 1e-9 else float('nan')
                ratios.append(r)
            valid = [r for r in ratios if not np.isnan(r)]
            if len(valid) < 2:
                continue
            mean_r = float(np.mean(valid))
            std_r  = float(np.std(valid, ddof=0))
            cv_pct = (std_r / abs(mean_r) * 100.0) if abs(mean_r) > 1e-9 else float('inf')
            if cv_pct < best_cv:
                best_cv   = cv_pct
                best_tref = T_ref_C
                best_C    = mean_r
        rows.append({
            't_hr': t,
            'best_T_ref_C': round(best_tref, 2),
            'best_CV_pct': round(best_cv, 3),
            'C_t_F_per_alpha': round(best_C, 4) if best_C else None,
        })
    return rows


# ---------------------------------------------------------------------------
# §D.3 + §D.4  CW effective α and headline table
# ---------------------------------------------------------------------------

def build_comparison_table(bulk_shifts: dict, tref_fits: list) -> list:
    """Build 15-row direct comparison table.

    For each (dataset, t): engine α at 23°C, engine α at 21°C, CW α(effective),
    Δ at 23°C (signed), Δ at 21°C (signed).
    """
    T_ref_23_K = 296.15
    T_ref_21_K = 294.15

    fit_by_t = {r['t_hr']: r for r in tref_fits}

    rows = []
    for ds, p in DATASETS.items():
        arr_23 = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_23_K)
        arr_21 = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_21_K)
        for t in T_DIAG_HR:
            fit = fit_by_t[t]
            T_ref_best_K = fit['best_T_ref_C'] + 273.15
            C_t = fit['C_t_F_per_alpha']
            arr_best = isothermal_alpha(T_PL_K, p['Ea'], p['tau'], p['beta'], p['au'], T_ref_best_K)

            a23   = alpha_at_t(arr_23, t)
            a21   = alpha_at_t(arr_21, t)
            a_cw  = alpha_at_t(arr_best, t)   # CW effective α (Mode 1 — per-timestep T_ref)

            # Cross-check: CW α via direct bulk-shift / C(t)
            dT_bulk = bulk_shifts[ds][t]
            a_cw_direct = dT_bulk / C_t if C_t and C_t > 0 else float('nan')

            delta_23 = a_cw - a23
            delta_21 = a_cw - a21

            rows.append({
                'dataset': ds,
                't_hr': t,
                'engine_alpha_23C': round(a23, 4),
                'engine_alpha_21C': round(a21, 4),
                'cw_alpha_effective': round(a_cw, 4),
                'cw_alpha_direct': round(a_cw_direct, 4),
                'delta_at_23C': round(delta_23, 4),
                'delta_at_21C': round(delta_21, 4),
                'delta_at_23C_pct': round(delta_23 / a23 * 100 if a23 > 0 else float('nan'), 2),
                'delta_at_21C_pct': round(delta_21 / a21 * 100 if a21 > 0 else float('nan'), 2),
            })
    return rows


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_linearity_md(lin: dict, path: str) -> None:
    lines = [
        '# Sprint 8 Stage 2-prep v2 — §D.1 Bulk-Shift Linearity Sanity\n',
        '**Gate:** CV(ΔT_bulk / α) < 5% at all 5 timesteps confirms the linear inversion '
        'model ΔT_bulk = C(t)·α is valid for Path D.\n',
        '| t (hr) | R_lowhyd | R_highhyd | R_highhyd_b010 | Mean C | Std | CV% | Gate |',
        '|---|---|---|---|---|---|---|---|',
    ]
    all_pass = True
    for t in T_DIAG_HR:
        r = lin[t]
        cv = r['cv_pct']
        gate = 'PASS' if cv < LINEARITY_CV_GATE else 'FAIL'
        if gate == 'FAIL':
            all_pass = False
        rs = r['ratios']
        pad = rs + [float('nan')] * (3 - len(rs))
        lines.append(
            f'| {t} | {pad[0]:.3f} | {pad[1]:.3f} | {pad[2]:.3f} | '
            f'{r["mean"]:.3f} | {r["std"]:.3f} | {cv:.1f}% | **{gate}** |'
        )
    lines.append('')
    lines.append(f'**Overall gate: {"PASS — linearity assumption valid" if all_pass else "FAIL — Path D linearity not confirmed"}**')
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def write_tref_csv(rows: list, path: str) -> None:
    fields = ['t_hr', 'best_T_ref_C', 'best_CV_pct', 'C_t_F_per_alpha']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_comparison_md(rows: list, tref_fits: list, path: str) -> None:
    lines = [
        '# Sprint 8 Stage 2-prep v2 — §D.4 Direct α(t) Comparison Table\n',
        '## Method',
        'CW α(effective) inverted from CW di=24 wi=0 bulk-shift signal via:',
        '  CW α(eff) = ΔT_bulk[d, t] / C(t)',
        'where C(t) is calibrated by fitting the T_ref that makes ΔT_bulk/α constant across',
        'all 3 datasets at each timestep (per-timestep best-fit T_ref from §D.2).\n',
        '## Per-timestep best-fit T_ref\n',
        '| t (hr) | Best T_ref (°C) | CV% |',
        '|---|---|---|',
    ]
    for r in tref_fits:
        lines.append(f'| {r["t_hr"]} | {r["best_T_ref_C"]} | {r["best_CV_pct"]:.2f}% |')

    lines += [
        '',
        '## Headline Table: Engine α vs CW α (effective)\n',
        '| Dataset | t (hr) | Engine α (T_ref=23°C) | Engine α (T_ref=21°C) | CW α (eff) | Δ at 23°C | Δ at 21°C | Δ at 23°C % | Δ at 21°C % |',
        '|---|---|---|---|---|---|---|---|---|',
    ]
    for r in rows:
        sign23 = '+' if r['delta_at_23C'] >= 0 else ''
        sign21 = '+' if r['delta_at_21C'] >= 0 else ''
        lines.append(
            f'| {r["dataset"]} | {r["t_hr"]} | {r["engine_alpha_23C"]:.4f} | '
            f'{r["engine_alpha_21C"]:.4f} | {r["cw_alpha_effective"]:.4f} | '
            f'{sign23}{r["delta_at_23C"]:.4f} | {sign21}{r["delta_at_21C"]:.4f} | '
            f'{sign23}{r["delta_at_23C_pct"]:.1f}% | {sign21}{r["delta_at_21C_pct"]:.1f}% |'
        )

    # Summary statistics
    d23 = [abs(r['delta_at_23C_pct']) for r in rows if not np.isnan(r['delta_at_23C_pct'])]
    d21 = [abs(r['delta_at_21C_pct']) for r in rows if not np.isnan(r['delta_at_21C_pct'])]
    lines += [
        '',
        '## Summary\n',
        '| Metric | Δ at T_ref=23°C | Δ at T_ref=21°C |',
        '|---|---|---|',
        f'| max|Δ| (%) | {max(d23):.1f}% | {max(d21):.1f}% |',
        f'| mean|Δ| (%) | {np.mean(d23):.1f}% | {np.mean(d21):.1f}% |',
        f'| RMS|Δ| (%) | {float(np.sqrt(np.mean([x**2 for x in d23]))):.1f}% | '
        f'{float(np.sqrt(np.mean([x**2 for x in d21]))):.1f}% |',
        '',
        '_Δ is signed: (CW α − Engine α). Positive = CW α > Engine α (CW hydrates faster)._',
    ]
    with open(path, 'w') as f:
        f.write('\n'.join(lines))


def write_comparison_csv(rows: list, path: str) -> None:
    fields = ['dataset', 't_hr', 'engine_alpha_23C', 'engine_alpha_21C',
              'cw_alpha_effective', 'delta_at_23C', 'delta_at_21C',
              'delta_at_23C_pct', 'delta_at_21C_pct']
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(FINDINGS, exist_ok=True)

    cache = pickle.load(open(os.path.join(STAGE1, 'cw_fields_cache.pkl'), 'rb'))
    cw_fields = cache['fields']

    print('\n=== Sprint 8 Stage 2-prep v2 — Path D: Bulk-Shift Inversion ===')
    print(f'T_placement = {T_PLACEMENT}°F  |  datasets: {list(DATASETS.keys())}')
    print(f'Timesteps: {T_DIAG_HR} hr\n')

    # Extract bulk shifts
    bulk = load_cw_bulk_shifts(cw_fields)

    print('--- CW di=24 wi=0 bulk shifts (ΔT from 73°F placement) ---')
    print(f'{"Dataset":<18} ', end='')
    for t in T_DIAG_HR:
        print(f'  t={t:3d}hr', end='')
    print()
    for ds in DATASETS:
        print(f'{ds:<18} ', end='')
        for t in T_DIAG_HR:
            print(f'  {bulk[ds][t]:+7.4f}', end='')
        print()

    # §D.1 Linearity check
    print('\n--- §D.1 Linearity sanity (CV of ΔT_bulk/α at T_ref=21°C) ---')
    lin = linearity_check(bulk)
    gate_overall = True
    for t in T_DIAG_HR:
        r = lin[t]
        cv = r['cv_pct']
        gate = 'PASS' if cv < LINEARITY_CV_GATE else 'FAIL'
        if gate == 'FAIL':
            gate_overall = False
        print(f't={t:3d}hr: ratios={[f"{x:.3f}" for x in r["ratios"]]}  '
              f'mean={r["mean"]:.3f}  CV={cv:.1f}%  [{gate}]')

    lin_path = os.path.join(FINDINGS, 'bulk_shift_linearity.md')
    write_linearity_md(lin, lin_path)
    print(f'\nWrote: {lin_path}')

    if not gate_overall:
        print('\nLINEARITY GATE FAIL — Path D inversion not valid. Stopping.')
        return

    print('\nLinearity gate PASS. Proceeding with per-timestep T_ref fit.')

    # §D.2 Per-timestep T_ref fit
    print('\n--- §D.2 Per-timestep best-fit T_ref ---')
    tref_fits = fit_tref_per_timestep(bulk)
    for r in tref_fits:
        print(f't={r["t_hr"]:3d}hr: best T_ref={r["best_T_ref_C"]:.2f}°C  '
              f'CV={r["best_CV_pct"]:.2f}%  C(t)={r["C_t_F_per_alpha"]:.4f} °F/α')

    tref_csv = os.path.join(FINDINGS, 't_ref_per_timestep.csv')
    write_tref_csv(tref_fits, tref_csv)
    print(f'\nWrote: {tref_csv}')

    # §D.3 + §D.4 Headline comparison table
    print('\n--- §D.4 Direct α(t) comparison table ---')
    table = build_comparison_table(bulk, tref_fits)

    header = f'{"Dataset":<18}  {"t":>4}  {"Eng α 23°C":>11}  {"Eng α 21°C":>11}  {"CW α eff":>10}  {"Δ@23°C":>8}  {"Δ@21°C":>8}  {"Δ@23%":>7}  {"Δ@21%":>7}'
    print(header)
    print('-' * len(header))
    for r in table:
        print(f'{r["dataset"]:<18}  {r["t_hr"]:>4}  {r["engine_alpha_23C"]:>11.4f}  '
              f'{r["engine_alpha_21C"]:>11.4f}  {r["cw_alpha_effective"]:>10.4f}  '
              f'{r["delta_at_23C"]:>+8.4f}  {r["delta_at_21C"]:>+8.4f}  '
              f'{r["delta_at_23C_pct"]:>+7.1f}%  {r["delta_at_21C_pct"]:>+7.1f}%')

    d21 = [abs(r['delta_at_21C_pct']) for r in table]
    d23 = [abs(r['delta_at_23C_pct']) for r in table]
    print(f'\nΔ@23°C: max={max(d23):.1f}%  mean={np.mean(d23):.1f}%  RMS={float(np.sqrt(np.mean([x**2 for x in d23]))):.1f}%')
    print(f'Δ@21°C: max={max(d21):.1f}%  mean={np.mean(d21):.1f}%  RMS={float(np.sqrt(np.mean([x**2 for x in d21]))):.1f}%')

    md_path = os.path.join(FINDINGS, 't_ref_validation.md')
    write_comparison_md(table, tref_fits, md_path)
    print(f'\nWrote: {md_path}')

    csv_path = os.path.join(FINDINGS, 't_ref_validation.csv')
    write_comparison_csv(table, csv_path)
    print(f'Wrote: {csv_path}')

    print('\n=== Path D complete ===')
    print('Headline: §D.4 delivers the 15-row α(t) comparison table.')
    print('CW α(effective) derived from CW di=24 bulk shift via per-timestep T_ref inversion.')


if __name__ == '__main__':
    main()
