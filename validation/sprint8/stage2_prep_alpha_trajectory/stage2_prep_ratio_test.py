"""Sprint 8 Stage 2-prep §3.2+§3.3 — α(t) trajectory consistency via pairwise ΔT ratio test.

CW does not output α(t). Indirect test:
  In a linear thermal diffusion model with Van Breugel k(α),
  ΔT ∝ Δk_c = k_uc·0.33·Δα  (same BCs, geometry, k_uc for all pairs).
  Therefore: ΔT_A/ΔT_B ≈ Δα_A/Δα_B at each timestep.

If engine and CW compute matching α(t), this ratio should equal 1.0 (within linearization
error). Systematic deviation across all timesteps and stats indicates an α trajectory
discrepancy.

Uses Stage 1 pairwise_dT_stats.csv and engine_alpha_reference.csv (no new runs needed).
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import csv

BASE = os.path.dirname(__file__)
STAGE1 = os.path.join(BASE, '../stage1_cw_k_alpha_empirical/findings')
FINDINGS = os.path.join(BASE, 'findings')

T_DIAG_HR = [24, 84, 168]
STATS = ['max_abs_dT_F', 'mean_abs_dT_F', 'rms_dT_F']
STAT_LABELS = {'max_abs_dT_F': 'max|ΔT|', 'mean_abs_dT_F': 'mean|ΔT|', 'rms_dT_F': 'RMS ΔT'}


def load_dT_stats():
    rows = {}
    with open(os.path.join(STAGE1, 'pairwise_dT_stats.csv')) as f:
        for r in csv.DictReader(f):
            key = (r['pair'], int(float(r['t_hr'])))
            rows[key] = {s: float(r[s]) for s in STATS}
    return rows


def load_alpha_ref():
    alpha = {}
    with open(os.path.join(STAGE1, 'engine_alpha_reference.csv')) as f:
        for r in csv.DictReader(f):
            if r['regime'] == 'isothermal':
                ds = r['dataset']
                alpha[ds] = {int(t): float(r[f'alpha_t{t}hr']) for t in T_DIAG_HR}
    return alpha


def delta_alpha(alpha, ds_a, ds_b, t):
    return alpha[ds_a][t] - alpha[ds_b][t]


def main():
    os.makedirs(FINDINGS, exist_ok=True)

    dT = load_dT_stats()
    alpha = load_alpha_ref()

    # Ratio test definitions:
    #   Test 1 (main): Pair A / Pair B — highhyd vs lowhyd divided by highhyd_b010 vs lowhyd
    #   Test 2 (cross-check): Pair A / Pair C — highhyd vs lowhyd divided by highhyd vs highhyd_b010
    ratio_tests = [
        ('A/B', 'A', 'B', 'highhyd', 'lowhyd', 'highhyd_b010', 'lowhyd'),
        ('A/C', 'A', 'C', 'highhyd', 'lowhyd', 'highhyd',      'highhyd_b010'),
    ]

    all_discrep = []
    report_rows = []

    for test_name, pairA, pairB, da_num, db_num, da_den, db_den in ratio_tests:
        for t in T_DIAG_HR:
            d_alpha_num = delta_alpha(alpha, da_num, db_num, t)  # Δα numerator
            d_alpha_den = delta_alpha(alpha, da_den, db_den, t)  # Δα denominator
            expected_ratio = d_alpha_num / d_alpha_den if d_alpha_den != 0 else float('nan')

            for stat in STATS:
                val_num = dT[(pairA, t)][stat]
                val_den = dT[(pairB, t)][stat]
                observed_ratio = val_num / val_den if val_den != 0 else float('nan')
                discrepancy_pct = (observed_ratio / expected_ratio - 1.0) * 100.0

                row = {
                    'test': test_name,
                    'stat': STAT_LABELS[stat],
                    't_hr': t,
                    'observed_ratio': round(observed_ratio, 4),
                    'expected_ratio': round(expected_ratio, 4),
                    'discrepancy_pct': round(discrepancy_pct, 2),
                    'Dalpha_num': round(d_alpha_num, 4),
                    'Dalpha_den': round(d_alpha_den, 4),
                }
                report_rows.append(row)
                all_discrep.append(abs(discrepancy_pct))

    max_discrep = max(all_discrep)
    if max_discrep < 1.0:
        verdict = 'PASS (<1%) — trajectories agree, Sprint 8 Schindler-param spec proceeds without modification'
    elif max_discrep <= 5.0:
        verdict = f'CAUTION ({max_discrep:.1f}% in 1–5% band) — noticeable but probably tolerable; quantify contamination'
    else:
        verdict = f'FAIL (>{max_discrep:.1f}%) — significant discrepancy; Sprint 8 dataset design needs compensation'

    # Print to stdout
    print(f'\n=== Stage 2-prep Ratio Test ===')
    print(f'{"Test":<6} {"Stat":<12} {"t(hr)":<7} {"Obs ratio":>10} {"Exp ratio":>10} {"Discrepancy%":>13}')
    print('-' * 62)
    for r in report_rows:
        print(f"{r['test']:<6} {r['stat']:<12} {r['t_hr']:<7} "
              f"{r['observed_ratio']:>10.4f} {r['expected_ratio']:>10.4f} "
              f"{r['discrepancy_pct']:>+12.2f}%")

    print(f'\nMax |discrepancy|: {max_discrep:.2f}%')
    print(f'Verdict: {verdict}')

    # Write CSV
    csv_path = os.path.join(FINDINGS, 'ratio_test_results.csv')
    fieldnames = ['test', 'stat', 't_hr', 'observed_ratio', 'expected_ratio',
                  'discrepancy_pct', 'Dalpha_num', 'Dalpha_den']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader(); w.writerows(report_rows)
    print(f'\nWrote: {csv_path}')

    # Write Markdown
    md_lines = [
        '# Sprint 8 Stage 2-prep — §3.2 + §3.3 Ratio Test Results\n',
        '## Method\n',
        'In a linear thermal model with Van Breugel k(α), ΔT ∝ Δk_c = k_uc·0.33·Δα.',
        'Therefore ΔT_A/ΔT_B ≈ Δα_A/Δα_B at each timestep (same BCs, geometry, k_uc).',
        'Discrepancy = (observed_ratio/expected_ratio − 1) × 100%.',
        'Computed for both ratio tests (A/B and A/C) and all three ΔT statistics.\n',
        '## Results\n',
        '| Test | Stat | t (hr) | Observed ratio | Expected ratio | Discrepancy% |',
        '|---|---|---|---|---|---|',
    ]
    for r in report_rows:
        md_lines.append(
            f"| {r['test']} | {r['stat']} | {r['t_hr']} "
            f"| {r['observed_ratio']:.4f} | {r['expected_ratio']:.4f} "
            f"| {r['discrepancy_pct']:+.2f}% |"
        )
    md_lines.append(f'\n## §3.3 Magnitude Verdict\n')
    md_lines.append(f'**Max |discrepancy| across all stats, timesteps, and ratio tests: {max_discrep:.2f}%**\n')
    md_lines.append(f'**{verdict}**\n')
    md_lines.append('Thresholds: <1% = PASS, 1–5% = CAUTION, >5% = FAIL\n')

    md_path = os.path.join(FINDINGS, 'ratio_test_results.md')
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    print(f'Wrote: {md_path}')


if __name__ == '__main__':
    main()
