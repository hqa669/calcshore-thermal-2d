"""Sprint 8 Stage 1 §3.3+§3.4 — load CW T-fields and compute pairwise ΔT statistics."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import numpy as np
import csv
import pickle

from cw_scenario_loader import parse_cw_temp_output

BASE = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE, 'cw_data')
FINDINGS_DIR = os.path.join(BASE, 'findings')

DATASET_FOLDERS = {
    'lowhyd':       'thermal_a_cal_73_100_lowhyd',
    'highhyd':      'thermal_a_cal_73_100_highhyd',
    'highhyd_b010': 'thermal_a_cal_73_100_highhyd_beta010',
}

T_DIAG_HR = [24, 84, 168]

# Mask per brief: di in [4..47] inclusive, all wi
DI_MIN, DI_MAX = 4, 47

PAIRS = [
    ('A', 'highhyd',      'lowhyd',       'highhyd − lowhyd'),
    ('B', 'highhyd_b010', 'lowhyd',       'highhyd_b010 − lowhyd'),
    ('C', 'highhyd',      'highhyd_b010', 'highhyd − highhyd_b010'),
]


def nearest_t_idx(time_hrs, t_target):
    return int(np.argmin(np.abs(time_hrs - t_target)))


def main():
    os.makedirs(FINDINGS_DIR, exist_ok=True)

    # Load all three datasets
    fields = {}
    for label, folder in DATASET_FOLDERS.items():
        out_path = os.path.join(DATA_DIR, folder, 'output.txt')
        print(f'Loading {label} from {out_path}...')
        series = parse_cw_temp_output(out_path)
        fields[label] = series
        nD, nW = series.T_field_F.shape[1], series.T_field_F.shape[2]
        print(f'  shape={series.T_field_F.shape}, t=[{series.time_hrs[0]:.1f}, {series.time_hrs[-1]:.1f}] hr, nD={nD}, nW={nW}')

    nD = fields['lowhyd'].T_field_F.shape[1]
    nW = fields['lowhyd'].T_field_F.shape[2]
    print(f'\nGrid: nD={nD}, nW={nW}')

    # Build mask
    mask = np.zeros((nD, nW), dtype=bool)
    mask[DI_MIN:DI_MAX + 1, :] = True
    n_cells = int(mask.sum())
    print(f'Mask: di=[{DI_MIN},{DI_MAX}], all wi → {n_cells} cells (expected {(DI_MAX-DI_MIN+1)*nW})')

    rows = []
    for pair_id, lab_a, lab_b, pair_label in PAIRS:
        for t_hr in T_DIAG_HR:
            idx_a = nearest_t_idx(fields[lab_a].time_hrs, t_hr)
            idx_b = nearest_t_idx(fields[lab_b].time_hrs, t_hr)
            t_a = float(fields[lab_a].time_hrs[idx_a])
            t_b = float(fields[lab_b].time_hrs[idx_b])

            T_a = fields[lab_a].T_field_F[idx_a]   # (nD, nW)
            T_b = fields[lab_b].T_field_F[idx_b]   # (nD, nW)
            dT  = T_a - T_b                         # (nD, nW)

            dT_masked = np.where(mask, dT, np.nan)
            abs_dT    = np.abs(dT_masked)

            max_abs   = float(np.nanmax(abs_dT))
            mean_abs  = float(np.nanmean(abs_dT))
            rms       = float(np.sqrt(np.nanmean(dT_masked ** 2)))

            flat_idx  = int(np.nanargmax(abs_dT))
            di_max, wi_max = divmod(flat_idx, nW)

            rows.append({
                'pair':         pair_id,
                'pair_label':   pair_label,
                't_hr':         t_hr,
                't_actual_a':   round(t_a, 2),
                't_actual_b':   round(t_b, 2),
                'max_abs_dT_F': round(max_abs, 4),
                'mean_abs_dT_F': round(mean_abs, 4),
                'rms_dT_F':     round(rms, 4),
                'loc_di':       di_max,
                'loc_wi':       wi_max,
            })

    # Write CSV
    csv_path = os.path.join(FINDINGS_DIR, 'pairwise_dT_stats.csv')
    fieldnames = ['pair', 'pair_label', 't_hr', 't_actual_a', 't_actual_b',
                  'max_abs_dT_F', 'mean_abs_dT_F', 'rms_dT_F', 'loc_di', 'loc_wi']
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f'\nWrote: {csv_path}')

    # Write Markdown
    md_lines = [
        '# Sprint 8 Stage 1 — Pairwise ΔT Statistics\n',
        f'Mask: di={DI_MIN}..{DI_MAX} (inclusive), all wi={0}..{nW-1}. Temperature in °F.\n',
        '| Pair | t (hr) | max\\|ΔT\\| (°F) | mean\\|ΔT\\| (°F) | RMS ΔT (°F) | location (di,wi) |',
        '|---|---|---|---|---|---|',
    ]
    for r in rows:
        md_lines.append(
            f"| {r['pair']}: {r['pair_label']} | {r['t_hr']} "
            f"| {r['max_abs_dT_F']:.4f} | {r['mean_abs_dT_F']:.4f} "
            f"| {r['rms_dT_F']:.4f} | ({r['loc_di']},{r['loc_wi']}) |"
        )
    md_path = os.path.join(FINDINGS_DIR, 'pairwise_dT_stats.md')
    with open(md_path, 'w') as f:
        f.write('\n'.join(md_lines))
    print(f'Wrote: {md_path}')

    print('\n' + '\n'.join(md_lines))

    # Cache fields for figures script
    pkl_path = os.path.join(FINDINGS_DIR, 'cw_fields_cache.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump({'fields': fields, 'mask': mask, 'pairs': PAIRS}, f)
    print(f'\nCached fields to {pkl_path} (used by stage1_figures.py)')


if __name__ == '__main__':
    main()
