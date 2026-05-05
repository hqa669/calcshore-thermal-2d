"""Sprint 8 Stage 1 §3.5+§3.6+§3.7 — three diagnostic figures."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle

from cw_scenario_loader import parse_cw_temp_output

BASE = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE, 'cw_data')
FINDINGS_DIR = os.path.join(BASE, 'findings')
FIGURES_DIR = os.path.join(BASE, 'figures')

DATASET_FOLDERS = {
    'lowhyd':       'thermal_a_cal_73_100_lowhyd',
    'highhyd':      'thermal_a_cal_73_100_highhyd',
    'highhyd_b010': 'thermal_a_cal_73_100_highhyd_beta010',
}

T_DIAG_HR = [24, 84, 168]
DI_MIN, DI_MAX = 4, 47

PAIRS = [
    ('A', 'highhyd',      'lowhyd',       'highhyd − lowhyd'),
    ('B', 'highhyd_b010', 'lowhyd',       'highhyd_b010 − lowhyd'),
    ('C', 'highhyd',      'highhyd_b010', 'highhyd − highhyd_b010'),
]

COLORS = {
    'lowhyd':       'tab:blue',
    'highhyd':      'tab:red',
    'highhyd_b010': 'tab:green',
}
LABELS = {
    'lowhyd':       'low-hyd  (α_u=0.10, τ=200, β=0.10)',
    'highhyd':      'high-hyd (α_u=0.70, τ=10,  β=0.85)',
    'highhyd_b010': 'high-hyd β=0.10 (α_u=0.70, τ=10,  β=0.10)',
}


def nearest_t_idx(time_hrs, t_target):
    return int(np.argmin(np.abs(time_hrs - t_target)))


def load_data():
    pkl_path = os.path.join(FINDINGS_DIR, 'cw_fields_cache.pkl')
    if os.path.exists(pkl_path):
        print(f'Loading cached fields from {pkl_path}')
        with open(pkl_path, 'rb') as f:
            return pickle.load(f)
    print('Cache not found — loading output.txt files directly...')
    fields = {}
    for label, folder in DATASET_FOLDERS.items():
        out_path = os.path.join(DATA_DIR, folder, 'output.txt')
        print(f'  {label}...')
        fields[label] = parse_cw_temp_output(out_path)
    nD = fields['lowhyd'].T_field_F.shape[1]
    nW = fields['lowhyd'].T_field_F.shape[2]
    mask = np.zeros((nD, nW), dtype=bool)
    mask[DI_MIN:DI_MAX + 1, :] = True
    return {'fields': fields, 'mask': mask, 'pairs': PAIRS}


def figure1_heatmaps(fields, pairs):
    nD = fields['lowhyd'].T_field_F.shape[1]
    nW = fields['lowhyd'].T_field_F.shape[2]

    # Pre-compute all ΔT fields
    all_dT = {}
    for pair_id, lab_a, lab_b, _ in pairs:
        all_dT[pair_id] = {}
        for t_hr in T_DIAG_HR:
            idx_a = nearest_t_idx(fields[lab_a].time_hrs, t_hr)
            idx_b = nearest_t_idx(fields[lab_b].time_hrs, t_hr)
            all_dT[pair_id][t_hr] = (fields[lab_a].T_field_F[idx_a]
                                     - fields[lab_b].T_field_F[idx_b])

    # Symmetric color scale: global max |ΔT| over masked region
    vmax = float(max(
        np.abs(all_dT[pid][t][DI_MIN:DI_MAX + 1, :]).max()
        for pid, *_ in pairs
        for t in T_DIAG_HR
    ))
    # Avoid zero vmax (would happen if all fields are identical)
    if vmax < 1e-6:
        vmax = 0.1

    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    im = None
    for row, (pair_id, lab_a, lab_b, pair_label) in enumerate(pairs):
        for col, t_hr in enumerate(T_DIAG_HR):
            ax = axes[row, col]
            dT = all_dT[pair_id][t_hr]
            im = ax.imshow(
                dT, aspect='auto', cmap='RdBu_r',
                vmin=-vmax, vmax=vmax,
                origin='upper',
                extent=[-0.5, nW - 0.5, nD - 0.5, -0.5],
            )
            ax.axhline(DI_MIN - 0.5, color='gray', lw=0.8, ls='--', alpha=0.6)
            ax.axhline(DI_MAX + 0.5, color='gray', lw=0.8, ls='--', alpha=0.6)
            if row == 0:
                ax.set_title(f't = {t_hr} hr', fontsize=11)
            if col == 0:
                ax.set_ylabel(f'Pair {pair_id}\n{pair_label}', fontsize=8)
            if row == 2:
                ax.set_xlabel('wi (0=CL → 12=form face)', fontsize=8)

    fig.subplots_adjust(right=0.86, hspace=0.35, wspace=0.3)
    cbar_ax = fig.add_axes([0.88, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label('ΔT (°F)', fontsize=11)

    fig.suptitle(
        f'CW Pairwise ΔT Fields (°F) — symmetric scale ±{vmax:.3f}°F\n'
        'di=0 top → di=48 bottom; dashed = mask boundary (di=4/47)',
        fontsize=10,
    )
    out = os.path.join(FIGURES_DIR, 'stage1_cw_pairwise_dT_fields.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Figure 1 saved: {out}  (scale: ±{vmax:.4f}°F)')
    return vmax


def figure2_midprofile(fields):
    di_mid = 24
    t_target = 168

    # widths_m stored decreasing: widths_m[0]=6.1 (CL), widths_m[-1]=0 (form face)
    # Distance from CL = widths_m[0] - widths_m
    widths_m = fields['lowhyd'].widths_m
    dist_from_CL = widths_m[0] - widths_m  # 0 at CL → ~6.1 at form face

    fig, ax = plt.subplots(figsize=(8, 5))
    for label in ['lowhyd', 'highhyd', 'highhyd_b010']:
        series = fields[label]
        idx = nearest_t_idx(series.time_hrs, t_target)
        T_profile = series.T_field_F[idx, di_mid, :]
        ax.plot(dist_from_CL, T_profile, color=COLORS[label], lw=2, label=LABELS[label])

    ax.set_xlabel('Lateral distance from CL (m)', fontsize=11)
    ax.set_ylabel('Temperature (°F)', fontsize=11)
    ax.set_title(f'CW mid-depth (di=24) lateral profile at t={t_target} hr\n'
                 'Placement=73°F, Soil=100°F', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    out = os.path.join(FIGURES_DIR, 'stage1_cw_midprofile_t168.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Figure 2 saved: {out}')


def figure3_cl_evolution(fields):
    di_mid = 24
    wi_cl  = 0   # centerline: wi=0 in CW loader convention

    fig, ax = plt.subplots(figsize=(9, 5))
    for label in ['lowhyd', 'highhyd', 'highhyd_b010']:
        series = fields[label]
        T_cl = series.T_field_F[:, di_mid, wi_cl]
        ax.plot(series.time_hrs, T_cl, color=COLORS[label], lw=2, label=LABELS[label])

    ax.set_xlabel('Time (hr)', fontsize=11)
    ax.set_ylabel('Temperature (°F)', fontsize=11)
    ax.set_title('CW centerline T(t) — di=24, wi=0\n'
                 'Placement=73°F, Soil=100°F — critical k(α) signal', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xticks(range(0, 169, 24))
    ax.set_xlim(0, 168)

    out = os.path.join(FIGURES_DIR, 'stage1_cw_CL_T_evolution.png')
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Figure 3 saved: {out}')


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    data = load_data()
    fields = data['fields']
    pairs  = data['pairs']

    vmax = figure1_heatmaps(fields, pairs)
    figure2_midprofile(fields)
    figure3_cl_evolution(fields)
    print(f'\nAll three figures written to {FIGURES_DIR}/')


if __name__ == '__main__':
    main()
