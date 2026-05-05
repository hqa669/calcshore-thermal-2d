"""Sprint 8 Stage 1 §3.2 — input.dat consistency check across 3 CW datasets."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

from cw_scenario_loader import parse_cw_dat

BASE = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE, 'cw_data')
FINDINGS_DIR = os.path.join(BASE, 'findings')

DATASETS = [
    ('lowhyd',       'thermal_a_cal_73_100_lowhyd'),
    ('highhyd',      'thermal_a_cal_73_100_highhyd'),
    ('highhyd_b010', 'thermal_a_cal_73_100_highhyd_beta010'),
]

# Only these are expected to differ
EXPECTED_DIFF_KEYS = {'alpha_u', 'tau_hrs', 'beta', 'Hu_J_kg'}


def extract_params(mix, geom, constr):
    return {
        # Mix design
        'cement_lb_yd3':      mix.cement_type_I_II_lb_yd3,
        'water_lb_yd3':       mix.water_lb_yd3,
        'coarse_agg_lb_yd3':  mix.coarse_agg_lb_yd3,
        'fine_agg_lb_yd3':    mix.fine_agg_lb_yd3,
        'fly_ash_F_lb_yd3':   mix.fly_ash_F_lb_yd3,
        'fly_ash_C_lb_yd3':   mix.fly_ash_C_lb_yd3,
        'ggbfs_lb_yd3':       mix.ggbfs_lb_yd3,
        'silica_fume_lb_yd3': mix.silica_fume_lb_yd3,
        'air_content_pct':    mix.air_content_pct,
        'cement_type':        mix.cement_type,
        'coarse_agg_type':    mix.coarse_agg_type,
        'fine_agg_type':      mix.fine_agg_type,
        # Hydration kinetics
        'activation_energy_J_mol': mix.activation_energy_J_mol,
        'tau_hrs':   mix.tau_hrs,
        'beta':      mix.beta,
        'alpha_u':   mix.alpha_u,
        'Hu_J_kg':   mix.Hu_J_kg,
        # Derived thermal
        'k_BTU_hr_ft_F':   mix.thermal_conductivity_BTU_hr_ft_F,
        'agg_Cp_BTU_lb_F': mix.aggregate_Cp_BTU_lb_F,
        'CTE_microstrain_F': mix.CTE_microstrain_F,
        # Geometry
        'depth_ft':  geom.depth_ft,
        'width_ft':  geom.width_ft,
        'length_ft': geom.length_ft,
        # Construction
        'placement_temp_F':  constr.placement_temp_F,
        'soil_temp_F':       constr.soil_temp_F,
        'form_type':         constr.form_type,
        'form_removal_hrs':  constr.form_removal_hrs,
        'blanket_R_value':   constr.blanket_R_value,
        'footing_subbase':   constr.footing_subbase,
    }


def main():
    os.makedirs(FINDINGS_DIR, exist_ok=True)

    all_params = {}
    for label, folder in DATASETS:
        dat_path = os.path.join(DATA_DIR, folder, 'input.dat')
        mix, geom, constr, _ = parse_cw_dat(dat_path)
        all_params[label] = extract_params(mix, geom, constr)
        print(f'{label}: tau_hrs={mix.tau_hrs}, beta={mix.beta}, alpha_u={mix.alpha_u}, Hu_J_kg={mix.Hu_J_kg}')

    labels = [d[0] for d in DATASETS]
    all_keys = list(next(iter(all_params.values())).keys())

    differing = [k for k in all_keys if len({str(all_params[lab][k]) for lab in labels}) > 1]
    unexpected = [k for k in differing if k not in EXPECTED_DIFF_KEYS]

    # Build markdown
    lines = ['# Sprint 8 Stage 1 — input.dat Consistency Check\n']
    lines.append(f'Datasets compared: {", ".join(labels)}\n')
    lines.append('| Parameter | lowhyd | highhyd | highhyd_b010 | Status |')
    lines.append('|---|---|---|---|---|')
    for key in all_keys:
        vals = [str(all_params[lab][key]) for lab in labels]
        differs = len(set(vals)) > 1
        if key in EXPECTED_DIFF_KEYS:
            status = 'EXPECTED DIFF' if differs else 'same (unexpected)'
        else:
            status = '**UNEXPECTED DIFF — INVESTIGATE**' if differs else 'same'
        lines.append(f'| {key} | {" | ".join(vals)} | {status} |')

    lines.append('\n## Verdict\n')
    if unexpected:
        lines.append(f'**FAIL** — unexpected differences in: {unexpected}')
        lines.append('These differences may invalidate the k(α) diagnostic.')
    else:
        lines.append('**PASS** — only hydration kinetics parameters differ (tau_hrs, beta, alpha_u, Hu_J_kg). '
                     'Diagnostic is valid.')

    md = '\n'.join(lines)
    out_path = os.path.join(FINDINGS_DIR, 'input_consistency.md')
    with open(out_path, 'w') as f:
        f.write(md)
    print(md)
    print(f'\nWrote: {out_path}')

    if unexpected:
        raise SystemExit(f'FAIL: unexpected parameter diffs in {unexpected}')
    print('\nconsistency PASS')


if __name__ == '__main__':
    main()
