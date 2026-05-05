"""Sprint 8 Stage 0 — detailed α(t) study over 0-168 hr.

Computes α(t) under two thermal regimes for all 15 mixes:
  1. Isothermal at placement temperature   (lower bound on how fast α grows)
  2. Adiabatic from placement temperature  (upper bound — self-heating accelerates Arrhenius)

Outputs (findings/):
  alpha_study.png   - 3-panel: α(t), k_c(t), and k vs α scatter coloured by time
  alpha_milestones.csv  - α and k at t = 0, 6, 12, 24, 48, 72, 120, 168 hr (all mixes)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv

from cw_scenario_loader import parse_cw_dat
from thermal_engine_2d import (
    K_UC_CALIBRATION_FACTOR_SPRINT7, BTU_HR_FT_F_TO_W_M_K,
    R_GAS, T_REF_K,
)

CW_EXPORTS_DIR = os.path.join(os.path.dirname(__file__), '../../../validation/cw_exports')
FINDINGS_DIR   = os.path.join(os.path.dirname(__file__), 'findings')
MILESTONES_HR  = [0, 6, 12, 24, 48, 72, 120, 168]
DT_HR          = 0.25   # integration step for adiabatic (hr)
T_HOURS_FINE   = np.arange(0, 168 + DT_HR, DT_HR)

# Approximate concrete density × Cp for adiabatic rise (J/m³·K)
RHO_CP_CONCRETE = 2350 * 900   # kg/m³ × J/(kg·K)

def arrhenius(T_K: np.ndarray, Ea: float) -> np.ndarray:
    return np.exp(-Ea / R_GAS * (1.0 / T_K - 1.0 / T_REF_K))

def alpha_from_te(te: np.ndarray, tau: float, beta: float, au: float) -> np.ndarray:
    te_s = np.maximum(te, 1e-6)
    return au * np.exp(-(tau / te_s) ** beta)

def k_vb(k_uc: float, alpha: np.ndarray) -> np.ndarray:
    return k_uc * (1.33 - 0.33 * alpha)

def isothermal_trajectory(T_pl_F, Ea, tau, beta, au):
    T_K = (T_pl_F - 32) * 5/9 + 273.15
    af  = float(arrhenius(np.array([T_K]), Ea)[0])
    te  = af * T_HOURS_FINE
    return alpha_from_te(te, tau, beta, au)

def adiabatic_trajectory(T_pl_F, Ea, tau, beta, au, Hu_J_kg, Cc_kg_m3):
    """Forward-Euler adiabatic integration."""
    T_K  = (T_pl_F - 32) * 5/9 + 273.15
    te   = 0.0
    alph = 0.0
    dt   = DT_HR
    alpha_arr = np.zeros_like(T_HOURS_FINE)
    T_arr     = np.full_like(T_HOURS_FINE, T_K)

    for i in range(1, len(T_HOURS_FINE)):
        af      = float(arrhenius(np.array([T_K]), Ea)[0])
        te     += af * dt
        alph    = float(alpha_from_te(np.array([te]), tau, beta, au)[0])
        alpha_arr[i] = alph
        # adiabatic temperature rise
        dT = Hu_J_kg * Cc_kg_m3 / RHO_CP_CONCRETE * (alph - alpha_arr[i-1])
        T_K += dT
        T_arr[i] = T_K

    return alpha_arr, T_arr

def Cc_from_mix(mix) -> float:
    """Total cementitious content kg/m³ from lb/yd³."""
    LB_YD3_TO_KG_M3 = 0.5933
    return ((mix.cement_type_I_II_lb_yd3 + mix.fly_ash_F_lb_yd3 + mix.fly_ash_C_lb_yd3
             + mix.ggbfs_lb_yd3 + mix.silica_fume_lb_yd3) * LB_YD3_TO_KG_M3)

def main():
    os.makedirs(FINDINGS_DIR, exist_ok=True)

    mixes = sorted(
        d for d in os.listdir(CW_EXPORTS_DIR)
        if d.startswith('MIX-') and os.path.isdir(os.path.join(CW_EXPORTS_DIR, d))
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    cmap_vals = cm.viridis(np.linspace(0, 1, len(mixes)))

    rows = []   # for milestone CSV

    sprint7_alpha = 0.036   # operating point
    sprint7_k_uc  = None    # set from first mix (same for all)

    for i, mix_dir in enumerate(mixes):
        mix, geom, constr, _ = parse_cw_dat(
            os.path.join(CW_EXPORTS_DIR, mix_dir, 'input.dat')
        )
        k_uc = (mix.thermal_conductivity_BTU_hr_ft_F
                * BTU_HR_FT_F_TO_W_M_K
                * K_UC_CALIBRATION_FACTOR_SPRINT7)
        if sprint7_k_uc is None:
            sprint7_k_uc = k_uc

        Cc = Cc_from_mix(mix)
        T_pl = constr.placement_temp_F

        alpha_iso  = isothermal_trajectory(T_pl, mix.activation_energy_J_mol,
                                           mix.tau_hrs, mix.beta, mix.alpha_u)
        alpha_adi, T_adi = adiabatic_trajectory(T_pl, mix.activation_energy_J_mol,
                                                 mix.tau_hrs, mix.beta, mix.alpha_u,
                                                 mix.Hu_J_kg, Cc)
        k_iso = k_vb(k_uc, alpha_iso)
        k_adi = k_vb(k_uc, alpha_adi)

        c = cmap_vals[i]
        label = mix_dir.replace('MIX-', '')

        axes[0].plot(T_HOURS_FINE, alpha_iso, color=c, lw=1.2, ls='-',
                     label=f'{label} iso')
        axes[0].plot(T_HOURS_FINE, alpha_adi, color=c, lw=1.2, ls='--')

        axes[1].plot(T_HOURS_FINE, k_iso, color=c, lw=1.2, ls='-')
        axes[1].plot(T_HOURS_FINE, k_adi, color=c, lw=1.2, ls='--')

        # panel 3: k vs α (parametric in time), colour = time
        sc = axes[2].scatter(alpha_adi, k_adi,
                             c=T_HOURS_FINE, cmap='plasma',
                             s=3, alpha=0.4, vmin=0, vmax=168)

        # milestones
        for t_m in MILESTONES_HR:
            idx = int(t_m / DT_HR)
            idx = min(idx, len(T_HOURS_FINE) - 1)
            rows.append({
                'mix': mix_dir, 'regime': 'isothermal',
                't_hr': t_m,
                'alpha': round(float(alpha_iso[idx]), 4),
                'k_W_mK': round(float(k_iso[idx]), 4),
                'T_K': round(float((T_pl - 32) * 5/9 + 273.15), 2),
            })
            rows.append({
                'mix': mix_dir, 'regime': 'adiabatic',
                't_hr': t_m,
                'alpha': round(float(alpha_adi[idx]), 4),
                'k_W_mK': round(float(k_adi[idx]), 4),
                'T_K': round(float(T_adi[idx]), 2),
            })

    # Sprint 7 operating-point reference lines
    k_sprint7 = float(k_vb(sprint7_k_uc, np.array([sprint7_alpha]))[0])
    for ax, val, lbl, yl in [
        (axes[0], sprint7_alpha,  f'Sprint 7 α={sprint7_alpha}', 'α'),
        (axes[1], k_sprint7,      f'Sprint 7 k={k_sprint7:.3f}', 'k_c (W/m·K)'),
    ]:
        ax.axhline(val, color='red', lw=1.2, ls=':', alpha=0.8)
        ax.text(5, val + (val * 0.01), lbl, color='red', fontsize=8, va='bottom')

    # Van Breugel reference line on k-vs-α panel
    alpha_ref = np.linspace(0, 1, 200)
    axes[2].plot(alpha_ref, k_vb(sprint7_k_uc, alpha_ref), 'r-', lw=1.5,
                 label=f'Van Breugel k_uc={sprint7_k_uc:.3f}')
    axes[2].axvline(sprint7_alpha, color='red', lw=1, ls=':', alpha=0.6)
    axes[2].text(sprint7_alpha + 0.01, k_vb(sprint7_k_uc, np.array([sprint7_alpha]))[0] + 0.02,
                 f'Sprint 7\nα={sprint7_alpha}', color='red', fontsize=7)
    plt.colorbar(sc, ax=axes[2], label='Time (hr)')

    # Panel formatting
    axes[0].set_xlabel('Time (hr)')
    axes[0].set_ylabel('α')
    axes[0].set_title('α(t): solid = isothermal, dashed = adiabatic')
    axes[0].set_xlim(0, 168); axes[0].set_ylim(0, None)
    axes[0].legend(title='Mix', ncol=3, fontsize=6, loc='upper left')

    axes[1].set_xlabel('Time (hr)')
    axes[1].set_ylabel('k_c  (W/m·K)')
    axes[1].set_title('k_c(t): solid = isothermal, dashed = adiabatic')
    axes[1].set_xlim(0, 168)

    axes[2].set_xlabel('α')
    axes[2].set_ylabel('k_c  (W/m·K)')
    axes[2].set_title('k_c vs α (adiabatic, all mixes, colour = time hr)')
    axes[2].legend(fontsize=8)

    plt.tight_layout()
    out_png = os.path.join(FINDINGS_DIR, 'alpha_study.png')
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Plot: {out_png}')

    # CSV
    out_csv = os.path.join(FINDINGS_DIR, 'alpha_milestones.csv')
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['mix','regime','t_hr','alpha','k_W_mK','T_K'])
        w.writeheader(); w.writerows(rows)
    print(f'CSV:  {out_csv}')

    # Terminal table — adiabatic only, MIX-01 as representative
    print('\n--- MIX-01 adiabatic milestones ---')
    print(f"{'t (hr)':>8} {'α':>8} {'k_c (W/m·K)':>13} {'T (°C)':>9}")
    print('-' * 44)
    mx01_rows = [r for r in rows if r['mix'] == 'MIX-01' and r['regime'] == 'adiabatic']
    for r in mx01_rows:
        print(f"{r['t_hr']:>8} {r['alpha']:>8.4f} {r['k_W_mK']:>13.4f} {r['T_K']-273.15:>9.2f}")

    print('\n--- All-mix α(168 hr) adiabatic summary ---')
    print(f"{'Mix':<10} {'α_iso_168':>10} {'k_iso_168':>11} {'α_adi_168':>10} {'k_adi_168':>11}")
    print('-' * 56)
    seen = set()
    for mix_dir in mixes:
        iso_row = next(r for r in rows
                       if r['mix']==mix_dir and r['regime']=='isothermal' and r['t_hr']==168)
        adi_row = next(r for r in rows
                       if r['mix']==mix_dir and r['regime']=='adiabatic' and r['t_hr']==168)
        print(f"{mix_dir:<10} {iso_row['alpha']:>10.4f} {iso_row['k_W_mK']:>11.4f} "
              f"{adi_row['alpha']:>10.4f} {adi_row['k_W_mK']:>11.4f}")

if __name__ == '__main__':
    main()
