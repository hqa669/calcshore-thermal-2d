"""Stage 0 diagnostic: k_c(alpha(t)) trajectory for all 15 CW mixes, 0-168 hr.

Each mix is run isothermally at its placement temperature to produce alpha(te(t))
via the Schindler formula. k_c(t) = k_uc * (1.33 - 0.33 * alpha(t)) is then plotted.

k_uc used = thermal_conductivity_BTU_hr_ft_F * BTU_HR_FT_F_TO_W_M_K * 0.96 (Sprint 7 calibration).

Usage:
    python stage0_k_over_time.py

Outputs (written to findings/):
    k_over_time.png     - k_c(t) for all 15 mixes, two panels
    k_alpha_range.csv   - table: mix, k_uc, k_t0, k_t168, delta_k_pct, alpha_t168
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../..'))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import csv

from cw_scenario_loader import parse_cw_dat, CW_DAT_INDEX
from thermal_engine_2d import (
    K_UC_CALIBRATION_FACTOR_SPRINT7,
    BTU_HR_FT_F_TO_W_M_K,
    R_GAS,
    T_REF_K,
)

# --- config ---
CW_EXPORTS_DIR = os.path.join(os.path.dirname(__file__), '../../../validation/cw_exports')
FINDINGS_DIR   = os.path.join(os.path.dirname(__file__), 'findings')
T_HOURS        = np.linspace(0.01, 168, 2000)   # avoid te=0 at t=0

def arrhenius_factor(T_placement_F: float, Ea_J_mol: float) -> float:
    T_K = (T_placement_F - 32) * 5/9 + 273.15
    return float(np.exp(-Ea_J_mol / R_GAS * (1.0 / T_K - 1.0 / T_REF_K)))

def alpha_isothermal(t_hr: np.ndarray, af: float, tau: float, beta: float, au: float) -> np.ndarray:
    te = af * t_hr
    te = np.maximum(te, 1e-6)
    return au * np.exp(-(tau / te) ** beta)

def k_van_breugel(k_uc: float, alpha: np.ndarray) -> np.ndarray:
    return k_uc * (1.33 - 0.33 * alpha)

def load_mix(mix_dir: str):
    dat_path = os.path.join(CW_EXPORTS_DIR, mix_dir, 'input.dat')
    mix, geom, constr, _ = parse_cw_dat(dat_path)
    k_uc = (mix.thermal_conductivity_BTU_hr_ft_F
            * BTU_HR_FT_F_TO_W_M_K
            * K_UC_CALIBRATION_FACTOR_SPRINT7)
    return mix, constr, k_uc

def main():
    os.makedirs(FINDINGS_DIR, exist_ok=True)

    mixes = sorted(
        d for d in os.listdir(CW_EXPORTS_DIR)
        if d.startswith('MIX-') and os.path.isdir(os.path.join(CW_EXPORTS_DIR, d))
    )

    results = []
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    cmap = plt.cm.viridis(np.linspace(0, 1, len(mixes)))

    for i, mix_dir in enumerate(mixes):
        mix, constr, k_uc = load_mix(mix_dir)
        af = arrhenius_factor(constr.placement_temp_F, mix.activation_energy_J_mol)
        alpha = alpha_isothermal(T_HOURS, af, mix.tau_hrs, mix.beta, mix.alpha_u)
        k_t   = k_van_breugel(k_uc, alpha)

        k_t0   = k_uc * 1.33                          # alpha -> 0 at t=0
        k_t168 = k_t[-1]
        alpha_t168 = alpha[-1]
        delta_pct = 100 * (k_t0 - k_t168) / k_t0

        results.append({
            'mix':         mix_dir,
            'k_uc_W_mK':  round(k_uc, 4),
            'k_t0_W_mK':  round(k_t0, 4),
            'k_t168_W_mK':round(k_t168, 4),
            'delta_k_pct':round(delta_pct, 2),
            'alpha_t168': round(alpha_t168, 4),
            'T_pl_F':     round(constr.placement_temp_F, 1),
            'tau_hr':     round(mix.tau_hrs, 2),
            'beta':       round(mix.beta, 3),
            'alpha_u':    round(mix.alpha_u, 4),
        })

        label = mix_dir.replace('MIX-', '')
        axes[0].plot(T_HOURS, k_t, color=cmap[i], lw=1.2, label=label)
        axes[1].plot(T_HOURS, alpha, color=cmap[i], lw=1.2, label=label)

    # --- panel 0: k_c(t) ---
    axes[0].set_xlabel('Time (hr)')
    axes[0].set_ylabel('k_c  (W/m·K)')
    axes[0].set_title('Concrete thermal conductivity k_c(α(t))\n'
                      'Isothermal at placement temp, Van Breugel formula, k_uc × 0.96')
    axes[0].legend(title='Mix', ncol=3, fontsize=7, loc='upper right')
    axes[0].set_xlim(0, 168)

    k_all_t0   = np.array([r['k_t0_W_mK']  for r in results])
    k_all_t168 = np.array([r['k_t168_W_mK'] for r in results])
    axes[0].axhline(k_all_t0.min(),   color='grey', lw=0.6, ls='--')
    axes[0].axhline(k_all_t168.min(), color='grey', lw=0.6, ls=':')
    axes[0].text(170, k_all_t0.min(),   f'k_t0 min={k_all_t0.min():.3f}', va='center', fontsize=7)
    axes[0].text(170, k_all_t168.min(), f'k_t168 min={k_all_t168.min():.3f}', va='center', fontsize=7)

    # --- panel 1: alpha(t) ---
    axes[1].set_xlabel('Time (hr)')
    axes[1].set_ylabel('α (degree of hydration)')
    axes[1].set_title('Degree of hydration α(t)\nIsothermal at placement temp, Schindler formula')
    axes[1].legend(title='Mix', ncol=3, fontsize=7, loc='lower right')
    axes[1].set_xlim(0, 168)

    plt.tight_layout()
    out_png = os.path.join(FINDINGS_DIR, 'k_over_time.png')
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Plot saved: {out_png}')

    # --- CSV summary ---
    out_csv = os.path.join(FINDINGS_DIR, 'k_alpha_range.csv')
    fields = ['mix', 'k_uc_W_mK', 'k_t0_W_mK', 'k_t168_W_mK', 'delta_k_pct',
              'alpha_t168', 'T_pl_F', 'tau_hr', 'beta', 'alpha_u']
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(results)
    print(f'CSV saved:  {out_csv}')

    # --- terminal summary ---
    print()
    print(f"{'Mix':<10} {'k_uc':>8} {'k(t=0)':>8} {'k(168hr)':>9} {'Δk%':>7} {'α(168hr)':>9} {'T_pl°F':>7}")
    print('-' * 70)
    for r in results:
        print(f"{r['mix']:<10} {r['k_uc_W_mK']:>8.4f} {r['k_t0_W_mK']:>8.4f} "
              f"{r['k_t168_W_mK']:>9.4f} {r['delta_k_pct']:>7.2f} "
              f"{r['alpha_t168']:>9.4f} {r['T_pl_F']:>7.1f}")

    delta_vals = np.array([r['delta_k_pct'] for r in results])
    print(f"\nk drop over 168 hr:  min {delta_vals.min():.1f}%  max {delta_vals.max():.1f}%  "
          f"mean {delta_vals.mean():.1f}%")

if __name__ == '__main__':
    main()
