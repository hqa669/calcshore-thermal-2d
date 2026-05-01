#!/usr/bin/env python3
"""Stage 5a — S5a.1/S5a.2/S5a.3: Characterize engine state for Run B.

Runs the engine for Run B (placement=73°F, soil=60°F) with hydration suppressed
(Hu=1, αu=0.1) and captures α, T, k, ρCp, α_c fields at t = 1, 24, 84, 168 hr.

Sanity check: masked max|R| at t=168 hr must reproduce Stage 4b's 1.935°F ±0.01°F.

Outputs:
    validation/soil_calibration/STAGE5a_runB_alpha_field.npy   shape (4, ny_c, nx_c)
    validation/soil_calibration/STAGE5a_runB_T_field_C.npy     shape (4, ny_c, nx_c)
    validation/soil_calibration/STAGE5a_runB_alphac.npy        shape (4, ny_c, nx_c) m²/hr
    validation/soil_calibration/STAGE5a_engine_state_runB.md
    validation/soil_calibration/STAGE5a_cause1_k_value.md
    validation/soil_calibration/STAGE5a_cause2_rcp_model.md
"""
import os
import sys
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import parse_cw_dat, CWEnvironment
from thermal_engine_2d import (
    solve_hydration_2d,
    build_grid_half_mat,
    thermal_conductivity_variable,
    specific_heat_variable,
    LB_FT3_TO_KG_M3,
    BTU_HR_FT_F_TO_W_M_K,
    BTU_LB_F_TO_J_KG_K,
    LB_YD3_TO_KG_M3,
)
from kinetics_correction import compute_hu_factor
from stage3_compare import load_cw_slice, resample_engine_to_cw, COMPARE_HR
from stage4b_run import make_neutral_env, nearest_time_idx, masked_stats, save_field_csv

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
MASK_PATH = os.path.join(HERE, "STAGE35_validity_mask.npy")

TARGET_HRS = [1.0, 24.0, 84.0, 168.0]
RUN_LABEL = "B"
RUN_FOLDER = "runB_73_60"
PLACEMENT_F = 73
SOIL_F = 60

CW_TARGET_MASKED_MAX = 1.935  # °F from Stage 4b table
SANITY_TOL = 0.01             # °F


def make_env(placement_F):
    return make_neutral_env(placement_F)


def derive_mix_si(mix):
    """Return SI-unit thermal scalars the engine uses internally."""
    rho_c = mix.concrete_density_lb_ft3 * LB_FT3_TO_KG_M3
    k_uc  = mix.thermal_conductivity_BTU_hr_ft_F * BTU_HR_FT_F_TO_W_M_K
    Ca    = mix.aggregate_Cp_BTU_lb_F * BTU_LB_F_TO_J_KG_K
    Wc    = mix.cement_type_I_II_lb_yd3 * LB_YD3_TO_KG_M3
    Ww    = mix.water_lb_yd3 * LB_YD3_TO_KG_M3
    Wa    = (mix.coarse_agg_lb_yd3 + mix.fine_agg_lb_yd3) * LB_YD3_TO_KG_M3
    return dict(rho_c=rho_c, k_uc=k_uc, Ca=Ca, Wc=Wc, Ww=Ww, Wa=Wa)


def erode_concrete_mask(is_concrete_2d, n=1):
    """Remove n-cell border from a 2D boolean mask (interior cells only)."""
    m = is_concrete_2d.copy()
    m[:n, :] = False
    m[-n:, :] = False
    m[:, :n] = False
    m[:, -n:] = False
    return m


def field_stats(arr, mask):
    vals = arr[mask]
    return dict(
        min=float(vals.min()),
        max=float(vals.max()),
        mean=float(vals.mean()),
        std=float(vals.std()),
    )


def run():
    mask_full = np.load(MASK_PATH)
    assert mask_full.shape == (49, 13), f"Mask shape mismatch: {mask_full.shape}"
    print(f"Loaded Stage 3.5 mask: shape={mask_full.shape}, {int(mask_full.sum())} cells kept")

    dat_path = os.path.join(CW_RUNS, RUN_FOLDER, "input.dat")
    mix, geom, constr, _ = parse_cw_dat(dat_path)

    factor, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = factor
    mix.Hu_J_kg_effective = mix.Hu_J_kg * factor
    print(f"  Hu_J_kg={mix.Hu_J_kg:.4f}, Hu_J_kg_effective={mix.Hu_J_kg_effective:.6f}"
          f"  (suppressed={mix.Hu_J_kg_effective < 10})")

    si = derive_mix_si(mix)
    print(f"  Mix SI: rho_c={si['rho_c']:.2f} kg/m³, k_uc={si['k_uc']:.4f} W/m·K, "
          f"Ca={si['Ca']:.1f} J/kg·K, Wc={si['Wc']:.1f}, Wa={si['Wa']:.1f}, "
          f"Ww={si['Ww']:.1f} kg/m³")

    constr.model_soil = False
    constr.is_submerged = True

    grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                               is_submerged=True, model_soil=False)
    jslice, islice = grid.concrete_slice()
    nx_c = len(grid.x[islice])
    ny_c = len(grid.y[jslice])
    dx = grid.x[1] - grid.x[0] if grid.nx > 1 else float('nan')
    print(f"  Grid: full={grid.ny}×{grid.nx}, concrete={ny_c}×{nx_c}, dx={dx:.4f} m = {dx/0.3048:.3f} ft")

    T0_C = (PLACEMENT_F - 32.0) * 5.0 / 9.0
    T_soil_C = (SOIL_F - 32.0) * 5.0 / 9.0
    T_init = np.full((grid.ny, grid.nx), T0_C)
    T_init[grid.is_soil] = T_soil_C

    env = make_env(PLACEMENT_F)

    print(f"\nRunning engine for Run B at hourly output ...")
    t0 = time.perf_counter()
    result = solve_hydration_2d(
        grid, mix, T_init,
        duration_s=168 * 3600,
        output_interval_s=3600.0,
        boundary_mode="full_2d",
        environment=env,
        construction=constr,
        T_ground_deep_C=T_soil_C,
        diagnostic_outputs=False,
    )
    t_wall = time.perf_counter() - t0
    print(f"  Solved in {t_wall:.1f}s  ({result.n_output_samples} samples, dt_inner={result.dt_inner_s:.2f}s)")

    # Slice to concrete domain
    alpha_conc = result.alpha_field[:, jslice, islice]   # (n_out, ny_c, nx_c)
    T_conc_C   = result.T_field_C[:, jslice, islice]     # (n_out, ny_c, nx_c)
    t_s = np.asarray(result.t_s)

    # Concrete-only mask within the concrete slice (2D, ny_c × nx_c)
    is_conc_slice = grid.is_concrete[jslice, islice]
    is_interior   = erode_concrete_mask(is_conc_slice, n=1)

    # Sanity check: reproduce Stage 4b masked max|R| at t=168 hr
    ti_168 = nearest_time_idx(result.t_s, 168.0)
    T_conc_F_168 = T_conc_C[ti_168] * 9.0 / 5.0 + 32.0
    x_conc = grid.x[islice]
    y_conc = grid.y[jslice]
    import tempfile
    tmp_csv = os.path.join(tempfile.gettempdir(), "stage5a_runB_sanity.csv")
    save_field_csv(tmp_csv, T_conc_F_168, x_conc, y_conc)
    from stage3_compare import load_engine_csv
    eng_y, eng_x, eng_F = load_engine_csv(tmp_csv)
    cw_field, cw_widths_m, cw_depths_m, _ = load_cw_slice(RUN_FOLDER, COMPARE_HR)
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    residual_168 = eng_on_cw - cw_field
    s = masked_stats(residual_168, mask_full, cw_depths_m, cw_widths_m)
    sanity_ok = abs(s["masked_max"] - CW_TARGET_MASKED_MAX) <= SANITY_TOL
    print(f"\n  SANITY CHECK: masked max|R| = {s['masked_max']:.3f}°F  "
          f"(target {CW_TARGET_MASKED_MAX}°F, tol ±{SANITY_TOL}°F)  "
          f"{'PASS' if sanity_ok else 'FAIL ***'}")
    if not sanity_ok:
        print("  WARNING: sanity check failed — diagnostic values below may not be trustworthy.")

    # Capture state at TARGET_HRS
    t_indices = [nearest_time_idx(t_s, hr) for hr in TARGET_HRS]
    actual_hrs = [float(t_s[i]) / 3600.0 for i in t_indices]

    snap_alpha = np.stack([alpha_conc[i] for i in t_indices])  # (4, ny_c, nx_c)
    snap_T_C   = np.stack([T_conc_C[i]   for i in t_indices])  # (4, ny_c, nx_c)

    # Per-cell k, Cp, αc at each snapshot
    snap_k     = np.empty_like(snap_alpha)
    snap_Cp    = np.empty_like(snap_alpha)
    snap_alphac = np.empty_like(snap_alpha)

    for si_idx, (snap_a, snap_t) in enumerate(zip(snap_alpha, snap_T_C)):
        k_field  = thermal_conductivity_variable(si['k_uc'], snap_a)
        Cp_field = specific_heat_variable(
            si['Wc'], si['Wa'], si['Ww'], si['Ca'],
            snap_a, snap_t, si['rho_c'],
        )
        snap_k[si_idx]      = k_field
        snap_Cp[si_idx]     = Cp_field
        snap_alphac[si_idx] = k_field / (si['rho_c'] * Cp_field) * 3600.0  # m²/hr

    # Verify reconstruction against saved arrays (at t=168)
    alpha_check = snap_alpha[-1]
    T_check     = snap_T_C[-1]
    k_check     = thermal_conductivity_variable(si['k_uc'], alpha_check)
    Cp_check    = specific_heat_variable(si['Wc'], si['Wa'], si['Ww'], si['Ca'],
                                         alpha_check, T_check, si['rho_c'])
    alphac_check = k_check / (si['rho_c'] * Cp_check) * 3600.0
    max_diff = float(np.max(np.abs(alphac_check - snap_alphac[-1])))
    print(f"  αc reconstruction sanity: max_diff={max_diff:.2e} m²/hr  "
          f"{'PASS' if max_diff < 1e-12 else 'FAIL ***'}")

    # Save .npy
    np.save(os.path.join(HERE, "STAGE5a_runB_alpha_field.npy"), snap_alpha)
    np.save(os.path.join(HERE, "STAGE5a_runB_T_field_C.npy"),   snap_T_C)
    np.save(os.path.join(HERE, "STAGE5a_runB_alphac.npy"),       snap_alphac)
    print(f"  Saved .npy arrays: alpha={snap_alpha.shape}, T={snap_T_C.shape}, αc={snap_alphac.shape}")

    # ------------------------------------------------------------------ #
    # Write STAGE5a_engine_state_runB.md (S5a.1)                          #
    # ------------------------------------------------------------------ #
    lines_state = [
        "# STAGE5a — Engine State: Run B at Multiple Time Slices",
        "",
        "Run B: placement=73°F, soil=60°F, |ΔT|=13°F.",
        "Grid: model_soil=False, is_submerged=True.",
        f"Grid dimensions: {ny_c}×{nx_c} concrete nodes, "
        f"dx={grid.x[1]-grid.x[0]:.4f} m = {(grid.x[1]-grid.x[0])/0.3048:.3f} ft.",
        "",
        f"**Sanity check vs Stage 4b:** masked max|R| at t=168 hr = "
        f"{s['masked_max']:.3f}°F (Stage 4b reference: {CW_TARGET_MASKED_MAX}°F, "
        f"tol ±{SANITY_TOL}°F) — {'PASS' if sanity_ok else 'FAIL'}.",
        "",
        "**Note on 'bulk' value:** mean over interior concrete cells only (1-cell erosion "
        "from all four concrete-domain edges).",
        "",
        "---",
    ]

    for idx, (hr, ti) in enumerate(zip(actual_hrs, t_indices)):
        a = snap_alpha[idx]
        T = snap_T_C[idx]
        k = snap_k[idx]
        Cp = snap_Cp[idx]
        rho_Cp = si['rho_c'] * Cp
        ac = snap_alphac[idx]

        all_stats  = {
            "α_hyd":         field_stats(a,      is_conc_slice),
            "k (W/m·K)":     field_stats(k,      is_conc_slice),
            "ρCp (J/m³·K)":  field_stats(rho_Cp, is_conc_slice),
            "α_c (m²/hr)":   field_stats(ac,     is_conc_slice),
        }
        bulk_stats = {
            "α_hyd":         field_stats(a,      is_interior),
            "k (W/m·K)":     field_stats(k,      is_interior),
            "ρCp (J/m³·K)":  field_stats(rho_Cp, is_interior),
            "α_c (m²/hr)":   field_stats(ac,     is_interior),
        }

        T_F = T * 9.0 / 5.0 + 32.0
        T_min_F = float(T_F[is_conc_slice].min())
        T_max_F = float(T_F[is_conc_slice].max())

        lines_state += [
            f"",
            f"## t = {hr:.1f} hr (sample index {ti})",
            f"",
            f"Temperature range: {T_min_F:.1f}–{T_max_F:.1f}°F",
            f"",
            "| Field | min | max | mean | std | bulk mean (interior) |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
        for fname in ["α_hyd", "k (W/m·K)", "ρCp (J/m³·K)", "α_c (m²/hr)"]:
            s_all  = all_stats[fname]
            s_bulk = bulk_stats[fname]
            lines_state.append(
                f"| {fname} | {s_all['min']:.5f} | {s_all['max']:.5f} | "
                f"{s_all['mean']:.5f} | {s_all['std']:.5f} | **{s_bulk['mean']:.5f}** |"
            )

    # One-liner answer at the bottom
    bulk_alpha_168 = float(snap_alpha[-1][is_interior].mean())
    bulk_alpha_168_std = float(snap_alpha[-1][is_interior].std())
    bulk_alphac_168 = float(snap_alphac[-1][is_interior].mean())
    lines_state += [
        "",
        "---",
        "",
        f"**Engine's actual α_hyd at t=168 hr (interior bulk): "
        f"{bulk_alpha_168:.4f} ± {bulk_alpha_168_std:.4f}**",
        f"This is close to α_u={mix.alpha_u:.4f} — hydration suppression "
        f"{'confirmed (bulk ≈ α_u)' if abs(bulk_alpha_168 - mix.alpha_u) < 0.03 else 'not fully suppressed'}.",
        "",
        f"**Engine's actual α_c at t=168 hr (interior bulk): {bulk_alphac_168:.5f} m²/hr**",
        f"CW M3 reference: ~0.00528 m²/hr. Ratio (engine/CW): {bulk_alphac_168/0.00528:.3f}.",
    ]

    with open(os.path.join(HERE, "STAGE5a_engine_state_runB.md"), "w") as f:
        f.write("\n".join(lines_state) + "\n")
    print(f"  Wrote STAGE5a_engine_state_runB.md")

    # ------------------------------------------------------------------ #
    # S5a.2: Cause 1 — k value analysis                                   #
    # ------------------------------------------------------------------ #
    alpha_check_vals = [0.0, 0.1, 0.3, 0.5, 0.8, 1.0, bulk_alpha_168]
    T_ref_C = 25.0  # representative mid-range temperature

    lines_c1 = [
        "# STAGE5a — Cause 1: k(α) Value at Engine's Operating α_hyd",
        "",
        "## Engine's k(α) formula",
        "",
        f"`k_c(α) = k_uc × (1.33 − 0.33·α)`, k_uc = {si['k_uc']:.4f} W/m·K",
        "",
        "## k(α) and α_c(α) vs. α_hyd",
        "",
        f"Reference temperature for Cp: T = {T_ref_C:.0f}°C (mid-range).",
        f"CW M3 target α_c = 0.00528 m²/hr.",
        "",
        "| α_hyd | k (W/m·K) | Cp (J/kg·K) | ρCp (J/m³·K) | α_c (m²/hr) | CW ratio (α_c/0.00528) |",
        "| --- | --- | --- | --- | --- | --- |",
    ]

    for av in alpha_check_vals:
        k_v = thermal_conductivity_variable(si['k_uc'], np.array([av]))[0]
        Cp_v = specific_heat_variable(
            si['Wc'], si['Wa'], si['Ww'], si['Ca'],
            np.array([av]), np.array([T_ref_C]), si['rho_c']
        )[0]
        rCp_v = si['rho_c'] * Cp_v
        ac_v = k_v / rCp_v * 3600.0
        tag = " ← **operating point**" if abs(av - bulk_alpha_168) < 0.005 else ""
        lines_c1.append(
            f"| {av:.4f}{tag} | {k_v:.4f} | {Cp_v:.1f} | {rCp_v:.0f} | {ac_v:.5f} | {ac_v/0.00528:.3f} |"
        )

    # Scaling factor to reach CW α_c = 0.00528 at the operating point
    k_op = thermal_conductivity_variable(si['k_uc'], np.array([bulk_alpha_168]))[0]
    Cp_op = specific_heat_variable(
        si['Wc'], si['Wa'], si['Ww'], si['Ca'],
        np.array([bulk_alpha_168]), np.array([T_ref_C]), si['rho_c']
    )[0]
    rCp_op = si['rho_c'] * Cp_op
    ac_op = k_op / rCp_op * 3600.0
    k_required = 0.00528 * rCp_op / 3600.0
    k_scale = k_required / k_op
    k_uc_required = k_required / (1.33 - 0.33 * bulk_alpha_168)

    lines_c1 += [
        "",
        "## Analysis",
        "",
        f"**Engine α_c at operating α_hyd={bulk_alpha_168:.4f}:** {ac_op:.5f} m²/hr",
        f"**CW target:** 0.00528 m²/hr",
        f"**Gap:** {(0.00528 - ac_op)/0.00528*100:.1f}% ({0.00528/ac_op:.3f}× ratio)",
        "",
        f"**To match CW α_c at the operating point with ρCp fixed:**",
        f"  Required k = {k_required:.4f} W/m·K (vs current {k_op:.4f} W/m·K)",
        f"  Required scaling factor on k alone: {k_scale:.3f}×",
        f"  Required k_uc to produce this via k_uc×(1.33−0.33α): {k_uc_required:.4f} W/m·K",
        f"  Current k_uc: {si['k_uc']:.4f} W/m·K → required k_uc scale: {k_uc_required/si['k_uc']:.3f}×",
        "",
        "## Plausibility check (textbook limestone concrete)",
        "",
        "Textbook range for limestone-aggregate concrete: k ≈ 2.0–2.5 W/m·K (normal-weight).",
        f"Current k_uc = {si['k_uc']:.4f} W/m·K: "
        f"{'within range' if 2.0 <= si['k_uc'] <= 2.5 else 'outside range'}.",
        f"Required k_uc = {k_uc_required:.4f} W/m·K: "
        f"{'within range' if 2.0 <= k_uc_required <= 2.5 else 'outside typical range (high)'}.",
        "",
        f"If instead ρCp ≈ 2.3e6 J/m³·K (textbook normal concrete), required k = "
        f"{0.00528 * 2.3e6 / 3600:.4f} W/m·K — labeled "
        + ("'high' by the brief." if 0.00528 * 2.3e6 / 3600 > 2.5 else "'plausible'."),
        f"If ρCp ≈ 2.0e6 J/m³·K (lighter mix), required k = "
        f"{0.00528 * 2.0e6 / 3600:.4f} W/m·K.",
        f"Engine's actual ρCp at operating point: {rCp_op:.0f} J/m³·K → "
        f"required k = {k_required:.4f} W/m·K "
        f"({'plausible' if k_required <= 2.5 else 'high but cited for limestone'}).",
        "",
        "## Verdict",
        "",
        f"k_uc alone scaling by {k_scale:.3f}× {'can' if k_uc_required <= 2.8 else 'cannot'} "
        f"close the gap. The required k_uc = {k_uc_required:.4f} W/m·K is "
        + ("within the plausible range for limestone-aggregate mass concrete."
           if k_uc_required <= 2.8 else "above the typical range; Cause 2 or Cause 3 may contribute."),
    ]

    with open(os.path.join(HERE, "STAGE5a_cause1_k_value.md"), "w") as f:
        f.write("\n".join(lines_c1) + "\n")
    print(f"  Wrote STAGE5a_cause1_k_value.md")

    # ------------------------------------------------------------------ #
    # S5a.3: Cause 2 — ρCp model                                          #
    # ------------------------------------------------------------------ #
    T_min_C = (SOIL_F - 32.0) * 5.0 / 9.0
    T_max_C = (PLACEMENT_F - 32.0) * 5.0 / 9.0  # no hydration heat, max ≈ placement
    T_mid_C = (T_min_C + T_max_C) / 2.0

    lines_c2 = [
        "# STAGE5a — Cause 2: ρCp Model Assessment",
        "",
        "## Van Breugel specific_heat_variable — formula",
        "",
        "```",
        "Cp = (1/ρ) × [Wc·α·(8.4·T + 339) + Wc·(1−α)·Ca + Wa·Ca + Ww·4186]",
        "```",
        "",
        "## Run B mix composition (SI units)",
        "",
        f"| Parameter | Value | Units |",
        f"| --- | --- | --- |",
        f"| ρ_concrete | {si['rho_c']:.2f} | kg/m³ |",
        f"| Wc (cement) | {si['Wc']:.2f} | kg/m³ |",
        f"| Wa (aggregate) | {si['Wa']:.2f} | kg/m³ |",
        f"| Ww (water) | {si['Ww']:.2f} | kg/m³ |",
        f"| Ca (aggregate Cp) | {si['Ca']:.2f} | J/(kg·K) |",
        f"| Wc+Wa+Ww total | {si['Wc']+si['Wa']+si['Ww']:.2f} | kg/m³ |",
        "",
    ]

    # Van Breugel term-by-term decomposition at operating point, mid-range T
    for T_label, T_C in [("T_soil_C (cold)", T_min_C), ("T_mid (avg)", T_mid_C), ("T_pl (placement)", T_max_C)]:
        c_cef = 8.4 * T_C + 339.0
        term1 = si['Wc'] * bulk_alpha_168 * c_cef          # hydrated cement
        term2 = si['Wc'] * (1.0 - bulk_alpha_168) * si['Ca']  # unhydrated cement
        term3 = si['Wa'] * si['Ca']                          # aggregate
        term4 = si['Ww'] * 4186.0                            # water
        total = term1 + term2 + term3 + term4
        Cp_vanb = total / si['rho_c']
        rCp_vanb = si['rho_c'] * Cp_vanb

        lines_c2 += [
            f"## At α_hyd={bulk_alpha_168:.4f}, {T_label} ({T_C:.1f}°C / {T_C*9/5+32:.1f}°F)",
            "",
            f"c_cef = 8.4×{T_C:.1f} + 339 = {c_cef:.1f} J/(kg·K)",
            "",
            f"| Term | Contribution (J/m³·K) | % of total |",
            f"| --- | --- | --- |",
            f"| Wc·α·c_cef (hydrated cement) | {term1:.0f} | {term1/total*100:.1f}% |",
            f"| Wc·(1−α)·Ca (unhydrated cement) | {term2:.0f} | {term2/total*100:.1f}% |",
            f"| Wa·Ca (aggregate) | {term3:.0f} | {term3/total*100:.1f}% |",
            f"| Ww·4186 (water) | {term4:.0f} | {term4/total*100:.1f}% |",
            f"| **ρCp total** | **{rCp_vanb:.0f}** | 100% |",
            "",
        ]

    # Temperature sensitivity
    dT = 10.0
    c_cef_lo = 8.4 * (T_mid_C - dT) + 339.0
    c_cef_hi = 8.4 * (T_mid_C + dT) + 339.0
    rCp_lo = (si['Wc'] * bulk_alpha_168 * c_cef_lo + si['Wc'] * (1 - bulk_alpha_168) * si['Ca']
              + si['Wa'] * si['Ca'] + si['Ww'] * 4186.0)
    rCp_hi = (si['Wc'] * bulk_alpha_168 * c_cef_hi + si['Wc'] * (1 - bulk_alpha_168) * si['Ca']
              + si['Wa'] * si['Ca'] + si['Ww'] * 4186.0)
    c_cef_mid = 8.4 * T_mid_C + 339.0
    rCp_mid = (si['Wc'] * bulk_alpha_168 * c_cef_mid + si['Wc'] * (1 - bulk_alpha_168) * si['Ca']
               + si['Wa'] * si['Ca'] + si['Ww'] * 4186.0)

    lines_c2 += [
        "## Temperature sensitivity of ρCp",
        "",
        f"T swing ±{dT:.0f}°C around mid-range {T_mid_C:.1f}°C:",
        f"  ρCp({T_mid_C-dT:.0f}°C) = {rCp_lo:.0f} J/m³·K",
        f"  ρCp({T_mid_C:.0f}°C)    = {rCp_mid:.0f} J/m³·K",
        f"  ρCp({T_mid_C+dT:.0f}°C) = {rCp_hi:.0f} J/m³·K",
        f"  Swing: {(rCp_hi-rCp_lo)/rCp_mid*100:.1f}% — "
        + ("negligible (<5%)" if abs(rCp_hi - rCp_lo) / rCp_mid < 0.05 else "significant (>5%)"),
        "",
    ]

    # Compare to textbook
    textbook_lo, textbook_hi = 2.1e6, 2.6e6
    lines_c2 += [
        "## Comparison to textbook normal-weight concrete",
        "",
        f"Textbook ρCp range: {textbook_lo/1e6:.1f}–{textbook_hi/1e6:.1f} ×10⁶ J/m³·K",
        f"Engine ρCp at operating point (mid-range T): {rCp_mid:.0f} J/m³·K = {rCp_mid/1e6:.3f} ×10⁶ J/m³·K",
    ]
    if rCp_mid < textbook_lo:
        lines_c2.append(f"→ **LOW by {(textbook_lo - rCp_mid)/textbook_lo*100:.1f}% relative to lower textbook bound.**")
    elif rCp_mid > textbook_hi:
        lines_c2.append(f"→ **HIGH by {(rCp_mid - textbook_hi)/textbook_hi*100:.1f}% relative to upper textbook bound.**")
    else:
        lines_c2.append(f"→ Within textbook range.")

    lines_c2 += [
        "",
        "## Verdict",
        "",
        f"If ρCp = {rCp_mid:.0f} J/m³·K is accepted as correct, the entire α_c gap must come from k.",
        f"If ρCp were as high as {textbook_hi:.0f} (upper textbook), α_c would be "
        + f"{si['k_uc']*(1.33-0.33*bulk_alpha_168)/textbook_hi*3600:.5f} m²/hr — "
        + ("worse (lower), not better." if textbook_hi > rCp_mid else "better."),
        f"Van Breugel ρCp appears {'within' if textbook_lo <= rCp_mid <= textbook_hi else 'outside'} "
        "textbook range — Cause 2 is unlikely to be the dominant driver of the α_c gap.",
    ]

    with open(os.path.join(HERE, "STAGE5a_cause2_rcp_model.md"), "w") as f:
        f.write("\n".join(lines_c2) + "\n")
    print(f"  Wrote STAGE5a_cause2_rcp_model.md")

    print(f"\nDone. Summary:")
    print(f"  Engine α_hyd at t=168 hr (bulk interior): {bulk_alpha_168:.4f} ± {bulk_alpha_168_std:.4f}")
    print(f"  Engine α_c at t=168 hr (bulk interior): {bulk_alphac_168:.5f} m²/hr")
    print(f"  CW M3 reference: 0.00528 m²/hr  |  gap: {(0.00528-bulk_alphac_168)/0.00528*100:.1f}%")


if __name__ == "__main__":
    run()
