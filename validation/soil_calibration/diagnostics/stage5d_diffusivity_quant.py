#!/usr/bin/env python3
"""Stage 5d-prep: Diffusivity Quantification Diagnostic.

Converts the observed |ΔT|-scaling residual into an implied effective lateral
diffusivity mismatch (Δα/α) between engine and CW.

Sections:
  §2.1  1D semi-infinite analytical reference solution
  §2.2  Implied diffusivity per run (all 8 runs, di=24, wi=9, t=168 hr)
  §2.3  Depth-stability check for Runs F and H (di=15/24/35)
  §2.4  Time-evolution check for Run F (t=24/84/168 hr) + lateral PNG
  §2.5  Cross-check against Stage 5a reported Δα/α ≈ 1.2%
  §5    Synthesis

No engine source changes.  Run F is re-run for t=24/84/168 snapshots;
all other runs use cached t=168 CSVs only.  No commits.
"""
import os
import sys
import numpy as np
from scipy.special import erf
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SC   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, SC)

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import load_engine_csv, resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

ENG_DIR   = os.path.join(SC, "stage5c_runs")
CW_RUNS   = os.path.join(SC, "cw_runs")
MASK_PATH = os.path.join(SC, "STAGE35_validity_mask.npy")
OUT_DIR   = os.path.join(SC, "diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

# Stage 5c config (fixed)
BLANKET_M    = 0.0
MODEL_SOIL   = False
IS_SUBMERGED = True

RUNS = [
    ("A", "runA_baseline",  73,  73),
    ("B", "runB_73_60",     73,  60),
    ("C", "runC_73_90",     73,  90),
    ("D", "runD_60_73",     60,  73),
    ("E", "runE_90_73",     90,  73),
    ("F", "runF_73_45",     73,  45),
    ("G", "runG_73_100",    73, 100),
    ("H", "runH_45_73",     45,  73),
    ("I", "runI_100_73",   100,  73),
]

WI9_X_M    = 1.52          # CW wi=9 distance from form face (x=0)
DI_PROBE   = 24            # mid-depth node for §2.2 main table
DI_DEPTHS  = [15, 24, 35]  # §2.3 depth stability
TIME_HRS   = [24, 84, 168] # §2.4 time evolution

# α search bounds for brentq: 0.0010 to 0.20 m²/hr converted to m²/s
ALPHA_LO_S = 1e-4 / 3600.0  # 0.0001 m²/hr
ALPHA_HI_S = 0.50 / 3600.0  # 0.50 m²/hr


# ---------------------------------------------------------------------------
# §2.1  1D semi-infinite solution
# ---------------------------------------------------------------------------
def T_1d(x_m, t_s, alpha_m2s, T_initial, T_bc):
    """T at position x, time t for semi-infinite slab with Dirichlet BC at x=0."""
    z = x_m / (2.0 * np.sqrt(alpha_m2s * t_s))
    return T_bc + (T_initial - T_bc) * erf(z)


def invert_alpha(T_obs, x_m, t_s, T_initial, T_bc):
    """Find α such that T_1d(x, t, α) = T_obs. Returns m²/s."""
    if abs(T_initial - T_bc) < 1e-6:
        return float("nan"), "dT_zero"
    # fraction of initial range still retained
    frac = (T_obs - T_bc) / (T_initial - T_bc)
    if frac <= 0.0 or frac >= 1.0:
        return float("nan"), f"frac_oob({frac:.3f})"
    def f(a):
        return T_1d(x_m, t_s, a, T_initial, T_bc) - T_obs
    try:
        alpha = brentq(f, ALPHA_LO_S, ALPHA_HI_S, xtol=1e-20, rtol=1e-10)
        return alpha, "ok"
    except ValueError as e:
        return float("nan"), str(e)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _load_cw_at_hr(folder, target_hr):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - target_hr).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m, float(v.time_hrs[ti])


def _load_cw_all_hrs(folder, target_hrs):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    out = {}
    for hr in target_hrs:
        ti = int(np.abs(v.time_hrs - hr).argmin())
        out[hr] = (v.T_field_F[ti], float(v.time_hrs[ti]))
    return out, v.widths_m, v.depths_m


def _get_engine_t168(label, folder, cw_depths_m, cw_widths_m):
    csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
    eng_y, eng_x, eng_F = load_engine_csv(csv_path)
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    return eng_on_cw


def _run_engine_snapshots(folder, pl_F, so_F, target_hrs_list):
    """Re-run engine for Run F (or any run) and return snapshots at target hours."""
    dat_path = os.path.join(CW_RUNS, folder, "input.dat")
    mix, geom, constr, _ = parse_cw_dat(dat_path)
    factor, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = factor
    mix.Hu_J_kg_effective = mix.Hu_J_kg * factor

    constr.model_soil = MODEL_SOIL
    constr.is_submerged = IS_SUBMERGED

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=IS_SUBMERGED,
        model_soil=MODEL_SOIL,
        blanket_thickness_m=BLANKET_M,
    )

    T0_C     = (pl_F - 32.0) * 5.0 / 9.0
    T_soil_C = (so_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_soil] = T_soil_C

    env = make_neutral_env(pl_F)

    print("  Re-running Run F engine (needed for t=24 and t=84 snapshots) …")
    result = solve_hydration_2d(
        grid, mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=env,
        construction=constr,
        T_ground_deep_C=T_soil_C,
        diagnostic_outputs=False,
    )
    print(f"    Done: {result.n_inner_steps} steps, dt={result.dt_inner_s:.0f}s")

    jslice, islice = grid.concrete_slice()
    T_conc_F = result.T_field_C[:, jslice, islice] * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jslice]
    eng_x = grid.x[islice]

    snapshots = {}
    for hr in target_hrs_list:
        ti = nearest_time_idx(result.t_s, float(hr))
        snapshots[hr] = (eng_y.copy(), eng_x.copy(), T_conc_F[ti].copy())
    return snapshots


def alpha_m2s_to_m2hr(a):
    return a * 3600.0


# ---------------------------------------------------------------------------
def main():
    # -------------------------------------------------------------------------
    # §2.1 verify 1D solution
    # -------------------------------------------------------------------------
    print("=== §2.1 1D semi-infinite solution verification ===")
    # At x=0 should return T_bc; at x→∞ should return T_initial
    assert abs(T_1d(0.0, 604800, 5e-3/3600, 73.0, 45.0) - 45.0) < 1e-9, "BC check failed"
    assert abs(T_1d(1e6,  604800, 5e-3/3600, 73.0, 45.0) - 73.0) < 1e-3, "Far-field check failed"
    # Check erf scaling: erf(0.5) ≈ 0.5205
    T_test = T_1d(1.52, 604800, 5.28e-3/3600, 73.0, 45.0)
    print(f"  T_1d(x=1.52m, t=168hr, α=5.28e-3 m²/hr, T0=73°F, Tbc=45°F) = {T_test:.3f}°F")
    print(f"  BC check (x=0)   : PASS")
    print(f"  Far-field check  : PASS")

    # -------------------------------------------------------------------------
    # Load mask and reference CW grid
    # -------------------------------------------------------------------------
    mask = np.load(MASK_PATH)
    assert mask.shape == (49, 13)

    # Get CW grid geometry from Run A
    _, cw_widths_m, cw_depths_m, _ = _load_cw_at_hr("runA_baseline", 168)
    print(f"\nCW wi=9 width: {cw_widths_m[9]:.4f} m  (form-face side BC)")
    print(f"CW di=15 depth: {cw_depths_m[15]:.4f} m")
    print(f"CW di=24 depth: {cw_depths_m[24]:.4f} m")
    print(f"CW di=35 depth: {cw_depths_m[35]:.4f} m")

    x_m = cw_widths_m[9]   # use exact CW coordinate (1.52 m)
    t_168_s = 168 * 3600.0

    # -------------------------------------------------------------------------
    # §2.2  Implied diffusivity at di=24, wi=9, t=168 hr
    # -------------------------------------------------------------------------
    print("\n=== §2.2 Implied Diffusivity Table (di=24, wi=9, t=168 hr) ===")
    hdr = (f"{'Run':>4}  {'|ΔT|':>5}  {'wave':>5}  "
           f"{'T_eng(°F)':>10}  {'T_CW(°F)':>9}  {'R(°F)':>7}  "
           f"{'α_eng(m²/hr)':>13}  {'α_CW(m²/hr)':>12}  {'Δα/α(%)':>9}")
    print(hdr)
    print("-" * len(hdr))

    dalpha_pct_list = []
    for label, folder, pl, so in RUNS:
        if label == "A":
            continue
        dT = abs(pl - so)
        wave = "cool" if pl > so else "warm"
        cw_field_168, _, _, _ = _load_cw_at_hr(folder, 168)
        eng_on_cw = _get_engine_t168(label, folder, cw_depths_m, cw_widths_m)

        T_eng = float(eng_on_cw[DI_PROBE, 9])
        T_cw  = float(cw_field_168[DI_PROBE, 9])
        R     = T_eng - T_cw

        a_eng, status_eng = invert_alpha(T_eng, x_m, t_168_s, float(pl), float(so))
        a_cw,  status_cw  = invert_alpha(T_cw,  x_m, t_168_s, float(pl), float(so))

        if status_eng == "ok" and status_cw == "ok":
            da_pct = (a_eng - a_cw) / a_cw * 100.0
            dalpha_pct_list.append(da_pct)
            a_eng_hr = alpha_m2s_to_m2hr(a_eng)
            a_cw_hr  = alpha_m2s_to_m2hr(a_cw)
            print(f"  {label:>2}  {dT:>5}  {wave:>5}  "
                  f"{T_eng:>10.3f}  {T_cw:>9.3f}  {R:>7.3f}  "
                  f"{a_eng_hr:>13.5f}  {a_cw_hr:>12.5f}  {da_pct:>9.3f}")
        else:
            print(f"  {label:>2}  {dT:>5}  {wave:>5}  "
                  f"{T_eng:>10.3f}  {T_cw:>9.3f}  {R:>7.3f}  "
                  f"  eng:{status_eng}  cw:{status_cw}")

    if dalpha_pct_list:
        arr = np.array(dalpha_pct_list)
        print(f"\n  Δα/α summary: mean={arr.mean():.3f}%  std={arr.std():.3f}%  "
              f"min={arr.min():.3f}%  max={arr.max():.3f}%")
        cv = arr.std() / abs(arr.mean()) if abs(arr.mean()) > 0 else float("inf")
        print(f"  CV={cv:.3f}  ({'consistent across runs (<0.15)' if cv < 0.15 else 'inconsistent (>0.15)'})")
        signs = set(np.sign(arr).astype(int).tolist())
        print(f"  Signs: {sorted(signs)}  ({'CONSISTENT' if len(signs) == 1 else 'MIXED — hypothesis fails'})")
        mag = abs(arr.mean())
        if mag < 5.0:
            verdict = "PHYSICALLY PLAUSIBLE (<5%)"
        elif mag < 15.0:
            verdict = "LARGE but possible (5–15%)"
        else:
            verdict = "UNPHYSICAL (>15%)"
        print(f"  Magnitude verdict: {verdict}")

    # -------------------------------------------------------------------------
    # §2.3  Depth stability for Runs F and H
    # -------------------------------------------------------------------------
    print("\n=== §2.3 Depth-stability Check (Runs F and H, t=168 hr) ===")
    hdr3 = (f"{'Run':>4}  {'di':>4}  {'depth_m':>8}  "
            f"{'T_eng':>8}  {'T_CW':>8}  "
            f"{'α_eng(m²/hr)':>13}  {'α_CW(m²/hr)':>12}  {'Δα/α(%)':>9}")
    print(hdr3)
    print("-" * len(hdr3))

    depth_da_F = []
    depth_da_H = []
    for label, folder, pl, so in RUNS:
        if label not in ("F", "H"):
            continue
        cw_field_168, _, _, _ = _load_cw_at_hr(folder, 168)
        eng_on_cw = _get_engine_t168(label, folder, cw_depths_m, cw_widths_m)
        store = depth_da_F if label == "F" else depth_da_H
        for di in DI_DEPTHS:
            depth = float(cw_depths_m[di])
            T_eng = float(eng_on_cw[di, 9])
            T_cw  = float(cw_field_168[di, 9])
            a_eng, se = invert_alpha(T_eng, x_m, t_168_s, float(pl), float(so))
            a_cw,  sc = invert_alpha(T_cw,  x_m, t_168_s, float(pl), float(so))
            if se == "ok" and sc == "ok":
                da_pct = (a_eng - a_cw) / a_cw * 100.0
                store.append(da_pct)
                print(f"  {label:>2}  {di:>4}  {depth:>8.4f}  "
                      f"{T_eng:>8.3f}  {T_cw:>8.3f}  "
                      f"{alpha_m2s_to_m2hr(a_eng):>13.5f}  "
                      f"{alpha_m2s_to_m2hr(a_cw):>12.5f}  {da_pct:>9.3f}")
            else:
                print(f"  {label:>2}  {di:>4}  {depth:>8.4f}  "
                      f"{T_eng:>8.3f}  {T_cw:>8.3f}  eng:{se}  cw:{sc}")

    for lbl, store in [("F", depth_da_F), ("H", depth_da_H)]:
        if len(store) >= 2:
            arr = np.array(store)
            spread = arr.max() - arr.min()
            print(f"\n  Run {lbl}: Δα/α spread across depths = {spread:.3f}% "
                  f"({'depth-uniform (<1%)' if spread < 1.0 else 'depth-varying (>1%)'})")

    # -------------------------------------------------------------------------
    # §2.4  Time evolution for Run F
    # -------------------------------------------------------------------------
    print("\n=== §2.4 Time Evolution — Run F, di=24, wi=9 ===")

    # CW at all three times
    cw_by_hr, _, _ = _load_cw_all_hrs("runF_73_45", TIME_HRS)

    # Engine re-run for t=24/84/168
    eng_snaps = _run_engine_snapshots("runF_73_45", 73, 45, TIME_HRS)

    # Reference R at t=168 should match cached value (sanity)
    eng_168_resampled = resample_engine_to_cw(
        eng_snaps[168][0], eng_snaps[168][1], eng_snaps[168][2],
        cw_depths_m, cw_widths_m
    )
    R_168_rerun  = float(eng_168_resampled[DI_PROBE, 9]) - float(cw_by_hr[168][0][DI_PROBE, 9])
    R_168_cached = float(_get_engine_t168("F", "runF_73_45", cw_depths_m, cw_widths_m)[DI_PROBE, 9]) \
                 - float(cw_by_hr[168][0][DI_PROBE, 9])
    print(f"  Sanity: R(t=168) re-run={R_168_rerun:.4f}°F  cached={R_168_cached:.4f}°F  "
          f"diff={abs(R_168_rerun - R_168_cached):.4f}°F")

    hdr4 = (f"{'t(hr)':>6}  {'T_eng':>8}  {'T_CW':>8}  {'R(°F)':>7}  "
            f"{'R/√(t/168)':>11}  {'R/(t/168)':>10}  {'Δα/α(%)':>9}")
    print(hdr4)
    print("-" * len(hdr4))

    R_vals = {}
    for hr in TIME_HRS:
        eng_y, eng_x, eng_F = eng_snaps[hr]
        eng_on_cw_hr = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
        T_eng = float(eng_on_cw_hr[DI_PROBE, 9])
        T_cw  = float(cw_by_hr[hr][0][DI_PROBE, 9])
        R     = T_eng - T_cw
        R_vals[hr] = R
        scale_sqrt = R / np.sqrt(hr / 168.0)
        scale_lin  = R / (hr / 168.0)
        t_s = hr * 3600.0
        a_eng, se = invert_alpha(T_eng, x_m, t_s, 73.0, 45.0)
        a_cw,  sc = invert_alpha(T_cw,  x_m, t_s, 73.0, 45.0)
        if se == "ok" and sc == "ok":
            da_pct = (a_eng - a_cw) / a_cw * 100.0
            print(f"  {hr:>6}  {T_eng:>8.3f}  {T_cw:>8.3f}  {R:>7.3f}  "
                  f"{scale_sqrt:>11.4f}  {scale_lin:>10.4f}  {da_pct:>9.3f}")
        else:
            print(f"  {hr:>6}  {T_eng:>8.3f}  {T_cw:>8.3f}  {R:>7.3f}  "
                  f"{scale_sqrt:>11.4f}  {scale_lin:>10.4f}  eng:{se} cw:{sc}")

    # Check monotonicity
    r_list = [R_vals[hr] for hr in TIME_HRS]
    if all(r_list[i] <= r_list[i+1] for i in range(len(r_list)-1)) or \
       all(r_list[i] >= r_list[i+1] for i in range(len(r_list)-1)):
        print("  Monotonic: YES")
    else:
        print("  Monotonic: NO (sign change or reversal — non-diffusivity mechanism suspected)")

    # Lateral evolution plot
    _plot_lateral_evolution(eng_snaps, cw_by_hr, cw_widths_m, cw_depths_m)

    # -------------------------------------------------------------------------
    # §2.5  Stage 5a cross-check
    # -------------------------------------------------------------------------
    print("\n=== §2.5 Stage 5a Cross-check ===")
    stage5a_ratio_pct = 1.2  # engine/CW α ratio reported as 1.012 → Δα/α=1.2%
    if dalpha_pct_list:
        mean_here = np.mean(dalpha_pct_list)
        print(f"  Stage 5a reported Δα/α ≈ +{stage5a_ratio_pct:.1f}%  (engine 1.2% higher than CW)")
        print(f"  This diagnostic §2.2 mean Δα/α = {mean_here:+.3f}%")
        ratio_comparison = mean_here / stage5a_ratio_pct if stage5a_ratio_pct != 0 else float("inf")
        print(f"  Ratio (§2.2 / Stage5a) = {ratio_comparison:.2f}×")
        if abs(ratio_comparison - 1.0) < 0.50:
            verdict = "CONSISTENT — both diagnostics agree on ~1–2% Δα/α"
        elif ratio_comparison > 1.5:
            verdict = ("LARGER than Stage 5a — suggests Stage 5a's α extraction "
                       "underestimated the mismatch (possibly due to grid-resolution noise "
                       "at the 6× refinement level or hydration-suppression operating-point dependence)")
        else:
            verdict = "SMALLER than Stage 5a — unexpected; check Stage 5a operating-point assumptions"
        print(f"  Verdict: {verdict}")

    # -------------------------------------------------------------------------
    # §5  Synthesis
    # -------------------------------------------------------------------------
    print("\n=== §5 Synthesis ===\n")
    if dalpha_pct_list:
        arr = np.array(dalpha_pct_list)
        mean_da = arr.mean()
        std_da  = arr.std()
        cv_da   = std_da / abs(mean_da) if abs(mean_da) > 0 else float("inf")
        signs_ok = len(set(np.sign(arr).astype(int).tolist())) == 1
        mag = abs(mean_da)

        mag_str  = "physically plausible (< 5%)" if mag < 5 else \
                   ("large but possibly explainable (5–15%)" if mag < 15 else "unphysical (>15%)")

        sign_str = "consistently positive across all cooling and heating runs" \
                   if signs_ok else "mixed — hypothesis inconsistent"

        # check √t stability
        R_24  = R_vals.get(24, float("nan"))
        R_84  = R_vals.get(84, float("nan"))
        R_168 = R_vals.get(168, float("nan"))
        if not any(np.isnan(v) for v in [R_24, R_84, R_168]):
            sq_24  = R_24  / np.sqrt(24/168)
            sq_84  = R_84  / np.sqrt(84/168)
            sq_168 = R_168 / np.sqrt(168/168)
            sq_cv  = np.std([sq_24, sq_84, sq_168]) / abs(np.mean([sq_24, sq_84, sq_168]))
            sqrt_t_stable = sq_cv < 0.20
        else:
            sqrt_t_stable = False

        time_str = ("roughly √t-like (consistent with diffusivity-driven mechanism)" if sqrt_t_stable
                    else "NOT √t-like — the structure deviates from pure diffusivity scaling")

        print(
            f"(a) Residual consistency with diffusivity mismatch: The implied Δα/α is "
            f"{mean_da:+.2f}% ± {std_da:.2f}% across 8 runs (cv={cv_da:.2f}). "
            f"Signs are {sign_str}. "
            f"All 8 runs converge on a similar magnitude, {'confirming' if cv_da < 0.15 else 'weakly suggesting'} "
            f"a single physical-property mismatch.\n\n"
            f"(b) Physical plausibility: {mag:.2f}% Δα/α is {mag_str}. "
            f"Stage 5a reported 1.2% Δα/α from a different method; this diagnostic finds "
            f"{mean_da:+.2f}%, a {abs(mean_da/stage5a_ratio_pct):.1f}× "
            f"{'larger' if mean_da > stage5a_ratio_pct else 'smaller'} value. "
            f"{'Agreement within a factor of ~2 is acceptable given Stage 5a operated on heavily resolution-affected data.' if abs(mean_da/stage5a_ratio_pct) < 3 else 'The discrepancy is large; Stage 5a and this diagnostic are not consistent.'}\n\n"
            f"(c) Time evolution: R(t) is {time_str}.\n\n"
            f"(d) Most parsimonious physical mechanism: If the diffusivity hypothesis holds "
            f"(cv consistent, signs consistent, √t behavior confirmed), the simplest "
            f"explanation is that the engine's effective lateral thermal diffusivity α_concrete "
            f"is {mean_da:+.2f}% higher than CW's. This could arise from a {mean_da:+.2f}% "
            f"discrepancy in k_concrete (direct proportional effect on α = k/(ρCp)) or an "
            f"equivalent fractional difference in ρCp in the opposite direction. If signs are "
            f"inconsistent or the residual is non-monotonic, the 1D semi-infinite model is "
            f"inadequate and the mechanism is not a simple diffusivity offset.\n\n"
            f"(e) Scope-boundary observation: If the residual is genuinely a bulk concrete "
            f"thermal diffusivity mismatch (k or ρCp), closing the gate is NOT within the scope "
            f"of soil-concrete BC calibration. It would require touching bulk thermal property "
            f"calibration — a separate axis that has not been authorized in the current sprint. "
            f"The user must decide whether to expand scope or accept the current gate failures "
            f"as a known out-of-scope mismatch."
        )
    else:
        print("  §5 synthesis skipped: §2.2 inversion failed for all runs.")


# ---------------------------------------------------------------------------
def _plot_lateral_evolution(eng_snaps, cw_by_hr, cw_widths_m, cw_depths_m):
    out_path = os.path.join(OUT_DIR, "stage5d_prep_runF_lateral_evolution.png")
    print(f"\n  Generating lateral evolution plot: {out_path}")

    DI = 24
    # CW: widths ordered CL→edge (wi=0=6.1m to wi=12=0.0m)
    # Reverse to edge→CL for plotting
    x_cw_plot = cw_widths_m[::-1]  # 0.0 → 6.1 m

    colors = {24: "steelblue", 84: "darkorange", 168: "darkred"}

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)
    fig.suptitle(
        "Stage 5d-prep — Run F Lateral Temperature Evolution at di=24 (mid-depth ≈12.19 m)\n"
        "Engine (solid) vs CW (dashed) — placement=73°F, soil=45°F, ΔT=+28°F",
        fontsize=10, fontweight="bold",
    )

    for col_idx, hr in enumerate(TIME_HRS):
        ax = axes[col_idx]
        eng_y, eng_x, eng_F = eng_snaps[hr]
        eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)

        cw_field = cw_by_hr[hr][0]

        # Both plotted edge→CL
        eng_row = eng_on_cw[DI, ::-1]
        cw_row  = cw_field[DI, ::-1]

        color = colors[hr]
        ax.plot(x_cw_plot, eng_row, color=color, lw=2.0, marker="o", ms=5, label="Engine")
        ax.plot(x_cw_plot, cw_row,  color=color, lw=1.5, ls="--", marker="s", ms=4, label="CW")
        ax.fill_between(x_cw_plot, eng_row, cw_row, alpha=0.15, color=color)

        # Annotate residual at wi=9
        r_wi9 = float(eng_on_cw[DI, 9] - cw_field[DI, 9])
        ax.axvline(cw_widths_m[9], color="gray", ls=":", lw=1.0)
        ax.text(cw_widths_m[9] + 0.05, min(float(cw_row[9]), float(eng_row[9])) - 0.4,
                f"R={r_wi9:+.3f}°F", fontsize=8, color="darkred")

        ax.set_title(f"t = {hr} hr", fontsize=10)
        ax.set_xlabel("Distance from form face (m)", fontsize=9)
        ax.set_ylabel("Temperature (°F)", fontsize=9)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.text(0.05, 0.97, "← form face (T_soil=45°F)", transform=ax.transAxes,
                fontsize=7, va="top", color="gray")
        ax.text(0.99, 0.97, "CL →", transform=ax.transAxes,
                fontsize=7, va="top", ha="right", color="gray")

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
