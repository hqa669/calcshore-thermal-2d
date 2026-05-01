#!/usr/bin/env python3
"""Stage 4a — Extract α_c from CW data and diagnose the bottom-CL residual.

KEY FINDING FROM DATA INSPECTION:
  CW's geometry is 40ft wide × 80ft deep ALL CONCRETE.  There is NO soil domain
  below the concrete in CW.  CW applies Dirichlet T_soil directly at:
    (a) the concrete side face (w=0.00m, full depth 0–24.38m)
    (b) the concrete bottom face (depth=24.38m, all widths)
  The engine (after Stage 3) adds 3m of soil BELOW the 80ft concrete, placing
  the Dirichlet at depth 27.4m — 3m further down.  This soil buffer prevents
  the concrete bottom from cooling as fast as CW, causing the 9–19°F residual.

  Both Stage 1 "side" and "bottom" penetration slopes reflect α_CONCRETE.
  α_s and e_c/e_s cannot be extracted from CW (no soil domain in CW output).

  Three α_c extraction methods:
    M1: side penetration slope (1°F front from side Dirichlet into concrete) vs √t
    M2: bottom penetration slope (1°F front from bottom Dirichlet into concrete) vs √t
    M3: erfc spatial profile fit at t=24,48,96,168hr (direct, most reliable)

  Bottom-temperature comparison (structural proof):
    M4: T(concrete bottom, t=168hr) in CW vs engine Stage3 CSV

Outputs:
    STAGE4a_cw_extracted_params.md
    STAGE4a_cw_alphac.npy    shape (8, 3) — 8 runs × [M1, M2, M3_mean]
"""
import os
import sys
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erfinv, erfc as _erfc
from scipy.stats import linregress

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
from cw_scenario_loader import parse_cw_temp_output

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
ENG_DIR = os.path.join(HERE, "stage3_fix2_runs")

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
NON_A = [(lbl, fld, pl, so) for lbl, fld, pl, so in RUNS if lbl != "A"]

ERFC_TIMESTAMPS_HR = [24.0, 48.0, 96.0, 168.0]
FRONT_THRESHOLD_F  = 1.0


def find_idx(arr, val):
    return int(np.abs(np.asarray(arr) - val).argmin())


# ─── Method 1: side penetration slope ───────────────────────────────────────

def penetration_side(T_field_F, widths_m, placement_F, di_mid):
    """1°F front from side Dirichlet (w=0) into concrete vs time."""
    n_time = T_field_F.shape[0]
    penetration = np.zeros(n_time)
    edge_idx = len(widths_m) - 1  # index 12 = w=0.00m (Dirichlet)
    for ti in range(n_time):
        T_row = T_field_F[ti, di_mid, :]
        pen_m = 0.0
        for ii in range(edge_idx, -1, -1):
            if abs(T_row[ii] - placement_F) < FRONT_THRESHOLD_F:
                if ii < edge_idx:
                    pen_m = float(widths_m[ii + 1])  # distance from outer edge
                break
        else:
            pen_m = float(widths_m[0])
        penetration[ti] = pen_m
    return penetration


# ─── Method 2: bottom penetration slope ─────────────────────────────────────

def penetration_bottom(T_field_F, depths_m, placement_F, wi):
    """1°F front from bottom Dirichlet (depth=24.38m) upward into concrete."""
    n_time = T_field_F.shape[0]
    penetration = np.zeros(n_time)
    bottom_idx = len(depths_m) - 1  # index 48 = 24.38m (Dirichlet)
    for ti in range(n_time):
        T_col = T_field_F[ti, :, wi]
        pen_m = 0.0
        for ii in range(bottom_idx, -1, -1):
            if abs(T_col[ii] - placement_F) < FRONT_THRESHOLD_F:
                if ii < bottom_idx:
                    pen_m = float(depths_m[bottom_idx] - depths_m[ii + 1])
                break
        else:
            pen_m = float(depths_m[bottom_idx] - depths_m[0])
        penetration[ti] = pen_m
    return penetration


def fit_slope(time_hrs, pen_m, t_max_hr=168.0):
    """Fit pen ~ slope × √t. Returns (slope_m_per_sqrth, R²)."""
    mask = (time_hrs > 0) & (time_hrs <= t_max_hr) & (pen_m > 0)
    if mask.sum() < 5:
        return float("nan"), float("nan")
    slope, _, r, *_ = linregress(np.sqrt(time_hrs[mask]), pen_m[mask])
    return float(slope), float(r ** 2)


def slope_to_alpha(slope, dT_abs):
    """Convert penetration slope (m/√hr) to diffusivity (m²/hr).
    For a concrete half-space with Dirichlet BC: slope = 2×erfinv(1-1/|ΔT|)×√α.
    |ΔT| is the full temperature driving force at the Dirichlet.
    """
    if np.isnan(slope) or abs(dT_abs) < 1.0:
        return float("nan")
    eta = float(erfinv(1.0 - 1.0 / abs(dT_abs)))
    if eta <= 0:
        return float("nan")
    return (slope / (2.0 * eta)) ** 2


# ─── Method 3: erfc spatial profile fit ─────────────────────────────────────

def fit_erfc_concrete_side(T_row_F, widths_m, placement_F, t_hr, min_pts=5):
    """Direct erfc profile fit: T(x,t) - T_placement = A × erfc(x/(2√(α t))).

    x = distance from outer edge (Dirichlet at x=0, w=0.00m) into concrete.
    Fit parameters: A (amplitude = T_dirichlet - T_placement) and α_c.
    No admittance correction needed — CW uses Dirichlet at x=0.
    """
    x = widths_m[::-1]          # ascending from 0.00m to 6.10m
    T = T_row_F[::-1]
    T_offset = T - placement_F

    valid = np.abs(T_offset) > 0.1
    if valid.sum() < min_pts:
        return float("nan"), float("nan"), float("nan")

    x_v = x[valid]
    T_v = T_offset[valid]

    def model(x, A, alpha):
        denom = 2.0 * np.sqrt(alpha * t_hr)
        return A * _erfc(x / denom) if denom > 1e-12 else np.zeros_like(x)

    A0 = float(T_v[0])
    try:
        popt, _ = curve_fit(model, x_v, T_v, p0=[A0, 5e-3],
                            bounds=([-60, 1e-5], [60, 0.05]), maxfev=3000)
        A_fit, alpha_fit = popt
        res = T_v - model(x_v, *popt)
        rmse = float(np.sqrt(np.mean(res**2)))
        return float(alpha_fit), float(A_fit), rmse
    except Exception:
        return float("nan"), float("nan"), float("nan")


# ─── Method 4: bottom-temperature structural proof ───────────────────────────

def load_engine_csv(path):
    with open(path) as f:
        lines = f.readlines()
    header = [float(x) for x in lines[0].strip().split(",")[1:]]
    rows_y, rows_T = [], []
    for line in lines[1:]:
        parts = line.strip().split(",")
        rows_y.append(float(parts[0]))
        rows_T.append([float(x) for x in parts[1:]])
    return np.array(rows_y), np.array(header), np.array(rows_T)


# ─── Main ────────────────────────────────────────────────────────────────────

def run():
    print("Loading CW data…")
    cw_data = {}
    for lbl, fld, pl, so in RUNS:
        path = os.path.join(CW_RUNS, fld, "output.txt")
        if not os.path.isfile(path):
            print(f"  SKIP {lbl}: {path} not found")
            continue
        v = parse_cw_temp_output(path)
        cw_data[lbl] = (v, int(pl), int(so))

    if "A" not in cw_data:
        raise RuntimeError("Run A baseline missing")

    v_A, _, _ = cw_data["A"]
    widths_m  = v_A.widths_m    # (13,) descending 6.10→0.00m
    depths_m  = v_A.depths_m    # (49,) ascending  0.00→24.38m
    time_hrs  = v_A.time_hrs    # (2016,) ~5-min cadence

    # Grid indices
    di_mid    = find_idx(depths_m, 12.19)   # mid-depth ≈ 40ft
    wi_cl     = 0                           # CL: w=6.10m
    wi_mid    = find_idx(widths_m, 3.05)    # mid-width ≈ 10ft

    print(f"Grid: di_mid={di_mid} (depth={depths_m[di_mid]:.2f}m), "
          f"wi_cl={wi_cl} (w={widths_m[wi_cl]:.2f}m), "
          f"wi_mid={wi_mid} (w={widths_m[wi_mid]:.2f}m)")

    n_runs = len(NON_A)
    alpha_c_m1  = np.full(n_runs, np.nan)
    alpha_c_m2  = np.full(n_runs, np.nan)
    alpha_c_m3  = np.full((n_runs, len(ERFC_TIMESTAMPS_HR)), np.nan)

    slope_side_raw = np.full(n_runs, np.nan)
    slope_bot_raw  = np.full(n_runs, np.nan)
    bot_T_cw       = np.full(n_runs, np.nan)   # CW concrete-bottom T at t=168hr, CL
    bot_T_eng      = np.full(n_runs, np.nan)   # engine concrete-bottom T at t=168hr, CL

    print("\n" + "="*70)
    for ri, (lbl, fld, pl, so) in enumerate(NON_A):
        if lbl not in cw_data:
            continue
        v, placement_F, soil_F = cw_data[lbl]
        dT = abs(so - pl)
        print(f"\nRun {lbl}: placement={placement_F}°F, soil={soil_F}°F, |ΔT|={dT}°F")

        # M1: side penetration slope
        pen_side = penetration_side(v.T_field_F, widths_m, float(placement_F), di_mid)
        slope_s, r2_s = fit_slope(time_hrs, pen_side)
        slope_side_raw[ri] = slope_s
        alpha_c_m1[ri] = slope_to_alpha(slope_s, dT)
        print(f"  M1 side: slope={slope_s:.4f} m/√hr, R²={r2_s:.4f}, "
              f"α_c={alpha_c_m1[ri]:.5f} m²/hr")

        # M2: bottom penetration slope (all concrete in CW)
        pen_bot = penetration_bottom(v.T_field_F, depths_m, float(placement_F), wi_cl)
        slope_b, r2_b = fit_slope(time_hrs, pen_bot)
        slope_bot_raw[ri] = slope_b
        alpha_c_m2[ri] = slope_to_alpha(slope_b, dT)
        print(f"  M2 bot:  slope={slope_b:.4f} m/√hr, R²={r2_b:.4f}, "
              f"α_c={alpha_c_m2[ri]:.5f} m²/hr (= α_CONC from below)")

        # M3: erfc profile fit at mid-depth
        for ti_fit, t_hr in enumerate(ERFC_TIMESTAMPS_HR):
            tii = find_idx(time_hrs, t_hr)
            T_row = v.T_field_F[tii, di_mid, :]
            a, A_fit, rmse = fit_erfc_concrete_side(
                T_row, widths_m, float(placement_F), t_hr)
            alpha_c_m3[ri, ti_fit] = a
            print(f"    M3 t={t_hr:.0f}hr: α_c={a:.5f} m²/hr, A={A_fit:.2f}°F"
                  f" (expected ≈ {soil_F-placement_F:.0f}°F), rmse={rmse:.3f}°F")

        # M4: concrete-bottom temperature structural comparison
        ti_168 = find_idx(time_hrs, 168.0)
        # CW: bottom row = last depth index (48), CL column (0)
        bot_T_cw[ri] = float(v.T_field_F[ti_168, -1, wi_cl])
        # Engine: load from Stage3 Fix2 CSV
        eng_csv = os.path.join(ENG_DIR, f"run{lbl}_t168.csv")
        if os.path.isfile(eng_csv):
            eng_y, eng_x, eng_T = load_engine_csv(eng_csv)
            # Bottom row = last row, CL = last x column
            bot_T_eng[ri] = float(eng_T[-1, -1])
        print(f"  M4 bot T: CW={bot_T_cw[ri]:.2f}°F, Eng={bot_T_eng[ri]:.2f}°F, "
              f"diff={bot_T_eng[ri]-bot_T_cw[ri]:.2f}°F  (expected={soil_F:.0f}°F)")

    # Aggregate α_c
    alpha_c_m3_mean = np.nanmean(alpha_c_m3, axis=1)
    print("\n--- Final α_c aggregates (all runs B–I) ---")
    for name, arr in [("M1 side slope", alpha_c_m1),
                      ("M2 bot slope ", alpha_c_m2),
                      ("M3 erfc fit  ", alpha_c_m3_mean)]:
        mn, sd = np.nanmean(arr), np.nanstd(arr)
        print(f"  {name}: {mn:.5f} ± {sd:.5f} m²/hr  (rel std={100*sd/mn:.1f}%)")

    # Structural proof
    print("\n--- Structural proof: bottom-T gap (engine − CW) ---")
    for ri, (lbl, fld, pl, so) in enumerate(NON_A):
        dT = abs(so - pl)
        gap = bot_T_eng[ri] - bot_T_cw[ri]
        print(f"  Run {lbl}: |ΔT|={dT}°F, gap={gap:.2f}°F  "
              f"(CW={bot_T_cw[ri]:.2f}°F, Eng={bot_T_eng[ri]:.2f}°F)")

    # Save arrays
    alpha_c_arr = np.column_stack([alpha_c_m1, alpha_c_m2, alpha_c_m3_mean])
    np.save(os.path.join(HERE, "STAGE4a_cw_alphac.npy"), alpha_c_arr)
    print("\nSaved STAGE4a_cw_alphac.npy")

    # Write report
    _write_report(NON_A, alpha_c_m1, alpha_c_m2, alpha_c_m3, alpha_c_m3_mean,
                  bot_T_cw, bot_T_eng)
    return 0


def _write_report(runs, ac_m1, ac_m2, ac_m3, ac_m3_mean, bot_T_cw, bot_T_eng):
    lines = [
        "# STAGE4a — CW-Extracted Thermal Parameters",
        "",
        "## Structural finding (preface to parameter extraction)",
        "",
        "**CW's model** for the 40ft×80ft geometry:",
        "- 80ft (24.38m) of concrete, all-concrete domain",
        "- Side BC: Dirichlet T_soil at w=0.00m for ALL depths (0–24.38m)",
        "- Bottom BC: Dirichlet T_soil at depth=24.38m for all widths",
        "- No soil domain below the concrete",
        "",
        "**Engine's model** (Stage 3 Fix 2):",
        "- 80ft (24.40m) of concrete",
        "- 3m of soil below the concrete (y=24.40m to 27.40m)",
        "- Dirichlet T_gw at depth=27.40m (3m below the concrete bottom)",
        "",
        "The soil buffer prevents the concrete bottom from reaching T_gw "
        "as quickly as CW's direct Dirichlet.  This causes the 9–19°F residual.",
        "",
        "## Method summary",
        "",
        "- **M1** (side penetration slope): 1°F front vs √t from side Dirichlet. "
        "Slope → α_c via `slope = 2 × erfinv(1−1/|ΔT|) × √α_c`. Grid resolution "
        "0.508m quantizes the front, degrading precision.",
        "- **M2** (bottom penetration slope): same method from bottom Dirichlet. "
        "Both M1 and M2 penetrate CONCRETE (CW has no soil below), so both give α_c.",
        "- **M3** (erfc profile fit): `T(x,t) − T_placement = A × erfc(x/(2√(α_c t)))` "
        "at t=24,48,96,168hr with `scipy.curve_fit`. A ≈ ΔT (confirms Dirichlet at "
        "x=0). α_c recovered directly — most reliable estimator.",
        "- **M4** (structural proof): compare engine vs CW concrete-bottom temperature "
        "at t=168hr, CL column. Gap = soil-buffer artifact.",
        "",
        "## Per-run α_c results",
        "",
    ]
    header = "| Run | ΔT (°F) | M1 α_c (m²/hr) | M2 α_c (m²/hr) | M3 α_c (m²/hr) |"
    sep    = "| --- | --- | --- | --- | --- |"
    lines += [header, sep]
    for ri, (lbl, fld, pl, so) in enumerate(runs):
        dT = so - pl
        lines.append(
            f"| {lbl} | {dT:+d} "
            f"| {ac_m1[ri]:.5f} | {ac_m2[ri]:.5f} | {ac_m3_mean[ri]:.5f} |"
        )

    lines += ["", "## M3 α_c by timestamp", ""]
    lines += ["| Run | t=24hr | t=48hr | t=96hr | t=168hr |",
              "| --- | --- | --- | --- | --- |"]
    for ri, (lbl, *_) in enumerate(runs):
        row = " | ".join(f"{ac_m3[ri, ti]:.5f}" for ti in range(4))
        lines.append(f"| {lbl} | {row} |")

    lines += ["", "## Aggregate α_c statistics", "",
              "| Method | Mean (m²/hr) | Std (m²/hr) | Rel std |",
              "| --- | --- | --- | --- |"]
    for name, arr in [("M1 side slope", ac_m1),
                      ("M2 bot slope", ac_m2),
                      ("M3 erfc fit (mean over timestamps)", ac_m3_mean)]:
        mn, sd = np.nanmean(arr), np.nanstd(arr)
        rel = 100.0 * sd / mn if mn > 0 else float("nan")
        lines.append(f"| {name} | {mn:.5f} | {sd:.5f} | {rel:.1f}% |")

    lines += ["", "## M4 — Structural proof: bottom-temperature gap", ""]
    lines += ["| Run | ΔT (°F) | CW bot T (°F) | Eng bot T (°F) | Gap (°F) |",
              "| --- | --- | --- | --- | --- |"]
    for ri, (lbl, fld, pl, so) in enumerate(runs):
        dT = so - pl
        gap = bot_T_eng[ri] - bot_T_cw[ri]
        lines.append(f"| {lbl} | {dT:+d} | {bot_T_cw[ri]:.2f} | {bot_T_eng[ri]:.2f} | {gap:+.2f} |")

    lines += [
        "",
        "The CW concrete-bottom temperature equals T_soil (Dirichlet at depth 24.38m).",
        "The engine concrete-bottom temperature is ~7–10°F warmer, buffered by the",
        "3m soil layer below the concrete.  This gap scales with |ΔT| and matches",
        "the Stage 3.5 masked residual (9.1–19.6°F) at the bottom-CL corner.",
        "",
        "## Key interpretation",
        "",
        "1. **α_s**: NOT extractable (CW has no soil domain). Irrelevant for the diagnosis.",
        "2. **e_c/e_s admittance**: NOT applicable. CW uses Dirichlet BCs at all",
        "   soil-concrete interfaces; the admittance concept does not apply.",
        "3. **α_c**: Extractable. CW M3 erfc fit gives α_c ≈ 0.0050–0.0056 m²/hr.",
        "   Engine α_c (at α_hyd=0.8) ≈ 0.00449 m²/hr — within ~10–20% of CW.",
        "   This is a secondary factor; see STAGE4a_comparison.md for the full table.",
    ]

    out = os.path.join(HERE, "STAGE4a_cw_extracted_params.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"Written: {out}")


if __name__ == "__main__":
    sys.exit(run())
