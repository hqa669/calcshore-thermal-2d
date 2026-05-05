#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-Hu §4.6–§4.8 — floor localization, scaling, α_u dependence.

Reads data/Hu_residual_table.csv.

Writes:
  data/floor_localization.csv
  data/scaling_fits.csv
  data/alpha_u_dependence.csv

Run from validation/sprint8/stage2_diag_Hu/.
"""
import csv
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"


def load_residual_table() -> list[dict]:
    rows = []
    with open(DATA / "Hu_residual_table.csv", newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "folder":   r["folder"],
                "alpha_u":  float(r["alpha_u"]),
                "Hu_J_kg":  float(r["Hu_J_kg"]),
                "dT_168_C": float(r["dT_168_C"]),
                "dT_168_F": float(r["dT_168_F"]),
                "dT_max_C": float(r["dT_max_C"]),
            })
    return rows


def r_squared(y: np.ndarray, y_fit: np.ndarray) -> float:
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    if ss_tot == 0:
        return 1.0
    return float(1.0 - ss_res / ss_tot)


def floor_localize(subset: list[dict], flat_tol_C: float = 0.05) -> dict:
    """Identify the Hu floor for a single α_u subset (sorted by Hu ascending)."""
    subset = sorted(subset, key=lambda r: r["Hu_J_kg"])
    baseline_dT = subset[0]["dT_168_C"]

    # Noise floor = std within the flat region
    # Start with all points, shrink from top until we find where it lifts
    flat_pts = [r for r in subset if abs(r["dT_168_C"] - baseline_dT) <= flat_tol_C]
    sigma = np.std([r["dT_168_C"] for r in flat_pts]) if len(flat_pts) > 1 else 0.0
    threshold = max(2 * sigma, flat_tol_C)

    # Hu_floor_lower_bound = highest Hu still in flat region
    hu_floor_lower = None
    for r in subset:
        if abs(r["dT_168_C"] - baseline_dT) <= threshold:
            hu_floor_lower = r["Hu_J_kg"]
        else:
            break

    # Hu_floor_upper_bound = lowest Hu clearly above flat
    hu_floor_upper = None
    for r in subset:
        if r["dT_168_C"] - baseline_dT > threshold:
            hu_floor_upper = r["Hu_J_kg"]
            break

    if hu_floor_lower is None and hu_floor_upper is None:
        # entire range is flat — no floor detected
        return {
            "hu_floor_lower": None,
            "hu_floor_upper": None,
            "hu_floor_best_est": None,
            "dT_flat_region_C": baseline_dT,
            "dT_upper_bound_C": None,
            "note": "response flat across all Hu — no floor detected",
        }
    if hu_floor_upper is None:
        return {
            "hu_floor_lower": hu_floor_lower,
            "hu_floor_upper": None,
            "hu_floor_best_est": None,
            "dT_flat_region_C": baseline_dT,
            "dT_upper_bound_C": None,
            "note": "response flat through highest Hu — no upper floor bound",
        }
    if hu_floor_lower is None:
        hu_floor_lower = subset[0]["Hu_J_kg"]

    best_est = math.sqrt(hu_floor_lower * hu_floor_upper)
    dT_at_upper = next(r["dT_168_C"] for r in subset if r["Hu_J_kg"] == hu_floor_upper)

    return {
        "hu_floor_lower": hu_floor_lower,
        "hu_floor_upper": hu_floor_upper,
        "hu_floor_best_est": best_est,
        "dT_flat_region_C": baseline_dT,
        "dT_upper_bound_C": dT_at_upper,
        "note": "",
    }


def fit_scaling(subset: list[dict], hu_min: float) -> list[dict]:
    """Fit linear, log-linear, and power-law to the above-floor points."""
    above = [r for r in subset if r["Hu_J_kg"] >= hu_min]
    if len(above) < 2:
        return [{"form": f, "a": None, "b": None, "R2": None,
                 "n_points": len(above), "note": "< 2 points above floor"}
                for f in ("linear", "log_linear", "power_law")]

    hu  = np.array([r["Hu_J_kg"]  for r in above])
    dT  = np.array([r["dT_168_C"] for r in above])

    fits = []

    # Linear: dT = a*Hu + b
    if len(above) >= 2:
        coeffs = np.polyfit(hu, dT, 1)
        dT_fit = np.polyval(coeffs, hu)
        fits.append({"form": "linear", "a": coeffs[0], "b": coeffs[1],
                     "R2": r_squared(dT, dT_fit), "n_points": len(above), "note": ""})

    # Log-linear: dT = a*log10(Hu) + b
    log_hu = np.log10(hu)
    coeffs = np.polyfit(log_hu, dT, 1)
    dT_fit = np.polyval(coeffs, log_hu)
    fits.append({"form": "log_linear", "a": coeffs[0], "b": coeffs[1],
                 "R2": r_squared(dT, dT_fit), "n_points": len(above), "note": ""})

    # Power law: dT = a*Hu^b  (fit in log-log space; skip if any dT ≤ 0)
    if np.all(dT > 0):
        log_dT = np.log10(dT)
        coeffs = np.polyfit(log_hu, log_dT, 1)
        b_pow = coeffs[0]
        a_pow = 10 ** coeffs[1]
        dT_fit = a_pow * hu ** b_pow
        fits.append({"form": "power_law", "a": a_pow, "b": b_pow,
                     "R2": r_squared(dT, dT_fit), "n_points": len(above), "note": ""})
    else:
        fits.append({"form": "power_law", "a": None, "b": None, "R2": None,
                     "n_points": len(above), "note": "dT ≤ 0 in above-floor region"})

    return fits


def alpha_u_dependence(rows: list[dict]) -> list[dict]:
    au02 = {r["Hu_J_kg"]: r["dT_168_C"] for r in rows if abs(r["alpha_u"] - 0.20) < 0.05}
    au08 = {r["Hu_J_kg"]: r["dT_168_C"] for r in rows if abs(r["alpha_u"] - 0.80) < 0.05}
    shared_hu = sorted(set(au02) & set(au08))

    result = []
    for hu in shared_hu:
        dT02 = au02[hu]
        dT08 = au08[hu]
        ratio = dT08 / dT02 if abs(dT02) > 1e-6 else float("nan")
        if math.isnan(ratio):
            comment = "dT02 ≈ 0"
        elif abs(ratio - 4.0) < 0.5:
            comment = "≈4x → linear in α_u"
        elif abs(ratio - 1.0) < 0.2:
            comment = "≈1x → independent of α_u"
        else:
            comment = f"{ratio:.2f}x → nonlinear"
        result.append({"Hu": hu, "dT_au02_C": dT02, "dT_au08_C": dT08,
                       "ratio": ratio, "comment": comment})
    return result


def main():
    rows = load_residual_table()

    # §4.6 Floor localization
    floor_rows = []
    for au in (0.20, 0.80):
        subset = [r for r in rows if abs(r["alpha_u"] - au) < 0.05]
        result = floor_localize(subset)
        print(f"\n§4.6 Floor localization α_u={au:.2f}:")
        for k, v in result.items():
            print(f"  {k}: {v}")
        floor_rows.append({"alpha_u": au, **result})

    fl_path = DATA / "floor_localization.csv"
    fl_fields = ["alpha_u", "hu_floor_lower", "hu_floor_upper", "hu_floor_best_est",
                 "dT_flat_region_C", "dT_upper_bound_C", "note"]
    with open(fl_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fl_fields)
        w.writeheader()
        w.writerows(floor_rows)
    print(f"\nWrote {fl_path}")

    # §4.7 Above-floor scaling
    sf_rows = []
    for au in (0.20, 0.80):
        subset = [r for r in rows if abs(r["alpha_u"] - au) < 0.05]
        fl_info = next(r for r in floor_rows if abs(r["alpha_u"] - au) < 0.05)
        hu_min = fl_info["hu_floor_upper"] if fl_info["hu_floor_upper"] else 0.0
        fits = fit_scaling(subset, hu_min)
        print(f"\n§4.7 Scaling fits α_u={au:.2f} (Hu ≥ {hu_min:.0f} J/kg):")
        for ft in fits:
            print(f"  {ft['form']:<12}  a={ft.get('a', '?')}  b={ft.get('b', '?')}  "
                  f"R²={ft.get('R2', '?')}  n={ft['n_points']}")
        for ft in fits:
            sf_rows.append({"alpha_u": au, **ft})

    sf_path = DATA / "scaling_fits.csv"
    sf_fields = ["alpha_u", "form", "a", "b", "R2", "n_points", "note"]
    with open(sf_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=sf_fields)
        w.writeheader()
        w.writerows(sf_rows)
    print(f"\nWrote {sf_path}")

    # §4.8 α_u dependence
    dep_rows = alpha_u_dependence(rows)
    print("\n§4.8 α_u dependence (Hu shared between α_u=0.20 and 0.80):")
    print(f"  {'Hu':>9}  {'dT_au02°C':>10}  {'dT_au08°C':>10}  {'ratio':>6}  comment")
    for r in dep_rows:
        print(f"  {r['Hu']:>9.2f}  {r['dT_au02_C']:>10.4f}  {r['dT_au08_C']:>10.4f}  "
              f"{r['ratio']:>6.2f}  {r['comment']}")

    dep_path = DATA / "alpha_u_dependence.csv"
    dep_fields = ["Hu", "dT_au02_C", "dT_au08_C", "ratio", "comment"]
    with open(dep_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=dep_fields)
        w.writeheader()
        w.writerows(dep_rows)
    print(f"\nWrote {dep_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
