#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.6 + §3.7 — calibration decision and curve fit.

Reads optimal factors from findings/kuc_sweep_optimal.txt (or _refined.txt
if present), adds the Sprint 7 anchor at α=0.036, then:

§3.6  Computes spread and classifies Outcome A / A-marginal / B.
§3.7  If Outcome B (spread > 0.04), fits linear, quadratic, and modified
      Van Breugel functions to the 5 (α, factor) points.

Output: findings/calibration_decision.md
"""

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).parent
FINDINGS = HERE / "findings"


# ── Load optimal factors ─────────────────────────────────────────────────

def load_optima() -> dict:
    # Prefer refined file
    for fname in ("kuc_sweep_optimal_refined.txt", "kuc_sweep_optimal.txt"):
        p = FINDINGS / fname
        if p.exists():
            print(f"Loading optima from: {p.name}")
            break
    else:
        print("ERROR: no kuc_sweep_optimal*.txt found. Run stage2_kuc_sweep.py first.")
        sys.exit(1)

    result = {}
    for line in p.read_text().splitlines():
        if line.startswith("#") or not line.strip() or "anchor" in line or "alpha_u" in line:
            continue
        parts = line.split()
        if len(parts) >= 2:
            try:
                au = float(parts[0])
                bf = float(parts[1])
                result[au] = bf
            except ValueError:
                pass
    return result


# ── Curve-fit models ─────────────────────────────────────────────────────

def linear_model(alpha, a, b):
    return a + b * alpha


def quadratic_model(alpha, a, b, c):
    return a + b * alpha + c * alpha ** 2


def van_breugel_model(alpha, a, b):
    # factor(α) = a * (1.33 − b·α) / (1.33 − 0.33·α)
    # Original k(α) = k_uc * (1.33 − 0.33·α); calibrated slope factor replaces 0.33 → b
    return a * (1.33 - b * alpha) / (1.33 - 0.33 * alpha)


def fit_model(model, alphas, factors, p0, name):
    try:
        popt, _ = curve_fit(model, alphas, factors, p0=p0, maxfev=5000)
        predicted = model(alphas, *popt)
        residuals = factors - predicted
        ss_res = np.sum(residuals ** 2)
        ss_tot = np.sum((factors - np.mean(factors)) ** 2)
        r2     = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else 1.0
        max_err = float(np.max(np.abs(residuals)))
        return {
            "name": name, "params": popt, "r2": r2,
            "max_err": max_err, "residuals": residuals, "ok": True,
        }
    except Exception as e:
        return {"name": name, "ok": False, "error": str(e)}


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    FINDINGS.mkdir(parents=True, exist_ok=True)

    optima = load_optima()
    # Add Sprint 7 anchor
    optima[0.036] = 0.96

    alphas_all  = sorted(optima)
    factors_all = np.array([optima[a] for a in alphas_all])

    print(f"\n§3.6 Optimal factor vs α — 5-point dataset")
    print("-" * 40)
    print(f"  {'α_u':>6}  {'factor':>8}  {'source'}")
    for a, f in zip(alphas_all, factors_all):
        src = "Sprint 7" if abs(a - 0.036) < 0.001 else "Sprint 8 Stage 2"
        print(f"  {a:>6.3f}  {f:>8.4f}  {src}")

    spread = float(np.max(factors_all) - np.min(factors_all))
    mean_f = float(np.mean(factors_all))
    std_f  = float(np.std(factors_all))
    print(f"\n  Range   = {spread:.4f}")
    print(f"  Mean    = {mean_f:.4f}")
    print(f"  Std dev = {std_f:.4f}")

    if spread <= 0.01:
        outcome = "A"
        label   = "Constant factor (variation ≤ 0.01 — within noise)"
    elif spread <= 0.04:
        outcome = "A-marginal"
        label   = "Weak α-dependent trend (0.01 < spread ≤ 0.04)"
    else:
        outcome = "B"
        label   = "α-dependent function required (spread > 0.04)"

    print(f"\n§3.6 DECISION: Outcome {outcome} — {label}")
    print(f"  Recommended constant:  {mean_f:.4f}  (or 0.96 if spread ≤ 0.01)")

    md_lines = [
        "# Sprint 8 Stage 2 §3.6 — Calibration Decision\n",
        "## Optimal factor table\n",
        "| α_u | Optimal factor | Source |",
        "|---|---|---|",
    ]
    for a, f in zip(alphas_all, factors_all):
        src = "Sprint 7" if abs(a - 0.036) < 0.001 else "Sprint 8 Stage 2"
        md_lines.append(f"| {a:.3f} | {f:.4f} | {src} |")

    md_lines += [
        f"\n## Spread analysis",
        f"- Range: {spread:.4f}",
        f"- Mean ± std: {mean_f:.4f} ± {std_f:.4f}",
        f"\n## Outcome",
        f"**{outcome}: {label}**",
    ]

    # §3.7 curve fit (only if Outcome B)
    if outcome == "B":
        alphas_arr  = np.array(alphas_all)
        factors_arr = factors_all

        print(f"\n§3.7 Curve fit (Outcome B)")
        fits = [
            fit_model(linear_model,     alphas_arr, factors_arr, [0.96, 0.0],     "Linear:        factor = a + b·α"),
            fit_model(quadratic_model,  alphas_arr, factors_arr, [0.96, 0.0, 0.0],"Quadratic:     factor = a + b·α + c·α²"),
            fit_model(van_breugel_model,alphas_arr, factors_arr, [0.96, 0.33],     "Van Breugel:   factor = a·(1.33−b·α)/(1.33−0.33·α)"),
        ]

        md_lines += ["\n## §3.7 Curve fit results\n"]

        for fit in fits:
            if not fit["ok"]:
                print(f"  {fit['name']}: FAILED — {fit['error']}")
                md_lines.append(f"- **{fit['name']}**: FAILED — {fit['error']}")
                continue
            r2      = fit["r2"]
            max_err = fit["max_err"]
            params  = fit["params"]
            param_str = ", ".join(f"{p:.5f}" for p in params)
            fits_ok = max_err <= 0.005
            flag    = "✓ (fits ±0.005)" if fits_ok else f"✗ (max err={max_err:.4f} > 0.005)"
            print(f"  {fit['name']}")
            print(f"    params=[{param_str}]  R²={r2:.5f}  max_err={max_err:.5f}  {flag}")
            for a, f_actual, f_pred in zip(alphas_all, factors_arr, fit["params"][0:0] or []):
                pass
            per_pt = "  ".join(
                f"α={a:.3f}→{e:+.4f}" for a, e in zip(alphas_all, fit["residuals"])
            )
            print(f"    per-point residuals: {per_pt}")
            md_lines.append(f"\n### {fit['name']}")
            md_lines.append(f"- params: [{param_str}]")
            md_lines.append(f"- R²: {r2:.5f}")
            md_lines.append(f"- max_err: {max_err:.5f}  {flag}")
            md_lines.append(f"- per-point: {per_pt}")

        # Recommend simplest adequate fit
        adequate = [ft for ft in fits if ft["ok"] and ft["max_err"] <= 0.005]
        if adequate:
            best_fit = adequate[0]
            print(f"\n  → Recommended: {best_fit['name']} (simplest adequate fit)")
            md_lines.append(f"\n**Recommended fit: {best_fit['name']}**")
        else:
            print(f"\n  → No fit satisfies ±0.005 at all 5 points — use best available")
            best_fit = min([ft for ft in fits if ft["ok"]], key=lambda x: x["max_err"]) if any(ft["ok"] for ft in fits) else None
            if best_fit:
                md_lines.append(f"\n**Best available fit: {best_fit['name']} (max_err={best_fit['max_err']:.5f})**")
    else:
        print(f"\n  § 3.7 skipped (Outcome {outcome} — no curve fit needed)")
        md_lines.append("\n## §3.7 Curve fit\nSkipped — Outcome is not B.")

    md_path = FINDINGS / "calibration_decision.md"
    md_path.write_text("\n".join(md_lines) + "\n")
    print(f"\nWrote: {md_path}")


if __name__ == "__main__":
    main()
