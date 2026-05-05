#!/usr/bin/env python3
"""Sprint 8 Stage 2 — figures 10, 11, 12.

Run after stage2_final_validation.py has produced findings/final_validation.csv
and findings/kuc_sweep.csv.

Figure 10: stage2_calibration_factor_vs_alpha.png
  Scatter of optimal k_uc factor vs α target (5 points).
  If calibration_decision.md names a curve fit, overlays it.

Figure 11: stage2_residuals_per_alpha.png
  Grouped bar chart: R1 max|R| and R2 max|R| at each α target,
  before (factor=0.96) and after (optimal factor) calibration.

Figure 12: stage2_R1_R2_field_at_optimal.png
  Residual field heatmaps for the single most-stressed scenario at
  each α target (the run with the highest R2 max|R| at optimal factor).

Output directory: /mnt/user-data/outputs/ (falls back to findings/ if not found)
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE     = Path(__file__).parent
ROOT     = (HERE / "../../..").resolve()
SC       = (HERE / "../../soil_calibration").resolve()
FINDINGS = HERE / "findings"

OUT_DIR  = Path("/mnt/user-data/outputs")
if not OUT_DIR.exists():
    OUT_DIR = FINDINGS
    print(f"Note: /mnt/user-data/outputs/ not found, saving figures to {FINDINGS}")

sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from engine_runner import run_and_compare, GATE_F
CW_DATA = HERE / "cw_data"


# ── Helpers ──────────────────────────────────────────────────────────────

def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def parse_float(v):
    return float(v)


def load_optimal_factors():
    """Returns dict {alpha_u: factor} from kuc_sweep_optimal*.txt."""
    for fname in ("kuc_sweep_optimal_refined.txt", "kuc_sweep_optimal.txt"):
        p = FINDINGS / fname
        if p.exists():
            result = {0.036: 0.96}
            for line in p.read_text().splitlines():
                if line.startswith("#") or "anchor" in line or "alpha_u" in line or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        result[float(parts[0])] = float(parts[1])
                    except ValueError:
                        pass
            return result
    return {0.036: 0.96}


# ── Figure 10 ─────────────────────────────────────────────────────────────

def make_fig10(optima):
    alphas  = sorted(optima.keys())
    factors = [optima[a] for a in alphas]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.scatter(alphas, factors, s=90, zorder=5, color="steelblue", label="Optimal k_uc factor")
    ax.axhline(0.96, color="gray", ls="--", lw=1.2, label="Sprint 7 factor (0.96)")
    ax.axhline(np.mean(factors), color="darkblue", ls=":", lw=1.2,
               label=f"Mean = {np.mean(factors):.4f}")

    spread = max(factors) - min(factors)
    outcome = "A" if spread <= 0.01 else ("A-marginal" if spread <= 0.04 else "B")
    ax.set_xlabel("α_u (ultimate degree of hydration)", fontsize=10)
    ax.set_ylabel("Optimal k_uc multiplier", fontsize=10)
    ax.set_title(
        f"Sprint 8 Stage 2 — Optimal k_uc factor vs α target\n"
        f"Outcome {outcome}: spread = {spread:.4f}",
        fontsize=10,
    )
    ax.set_xlim(-0.02, 0.90)
    ax.set_ylim(min(factors) - 0.02, max(factors) + 0.02)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.25)
    ax.tick_params(labelsize=8)

    path = OUT_DIR / "stage2_calibration_factor_vs_alpha.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {path}")


# ── Figure 11 ─────────────────────────────────────────────────────────────

def make_fig11(optima):
    baseline_data = {}   # alpha_u → {"R1": max, "R2": max}
    optimal_data  = {}

    sweep_path = FINDINGS / "kuc_sweep.csv"
    baseline_path = FINDINGS / "baseline_residuals.csv"
    final_path    = FINDINGS / "final_validation.csv"

    if not baseline_path.exists() or not final_path.exists():
        print("Skipping Figure 11 — baseline_residuals.csv or final_validation.csv not found")
        return

    # Baseline
    for row in load_csv(baseline_path):
        au = parse_float(row["alpha_u"])
        if au not in baseline_data:
            baseline_data[au] = {"R1": [], "R2": []}
        baseline_data[au]["R1"].append(parse_float(row["r1_max_F"]))
        baseline_data[au]["R2"].append(parse_float(row["r2_max_F"]))

    # Final (optimal factor, Sprint 8 runs only)
    for row in load_csv(final_path):
        au = parse_float(row["alpha_u"])
        if abs(au - 0.036) < 0.001:
            continue
        if au not in optimal_data:
            optimal_data[au] = {"R1": [], "R2": []}
        optimal_data[au]["R1"].append(parse_float(row["r1_max_F"]))
        optimal_data[au]["R2"].append(parse_float(row["r2_max_F"]))

    aus   = sorted(baseline_data.keys())
    x     = np.arange(len(aus))
    width = 0.18

    fig, ax = plt.subplots(figsize=(10, 5))
    for i, (au, col, prefix) in enumerate([
        (None, "steelblue",  "Baseline (0.96) R1"),
        (None, "royalblue",  "Baseline (0.96) R2"),
        (None, "darkorange", "Optimal  R1"),
        (None, "orangered",  "Optimal  R2"),
    ]):
        pass  # will draw below

    base_r1 = [max(baseline_data[a]["R1"]) if a in baseline_data else 0.0 for a in aus]
    base_r2 = [max(baseline_data[a]["R2"]) if a in baseline_data else 0.0 for a in aus]
    opt_r1  = [max(optimal_data[a]["R1"])  if a in optimal_data  else 0.0 for a in aus]
    opt_r2  = [max(optimal_data[a]["R2"])  if a in optimal_data  else 0.0 for a in aus]

    ax.bar(x - 1.5*width, base_r1, width, label="Baseline R1 (factor=0.96)", color="steelblue", alpha=0.85)
    ax.bar(x - 0.5*width, base_r2, width, label="Baseline R2 (factor=0.96)", color="royalblue",  alpha=0.85)
    ax.bar(x + 0.5*width, opt_r1,  width, label="Optimal R1",                 color="darkorange", alpha=0.85)
    ax.bar(x + 1.5*width, opt_r2,  width, label="Optimal R2",                 color="orangered",  alpha=0.85)

    ax.axhline(GATE_F, color="green", ls="--", lw=1.5, label=f"Gate: {GATE_F}°F")
    ax.set_xticks(x)
    ax.set_xticklabels([f"α_u={a:.2f}" for a in aus], fontsize=8)
    ax.set_ylabel("max|R| (°F)", fontsize=10)
    ax.set_title("Sprint 8 Stage 2 — R1/R2 residuals before and after calibration", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.2, axis="y")
    ax.tick_params(labelsize=8)

    path = OUT_DIR / "stage2_residuals_per_alpha.png"
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {path}")


# ── Figure 12 ─────────────────────────────────────────────────────────────

def make_fig12(optima):
    """Residual field heatmaps for the worst-case run at each α target."""
    final_path = FINDINGS / "final_validation.csv"
    if not final_path.exists():
        print("Skipping Figure 12 — final_validation.csv not found")
        return

    from cw_scenario_loader import parse_cw_temp_output

    rows = load_csv(final_path)
    # Find worst-case run per α_u (highest R2 max|R|)
    by_au = {}
    for row in rows:
        au = parse_float(row["alpha_u"])
        if abs(au - 0.036) < 0.001:
            continue
        if au not in by_au or parse_float(row["r2_max_F"]) > parse_float(by_au[au]["r2_max_F"]):
            by_au[au] = row

    aus_plot = sorted(by_au.keys())[:4]
    if not aus_plot:
        print("Skipping Figure 12 — no Sprint 8 rows in final_validation.csv")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, au in enumerate(aus_plot):
        row  = by_au[au]
        fld  = Path(row["folder"])
        T_pl = parse_float(row["T_pl_F"])
        T_so = parse_float(row["T_so_F"])
        fac  = parse_float(row["factor"])
        r1   = parse_float(row["r1_max_F"])
        r2   = parse_float(row["r2_max_F"])

        rpt  = run_and_compare(fld, T_pl, T_so, factor=fac)
        R    = rpt.R_field

        v    = parse_cw_temp_output(str(fld / "output.txt"))
        ti   = int(np.abs(v.time_hrs - 168.0).argmin())
        d_m  = v.depths_m
        w_m  = v.widths_m

        extent = [float(w_m[-1]), float(w_m[0]), float(d_m[-1]), 0.0]
        sym = max(0.01, float(np.abs(R).max()))

        ax = axes[idx]
        im = ax.imshow(R[:, ::-1], aspect="auto", origin="upper", extent=extent,
                       cmap="RdBu_r", vmin=-sym, vmax=sym, interpolation="nearest")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04).ax.tick_params(labelsize=6)
        verdict = "PASS" if rpt.passes else "FAIL"
        ax.set_title(
            f"α_u={au:.2f}  {row['label']}  (T_pl={T_pl:.0f}°F, T_so={T_so:.0f}°F)\n"
            f"factor={fac:.4f}  R1={r1:.4f}  R2={r2:.4f}  [{verdict}]",
            fontsize=8,
        )
        ax.set_xlabel("Width (m)", fontsize=7)
        ax.set_ylabel("Depth (m)", fontsize=7)
        ax.tick_params(labelsize=6)

    fig.suptitle("Sprint 8 Stage 2 — Residual fields at optimal calibration (t=168 hr)\n"
                 "Worst-case scenario per α target", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    path = OUT_DIR / "stage2_R1_R2_field_at_optimal.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote: {path}")


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINDINGS.mkdir(parents=True, exist_ok=True)

    optima = load_optimal_factors()
    print(f"\nGenerating Stage 2 figures...")
    make_fig10(optima)
    make_fig11(optima)
    make_fig12(optima)
    print("\nDone.")


if __name__ == "__main__":
    main()
