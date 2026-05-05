#!/usr/bin/env python3
"""Sprint 9 Stage 1-pilot sweep — residual surface visualization.

Reads sweep_residual_table.csv produced by run_engine_sweep.py and generates
8 heatmaps (2 mixes × 4 metrics) over the (Hu_factor_override, c_multiplier) grid.

Deliverables:
  figures/{mix}_heatmap_max_R_F.png    — max|R| over di∈[24,48]
  figures/{mix}_heatmap_R_di47_F.png   — R at dip cell di=47 (stencil-decoupling test)
  figures/{mix}_heatmap_R_di48_F.png   — R at Dirichlet-snap cell di=48
  figures/{mix}_heatmap_R_di36_F.png   — R at mid-section di=36

Usage:
    cd /Users/hqa668/calcshore-thermal-2d
    python validation/sprint9/stage1_pilot_sweep/scripts/analyze_sweep.py
"""

import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

HERE  = Path(__file__).resolve().parent
SWEEP = HERE.parent
ROOT  = (SWEEP / "../../..").resolve()

sys.path.insert(0, str(ROOT))

DATA_DIR = SWEEP / "data"
FIGS_DIR = SWEEP / "figures"
FIGS_DIR.mkdir(parents=True, exist_ok=True)

# Sweep grid — must match run_engine_sweep.py
HU_FACTOR_GRID = [0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95]
C_MULT_GRID    = [0.95, 0.97, 1.00, 1.03, 1.05]
GATE_F         = 0.5

# Nearest grid point to each mix's pilot-baseline Hu_factor
# MIX-01 pilot: Hu_factor=0.951403 → nearest grid point 0.95
# MIX-07 pilot: Hu_factor=0.894603 → nearest grid point 0.89
PILOT_HU_NEAREST = {"mix01": 0.95, "mix07": 0.89}
PILOT_C_MULT     = 1.00   # pilot used c_multiplier=1.00 (c_eff = 1×1.0532)

# (column label, colormap, divergent-about-zero)
METRICS = [
    ("max_R_F",  "max|R| over di∈[24,48] (°F)", "viridis", False),
    ("R_di47_F", "R(di=47, t=168 h) (°F)",       "RdBu_r",  True),
    ("R_di48_F", "R(di=48, t=168 h) (°F)",        "RdBu_r",  True),
    ("R_di36_F", "R(di=36, t=168 h) (°F)",        "RdBu_r",  True),
]


def _find_idx(lst, val):
    return int(np.argmin([abs(v - val) for v in lst]))


def load_sweep_table():
    path = DATA_DIR / "sweep_residual_table.csv"
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append({
                "mix":          r["mix"],
                "Hu_factor":    float(r["Hu_factor"]),
                "c_multiplier": float(r["c_multiplier"]),
                "max_R_F":      float(r["max_R_F"]),
                "di_at_max":    int(r["di_at_max"]),
                "R_di47_F":     float(r["R_di47_F"]),
                "R_di48_F":     float(r["R_di48_F"]),
                "R_di36_F":     float(r["R_di36_F"]),
                "pass":         r["pass"],
            })
    return rows


def pivot_surface(rows, mix_name, metric):
    """Build (n_Hu × n_c) grid; rows=Hu_factor, cols=c_multiplier."""
    grid = np.full((len(HU_FACTOR_GRID), len(C_MULT_GRID)), np.nan)
    for r in rows:
        if r["mix"] != mix_name:
            continue
        i = _find_idx(HU_FACTOR_GRID, r["Hu_factor"])
        j = _find_idx(C_MULT_GRID,    r["c_multiplier"])
        grid[i, j] = r[metric]
    return grid


def make_heatmap(mix_name, metric, label, cmap, divergent, rows):
    surface = pivot_surface(rows, mix_name, metric)

    fig, ax = plt.subplots(figsize=(7, 6))

    if divergent:
        finite = surface[np.isfinite(surface)]
        vmax   = max(abs(finite.min()), abs(finite.max())) if len(finite) else 1.0
        vmin   = -vmax
        im     = ax.imshow(surface, aspect="auto", cmap=cmap,
                           vmin=vmin, vmax=vmax, origin="lower")
    else:
        finite = surface[np.isfinite(surface)]
        im     = ax.imshow(surface, aspect="auto", cmap=cmap,
                           vmin=0, vmax=(finite.max() if len(finite) else 1.0),
                           origin="lower")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label, fontsize=9)

    # Cell annotations
    for i in range(len(HU_FACTOR_GRID)):
        for j in range(len(C_MULT_GRID)):
            val = surface[i, j]
            if not np.isfinite(val):
                ax.text(j, i, "N/A", ha="center", va="center", fontsize=7)
                continue
            bold = abs(val) >= GATE_F
            ax.text(j, i, f"{val:+.2f}",
                    ha="center", va="center", fontsize=7.5,
                    fontweight="bold" if bold else "normal")

    # Axis tick labels
    ax.set_xticks(range(len(C_MULT_GRID)))
    ax.set_xticklabels([f"{v:.2f}" for v in C_MULT_GRID], fontsize=8)
    ax.set_yticks(range(len(HU_FACTOR_GRID)))
    ax.set_yticklabels([f"{v:.2f}" for v in HU_FACTOR_GRID], fontsize=8)
    ax.set_xlabel("c_multiplier  (c_eff = c_multiplier × 1.0532)", fontsize=9)
    ax.set_ylabel("Hu_factor_override", fontsize=9)

    # Red border on the pilot-nearest grid cell
    pi = _find_idx(HU_FACTOR_GRID, PILOT_HU_NEAREST[mix_name])
    pj = _find_idx(C_MULT_GRID,    PILOT_C_MULT)
    ax.add_patch(plt.Rectangle(
        (pj - 0.5, pi - 0.5), 1.0, 1.0,
        fill=False, edgecolor="red", lw=2.0, zorder=5
    ))

    gate_note = f" | gate={GATE_F}°F (bold cells)" if not divergent else ""
    ax.set_title(
        f"{mix_name.upper()} — {label}\n"
        f"t=168 h, wi=0, di=[24,48]{gate_note}\n"
        f"red border = pilot-nearest grid point "
        f"(Hu_fac={PILOT_HU_NEAREST[mix_name]:.2f}, c_mult={PILOT_C_MULT:.2f})",
        fontsize=8.5,
    )

    fig.tight_layout()
    fname = f"{mix_name}_heatmap_{metric}.png"
    out   = FIGS_DIR / fname
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  Saved {out.name}")
    return out


def print_headline_stats(rows):
    print("\n" + "=" * 70)
    print("Headline diagnostic — R(di=47) across full 8×5 (Hu_factor, c_multiplier) grid:")
    for mix_name in ["mix01", "mix07"]:
        vals = np.array([r["R_di47_F"] for r in rows if r["mix"] == mix_name])
        print(f"  {mix_name.upper()}: mean={vals.mean():+.4f}°F  std={vals.std():.4f}°F  "
              f"range=[{vals.min():+.4f}, {vals.max():+.4f}]°F")
    print()
    print("Principal question: if std(R_di47) << gate (0.5°F) across the full grid,")
    print("R(di=47) is decoupled from (Hu_factor, c) — consistent with an INTERNAL-stencil")
    print("discretization mechanism inside the concrete domain near the bottom face.")
    print("(Both solvers apply T_soil=85°F as a direct Dirichlet at the bottom face under")
    print(" model_soil=False, so this is NOT a BC physics mismatch.)")
    print("=" * 70)


def main():
    print("Sprint 9 Stage 1-pilot sweep — residual surface analysis")

    rows = load_sweep_table()
    n_expected = len(HU_FACTOR_GRID) * len(C_MULT_GRID) * 2
    print(f"Loaded {len(rows)} rows from sweep_residual_table.csv  (expected {n_expected})")
    if len(rows) != n_expected:
        print(f"WARNING: row count mismatch — expected {n_expected}, got {len(rows)}")

    print_headline_stats(rows)

    print("\nGenerating heatmaps...")
    for mix_name in ["mix01", "mix07"]:
        print(f"\n  {mix_name.upper()}:")
        for metric, label, cmap, divergent in METRICS:
            make_heatmap(mix_name, metric, label, cmap, divergent, rows)

    print(f"\nAll figures written to {FIGS_DIR}/")


if __name__ == "__main__":
    main()
