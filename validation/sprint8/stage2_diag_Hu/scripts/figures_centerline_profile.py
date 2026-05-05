#!/usr/bin/env python3
"""Sprint 8 Stage 2-floor-test — centerline T(depth) profiles, A and I scenarios.

Plots T at the centerline (symmetry boundary) vs depth, from surface (top) to
bottom, at t = 0, 84, 168 hr, for α_u ∈ {0.20, 0.40, 0.60, 0.80}, A and I scenarios.

Engine: solid lines (Hu_J_kg_effective = 13,641.5 J/kg, k_uc=0.96).
CW: dashed lines (native dataset Hu = 1 J/kg).

Engine centerline = T_field[:, -1] (x=6.096m, symmetry boundary in half-mat).
CW centerline    = T_field[:, 0]  (widths_m[0]=6.1m).

Reads:
  data/T_fields_engine_Hu_residual/{alpha0X_A,alpha0X_I}.npz
  ../stage2_calibration/cw_data/thermal_alpha0X_{A_73_73,I_100_73}/output.txt

Writes:
  figures/centerline_T_vs_depth_A_I.png
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))
from cw_scenario_loader import parse_cw_temp_output

CW_DATA_BASE = REPO / "validation" / "sprint8" / "stage2_calibration" / "cw_data"
DATA  = HERE / "data"
TFLD  = DATA / "T_fields_engine_Hu_residual"
FIGS  = HERE / "figures"

ALPHAS = [(0.20, "alpha02"), (0.40, "alpha04"), (0.60, "alpha06"), (0.80, "alpha08")]
SCENARIOS = [
    ("A", "73_73",  73, 73),
    ("I", "100_73", 100, 73),
]
SNAPS = (0, 84, 168)
COLORS = {0: "#444444", 84: "#d6604d", 168: "#2166ac"}


def load_engine(label: str):
    d = np.load(TFLD / f"{label}.npz")
    depths = d["depths_m"]
    # Engine centerline is the LAST x index (x = 6.096m, symmetry boundary)
    T0   = None  # synthesize from initial condition (uniform T_pl)
    T84  = d["T_C_t84hr"][:, -1]
    T168 = d["T_C_t168hr"][:, -1]
    return depths, T0, T84, T168


def load_cw(folder_name: str):
    v = parse_cw_temp_output(str(CW_DATA_BASE / folder_name / "output.txt"))
    depths = v.depths_m
    # CW centerline is wi=0 (widths_m[0] = 6.1 m, the symmetry side)
    T_F = v.T_field_F[:, :, 0]
    T_C = (T_F - 32.0) * 5.0 / 9.0
    t   = v.time_hrs
    def at_hr(hr):
        if hr == 0:
            return T_C[0]                       # CW first sample is t≈0.083 hr; use as t=0 proxy
        ti = int(np.abs(t - hr).argmin())
        return T_C[ti]
    return depths, at_hr(0), at_hr(84), at_hr(168)


def main():
    FIGS.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 4, figsize=(16, 11), sharey=True)
    # Full slab depth: 0 (top surface) → 24.384 m (bottom surface, 80 ft).
    # Both surfaces have boundary effects; the bulk in between holds T_pl plus
    # internal hydration warming.

    for col, (au_val, au_tag) in enumerate(ALPHAS):
        for row, (sc_tag, sc_suffix, T_pl_F, T_so_F) in enumerate(SCENARIOS):
            ax = axes[row, col]
            label   = f"{au_tag}_{sc_tag}"
            cw_dir  = f"thermal_{au_tag}_{sc_tag}_{sc_suffix}"

            T_pl_C = (T_pl_F - 32.0) * 5.0 / 9.0

            eng_d, _, eng_84, eng_168  = load_engine(label)
            cw_d,  cw_0, cw_84, cw_168 = load_cw(cw_dir)

            # t=0 — engine IC is uniform T_pl across concrete grid
            eng_0 = np.full_like(eng_d, T_pl_C)

            # Plot t=0 (overlay; engine and CW agree by construction)
            ax.plot(eng_0,  eng_d, color=COLORS[0],   linestyle="-",  linewidth=1.4,
                    label="t=0 (T_pl, engine & CW)")
            ax.plot(cw_0,   cw_d,  color=COLORS[0],   linestyle=":",  linewidth=1.0,
                    alpha=0.7)

            # t=84
            ax.plot(eng_84, eng_d, color=COLORS[84],  linestyle="-",  linewidth=1.5,
                    label="t=84 hr engine")
            ax.plot(cw_84,  cw_d,  color=COLORS[84],  linestyle="--", linewidth=1.5,
                    label="t=84 hr CW")

            # t=168
            ax.plot(eng_168, eng_d, color=COLORS[168], linestyle="-",  linewidth=1.5,
                    label="t=168 hr engine")
            ax.plot(cw_168,  cw_d,  color=COLORS[168], linestyle="--", linewidth=1.5,
                    label="t=168 hr CW")

            ax.invert_yaxis()   # surface at top, bottom at bottom
            ax.grid(True, linestyle="--", alpha=0.4)
            ax.set_title(f"{sc_tag} scenario  α_u={au_val:.2f}\n"
                         f"T_pl={T_pl_F}°F  T_so={T_so_F}°F",
                         fontsize=10)

            if col == 0:
                ax.set_ylabel("Depth (m)  —  surface at top, bottom at bottom",
                              fontsize=9)
            if row == 1:
                ax.set_xlabel("T (°C)", fontsize=10)

            if row == 0 and col == 3:
                ax.legend(fontsize=7, loc="lower right")

    fig.suptitle("Centerline T(depth) profiles — A and I scenarios at Hu=Hu_residual "
                 "(13,641.5 J/kg, k_uc=0.96). Engine: solid. CW: dashed. "
                 "Full slab depth: surface (top) to bottom (24.4 m = 80 ft).",
                 fontsize=10)
    fig.tight_layout()
    out = FIGS / "centerline_T_vs_depth_A_I.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Wrote {out}")

    # Also dump a CSV of the depth profiles for reference
    out_csv = DATA / "centerline_T_vs_depth_A_I.csv"
    rows = []
    for au_val, au_tag in ALPHAS:
        for sc_tag, sc_suffix, T_pl_F, T_so_F in SCENARIOS:
            label  = f"{au_tag}_{sc_tag}"
            cw_dir = f"thermal_{au_tag}_{sc_tag}_{sc_suffix}"
            T_pl_C = (T_pl_F - 32.0) * 5.0 / 9.0
            eng_d, _, eng_84, eng_168  = load_engine(label)
            cw_d,  cw_0, cw_84, cw_168 = load_cw(cw_dir)
            for k, d in enumerate(eng_d):
                rows.append({
                    "alpha_u": au_val, "scenario": sc_tag,
                    "source": "engine", "depth_m": float(d),
                    "T_C_t0":   float(T_pl_C),
                    "T_C_t84":  float(eng_84[k]),
                    "T_C_t168": float(eng_168[k]),
                })
            for k, d in enumerate(cw_d):
                rows.append({
                    "alpha_u": au_val, "scenario": sc_tag,
                    "source": "cw", "depth_m": float(d),
                    "T_C_t0":   float(cw_0[k]),
                    "T_C_t84":  float(cw_84[k]),
                    "T_C_t168": float(cw_168[k]),
                })
    import csv as _csv
    with open(out_csv, "w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=["alpha_u","scenario","source",
                                            "depth_m","T_C_t0","T_C_t84","T_C_t168"])
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
