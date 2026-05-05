#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-pre — central side profile, F scenario.

Plots T(di=0..48, wi=0) at t=168hr for α_u ∈ {0.20, 0.40, 0.60, 0.80},
F scenario (T_pl=73°F, T_soil=45°F), engine vs CW overlaid.

Outputs:
  figures/central_side_profile_F_scenario.png
  data/central_side_profile_F_scenario.csv
"""

import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

HERE = Path(__file__).parent
SC2  = (HERE / "../stage2_calibration").resolve()
ROOT = (HERE / "../../..").resolve()
SC   = (HERE / "../../soil_calibration").resolve()

sys.path.insert(0, str(SC2))
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from engine_runner import run_and_compare
from cw_scenario_loader import parse_cw_temp_output

CW_DATA = SC2 / "cw_data"

RUNS = [
    ("alpha02_F", 0.20, CW_DATA / "thermal_alpha02_F_73_45"),
    ("alpha04_F", 0.40, CW_DATA / "thermal_alpha04_F_73_45"),
    ("alpha06_F", 0.60, CW_DATA / "thermal_alpha06_F_73_45"),
    ("alpha08_F", 0.80, CW_DATA / "thermal_alpha08_F_73_45"),
]

T_PL_F = 73.0
T_SO_F = 45.0
FACTOR = 0.96
WI     = 0   # side-face column

FIG_DIR  = HERE / "figures"
DATA_DIR = HERE / "data"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR.mkdir(parents=True, exist_ok=True)

COLORS = [cm.viridis(v) for v in [0.10, 0.37, 0.63, 0.90]]


def main():
    fig, (ax_main, ax_res) = plt.subplots(
        2, 1, figsize=(8, 9),
        gridspec_kw={"height_ratios": [3, 1]},
        sharex=True,
    )

    csv_rows = []

    print(f"\nRunning 4 engine solves (F scenario, factor={FACTOR})...")
    print(f"\n{'α_u':>5}  {'source':>8}  {'di=0':>8}  {'di=48':>8}")
    print("-" * 38)

    for (label, alpha_u, folder), color in zip(RUNS, COLORS):
        print(f"  {alpha_u:.2f}  {label} ...", end="", flush=True)
        rpt = run_and_compare(folder, T_PL_F, T_SO_F, factor=FACTOR)

        v   = parse_cw_temp_output(str(folder / "output.txt"))
        depths_m = v.depths_m   # shape (n_d,)
        n_di = len(depths_m)

        eng_col = rpt.eng_F[:n_di, WI]
        cw_col  = rpt.cw_F [:n_di, WI]
        res_col = rpt.R_field[:n_di, WI]

        ax_main.plot(depths_m, cw_col,  color=color, ls="-",  lw=1.8,
                     label=f"α_u={alpha_u:.2f} CW")
        ax_main.plot(depths_m, eng_col, color=color, ls="--", lw=1.2,
                     label=f"α_u={alpha_u:.2f} engine")
        ax_res.plot(depths_m, res_col, color=color, ls="-", lw=1.2)

        print(f"\r  {alpha_u:.2f}  CW      T(di=0)={cw_col[0]:7.2f}°F  T(di={n_di-1})={cw_col[-1]:7.2f}°F")
        print(f"       engine   T(di=0)={eng_col[0]:7.2f}°F  T(di={n_di-1})={eng_col[-1]:7.2f}°F")

        for di in range(n_di):
            d = float(depths_m[di])
            for source, col in [("cw", cw_col), ("engine", eng_col)]:
                csv_rows.append({
                    "di": di, "depth_m": f"{d:.4f}",
                    "alpha_u": alpha_u, "source": source,
                    "T_F": f"{col[di]:.4f}",
                })

    ax_main.set_ylabel("Temperature (°F)", fontsize=10)
    ax_main.set_title(
        "F scenario (T_pl=73°F, T_soil=45°F), t=168 hr — wi=0\n"
        "Solid=CW, dashed=engine",
        fontsize=10,
    )
    ax_main.legend(fontsize=7, ncol=2, loc="upper right")
    ax_main.grid(True, alpha=0.25)
    ax_main.tick_params(labelsize=8)

    ax_res.axhline(0.0, color="gray", lw=0.8, ls=":")
    ax_res.set_xlabel("Depth (m)", fontsize=10)
    ax_res.set_ylabel("engine − CW (°F)", fontsize=9)
    ax_res.grid(True, alpha=0.25)
    ax_res.tick_params(labelsize=8)

    fig.tight_layout()
    png_path = FIG_DIR / "central_side_profile_F_scenario.png"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nWrote: {png_path}")

    csv_path = DATA_DIR / "central_side_profile_F_scenario.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["di","depth_m","alpha_u","source","T_F"])
        w.writeheader()
        w.writerows(csv_rows)
    print(f"Wrote: {csv_path}  ({len(csv_rows)} rows)")


if __name__ == "__main__":
    main()
