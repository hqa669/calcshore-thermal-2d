#!/usr/bin/env python3
"""Sprint 8 follow-up — centerline T(depth) profiles at T_ref = 21°C.

Re-runs the 4 A and 4 I datasets with te2d.T_REF_K = 294.15 (21°C) instead
of 296.15 (23°C), saves T-field snapshots at 84 and 168 hr, and produces
a comparison plot with three curve sets per panel:

  - engine, T_ref=23°C (solid)        — the §4.10 baseline (loaded from existing npz files)
  - engine, T_ref=21°C (dotted)       — the hypothesis-test re-runs
  - CW reference                      (dashed)

Hu_J_kg_effective = 13,641.5 J/kg, k_uc multiplier = 0.96.
Full slab depth view: surface (top) → bottom (24.4 m).

Writes:
  data/T_fields_engine_Tref21/<label>.npz
  figures/centerline_T_vs_depth_A_I_Tref21.png
  data/centerline_T_vs_depth_A_I_Tref21.csv
"""
import csv as _csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE  = Path(__file__).resolve().parents[1]
REPO  = Path(__file__).resolve().parents[4]
STAGE = REPO / "validation" / "sprint8" / "stage2_calibration"
SC    = REPO / "validation" / "soil_calibration"

sys.path.insert(0, str(REPO))
sys.path.insert(0, str(STAGE))
sys.path.insert(0, str(SC))
sys.path.insert(0, str(HERE / "scripts"))

import thermal_engine_2d as te2d
from cw_scenario_loader import parse_cw_temp_output
from stage4b_run import nearest_time_idx
from run_engine_Hu_residual import run_one, HU_RESIDUAL_J_KG, K_UC_FACTOR

DATA          = HERE / "data"
TFLD_REF23    = DATA / "T_fields_engine_Hu_residual"
TFLD_REF21    = DATA / "T_fields_engine_Tref21"
FIGS          = HERE / "figures"
CW_DATA_BASE  = STAGE / "cw_data"

T_REF_K_NEW = 294.15
SNAP_HRS    = (24.0, 84.0, 168.0)

ALPHAS    = [(0.20, "alpha02"), (0.40, "alpha04"),
             (0.60, "alpha06"), (0.80, "alpha08")]
SCENARIOS = [
    ("A", "73_73",  73, 73),
    ("I", "100_73", 100, 73),
]
COLORS = {0: "#444444", 84: "#d6604d", 168: "#2166ac"}


def rerun_save_Tref21():
    TFLD_REF21.mkdir(parents=True, exist_ok=True)
    print(f"\nRe-running 8 A/I datasets at T_REF_K = {T_REF_K_NEW} (21°C)")
    print(f"  Hu_J_kg_effective = {HU_RESIDUAL_J_KG} J/kg, k_uc = {K_UC_FACTOR}")
    original = te2d.T_REF_K
    try:
        te2d.T_REF_K = T_REF_K_NEW
        for au_val, au_tag in ALPHAS:
            for sc_tag, sc_suffix, T_pl, T_so in SCENARIOS:
                label = f"{au_tag}_{sc_tag}"
                folder = f"thermal_{au_tag}_{sc_tag}_{sc_suffix}"
                out = TFLD_REF21 / f"{label}.npz"
                if out.exists():
                    print(f"  {label}: already saved, skipping")
                    continue
                print(f"  {label} ...", end="", flush=True)
                d = run_one(CW_DATA_BASE / folder, T_pl, T_so, factor=K_UC_FACTOR)
                jsl, isl = d["grid"].concrete_slice()
                T_field_C = d["T_field_C"]
                t_s       = d["t_s"]
                snaps = {}
                for hr in SNAP_HRS:
                    ti = nearest_time_idx(t_s, hr)
                    snaps[f"T_C_t{int(hr)}hr"] = T_field_C[ti, jsl, isl]
                snaps["depths_m"]     = d["grid"].y[jsl]
                snaps["widths_m"]     = d["grid"].x[isl]
                snaps["snapshot_hrs"] = np.array(SNAP_HRS)
                np.savez(out, **snaps)
                print(" saved")
    finally:
        te2d.T_REF_K = original


def load_engine_npz(npz_path: Path):
    d = np.load(npz_path)
    depths = d["depths_m"]
    T84  = d["T_C_t84hr"][:, -1]      # centerline = last x index
    T168 = d["T_C_t168hr"][:, -1]
    return depths, T84, T168


def load_cw(folder_name: str):
    v = parse_cw_temp_output(str(CW_DATA_BASE / folder_name / "output.txt"))
    depths = v.depths_m
    T_F = v.T_field_F[:, :, 0]         # centerline = wi=0
    T_C = (T_F - 32.0) * 5.0 / 9.0
    t   = v.time_hrs
    def at_hr(hr):
        if hr == 0:
            return T_C[0]
        ti = int(np.abs(t - hr).argmin())
        return T_C[ti]
    return depths, at_hr(0), at_hr(84), at_hr(168)


def build_figure():
    FIGS.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 4, figsize=(16, 11), sharey=True)
    rows_csv = []

    for col, (au_val, au_tag) in enumerate(ALPHAS):
        for row_i, (sc_tag, sc_suffix, T_pl_F, T_so_F) in enumerate(SCENARIOS):
            ax = axes[row_i, col]
            label  = f"{au_tag}_{sc_tag}"
            cw_dir = f"thermal_{au_tag}_{sc_tag}_{sc_suffix}"
            T_pl_C = (T_pl_F - 32.0) * 5.0 / 9.0

            d23, e23_84, e23_168 = load_engine_npz(TFLD_REF23 / f"{label}.npz")
            d21, e21_84, e21_168 = load_engine_npz(TFLD_REF21 / f"{label}.npz")
            cw_d, cw_0, cw_84, cw_168 = load_cw(cw_dir)

            # t=0 — uniform IC (engine + CW)
            ax.plot(np.full_like(d23, T_pl_C), d23,
                    color=COLORS[0], linestyle="-", linewidth=1.2,
                    label="t=0 (T_pl, engine & CW)")

            # t=84
            ax.plot(e23_84, d23, color=COLORS[84], linestyle="-",  linewidth=1.4,
                    label="t=84 hr engine (T_ref=23°C)")
            ax.plot(e21_84, d21, color=COLORS[84], linestyle=":",  linewidth=1.6,
                    label="t=84 hr engine (T_ref=21°C)")
            ax.plot(cw_84,  cw_d, color=COLORS[84], linestyle="--", linewidth=1.4,
                    label="t=84 hr CW")

            # t=168
            ax.plot(e23_168, d23, color=COLORS[168], linestyle="-",  linewidth=1.4,
                    label="t=168 hr engine (T_ref=23°C)")
            ax.plot(e21_168, d21, color=COLORS[168], linestyle=":",  linewidth=1.6,
                    label="t=168 hr engine (T_ref=21°C)")
            ax.plot(cw_168,  cw_d, color=COLORS[168], linestyle="--", linewidth=1.4,
                    label="t=168 hr CW")

            ax.invert_yaxis()
            ax.grid(True, linestyle="--", alpha=0.4)
            ax.set_title(f"{sc_tag} scenario  α_u={au_val:.2f}\n"
                         f"T_pl={T_pl_F}°F  T_so={T_so_F}°F",
                         fontsize=10)

            if col == 0:
                ax.set_ylabel("Depth (m)  —  surface at top, bottom at bottom",
                              fontsize=9)
            if row_i == 1:
                ax.set_xlabel("T (°C)", fontsize=10)
            if row_i == 0 and col == 3:
                ax.legend(fontsize=6, loc="lower right")

            for k, dpt in enumerate(d21):
                rows_csv.append({
                    "alpha_u": au_val, "scenario": sc_tag,
                    "source": "engine_Tref21", "depth_m": float(dpt),
                    "T_C_t0":   float(T_pl_C),
                    "T_C_t84":  float(e21_84[k]),
                    "T_C_t168": float(e21_168[k]),
                })

    fig.suptitle("Centerline T(depth) profiles — A and I scenarios at Hu=Hu_residual "
                 "(13,641.5 J/kg, k_uc=0.96).\n"
                 "Engine T_ref=23°C: solid.  Engine T_ref=21°C: dotted.  CW: dashed.  "
                 "Full slab depth: surface (top) to bottom (24.4 m = 80 ft).",
                 fontsize=10)
    fig.tight_layout()
    out_png = FIGS / "centerline_T_vs_depth_A_I_Tref21.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"\nWrote {out_png}")

    out_csv = DATA / "centerline_T_vs_depth_A_I_Tref21.csv"
    with open(out_csv, "w", newline="") as f:
        w = _csv.DictWriter(f, fieldnames=["alpha_u","scenario","source",
                                            "depth_m","T_C_t0","T_C_t84","T_C_t168"])
        w.writeheader()
        w.writerows(rows_csv)
    print(f"Wrote {out_csv}")


def main():
    rerun_save_Tref21()
    build_figure()


if __name__ == "__main__":
    main()
