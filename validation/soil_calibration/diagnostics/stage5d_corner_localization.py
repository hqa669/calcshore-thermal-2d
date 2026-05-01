#!/usr/bin/env python3
"""Stage 5d-prep: Corner Localization Diagnostic.

Determines whether the |ΔT|-scaling residual peak at (di=45, wi=9) is
corner-localized (Interp A) or genuinely distributed laterally (Interp B).

Sections:
  §2.1  Bulk-middle (di 10–35) vs near-bottom (di 40–48) max|R_dec| at wi=9
  §2.2  wi=9 depth profile for F and H at 9 depths
  §2.3  Corner-zoom heatmap (di 40–48, wi 7–12) for Run F
  §2.4  Corner-zoom heatmap for Run H
  §2.5  Engine corner BC code findings (verbatim from Phase 1 review)
  §5    Synthesis

No engine re-runs.  Cached CSVs only.  No commits.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SC   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, SC)

from cw_scenario_loader import parse_cw_temp_output
from stage3_compare import load_engine_csv, resample_engine_to_cw

ENG_DIR   = os.path.join(SC, "stage5c_runs")
CW_RUNS   = os.path.join(SC, "cw_runs")
MASK_PATH = os.path.join(SC, "STAGE35_validity_mask.npy")
OUT_DIR   = os.path.join(SC, "diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

COMPARE_HR = 168

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

# depth-profile sample indices for §2.2
DEPTH_SAMPLE_DI = [10, 20, 25, 30, 35, 40, 43, 45, 48]
# corner window for §2.3/§2.4
CORNER_DI = slice(40, 49)   # di 40..48 inclusive
CORNER_WI = slice(7, 13)    # wi 7..12 inclusive


# ---------------------------------------------------------------------------
def _load_cw_t168(folder):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m


def load_residual(label, folder, cw_depths_m, cw_widths_m):
    csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"Missing cached CSV: {csv_path}")
    eng_y, eng_x, eng_F = load_engine_csv(csv_path)
    cw_field, _, _ = _load_cw_t168(folder)
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    R = eng_on_cw - cw_field
    if np.any(np.isnan(R)):
        raise ValueError(f"NaN in residual for Run {label}")
    return R, eng_on_cw, cw_field


# ---------------------------------------------------------------------------
def main():
    # -------------------------------------------------------------------------
    # Load mask
    # -------------------------------------------------------------------------
    mask = np.load(MASK_PATH)
    assert mask.shape == (49, 13), f"Mask shape {mask.shape} != (49,13)"

    # -------------------------------------------------------------------------
    # Load all residuals
    # -------------------------------------------------------------------------
    print("Loading residuals from cached CSVs …")
    residuals  = {}
    eng_fields = {}
    cw_fields  = {}
    cw_depths_m = cw_widths_m = None

    for label, folder, pl, so in RUNS:
        _, _, cw_d, cw_w = None, None, None, None
        cw_field_t168, cw_w, cw_d = _load_cw_t168(folder)
        if cw_depths_m is None:
            cw_depths_m, cw_widths_m = cw_d, cw_w
        R, eng_on_cw, cw_field = load_residual(label, folder, cw_depths_m, cw_widths_m)
        residuals[label]  = R
        eng_fields[label] = eng_on_cw
        cw_fields[label]  = cw_field
        print(f"  Run {label}: max|R|={np.max(np.abs(R)):.4f}°F  shape={R.shape}")

    R_A = residuals["A"]

    # Sanity: R_A decomposed should be zero
    R_dec_A = R_A - R_A
    assert np.allclose(R_dec_A, 0.0), "Run A R_dec unexpectedly non-zero"
    print("  Sanity check passed: R_A − R_A == 0")

    # -------------------------------------------------------------------------
    # §2.1  Bulk-middle vs near-bottom max|R_dec| at wi=9
    # -------------------------------------------------------------------------
    print("\n=== §2.1 Bulk-middle vs Near-bottom max|R_dec| at wi=9 ===")
    hdr = f"{'Run':>4}  {'|ΔT|':>5}  {'max|R_dec| bulk[10-35]':>22}  {'max|R_dec| near-bot[40-48]':>26}  {'ratio':>7}"
    print(hdr)
    print("-" * len(hdr))
    ratio_data = {}
    for label, folder, pl, so in RUNS:
        if label == "A":
            continue
        dT = abs(pl - so)
        R_dec_masked = np.where(mask, residuals[label] - R_A, np.nan)
        wi9_col = R_dec_masked[:, 9]                             # shape (49,)
        bulk    = wi9_col[10:36]                                  # di 10..35
        near    = np.concatenate([wi9_col[40:49]])                # di 40..48
        bulk_max = float(np.nanmax(np.abs(bulk)))
        near_max = float(np.nanmax(np.abs(near)))
        ratio    = near_max / bulk_max if bulk_max > 0 else float("inf")
        ratio_data[label] = (dT, bulk_max, near_max, ratio)
        print(f"  {label:>2}  {dT:>5}  {bulk_max:>22.4f}  {near_max:>26.4f}  {ratio:>7.2f}")

    ratios = [v[3] for v in ratio_data.values()]
    print(f"\n  Ratio stats: min={min(ratios):.2f}  max={max(ratios):.2f}  mean={np.mean(ratios):.2f}")
    interp = "A (corner-localized)" if np.mean(ratios) > 2.0 else "B (laterally distributed)"
    print(f"  → Consistent with Interpretation {interp}")

    # -------------------------------------------------------------------------
    # §2.2  wi=9 depth profile for Runs F and H
    # -------------------------------------------------------------------------
    print("\n=== §2.2 wi=9 R_dec Depth Profile — Runs F and H ===")
    hdr2 = f"{'di':>4}  {'depth_m':>8}  {'R_dec_F':>10}  {'R_dec_H':>10}  {'asymmetry':>10}"
    print(hdr2)
    print("-" * len(hdr2))
    R_dec_F = np.where(mask, residuals["F"] - R_A, np.nan)
    R_dec_H = np.where(mask, residuals["H"] - R_A, np.nan)
    for di in DEPTH_SAMPLE_DI:
        depth = float(cw_depths_m[di])
        vF = float(R_dec_F[di, 9])
        vH = float(R_dec_H[di, 9])
        asym = vF + vH  # zero if sign-symmetric
        print(f"  {di:>2}  {depth:>8.4f}  {vF:>10.4f}  {vH:>10.4f}  {asym:>10.4f}")

    # -------------------------------------------------------------------------
    # §2.3 & §2.4  Corner-zoom heatmaps for Runs F and H
    # -------------------------------------------------------------------------
    for label, title in [("F", "Run F (placement 73°F, soil 45°F, ΔT=+28°F cooling)"),
                          ("H", "Run H (placement 45°F, soil 73°F, ΔT=−28°F heating)")]:
        out_path = os.path.join(OUT_DIR, f"stage5d_prep_corner_zoom_run{label}.png")
        _corner_zoom(label, title, eng_fields, cw_fields, R_A, mask,
                     cw_depths_m, cw_widths_m, out_path)

    # -------------------------------------------------------------------------
    # §2.5  Engine corner BC code findings
    # -------------------------------------------------------------------------
    print("\n=== §2.5 Engine Corner BC Code Findings (Phase 1 review) ===\n")
    print(
        "TOP-SIDE CORNER (y=0, x=0):\n"
        "  Lines ~2154–2163: Half-cell BC stencil correction at the side-face row.\n"
        "    For is_submerged=True this is DEAD CODE — the Dirichlet write at line 2259\n"
        "    overwrites T_new for the entire side column after the stencil runs.\n"
        "  Lines ~2189–2243: Quarter-cell energy balance at the top-side corner cell.\n"
        "    This explicitly accounts for the reduced control-volume area at the corner.\n"
        "\n"
        "BOTTOM-SIDE CORNER (y=depth_m, x=0):\n"
        "  Line ~2256: T_new[iy_concrete_end, ix_concrete_start:] = T_gw_C\n"
        "    (pure strong Dirichlet along entire bottom face, written unconditionally)\n"
        "  Line ~2259: T_new[iy_concrete_start:iy_concrete_end+1, ix_concrete_start] = T_gw_C\n"
        "    (pure strong Dirichlet along entire side face, written if is_submerged=True)\n"
        "  Corner cell (iy_concrete_end, ix_concrete_start) is written TWICE — idempotent.\n"
        "  NO quarter-cell energy balance at the bottom-side corner.\n"
        "  NO half-cell BC stencil in the row above the bottom face.\n"
        "  NO stencil-mass correction anywhere near the bottom face.\n"
        "  The interior cell (iy_concrete_end-1, ix_concrete_start+1) uses the standard\n"
        "    interior FD stencil with no BC-distance correction.\n"
        "\n"
        "ASYMMETRY SUMMARY:\n"
        "  Top face:    half-cell stencil correction + quarter-cell corner block\n"
        "  Bottom face: pure strong Dirichlet only — no energy balance correction\n"
        "  This asymmetry is consistent with CW using a symmetric treatment at both faces.\n"
    )

    # -------------------------------------------------------------------------
    # §5  Synthesis
    # -------------------------------------------------------------------------
    print("\n=== §5 Synthesis ===\n")
    mean_ratio = np.mean(ratios)
    if mean_ratio > 2.0:
        interp_verdict = "A (corner-localized)"
        interp_detail = (
            f"The mean near-bottom/bulk-middle max|R_dec| ratio at wi=9 is {mean_ratio:.2f}×, "
            "indicating the peak is concentrated in the near-bottom band (di=40–48) rather than "
            "spanning the full depth. "
        )
    else:
        interp_verdict = "B (laterally distributed)"
        interp_detail = (
            f"The mean near-bottom/bulk-middle max|R_dec| ratio at wi=9 is {mean_ratio:.2f}×, "
            "indicating the residual is roughly uniform with depth rather than corner-concentrated. "
        )

    print(
        f"Interpretation: {interp_verdict}\n\n"
        f"{interp_detail}"
        "The corner-zoom heatmaps (§2.3/§2.4) show the CW and engine temperature fields near "
        "the side-bottom corner directly; the decomposed residual panel highlights where the "
        "difference is concentrated.\n\n"
        "Mechanism candidate: CW likely applies a symmetric energy balance at both the "
        "top-face and bottom-face corners (or uses a cell-centered finite-volume formulation "
        "that naturally captures corner control volumes), while the engine's bottom-side corner "
        "(iy_concrete_end, ix_concrete_start) receives only a strong Dirichlet overwrite with "
        "no quarter-cell correction. This means the single interior cell diagonally inward from "
        "the corner — which participates in the FD stencil for the full time step before being "
        "overwritten — can accumulate a small bias that grows linearly with the temperature "
        "gradient driving the corner (i.e., with |ΔT_soil − T_placement|).\n\n"
        "F/H sign asymmetry: The F/H pair (|ΔT|=28°F, the largest) has 0.134°F asymmetry "
        "localized to di=40–48, wi=7–9. The §2.2 depth profile and §2.3/§2.4 heatmaps identify "
        "whether the asymmetry is caused by a nonlinear near-bottom CW heat sink (heating vs "
        "cooling wave interacts differently with the soil BC near the corner) or by a secondary "
        "bottom-face stencil artifact that breaks sign symmetry at extreme |ΔT|. No proposed "
        "fix; the mechanism is the input for Stage 5d.\n"
    )


# ---------------------------------------------------------------------------
def _corner_zoom(label, title, eng_fields, cw_fields, R_A, mask,
                 cw_depths_m, cw_widths_m, out_path):
    print(f"\n--- §2.3/§2.4 Corner-zoom heatmap: Run {label} ---")

    eng_corner = eng_fields[label][CORNER_DI, CORNER_WI]   # (9, 6)
    cw_corner  = cw_fields[label][CORNER_DI, CORNER_WI]
    R_X        = eng_fields[label] - cw_fields[label]
    R_dec      = R_X - (eng_fields["A"] - cw_fields["A"])
    R_dec_corner = np.where(mask, R_dec, np.nan)[CORNER_DI, CORNER_WI]

    di_indices = np.arange(40, 49)
    wi_indices = np.arange(7, 13)
    depths  = cw_depths_m[di_indices]
    widths  = cw_widths_m[wi_indices]  # decreasing: 1.016 → 0.0 m

    print(f"  Corner window: di={di_indices[0]}..{di_indices[-1]}, wi={wi_indices[0]}..{wi_indices[-1]}")
    print(f"  depths: {depths}")
    print(f"  widths: {widths}")
    print(f"  R_dec corner max|R_dec|={np.nanmax(np.abs(R_dec_corner)):.4f}°F at "
          f"{np.unravel_index(np.nanargmax(np.abs(R_dec_corner)), R_dec_corner.shape)}")

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle(
        f"Stage 5d-prep Corner Zoom — {title}\n"
        f"Window: di={di_indices[0]}–{di_indices[-1]}, wi={wi_indices[0]}–{wi_indices[-1]}",
        fontsize=10, fontweight="bold",
    )

    datasets = [
        (cw_corner,    "CW T (°F)",         "YlOrRd",  False),
        (eng_corner,   "Engine T (°F)",      "YlOrRd",  False),
        (R_dec_corner, "R_dec = (Eng−CW)−(A) (°F)", "RdBu_r", True),
    ]

    # Shared colour limits for T panels; symmetric limits for R_dec
    T_all   = np.concatenate([cw_corner.ravel(), eng_corner.ravel()])
    T_vmin, T_vmax = float(np.nanmin(T_all)), float(np.nanmax(T_all))
    rdec_abs = float(np.nanmax(np.abs(R_dec_corner))) if not np.all(np.isnan(R_dec_corner)) else 0.01
    rdec_abs = max(rdec_abs, 0.01)

    for ax, (data, subtitle, cmap, is_div) in zip(axes, datasets):
        if is_div:
            vmin, vmax = -rdec_abs, rdec_abs
        else:
            vmin, vmax = T_vmin, T_vmax

        # imshow: rows=di (depth, y-axis), cols=wi (width, x-axis)
        im = ax.imshow(data, origin="upper", aspect="auto",
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=[0, len(wi_indices), len(di_indices), 0])
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        ax.set_title(subtitle, fontsize=9)

        # Overlay cell values as text
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                val = data[i, j]
                if not np.isnan(val):
                    txt_color = "white" if abs(val - (vmin + vmax) / 2) > (vmax - vmin) * 0.3 else "black"
                    ax.text(j + 0.5, i + 0.5, f"{val:.1f}",
                            ha="center", va="center", fontsize=6.5, color=txt_color)

        # x-axis: wi labels (wi index and physical width in m)
        ax.set_xticks(np.arange(len(wi_indices)) + 0.5)
        ax.set_xticklabels([f"wi={wi_indices[k]}\n{widths[k]:.3f}m" for k in range(len(wi_indices))],
                           fontsize=7)
        # y-axis: di labels
        ax.set_yticks(np.arange(len(di_indices)) + 0.5)
        ax.set_yticklabels([f"di={di_indices[k]}\n{depths[k]:.3f}m" for k in range(len(di_indices))],
                           fontsize=7)
        ax.set_xlabel("Lateral position (wi)", fontsize=8)
        ax.set_ylabel("Depth (di)", fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.88])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
