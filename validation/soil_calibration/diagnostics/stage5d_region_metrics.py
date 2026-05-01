#!/usr/bin/env python3
"""Stage 5d-prep: Region-Based Metric Diagnostic.

Computes residuals over three physically-distinct regions and evaluates
three candidate threshold structures as replacements for the current
single-metric max|R| gate.

Regions:
  R1 — Side profile at mid-depth:   di=24, wi=0..12  (13 cells)
  R2 — Bottom profile at centerline: wi=0, di=24..48  (25 cells)
  R3 — Bottom-side corner:           di=43..48, wi=8..12 (30 cells)
  R4 — Bulk middle (reference):      everything else in mask

Threshold structures A/B/C evaluated for all 9 runs.
No engine re-runs — cached t=168 CSVs only.
No engine source changes, no commits.
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

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

# Region definitions (CW grid indices)
# R1: side profile at mid-depth
R1_DI = 24
R1_WI = slice(0, 13)       # wi=0..12 (all lateral columns)

# R2: bottom profile at centerline
R2_WI = 0
R2_DI = slice(24, 49)      # di=24..48 (25 cells)

# R3: bottom-side corner
R3_DI = slice(43, 49)      # di=43..48 (6 rows)
R3_WI = slice(8, 13)       # wi=8..12  (5 cols)

# Threshold structures
STRUCTURES = {
    "A": {"R1_max": 0.35, "R2_max": 0.35, "R3_rmse": 0.30, "label": "Tight physics, loose corner"},
    "B": {"R1_max": 0.40, "R2_max": 0.40, "R3_rmse": 0.35, "label": "Standard physics, looser corner"},
    "C": {"R1_max": 0.35, "R2_max": 0.35, "R3_rmse": None,  "label": "Tight physics, no corner gate"},
}


# ---------------------------------------------------------------------------
def _load_residual(label, folder, cw_depths_m, cw_widths_m):
    csv_path = os.path.join(ENG_DIR, f"run{label}_t168.csv")
    eng_y, eng_x, eng_F = load_engine_csv(csv_path)
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    cw_field = v.T_field_F[ti]
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, eng_F, cw_depths_m, cw_widths_m)
    return eng_on_cw - cw_field, eng_on_cw, cw_field


def _region_metrics(R, mask):
    """Compute per-region metrics from residual field (49,13)."""
    # R1: single row
    r1 = R[R1_DI, R1_WI]                             # shape (13,)
    r1_max = float(np.max(np.abs(r1)))

    # R2: single column
    r2_all = R[R2_DI, R2_WI]                         # shape (25,)
    r2_excl = R[R2_DI.start:R2_DI.stop - 1, R2_WI]  # shape (24,) excl di=48
    r2_max_incl = float(np.max(np.abs(r2_all)))
    r2_max_excl = float(np.max(np.abs(r2_excl)))
    r2_di48_val = float(R[48, R2_WI])

    # R3: corner block
    r3 = R[R3_DI, R3_WI]                             # shape (6,5)
    r3_rmse = float(np.sqrt(np.mean(r3**2)))
    r3_max  = float(np.max(np.abs(r3)))

    # R4: bulk middle = masked cells outside R1 row, R2 col, R3 block
    r4_mask = mask.copy()
    r4_mask[R1_DI, :] = False        # exclude R1 row
    r4_mask[:, R2_WI] = False        # exclude R2 column
    r4_mask[R3_DI, R3_WI] = False    # exclude R3 block
    r4_vals = R[r4_mask]
    r4_max  = float(np.max(np.abs(r4_vals))) if r4_vals.size else float("nan")
    r4_rmse = float(np.sqrt(np.mean(r4_vals**2))) if r4_vals.size else float("nan")

    # R1 wi=12 (form face — Dirichlet forced) sanity
    r1_wi12_val = float(R[R1_DI, 12])

    return {
        "r1_max": r1_max, "r1_wi12": r1_wi12_val,
        "r2_max_incl": r2_max_incl, "r2_max_excl": r2_max_excl, "r2_di48": r2_di48_val,
        "r3_rmse": r3_rmse, "r3_max": r3_max,
        "r4_max": r4_max, "r4_rmse": r4_rmse,
    }


def _full_field_max(R, mask):
    return float(np.max(np.abs(np.where(mask, R, 0.0))))


# ---------------------------------------------------------------------------
def main():
    mask = np.load(MASK_PATH)
    assert mask.shape == (49, 13)

    # Get CW grid from Run A
    from cw_scenario_loader import parse_cw_temp_output
    v0 = parse_cw_temp_output(os.path.join(CW_RUNS, "runA_baseline", "output.txt"))
    cw_depths_m, cw_widths_m = v0.depths_m, v0.widths_m

    print("Loading residuals …")
    metrics   = {}
    resid_all = {}
    eng_all   = {}
    cw_all    = {}
    for label, folder, pl, so in RUNS:
        R, eng_on_cw, cw_field = _load_residual(label, folder, cw_depths_m, cw_widths_m)
        metrics[label]   = _region_metrics(R, mask)
        metrics[label]["stage5c_max"] = _full_field_max(R, mask)
        resid_all[label] = R
        eng_all[label]   = eng_on_cw
        cw_all[label]    = cw_field
        print(f"  Run {label}: R1={metrics[label]['r1_max']:.4f}  "
              f"R2={metrics[label]['r2_max_excl']:.4f}  "
              f"R3_rmse={metrics[label]['r3_rmse']:.4f}  "
              f"full={metrics[label]['stage5c_max']:.4f}°F")

    # -------------------------------------------------------------------------
    # §3.1 Region-by-region table
    # -------------------------------------------------------------------------
    print("\n=== §3.1 Region-by-Region Residual Table (t=168 hr) ===")
    hdr = (f"{'Run':>4}  {'|ΔT|':>5}  "
           f"{'R1 max|R|':>10}  {'R2 max|R|':>10}  {'R2 excl48':>10}  "
           f"{'R3 RMSE':>8}  {'R3 max':>8}  "
           f"{'R4 max|R|':>10}  {'R4 RMSE':>8}  {'Full max|R|':>11}")
    print(hdr)
    print("-" * len(hdr))
    for label, folder, pl, so in RUNS:
        m = metrics[label]
        dT = abs(pl - so)
        print(f"  {label:>2}  {dT:>5}  "
              f"{m['r1_max']:>10.4f}  {m['r2_max_incl']:>10.4f}  {m['r2_max_excl']:>10.4f}  "
              f"{m['r3_rmse']:>8.4f}  {m['r3_max']:>8.4f}  "
              f"{m['r4_max']:>10.4f}  {m['r4_rmse']:>8.4f}  {m['stage5c_max']:>11.4f}")

    # -------------------------------------------------------------------------
    # §3.2 Sanity checks
    # -------------------------------------------------------------------------
    print("\n=== §3.2 Sanity Checks ===")
    print("\n(a) Region 2 di=48 |R| (Dirichlet-forced bottom face):")
    all_ok_a = True
    for label, *_ in RUNS:
        val = metrics[label]["r2_di48"]
        ok = abs(val) < 0.05
        if not ok: all_ok_a = False
        print(f"  Run {label}: R[di=48, wi=0] = {val:+.4f}°F  {'OK' if ok else 'WARN >0.05'}")
    print(f"  Result: {'ALL PASS' if all_ok_a else 'SOME FAIL'}")

    print("\n(b) Region 1 wi=12 |R| (Dirichlet-forced form face):")
    all_ok_b = True
    for label, *_ in RUNS:
        val = metrics[label]["r1_wi12"]
        ok = abs(val) < 0.05
        if not ok: all_ok_b = False
        print(f"  Run {label}: R[di=24, wi=12] = {val:+.4f}°F  {'OK' if ok else 'WARN >0.05'}")
    print(f"  Result: {'ALL PASS' if all_ok_b else 'SOME FAIL'}")

    print("\n(c) Run A (|ΔT|=0) residuals (hydration noise only):")
    m_A = metrics["A"]
    print(f"  Region 1 max|R| = {m_A['r1_max']:.4f}°F")
    print(f"  Region 2 max|R| = {m_A['r2_max_excl']:.4f}°F")
    print(f"  Region 3 RMSE   = {m_A['r3_rmse']:.4f}°F")
    print(f"  (all should be small — expected ~0.05–0.15°F from bulk hydration noise)")

    # -------------------------------------------------------------------------
    # §3.3 Candidate threshold evaluation
    # -------------------------------------------------------------------------
    for struct_name, struct in STRUCTURES.items():
        print(f"\n=== §3.3 Structure {struct_name}: {struct['label']} ===")
        print(f"  R1 max|R| ≤ {struct['R1_max']}°F  |  "
              f"R2 max|R| ≤ {struct['R2_max']}°F  |  "
              f"R3 RMSE ≤ {struct['R3_rmse']}°F"
              if struct['R3_rmse'] else
              f"  R1 max|R| ≤ {struct['R1_max']}°F  |  "
              f"R2 max|R| ≤ {struct['R2_max']}°F  |  "
              f"R3 RMSE: reported only")
        hdr_s = (f"{'Run':>4}  {'|ΔT|':>5}  {'R1':>5}  {'R2':>5}  "
                 f"{'R3':>5}  {'Overall':>8}  {'Fail region(s)':>20}")
        print(hdr_s)
        print("-" * len(hdr_s))
        pass_count = 0
        for label, folder, pl, so in RUNS:
            m = metrics[label]
            dT = abs(pl - so)
            p1 = m["r1_max"]    <= struct["R1_max"]
            p2 = m["r2_max_excl"] <= struct["R2_max"]
            p3 = (m["r3_rmse"]  <= struct["R3_rmse"]) if struct["R3_rmse"] else True
            overall = p1 and p2 and p3
            if overall: pass_count += 1
            fails = []
            if not p1: fails.append(f"R1={m['r1_max']:.3f}")
            if not p2: fails.append(f"R2={m['r2_max_excl']:.3f}")
            if not p3: fails.append(f"R3_rmse={m['r3_rmse']:.3f}")
            fail_str = ", ".join(fails) if fails else "—"
            print(f"  {label:>2}  {dT:>5}  "
                  f"{'✓' if p1 else '✗':>5}  {'✓' if p2 else '✗':>5}  "
                  f"{'✓' if p3 else '✗':>5}  "
                  f"{'PASS' if overall else 'FAIL':>8}  {fail_str:>20}")
        print(f"  → {pass_count}/9 runs pass under Structure {struct_name}")

    # -------------------------------------------------------------------------
    # §3.4 9×3 visual figure
    # -------------------------------------------------------------------------
    _plot_regions(resid_all, eng_all, cw_all, cw_depths_m, cw_widths_m, metrics, mask)

    # -------------------------------------------------------------------------
    # §3.5 Single-metric vs region-based comparison
    # -------------------------------------------------------------------------
    print("\n=== §3.5 Single-Metric vs Region-Based Gate Comparison ===")
    print(f"{'Run':>4}  {'|ΔT|':>5}  {'Stage5c max|R|':>15}  {'5c gate':>8}  "
          f"{'Struct A':>9}  {'Struct B':>9}  {'Struct C':>9}")
    print("-" * 75)
    for label, folder, pl, so in RUNS:
        m  = metrics[label]
        dT = abs(pl - so)
        # Stage 5c gate: max|R| ≤ 0.35
        g5c = "PASS" if m["stage5c_max"] <= 0.35 else "FAIL"
        results = {}
        for sn, struct in STRUCTURES.items():
            p1 = m["r1_max"]      <= struct["R1_max"]
            p2 = m["r2_max_excl"] <= struct["R2_max"]
            p3 = (m["r3_rmse"]    <= struct["R3_rmse"]) if struct["R3_rmse"] else True
            results[sn] = "PASS" if (p1 and p2 and p3) else "FAIL"
        print(f"  {label:>2}  {dT:>5}  {m['stage5c_max']:>15.4f}  {g5c:>8}  "
              f"{results['A']:>9}  {results['B']:>9}  {results['C']:>9}")

    # -------------------------------------------------------------------------
    # §5 Synthesis
    # -------------------------------------------------------------------------
    print("\n=== §5 Synthesis ===\n")
    _synthesis(metrics, cw_widths_m, cw_depths_m)


# ---------------------------------------------------------------------------
def _synthesis(metrics, cw_widths_m, cw_depths_m):
    r1_vals = {lbl: metrics[lbl]["r1_max"] for lbl, *_ in RUNS}
    r2_vals = {lbl: metrics[lbl]["r2_max_excl"] for lbl, *_ in RUNS}
    r3_vals = {lbl: metrics[lbl]["r3_rmse"] for lbl, *_ in RUNS}

    # Count passes per structure
    pass_counts = {}
    fail_info   = {}
    for sn, struct in STRUCTURES.items():
        passes, fails = 0, []
        for label, folder, pl, so in RUNS:
            m = metrics[label]
            p1 = m["r1_max"]      <= struct["R1_max"]
            p2 = m["r2_max_excl"] <= struct["R2_max"]
            p3 = (m["r3_rmse"]    <= struct["R3_rmse"]) if struct["R3_rmse"] else True
            if p1 and p2 and p3:
                passes += 1
            else:
                fail_regions = []
                if not p1: fail_regions.append("R1")
                if not p2: fail_regions.append("R2")
                if not p3: fail_regions.append("R3")
                fails.append((label, fail_regions))
        pass_counts[sn] = passes
        fail_info[sn]   = fails

    r1_list = list(r1_vals.values())
    r2_list = list(r2_vals.values())
    r3_list = list(r3_vals.values())

    r1_spread = max(r1_list) - min(r1_list)
    r2_spread = max(r2_list) - min(r2_list)

    r1_cluster = r1_spread < 0.15
    r2_cluster = r2_spread < 0.15

    print(
        f"(a) Pass counts: Structure A={pass_counts['A']}/9, B={pass_counts['B']}/9, C={pass_counts['C']}/9.\n"
        f"    Struct A failures: {[(lbl, regs) for lbl, regs in fail_info['A']] if fail_info['A'] else 'none'}.\n"
        f"    Struct B failures: {[(lbl, regs) for lbl, regs in fail_info['B']] if fail_info['B'] else 'none'}.\n"
        f"    Struct C failures: {[(lbl, regs) for lbl, regs in fail_info['C']] if fail_info['C'] else 'none'}.\n"
    )
    print(
        f"(b) Region 1 (side profile) max|R| range: [{min(r1_list):.3f}, {max(r1_list):.3f}]°F "
        f"(spread={r1_spread:.3f}°F). "
        f"{'CLUSTERED — lateral physics validates cleanly if threshold is above ~' + f'{max(r1_list):.2f}°F' if r1_cluster else 'SPREAD — |ΔT|-scaling residuals present in lateral physics'}.\n"
    )
    print(
        f"(c) Region 2 (bottom profile) max|R| range: [{min(r2_list):.3f}, {max(r2_list):.3f}]°F "
        f"(spread={r2_spread:.3f}°F). "
        f"{'CLUSTERED — vertical physics validates cleanly' if r2_cluster else 'SPREAD — |ΔT|-scaling residuals present in vertical physics'}.\n"
    )

    # R3 distribution check: compare max vs RMSE ratio
    r3_max_list  = [metrics[lbl]["r3_max"] for lbl, *_ in RUNS]
    ratio_max_rmse = [r3_max_list[i] / r3_list[i] for i in range(len(RUNS)) if r3_list[i] > 0]
    mean_ratio = np.mean(ratio_max_rmse)
    r3_outlier = mean_ratio > 2.0
    print(
        f"(d) Region 3 (corner) RMSE vs max|R| ratio: mean {mean_ratio:.2f}× across 9 runs. "
        f"{'A ratio > 2× suggests a single dominant outlier cell — RMSE is more robust than max|R| here. ' if r3_outlier else 'Ratio near 1× suggests roughly uniform distribution — RMSE and max|R| give similar information. '}"
        f"RMSE is the {'appropriate' if r3_outlier else 'conservative'} metric for Region 3.\n"
    )

    print(
        f"(e) Recommendation framing: "
        f"Region 1 captures lateral physics (side-BC diffusion) and Region 2 captures vertical "
        f"physics (bottom-BC diffusion) — both are far from the corner BC interaction and reflect "
        f"the bulk thermal material properties the sprint is calibrating. "
        f"Structure C (gate only R1 and R2, report R3 without threshold) most cleanly distinguishes "
        f"'engine physics validated' from 'engine has calibration gaps' because it excludes the "
        f"corner region whose residuals are amplified by a known BC-stencil asymmetry unrelated to "
        f"the bulk calibration objective. "
        f"Structure A is the correct choice if the corner stencil asymmetry is within scope to fix; "
        f"Structure C is correct if it is out of scope."
    )


# ---------------------------------------------------------------------------
def _plot_regions(resid_all, eng_all, cw_all, cw_depths_m, cw_widths_m, metrics, mask):
    out_path = os.path.join(OUT_DIR, "stage5d_prep_region_based_metrics.png")
    print(f"\n--- §3.4 Generating 9×3 figure ---")

    n_runs = len(RUNS)
    fig, axes = plt.subplots(n_runs, 3, figsize=(16, 3.2 * n_runs))
    fig.suptitle(
        "Stage 5d-prep Region-Based Metrics (t=168 hr)\n"
        "Col 1: Region 1 side profile (di=24)  |  "
        "Col 2: Region 2 bottom profile (wi=0)  |  "
        "Col 3: Region 3 corner heatmap (di=43–48, wi=8–12)",
        fontsize=10, fontweight="bold",
    )

    # Shared corner color range across all runs
    corner_vals = []
    for label, *_ in RUNS:
        corner_vals.extend(resid_all[label][R3_DI, R3_WI].ravel().tolist())
    corner_abs = max(float(np.max(np.abs(corner_vals))), 0.01)

    for row_idx, (label, folder, pl, so) in enumerate(RUNS):
        dT  = abs(pl - so)
        R   = resid_all[label]
        eng = eng_all[label]
        cw  = cw_all[label]
        m   = metrics[label]

        ax1, ax2, ax3 = axes[row_idx]

        # ------ Column 1: Region 1 side profile (di=24) ------
        # CW wi order: wi=0=CL (6.1m), wi=12=face (0.0m). Plot edge→CL
        x_plot = cw_widths_m[::-1]         # 0.0 → 6.1 m
        eng_r1 = eng[R1_DI, ::-1]
        cw_r1  = cw[R1_DI, ::-1]
        res_r1 = R[R1_DI, ::-1]

        ax1.plot(x_plot, eng_r1, "r-o", ms=4, lw=1.5, label="Engine")
        ax1.plot(x_plot, cw_r1,  "b-s", ms=4, lw=1.5, label="CW")
        ax1b = ax1.twinx()
        ax1b.plot(x_plot, res_r1, "k--", ms=3, lw=1.0, alpha=0.6)
        ax1b.axhline(0, color="gray", lw=0.5)
        ax1b.axhspan(-0.35, 0.35, color="green", alpha=0.05)
        ax1b.set_ylabel("Residual (°F)", fontsize=7, color="gray")
        ax1b.tick_params(labelsize=6, colors="gray")

        # Mark max|R| location
        max_wi_rev = int(np.argmax(np.abs(res_r1)))
        ax1b.axvline(x_plot[max_wi_rev], color="red", ls=":", lw=0.8, alpha=0.5)

        ax1.set_title(f"Run {label} (|ΔT|={dT}°F) — R1 side max|R|={m['r1_max']:.3f}°F",
                      fontsize=8)
        ax1.set_xlabel("Dist from face (m)", fontsize=7)
        ax1.set_ylabel("T (°F)", fontsize=7)
        ax1.tick_params(labelsize=6)
        if row_idx == 0:
            ax1.legend(fontsize=6, loc="lower right")

        # ------ Column 2: Region 2 bottom profile (wi=0) ------
        di_range = np.arange(24, 49)
        depths   = cw_depths_m[24:49]
        eng_r2   = eng[24:49, R2_WI]
        cw_r2    = cw[24:49, R2_WI]
        res_r2   = R[24:49, R2_WI]

        ax2.plot(depths, eng_r2, "r-o", ms=4, lw=1.5, label="Engine")
        ax2.plot(depths, cw_r2,  "b-s", ms=4, lw=1.5, label="CW")
        ax2b = ax2.twinx()
        ax2b.plot(depths, res_r2, "k--", ms=3, lw=1.0, alpha=0.6)
        ax2b.axhline(0, color="gray", lw=0.5)
        ax2b.axhspan(-0.35, 0.35, color="green", alpha=0.05)
        ax2b.set_ylabel("Residual (°F)", fontsize=7, color="gray")
        ax2b.tick_params(labelsize=6, colors="gray")

        max_di_idx = int(np.argmax(np.abs(res_r2)))
        ax2b.axvline(depths[max_di_idx], color="red", ls=":", lw=0.8, alpha=0.5)

        ax2.set_title(f"Run {label} — R2 bottom max|R|={m['r2_max_excl']:.3f}°F (excl di=48)",
                      fontsize=8)
        ax2.set_xlabel("Depth (m)", fontsize=7)
        ax2.set_ylabel("T (°F)", fontsize=7)
        ax2.tick_params(labelsize=6)

        # ------ Column 3: Region 3 corner heatmap ------
        corner = R[R3_DI, R3_WI]          # (6, 5)
        di_idx = np.arange(43, 49)
        wi_idx = np.arange(8, 13)

        im = ax3.imshow(corner, origin="upper", aspect="auto",
                        cmap="RdBu_r", vmin=-corner_abs, vmax=corner_abs,
                        extent=[0, len(wi_idx), len(di_idx), 0])
        plt.colorbar(im, ax=ax3, fraction=0.04, pad=0.02)

        for i in range(corner.shape[0]):
            for j in range(corner.shape[1]):
                val = corner[i, j]
                txt_c = "white" if abs(val) > corner_abs * 0.6 else "black"
                ax3.text(j + 0.5, i + 0.5, f"{val:.2f}",
                         ha="center", va="center", fontsize=6, color=txt_c)

        ax3.set_xticks(np.arange(len(wi_idx)) + 0.5)
        ax3.set_xticklabels([f"wi={wi_idx[k]}\n{cw_widths_m[wi_idx[k]]:.2f}m"
                              for k in range(len(wi_idx))], fontsize=6)
        ax3.set_yticks(np.arange(len(di_idx)) + 0.5)
        ax3.set_yticklabels([f"di={di_idx[k]}" for k in range(len(di_idx))], fontsize=6)
        ax3.set_title(f"Run {label} — R3 corner RMSE={m['r3_rmse']:.3f} max={m['r3_max']:.3f}°F",
                      fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
