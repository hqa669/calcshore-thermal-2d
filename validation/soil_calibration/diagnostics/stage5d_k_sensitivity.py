#!/usr/bin/env python3
"""Stage 5d-prep: k_concrete Bidirectional Sensitivity Test (±2%).

Probes whether multiplying k_uc by 0.98 or 1.02 reduces R2 (bottom profile
at centerline) residuals for failing runs F and I, and checks symmetry across
cooling/warming pairs.

No engine source changes — wrapper multiplies mix.thermal_conductivity_BTU_hr_ft_F
before passing to solve_hydration_2d.  No commits.  t=168 hr only.
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

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import load_engine_csv, resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

CW_RUNS   = os.path.join(SC, "cw_runs")
ENG_DIR   = os.path.join(SC, "stage5c_runs")
OUT_DIR   = os.path.join(SC, "diagnostics")
os.makedirs(OUT_DIR, exist_ok=True)

BLANKET_M    = 0.0
MODEL_SOIL   = False
IS_SUBMERGED = True
COMPARE_HR   = 168
GATE_TOL_F   = 0.35

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

K_FACTORS = [1.00, 0.98, 1.02]

R1_DI = 24        # side profile row
R2_WI = 0         # bottom profile column
R2_DI = slice(24, 49)   # di=24..48


# ---------------------------------------------------------------------------
def _load_cw_t168(folder, cw_depths_m, cw_widths_m):
    path = os.path.join(CW_RUNS, folder, "output.txt")
    v = parse_cw_temp_output(path)
    ti = int(np.abs(v.time_hrs - COMPARE_HR).argmin())
    return v.T_field_F[ti]


def _run_engine(folder, pl_F, so_F, k_factor, cw_depths_m, cw_widths_m):
    dat_path = os.path.join(CW_RUNS, folder, "input.dat")
    mix, geom, constr, _ = parse_cw_dat(dat_path)
    factor, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = factor
    mix.Hu_J_kg_effective = mix.Hu_J_kg * factor

    # ── k override (wrapper-only, no engine source change) ──
    mix.thermal_conductivity_BTU_hr_ft_F *= k_factor

    constr.model_soil = MODEL_SOIL
    constr.is_submerged = IS_SUBMERGED

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=IS_SUBMERGED,
        model_soil=MODEL_SOIL,
        blanket_thickness_m=BLANKET_M,
    )
    T0_C     = (pl_F - 32.0) * 5.0 / 9.0
    T_soil_C = (so_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    T_initial[grid.is_soil] = T_soil_C

    env = make_neutral_env(pl_F)
    result = solve_hydration_2d(
        grid, mix, T_initial,
        duration_s=COMPARE_HR * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=env,
        construction=constr,
        T_ground_deep_C=T_soil_C,
        diagnostic_outputs=False,
    )
    jslice, islice = grid.concrete_slice()
    ti = nearest_time_idx(result.t_s, float(COMPARE_HR))
    T_conc_F = result.T_field_C[ti, jslice, islice] * 9.0 / 5.0 + 32.0
    eng_y = grid.y[jslice]
    eng_x = grid.x[islice]
    eng_on_cw = resample_engine_to_cw(eng_y, eng_x, T_conc_F,
                                       cw_depths_m, cw_widths_m)
    return eng_on_cw


def _region_maxR(R, r1_di, r2_wi, r2_di):
    r1 = R[r1_di, :]
    r2 = R[r2_di, r2_wi]
    return float(np.max(np.abs(r1))), float(np.max(np.abs(r2)))


# ---------------------------------------------------------------------------
def main():
    # CW reference grid
    v0 = parse_cw_temp_output(os.path.join(CW_RUNS, "runA_baseline", "output.txt"))
    cw_depths_m, cw_widths_m = v0.depths_m, v0.widths_m
    di_r2 = np.arange(24, 49)

    # -------------------------------------------------------------------------
    # Pre-load CW fields and cached baseline
    # -------------------------------------------------------------------------
    cw_fields = {}
    baseline_cached_r2 = {}
    baseline_cached_r1 = {}
    for label, folder, pl, so in RUNS:
        cw_fields[label] = _load_cw_t168(folder, cw_depths_m, cw_widths_m)
        csv = os.path.join(ENG_DIR, f"run{label}_t168.csv")
        ey, ex, eT = load_engine_csv(csv)
        eng_cache = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)
        R_cache = eng_cache - cw_fields[label]
        r1, r2 = _region_maxR(R_cache, R1_DI, R2_WI,
                               slice(24, 48))   # excl di=48 for R2
        baseline_cached_r1[label] = r1
        baseline_cached_r2[label] = r2

    # -------------------------------------------------------------------------
    # Run all 9 × 3 configurations
    # -------------------------------------------------------------------------
    results = {}   # (label, k_factor) -> (R field, eng_on_cw)
    for k_factor in K_FACTORS:
        label_str = f"k×{k_factor:.2f}"
        print(f"\n--- {label_str} ---")
        for label, folder, pl, so in RUNS:
            print(f"  Run {label} ...", end=" ", flush=True)
            eng_on_cw = _run_engine(folder, pl, so, k_factor,
                                     cw_depths_m, cw_widths_m)
            R = eng_on_cw - cw_fields[label]
            results[(label, k_factor)] = (R, eng_on_cw)
            r1_max, r2_max = _region_maxR(R, R1_DI, R2_WI, slice(24, 48))
            print(f"R1={r1_max:.4f}  R2={r2_max:.4f}")

    # -------------------------------------------------------------------------
    # §4.1 Sanity check: factor=1.00 reproduces cached baseline
    # -------------------------------------------------------------------------
    print("\n=== §4.1 Sanity Check: k×1.00 vs Cached Baseline ===")
    hdr = f"{'Run':>4}  {'Cached R2':>10}  {'k×1.00 R2':>10}  {'Diff':>8}"
    print(hdr)
    print("-" * len(hdr))
    sanity_ok = True
    for label, *_ in RUNS:
        R_10, _ = results[(label, 1.00)]
        _, r2_10 = _region_maxR(R_10, R1_DI, R2_WI, slice(24, 48))
        r2_base  = baseline_cached_r2[label]
        diff = abs(r2_10 - r2_base)
        if diff > 0.001:
            sanity_ok = False
        print(f"  {label:>2}  {r2_base:>10.4f}  {r2_10:>10.4f}  "
              f"{diff:>8.5f}  {'OK' if diff <= 0.001 else 'FAIL'}")
    print(f"\n  Overall sanity: {'PASS' if sanity_ok else '*** FAIL — check wrapper ***'}")
    if not sanity_ok:
        print("  Stopping — sanity check failed.")
        return

    # -------------------------------------------------------------------------
    # §4.2 Bidirectional sensitivity tables
    # -------------------------------------------------------------------------
    for region, r_idx, cached_map, excl_di48 in [
        ("R2 (bottom profile, wi=0)", 1, baseline_cached_r2, True),
        ("R1 (side profile, di=24)",  0, baseline_cached_r1, False),
    ]:
        print(f"\n=== §4.2 Bidirectional Sensitivity — {region} ===")
        r2_slice = slice(24, 48) if excl_di48 else None
        hdr2 = (f"{'Run':>4}  {'|ΔT|':>5}  {'dir':>5}  "
                f"{'base':>8}  {'k×0.98':>8}  {'Δ0.98':>7}  "
                f"{'k×1.02':>8}  {'Δ1.02':>7}")
        print(hdr2)
        print("-" * len(hdr2))
        for label, folder, pl, so in RUNS:
            dT = abs(pl - so)
            direction = "cool" if pl > so else "warm"
            base_val = cached_map[label]
            R_098, _ = results[(label, 0.98)]
            R_102, _ = results[(label, 1.02)]
            if excl_di48:
                _, v098 = _region_maxR(R_098, R1_DI, R2_WI, r2_slice)
                _, v102 = _region_maxR(R_102, R1_DI, R2_WI, r2_slice)
            else:
                v098, _ = _region_maxR(R_098, R1_DI, R2_WI, r2_slice or slice(24,48))
                v102, _ = _region_maxR(R_102, R1_DI, R2_WI, r2_slice or slice(24,48))
                v098 = float(np.max(np.abs(R_098[R1_DI, :])))
                v102 = float(np.max(np.abs(R_102[R1_DI, :])))
            d098 = v098 - base_val
            d102 = v102 - base_val
            print(f"  {label:>2}  {dT:>5}  {direction:>5}  "
                  f"{base_val:>8.4f}  {v098:>8.4f}  {d098:>+7.4f}  "
                  f"{v102:>8.4f}  {d102:>+7.4f}")

    # -------------------------------------------------------------------------
    # §4.3 Detailed R2 profiles for Runs F and I
    # -------------------------------------------------------------------------
    for label in ("F", "I"):
        folder = next(f for lbl, f, *_ in RUNS if lbl == label)
        pl, so = next((pl, so) for lbl, _, pl, so in RUNS if lbl == label)
        csv = os.path.join(ENG_DIR, f"run{label}_t168.csv")
        ey, ex, eT = load_engine_csv(csv)
        eng_base = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)

        print(f"\n=== §4.3 Run {label} R2 Detailed Profile ===")
        hdr3 = (f"{'di':>4}  {'z(ft)':>7}  "
                f"{'T_eng_base':>12}  {'T_CW':>8}  "
                f"{'R_base':>8}  {'R_0.98':>8}  {'R_1.02':>8}")
        print(hdr3)
        print("-" * len(hdr3))
        for di in di_r2:
            z_ft = -float(cw_depths_m[di]) / 0.3048
            T_eng_b = float(eng_base[di, R2_WI])
            T_cw    = float(cw_fields[label][di, R2_WI])
            R_b     = T_eng_b - T_cw
            R_098   = float(results[(label, 0.98)][0][di, R2_WI])
            R_102   = float(results[(label, 1.02)][0][di, R2_WI])
            print(f"  {di:>4}  {z_ft:>7.1f}  "
                  f"{T_eng_b:>12.4f}  {T_cw:>8.4f}  "
                  f"{R_b:>+8.4f}  {R_098:>+8.4f}  {R_102:>+8.4f}")

    # -------------------------------------------------------------------------
    # §4.4 Visual figure
    # -------------------------------------------------------------------------
    _plot_fi(results, cw_fields, cw_depths_m, cw_widths_m, di_r2)

    # -------------------------------------------------------------------------
    # §4.5 Symmetry check (cooling vs warming pairs)
    # -------------------------------------------------------------------------
    print("\n=== §4.5 Cooling/Warming Symmetry Check (R2) ===")
    PAIRS = [("F", "H", 28), ("I", "G", 27), ("E", "C", 17), ("B", "D", 13)]
    hdr5 = (f"{'Pair':>8}  {'|ΔT|':>5}  "
            f"{'cool Δ(0.98)':>13}  {'warm Δ(1.02)':>13}  {'asymm':>8}")
    print(hdr5)
    print("-" * len(hdr5))
    for lbl_c, lbl_w, dT in PAIRS:
        R_c_098, _ = results[(lbl_c, 0.98)]
        R_c_base   = baseline_cached_r2[lbl_c]
        _, v_c_098 = _region_maxR(R_c_098, R1_DI, R2_WI, slice(24, 48))
        d_cool = v_c_098 - R_c_base

        R_w_102, _ = results[(lbl_w, 1.02)]
        R_w_base   = baseline_cached_r2[lbl_w]
        _, v_w_102 = _region_maxR(R_w_102, R1_DI, R2_WI, slice(24, 48))
        d_warm = v_w_102 - R_w_base

        asym = abs(d_cool - d_warm)
        print(f"  {lbl_c+'/'+lbl_w:>8}  {dT:>5}  "
              f"{d_cool:>+13.4f}  {d_warm:>+13.4f}  {asym:>8.4f}  "
              f"{'clean (<0.05)' if asym < 0.05 else ('borderline' if asym < 0.10 else 'ASYMMETRIC (>0.10)')}")

    # -------------------------------------------------------------------------
    # §7 Synthesis
    # -------------------------------------------------------------------------
    print("\n=== §7 Synthesis ===\n")
    _synthesis(results, baseline_cached_r2, baseline_cached_r1, PAIRS)


# ---------------------------------------------------------------------------
def _synthesis(results, cached_r2, cached_r1, pairs):
    def r2_max(label, kf):
        R, _ = results[(label, kf)]
        _, v = _region_maxR(R, R1_DI, R2_WI, slice(24, 48))
        return v

    f_base  = cached_r2["F"];  f_098 = r2_max("F", 0.98); f_102 = r2_max("F", 1.02)
    i_base  = cached_r2["I"];  i_098 = r2_max("I", 0.98); i_102 = r2_max("I", 1.02)

    f_improves = f_098 < f_base
    i_improves = i_098 < i_base
    f_closes   = f_098 <= GATE_TOL_F
    i_closes   = i_098 <= GATE_TOL_F

    # Check heating runs improve under k×1.02
    heating_labels = ["C", "D", "G", "H"]
    heat_improve_102 = all(r2_max(lbl, 1.02) < cached_r2[lbl] for lbl in heating_labels)
    heat_worsen_098  = all(r2_max(lbl, 0.98) > cached_r2[lbl] for lbl in heating_labels)

    # Symmetry asym for F/H and I/G
    d_f_098 = f_098 - f_base
    d_h_102 = r2_max("H", 1.02) - cached_r2["H"]
    asym_fh = abs(d_f_098 - d_h_102)

    d_i_098 = i_098 - i_base
    d_g_102 = r2_max("G", 1.02) - cached_r2["G"]
    asym_ig = abs(d_i_098 - d_g_102)

    print(
        f"(a) k×0.98 effect on F and I (R2 max|R|):\n"
        f"    Run F: {f_base:.4f} → {f_098:.4f} (Δ={f_098-f_base:+.4f}°F, "
        f"{'IMPROVES' if f_improves else 'WORSENS'}, "
        f"{'closes gate' if f_closes else 'does not close gate'}).\n"
        f"    Run I: {i_base:.4f} → {i_098:.4f} (Δ={i_098-i_base:+.4f}°F, "
        f"{'IMPROVES' if i_improves else 'WORSENS'}, "
        f"{'closes gate' if i_closes else 'does not close gate'}).\n\n"
        f"(b) k×1.02 worsens F and I as predicted: "
        f"F {f_base:.4f}→{f_102:.4f}, I {i_base:.4f}→{i_102:.4f}. "
        f"{'Confirmed — both worsen under k×1.02.' if (f_102 > f_base and i_102 > i_base) else 'UNEXPECTED — one or both improve under k×1.02.'}\n\n"
        f"(c) Heating run response to k perturbation:\n"
        f"    C/D/G/H improve under k×1.02: {'YES' if heat_improve_102 else 'NO'}.\n"
        f"    C/D/G/H worsen under k×0.98:  {'YES' if heat_worsen_098 else 'NO'}.\n"
        f"    {'Heating runs respond in the OPPOSITE direction to cooling runs — bulk α hypothesis supported.' if (heat_improve_102 and heat_worsen_098) else 'Heating runs do NOT consistently reverse — hypothesis weakened.'}\n\n"
        f"(d) Cooling/warming symmetry:\n"
        f"    F/H (|ΔT|=28°F): cooling Δ under k×0.98 = {d_f_098:+.4f}, "
        f"warming Δ under k×1.02 = {d_h_102:+.4f}, asymmetry = {asym_fh:.4f}°F "
        f"({'clean' if asym_fh < 0.05 else 'borderline' if asym_fh < 0.10 else 'ASYMMETRIC'}).\n"
        f"    I/G (|ΔT|=27°F): cooling Δ under k×0.98 = {d_i_098:+.4f}, "
        f"warming Δ under k×1.02 = {d_g_102:+.4f}, asymmetry = {asym_ig:.4f}°F "
        f"({'clean' if asym_ig < 0.05 else 'borderline' if asym_ig < 0.10 else 'ASYMMETRIC'}).\n\n"
        f"(e) The Run F and I R2 detailed profiles (§4.3) show whether the k perturbation\n"
        f"    shifts the residual curve uniformly (bulk α effect) or changes shape (localized BC).\n"
        f"    Inspect the di=24..36 (flat bulk) vs di=40..46 (near-BC ramp) segments separately."
    )


# ---------------------------------------------------------------------------
def _plot_fi(results, cw_fields, cw_depths_m, cw_widths_m, di_r2):
    out_path = os.path.join(OUT_DIR, "stage5d_k_sensitivity_F_I.png")
    print(f"\n--- §4.4 Generating figure ---")

    z_ft = -cw_depths_m[di_r2] / 0.3048   # −40 to −80 ft

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
    fig.suptitle(
        "Stage 5d-prep — k_concrete ±2% Sensitivity: R2 Bottom Profile (wi=0, t=168 hr)\n"
        "Residual = Engine − CW (°F)",
        fontsize=10, fontweight="bold",
    )

    for ax, label in zip(axes, ("F", "I")):
        pl, so = next((pl, so) for lbl, _, pl, so in RUNS if lbl == label)
        # Baseline from cached CSV
        from stage3_compare import load_engine_csv
        csv = os.path.join(ENG_DIR, f"run{label}_t168.csv")
        ey, ex, eT = load_engine_csv(csv)
        eng_base = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)
        R_base = (eng_base - cw_fields[label])[di_r2, R2_WI]

        R_098 = results[(label, 0.98)][0][di_r2, R2_WI]
        R_102 = results[(label, 1.02)][0][di_r2, R2_WI]

        ax.plot(z_ft, R_base, "r-o",  ms=5, lw=2.0, label="Baseline (k×1.00)")
        ax.plot(z_ft, R_098,  "r--s", ms=5, lw=1.8, label="k×0.98 (−2%)")
        ax.plot(z_ft, R_102,  "r:^",  ms=5, lw=1.8, label="k×1.02 (+2%)")

        ax.axhline(0,           color="black", lw=0.8)
        ax.axhline( GATE_TOL_F, color="green", ls="--", lw=1.2)
        ax.axhline(-GATE_TOL_F, color="green", ls="--", lw=1.2)
        ax.axhspan(-GATE_TOL_F, GATE_TOL_F, color="green", alpha=0.06)

        # di=45 marker
        di45_z = -float(cw_depths_m[45]) / 0.3048
        ax.axvline(di45_z, color="gray", ls=":", lw=1.0)
        ax.text(di45_z + 0.3, -0.48, "di=45", fontsize=7, color="gray")

        ax.set_xlim(-82, -38)
        ax.set_ylim(-0.55, 0.55)
        ax.set_xlabel("Depth z (ft)", fontsize=9)
        ax.set_ylabel("Engine − CW (°F)", fontsize=9)
        ax.set_title(
            f"Run {label}  (placement={pl}°F, soil={so}°F, |ΔT|={abs(pl-so)}°F)",
            fontsize=9,
        )
        ax.legend(fontsize=8, loc="lower left")
        ax.grid(True, alpha=0.25)
        ax.tick_params(labelsize=8)
        ax.text(-81, 0.38, "±0.35°F gate", fontsize=7, color="green")

    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


if __name__ == "__main__":
    main()
