#!/usr/bin/env python3
"""Sprint 9 Stage 1-pilot — engine runner with B1 wrapper.

Applies the Sprint 8/9 corrections to both pilot mixes:
  - compute_hu_factor(mix): composition-based Hu scaling (Sprint 7/8 standard)
  - Hu_residual conditional: passthrough for realistic Hu (>> 10,000 J/kg)
  - c(T_pl=73) shape correction to alpha_u; Hu inverse-compensated (B1)
  - k_uc × 0.96 (Sprint 7, already in engine source — no override)
  - model_soil=False, is_submerged=True, blanket=0.0

Pipeline order: Hu_raw → Hu_residual conditional → Hu_factored (×fac) → Hu_eff (÷c)

Saves per-mix hourly T(di, wi) and alpha(di, wi) on the CW grid as .npz files
for downstream analysis by analyze_pilot.py.

Also writes pilot_dataset_inventory.csv and pilot_wrapper_log.csv.

Usage:
    cd /Users/hqa668/calcshore-thermal-2d
    python validation/sprint9/stage1_pilot/scripts/run_engine_pilot.py
"""

import csv
import sys
from pathlib import Path

import numpy as np

# ---- Path setup ----
HERE  = Path(__file__).resolve().parent
S9    = HERE.parent                              # sprint9/stage1_pilot/
ROOT  = (S9 / "../../..").resolve()             # calcshore-thermal-2d/
SC    = ROOT / "validation" / "soil_calibration"

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

# ---- Sprint 9 constants ----
HU_FLOOR_THRESHOLD = 10_000.0   # J/kg
HU_RESIDUAL        = 12_937.0   # J/kg
T_PL_F             = 73.0       # °F  (both mixes)
T_SO_F             = 85.0       # °F  (both mixes)
DURATION_HR        = 168.0
OUTPUT_INTERVAL_S  = 3600.0     # hourly

# c(T_pl=73) — computed once; must equal ≈1.0532
C_VAL = -1.3025 + 0.04746 * T_PL_F - 2.081e-4 * T_PL_F**2

# §2.2 expected values for 0.1% tolerance check
# Note: Hu_eff is NOT checked here — it depends on compute_hu_factor which is
# composition-derived and verified by the wrapper log output.
EXPECTED_22 = {
    "mix01": {"alpha_u_raw": 0.7585, "Hu_raw": 424143.0, "alpha_u_eff": 0.7989},
    "mix07": {"alpha_u_raw": 0.8935, "Hu_raw": 463076.0, "alpha_u_eff": 0.9410},
}

# Post-run sanity check ranges
# Corrected bounds for 80ft deep mat with realistic Hu (424–463 kJ/kg).
# Original brief bounds ([115–135]°F, [110–130]°F) were calibrated for thin mat
# or suppressed Hu. CW itself produces 149°F core for these inputs. (§0 of report.)
SANITY = {
    "mix01": {"T_core_lo": 140.0, "T_core_hi": 165.0, "alpha_lo": 0.60, "alpha_hi": 0.85},
    "mix07": {"T_core_lo": 135.0, "T_core_hi": 170.0, "alpha_lo": 0.55, "alpha_hi": 0.85},
}

DATA_DIR = S9 / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def pct_delta(computed, expected):
    if expected == 0:
        return abs(computed)
    return abs(computed - expected) / abs(expected)


def check_tol(name, computed, expected, tol=0.001):
    d = pct_delta(computed, expected)
    ok = d <= tol
    tag = "OK" if ok else f"FAIL (delta={d*100:.4f}%, tol={tol*100:.1f}%)"
    print(f"    {name}: computed={computed:.6g}, expected={expected:.6g} → {tag}")
    if not ok:
        raise RuntimeError(
            f"HALT §5 pre-run: {name} tolerance failed: {computed:.6g} vs {expected:.6g} "
            f"(delta={d*100:.4f}%, tol={tol*100:.1f}%)"
        )


def run_mix(mix_name: str, cw_folder: Path) -> dict:
    print(f"\n{'='*70}")
    print(f"Mix: {mix_name.upper()}")

    # ---- §3.1 Parse raw CW input ----
    mix, geom, constr, _ = parse_cw_dat(str(cw_folder / "input.dat"))

    alpha_u_raw = mix.alpha_u
    Hu_raw      = mix.Hu_J_kg
    tau_hrs     = mix.tau_hrs
    beta        = mix.beta
    Ea          = mix.activation_energy_J_mol

    # ---- §2.2 Step 1: composition-based Hu scaling (Sprint 7/8 standard) ----
    fac, fac_note = compute_hu_factor(mix)

    # ---- §2.2 Step 2: Hu_residual conditional (applied to raw CW Hu before factor) ----
    if Hu_raw < HU_FLOOR_THRESHOLD:
        Hu_base  = HU_RESIDUAL
        cond_txt = f"APPLIED (Hu_raw={Hu_raw:.0f} < {HU_FLOOR_THRESHOLD:.0f})"
    else:
        Hu_base  = Hu_raw
        cond_txt = f"passthrough (Hu_raw={Hu_raw:.0f} >= {HU_FLOOR_THRESHOLD:.0f})"

    # ---- Apply composition factor to produce factored Hu ----
    Hu_factored = Hu_base * fac

    # ---- §2.2 Step 3/4: B1 corrections ----
    alpha_u_eff = C_VAL * alpha_u_raw
    Hu_eff      = Hu_factored / C_VAL

    # ---- Pre-run wrapper log (§3.3 required format) ----
    print(f"Raw from input.dat: alpha_u={alpha_u_raw}, Hu={Hu_raw:.0f} J/kg, "
          f"tau={tau_hrs:.1f} hr, beta={beta:.3f}, Ea={Ea:.0f} J/mol")
    print(f"Hu_residual: {cond_txt}")
    print(f"Hu_factor (compute_hu_factor): {fac:.6f}  ({fac_note or 'composition-derived'})")
    print(f"Hu_factored = {Hu_base:.0f} x {fac:.6f} = {Hu_factored:.1f} J/kg")
    print(f"c(T_pl={T_PL_F:.0f}) = {C_VAL:.4f}")
    print(f"alpha_u_effective = {C_VAL:.4f} x {alpha_u_raw} = {alpha_u_eff:.4f}")
    print(f"Hu_J_kg_effective = {Hu_factored:.1f} / {C_VAL:.4f} = {Hu_eff:.1f}")
    print(f"Engine settings: model_soil=False, is_submerged=True, blanket=0.0, k_uc×0.96")

    # ---- §5 Pre-run tolerance checks ----
    exp = EXPECTED_22[mix_name]
    print("  Pre-run tolerance checks (0.1% tol):")
    check_tol("c(73)", C_VAL, 1.0532)
    check_tol("alpha_u_eff", alpha_u_eff, exp["alpha_u_eff"])
    # Confirm Hu_residual did NOT fire (passthrough for realistic mixes)
    if Hu_raw < HU_FLOOR_THRESHOLD:
        raise RuntimeError(
            f"HALT §5: Hu_residual conditional FIRED for {mix_name} "
            f"(Hu_raw={Hu_raw}). This should not happen for realistic mixes."
        )
    print(f"    Hu_residual: passthrough confirmed (Hu_raw={Hu_raw:.0f} >> 10,000) → OK")
    print(f"    Hu_factor={fac:.6f}, Hu_factored={Hu_factored:.1f}, Hu_eff={Hu_eff:.1f} → (logged, not tolerance-checked)")

    # ---- Apply B1 to mix object ----
    mix.alpha_u           = alpha_u_eff
    mix.Hu_J_kg_effective = Hu_eff
    # τ, β, Ea pass through unmodified (already on mix object from parser)

    # ---- Override engine settings ----
    constr.model_soil   = False
    constr.is_submerged = True

    # ---- Build grid (6× refinement, default) ----
    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )

    T0_C = (T_PL_F - 32.0) * 5.0 / 9.0
    Ts_C = (T_SO_F - 32.0) * 5.0 / 9.0
    Ti   = np.full((grid.ny, grid.nx), T0_C)
    if hasattr(grid, "is_soil") and grid.is_soil is not None:
        Ti[grid.is_soil] = Ts_C

    # ---- Engine run ----
    print(f"\n  Running engine ({DURATION_HR:.0f} hr, hourly outputs)...")
    res = solve_hydration_2d(
        grid, mix, Ti,
        duration_s=DURATION_HR * 3600.0,
        output_interval_s=OUTPUT_INTERVAL_S,
        boundary_mode="full_2d",
        environment=make_neutral_env(T_PL_F),
        construction=constr,
        T_ground_deep_C=Ts_C,
        diagnostic_outputs=False,
    )
    t_hrs = np.asarray(res.t_s) / 3600.0
    print(f"  Engine done. n_output_times={len(t_hrs)}, "
          f"t_range=[{t_hrs[0]:.1f}, {t_hrs[-1]:.1f}] hr")

    # ---- Load CW output for grid coordinates ----
    cw_v = parse_cw_temp_output(str(cw_folder / "output.txt"))
    cw_depths_m = cw_v.depths_m
    cw_widths_m = cw_v.widths_m
    nD = len(cw_depths_m)
    nW = len(cw_widths_m)
    print(f"  CW grid: {nD} depth pts × {nW} width pts")

    # ---- Resample engine trajectory to CW grid (hourly) ----
    jsl, isl = grid.concrete_slice()
    eng_y = grid.y[jsl]
    eng_x = grid.x[isl]

    T_engine_F   = np.zeros((len(t_hrs), nD, nW), dtype=np.float32)
    alpha_engine = np.zeros((len(t_hrs), nD, nW), dtype=np.float32)

    for ti in range(len(t_hrs)):
        slab_T_F = res.T_field_C[ti, jsl, isl] * 9.0 / 5.0 + 32.0
        T_engine_F[ti]   = resample_engine_to_cw(eng_y, eng_x, slab_T_F,
                                                  cw_depths_m, cw_widths_m)
        slab_alpha = res.alpha_field[ti, jsl, isl]
        alpha_engine[ti] = resample_engine_to_cw(eng_y, eng_x, slab_alpha,
                                                  cw_depths_m, cw_widths_m)

    # ---- §5 Post-run sanity checks ----
    def ti_at(hr):
        return int(np.abs(t_hrs - hr).argmin())

    ti0   = ti_at(0.0)
    ti168 = ti_at(168.0)

    # Check 1: IC — bottom-half centerline at t=0 must equal T_pl
    T_init_bottom = T_engine_F[ti0, 24:49, 0]  # di ∈ [24,48], wi=0
    IC_max_dev = float(np.max(np.abs(T_init_bottom - T_PL_F)))

    # Check 2: bottom Dirichlet at t=168
    T_bot_168 = float(T_engine_F[ti168, 48, 0])

    # Check 3: core temperature at t=168
    T_core_168 = float(T_engine_F[ti168, 24, 0])

    # Check 4: core alpha at t=168
    alpha_core_168 = float(alpha_engine[ti168, 24, 0])

    san = SANITY[mix_name]
    print(f"\n  Post-run sanity checks:")
    IC_ok   = IC_max_dev  <= 0.05
    bot_ok  = abs(T_bot_168 - T_SO_F) <= 0.5
    core_ok = san["T_core_lo"] <= T_core_168 <= san["T_core_hi"]
    alp_ok  = san["alpha_lo"]  <= alpha_core_168 <= san["alpha_hi"]

    print(f"    T(di=24..48,t=0,wi=0) max dev from {T_PL_F}°F: "
          f"{IC_max_dev:.4f}°F (tol ±0.05°F) → {'OK' if IC_ok else 'FAIL'}")
    print(f"    T(di=48,t=168,wi=0): {T_bot_168:.4f}°F "
          f"(expect {T_SO_F}±0.5°F) → {'OK' if bot_ok else 'FAIL'}")
    print(f"    T_core(di=24,t=168,wi=0): {T_core_168:.4f}°F "
          f"(expect [{san['T_core_lo']},{san['T_core_hi']}]°F) → {'OK' if core_ok else 'FAIL'}")
    print(f"    α(di=24,t=168,wi=0): {alpha_core_168:.4f} "
          f"(expect [{san['alpha_lo']},{san['alpha_hi']}]) → {'OK' if alp_ok else 'FAIL'}")

    if not IC_ok:
        raise RuntimeError(
            f"HALT §5 post-run: IC check failed for {mix_name}: "
            f"max dev={IC_max_dev:.4f}°F > 0.05°F"
        )
    if not bot_ok:
        raise RuntimeError(
            f"HALT §5 post-run: Bottom Dirichlet check failed for {mix_name}: "
            f"T={T_bot_168:.4f}°F, expected {T_SO_F}±0.5°F"
        )
    if not core_ok:
        raise RuntimeError(
            f"HALT §5 post-run: Core temperature check failed for {mix_name}: "
            f"T_core={T_core_168:.4f}°F, expected [{san['T_core_lo']},{san['T_core_hi']}]°F. "
            f"If T_core ≈ 73°F, hydration heat path did not fire — wrapper broken."
        )
    if not alp_ok:
        raise RuntimeError(
            f"HALT §5 post-run: Alpha check failed for {mix_name}: "
            f"alpha={alpha_core_168:.4f}, expected [{san['alpha_lo']},{san['alpha_hi']}]. "
            f"If MIX-07 alpha > 0.80, kinetics too fast — check τ/β units."
        )
    print(f"\n  All sanity checks PASSED for {mix_name.upper()}.")

    # ---- Save engine trajectory .npz ----
    out_path = DATA_DIR / f"{mix_name}_engine_traj.npz"
    np.savez_compressed(
        out_path,
        t_hrs=t_hrs,
        T_engine_F=T_engine_F,
        alpha_field=alpha_engine,
        cw_depths_m=cw_depths_m,
        cw_widths_m=cw_widths_m,
    )
    print(f"  Saved trajectory → {out_path}")

    return {
        "mix":              mix_name,
        "alpha_u_raw":      alpha_u_raw,
        "Hu_raw":           Hu_raw,
        "tau_hrs":          tau_hrs,
        "beta":             beta,
        "Ea_J_mol":         Ea,
        "hu_factor":        fac,
        "Hu_factored":      Hu_factored,
        "c_val":            C_VAL,
        "conditional":      cond_txt,
        "alpha_u_eff":      alpha_u_eff,
        "Hu_eff":           Hu_eff,
        "T_IC_max_dev":     IC_max_dev,
        "T_bot_168":        T_bot_168,
        "T_core_168":       T_core_168,
        "alpha_core_168":   alpha_core_168,
        "IC_ok":            IC_ok,
        "bot_ok":           bot_ok,
        "core_ok":          core_ok,
        "alpha_ok":         alp_ok,
    }


def write_dataset_inventory(mix01_row, mix07_row):
    path = DATA_DIR / "pilot_dataset_inventory.csv"
    fields = ["mix", "T_pl_F", "T_soil_F", "alpha_u_raw", "Hu_raw_J_kg",
              "tau_hrs", "beta", "Ea_J_mol", "cement_lb_yd3", "wcm",
              "geometry_width_ft", "geometry_depth_ft", "T_ambient_line440"]

    rows = []
    for mix_name, folder in [("mix01", S9 / "cw_data/mix01"),
                              ("mix07", S9 / "cw_data/mix07")]:
        mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))
        # Read line 440 (T_ambient override) with latin-1
        with open(folder / "input.dat", encoding="latin-1") as f:
            raw_lines = f.readlines()
        t_ambient_line = raw_lines[440].strip() if len(raw_lines) > 440 else "N/A"
        rows.append({
            "mix":                  mix_name,
            "T_pl_F":               constr.placement_temp_F,
            "T_soil_F":             constr.soil_temp_F,
            "alpha_u_raw":          mix.alpha_u,
            "Hu_raw_J_kg":          mix.Hu_J_kg,
            "tau_hrs":              mix.tau_hrs,
            "beta":                 mix.beta,
            "Ea_J_mol":             mix.activation_energy_J_mol,
            "cement_lb_yd3":        mix.cement_type_I_II_lb_yd3,
            "wcm":                  f"{mix.wcm:.4f}",
            "geometry_width_ft":    geom.width_ft,
            "geometry_depth_ft":    geom.depth_ft,
            "T_ambient_line440":    t_ambient_line,
        })

    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {path}")


def write_wrapper_log(logs: list):
    path = DATA_DIR / "pilot_wrapper_log.csv"
    fields = [
        "mix", "alpha_u_raw", "Hu_raw", "tau_hrs", "beta", "Ea_J_mol",
        "hu_factor", "Hu_factored", "c_val", "conditional", "alpha_u_eff", "Hu_eff",
        "T_IC_max_dev", "T_bot_168", "T_core_168", "alpha_core_168",
        "IC_ok", "bot_ok", "core_ok", "alpha_ok",
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(logs)
    print(f"Wrote {path}")


def main():
    print("Sprint 9 Stage 1-pilot — Engine Runner (B1 wrapper)")
    print(f"c(T_pl=73) = {C_VAL:.6f}  (expect ≈1.0532)")
    print(f"HU_FLOOR_THRESHOLD = {HU_FLOOR_THRESHOLD}  HU_RESIDUAL = {HU_RESIDUAL}")

    mix01_log = run_mix("mix01", S9 / "cw_data/mix01")
    mix07_log = run_mix("mix07", S9 / "cw_data/mix07")

    write_dataset_inventory(mix01_log, mix07_log)
    write_wrapper_log([mix01_log, mix07_log])

    print("\n" + "="*70)
    print("Both engine runs complete. All sanity checks passed.")
    print(f"Trajectory files: {DATA_DIR}/mix0X_engine_traj.npz")
    print("Next: run analyze_pilot.py")


if __name__ == "__main__":
    main()
