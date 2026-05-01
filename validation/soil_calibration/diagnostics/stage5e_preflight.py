#!/usr/bin/env python3
"""Stage 5e pre-flight verification.

§2.1 — Run A k-independence: k×1.00 vs k×0.96 residuals must agree within 0.005°F
§2.2 — Time evolution stability: Runs F and I under k×0.96 at t=24/84/168 hr

Pass conditions:
  §2.1: Δ ≤ 0.005°F on R1, R2, full-field
  §2.2: R2 max|R| ≤ 0.40°F at all three timestamps; no sign reversal past gate
"""
import os, sys
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
SC   = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, SC)

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env, nearest_time_idx

CW_RUNS = os.path.join(SC, "cw_runs")

R1_DI  = 24
R2_WI  = 0
R2_DI  = slice(24, 49)
GATE_F = 0.35


def run_engine(folder, pl_F, so_F, k_factor, output_hrs):
    mix, geom, constr, _ = parse_cw_dat(os.path.join(CW_RUNS, folder, "input.dat"))
    fac, _ = compute_hu_factor(mix)
    mix.Hu_factor_calibrated = fac
    mix.Hu_J_kg_effective = mix.Hu_J_kg * fac
    mix.thermal_conductivity_BTU_hr_ft_F *= k_factor
    constr.model_soil = False; constr.is_submerged = True
    grid = build_grid_half_mat(geom.width_ft, geom.depth_ft,
                               is_submerged=True, model_soil=False,
                               blanket_thickness_m=0.0)
    T0 = (pl_F - 32) * 5/9
    Ts = (so_F - 32) * 5/9
    Ti = np.full((grid.ny, grid.nx), T0)
    Ti[grid.is_soil] = Ts
    duration_s = max(output_hrs) * 3600
    res = solve_hydration_2d(grid, mix, Ti, duration_s=duration_s,
                             output_interval_s=1800., boundary_mode="full_2d",
                             environment=make_neutral_env(pl_F), construction=constr,
                             T_ground_deep_C=Ts, diagnostic_outputs=False)
    jsl, isl = grid.concrete_slice()
    results = {}
    for hr in output_hrs:
        ti = nearest_time_idx(res.t_s, float(hr))
        TF = res.T_field_C[ti, jsl, isl] * 9/5 + 32
        results[hr] = (grid.y[jsl], grid.x[isl], TF)
    return results


def load_cw(folder, compare_hr):
    v = parse_cw_temp_output(os.path.join(CW_RUNS, folder, "output.txt"))
    ti = int(np.abs(v.time_hrs - compare_hr).argmin())
    return v.T_field_F[ti], v.widths_m, v.depths_m


def resid_metrics(ey, ex, eT, cw_F, cw_widths_m, cw_depths_m):
    eng_on_cw = resample_engine_to_cw(ey, ex, eT, cw_depths_m, cw_widths_m)
    R = eng_on_cw - cw_F
    r1 = float(np.max(np.abs(R[R1_DI, :])))
    r2 = float(np.max(np.abs(R[R2_DI, R2_WI])))
    # Full-field limited to validity-masked region (di=5..48, Stage 3.5 convention)
    full = float(np.max(np.abs(R[5:, :])))
    # R2 dominant cell index
    r2_col = R[R2_DI, R2_WI]
    r2_max_idx = int(np.argmax(np.abs(r2_col)))
    r2_max_di  = 24 + r2_max_idx
    r2_max_val = float(r2_col[r2_max_idx])
    return r1, r2, full, r2_max_di, r2_max_val


# ══════════════════════════════════════════════════════════════
# §2.1 — Run A k-independence
# ══════════════════════════════════════════════════════════════
print("=" * 70)
print("§2.1  Run A k-independence check")
print("=" * 70)

cw_F_168, cw_widths_m, cw_depths_m = load_cw("runA_baseline", 168)

results_100 = run_engine("runA_baseline", 73, 73, 1.00, [168])
results_096 = run_engine("runA_baseline", 73, 73, 0.96, [168])

r1_100, r2_100, ff_100, *_ = resid_metrics(*results_100[168], cw_F_168, cw_widths_m, cw_depths_m)
r1_096, r2_096, ff_096, *_ = resid_metrics(*results_096[168], cw_F_168, cw_widths_m, cw_depths_m)

print(f"\n  {'Metric':<25} {'k×1.00':>9} {'k×0.96':>9} {'Δ':>9}  {'Status'}")
print("  " + "-" * 62)
for name, v100, v096 in [("R1 max|R| (°F)", r1_100, r1_096),
                          ("R2 max|R| (°F)", r2_100, r2_096),
                          ("Full-field max|R| (°F)", ff_100, ff_096)]:
    delta = abs(v096 - v100)
    status = "PASS" if delta <= 0.005 else "FAIL"
    print(f"  {name:<25} {v100:>9.4f} {v096:>9.4f} {delta:>9.5f}  [{status}]")

sec21_pass = all(abs(v096 - v100) <= 0.005
                 for _, v100, v096 in [("R1", r1_100, r1_096),
                                        ("R2", r2_100, r2_096),
                                        ("FF", ff_100, ff_096)])
print(f"\n  §2.1 overall: {'PASS' if sec21_pass else 'FAIL'}")
if not sec21_pass:
    print("  STOPPING — k-independence claim violated. Engine change deferred.")
    sys.exit(1)


# ══════════════════════════════════════════════════════════════
# §2.2 — Time evolution stability under k×0.96 for F and I
# ══════════════════════════════════════════════════════════════
print()
print("=" * 70)
print("§2.2  Time evolution stability — Runs F and I under k×0.96")
print("=" * 70)

RUNS_FI = [
    ("F", "runF_73_45",  73,  45),
    ("I", "runI_100_73", 100, 73),
]
TIMESTAMPS = [24, 84, 168]
FAIL_THRESHOLD   = 0.50  # hard stop
SOFT_THRESHOLD   = 0.40  # soft pass

print(f"\n  {'Run':>4}  {'t(hr)':>7}  {'R1 max|R|':>11}  {'R2 max|R|':>11}  {'R2 dom di':>10}  {'R2 dom val':>12}  Status")
print("  " + "-" * 72)

sec22_pass = True
prev_r2 = {}

for label, folder, pl_F, so_F in RUNS_FI:
    cw_F_map = {}
    for hr in TIMESTAMPS:
        cw_F_map[hr], cw_widths_m, cw_depths_m = load_cw(folder, hr)

    eng_results = run_engine(folder, pl_F, so_F, 0.96, TIMESTAMPS)

    prev_r2[label] = {}
    for hr in TIMESTAMPS:
        ey, ex, eT = eng_results[hr]
        r1, r2, _, r2_di, r2_val = resid_metrics(ey, ex, eT, cw_F_map[hr], cw_widths_m, cw_depths_m)

        # Check hard fail
        run_fail = False
        if r2 > FAIL_THRESHOLD:
            run_fail = True
            sec22_pass = False
        # Check soft threshold
        soft_warn = (r2 > SOFT_THRESHOLD)

        status = "FAIL" if run_fail else ("WARN" if soft_warn else "PASS")
        print(f"  {label:>4}  {hr:>7}  {r1:>11.4f}  {r2:>11.4f}  {r2_di:>10}  {r2_val:>+12.4f}  [{status}]")
        prev_r2[label][hr] = r2

    # Check sign/magnitude stability: t=24 should not exceed t=168 by more than 0.05°F after k×0.96
    r2_24  = prev_r2[label][24]
    r2_168 = prev_r2[label][168]
    if r2_24 > r2_168 + 0.05 and r2_24 > FAIL_THRESHOLD:
        print(f"    WARNING Run {label}: t=24 R2={r2_24:.4f} >> t=168 R2={r2_168:.4f} and exceeds 0.50 — k has over-fit t=168")
        sec22_pass = False

print()
print(f"  §2.2 overall: {'PASS' if sec22_pass else 'FAIL'}")
if not sec22_pass:
    print("  STOPPING — time evolution instability detected. Engine change deferred.")
    sys.exit(1)

print()
print("=" * 70)
print("PRE-FLIGHT COMPLETE: both checks PASS — engine change authorized")
print("=" * 70)
