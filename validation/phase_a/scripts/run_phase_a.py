#!/usr/bin/env python3
"""Phase A — full 25-cell engineering-tolerance validation sweep.

Applies the Sprint 9 closing wrapper (unchanged) to each (mix, scenario) cell
and computes T_max, ΔT_max, max|R|, and supporting diagnostic metrics.

Usage:
    cd /Users/hqa668/calcshore-thermal-2d
    python validation/phase_a/scripts/run_phase_a.py
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np

HERE  = Path(__file__).resolve().parent
ROOT  = (HERE / "../../..").resolve()
SC    = ROOT / "validation" / "soil_calibration"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output
from thermal_engine_2d import solve_hydration_2d, build_grid_half_mat
from kinetics_correction import compute_hu_factor
from stage3_compare import resample_engine_to_cw
from stage4b_run import make_neutral_env

HU_FLOOR_THRESHOLD = 10_000.0
HU_RESIDUAL        = 12_937.0
DURATION_HR        = 168.0
OUTPUT_INTERVAL_S  = 3600.0

CW_DIR  = ROOT / "validation" / "phase_a" / "cw_data"
NPZ_DIR = ROOT / "validation" / "phase_a" / "data"
RES_DIR = ROOT / "validation" / "phase_a" / "results"

MIXES     = ["mix01", "mix04", "mix06", "mix07", "mix10"]
SCENARIOS = ["S1", "S2", "S3", "S4", "S5"]


def run_phase_a_cell(mix_id: str, scenario_id: str, cw_dir: Path, output_dir: Path) -> dict:
    folder = cw_dir / f"{mix_id}_{scenario_id}"
    mix, geom, constr, _ = parse_cw_dat(str(folder / "input.dat"))

    T_pl_F = constr.placement_temp_F
    T_so_F = constr.soil_temp_F

    # Wrapper — verbatim Sprint 9 pipeline (no modifications allowed)
    Hu_raw = mix.Hu_J_kg
    fac, _ = compute_hu_factor(mix)
    if Hu_raw < HU_FLOOR_THRESHOLD:
        Hu_base = HU_RESIDUAL
    else:
        Hu_base = Hu_raw
    Hu_factored = Hu_base * fac

    c = -1.3025 + 0.04746 * T_pl_F - 2.081e-4 * T_pl_F**2
    alpha_u_eff = c * mix.alpha_u
    Hu_eff      = Hu_factored / c

    mix.alpha_u           = alpha_u_eff
    mix.Hu_J_kg_effective = Hu_eff
    constr.model_soil     = False
    constr.is_submerged   = True

    grid = build_grid_half_mat(
        geom.width_ft, geom.depth_ft,
        is_submerged=True, model_soil=False, blanket_thickness_m=0.0,
    )
    T0_C = (T_pl_F - 32.0) * 5.0 / 9.0
    Ts_C = (T_so_F - 32.0) * 5.0 / 9.0
    Ti = np.full((grid.ny, grid.nx), T0_C)

    res = solve_hydration_2d(
        grid, mix, Ti,
        duration_s=DURATION_HR * 3600.0,
        output_interval_s=OUTPUT_INTERVAL_S,
        boundary_mode="full_2d",
        environment=make_neutral_env(T_pl_F),
        construction=constr,
        T_ground_deep_C=Ts_C,
        diagnostic_outputs=False,
    )
    t_hrs = np.asarray(res.t_s) / 3600.0

    # Load CW
    cw_v = parse_cw_temp_output(str(folder / "output.txt"))
    cw_depths_m = cw_v.depths_m
    cw_widths_m = cw_v.widths_m
    nD = len(cw_depths_m)
    nW = len(cw_widths_m)

    # Resample engine → CW grid (hourly)
    jsl, isl = grid.concrete_slice()
    eng_y = grid.y[jsl]
    eng_x = grid.x[isl]
    T_engine_F   = np.zeros((len(t_hrs), nD, nW), dtype=np.float32)
    alpha_engine = np.zeros((len(t_hrs), nD, nW), dtype=np.float32)
    for ti in range(len(t_hrs)):
        slab_T_F = res.T_field_C[ti, jsl, isl] * 9.0 / 5.0 + 32.0
        T_engine_F[ti]   = resample_engine_to_cw(eng_y, eng_x, slab_T_F, cw_depths_m, cw_widths_m)
        alpha_engine[ti] = resample_engine_to_cw(eng_y, eng_x, res.alpha_field[ti, jsl, isl], cw_depths_m, cw_widths_m)

    # CW time alignment
    cw_t_hrs = cw_v.time_hrs
    cw_T_F   = cw_v.T_field_F  # (n_time_cw, nD, nW)

    # Align CW time to engine times by nearest index (both are hourly, same duration)
    def align_cw_to_eng(cw_times, eng_times):
        idxs = []
        for t in eng_times:
            idxs.append(int(np.abs(cw_times - t).argmin()))
        return np.array(idxs)

    cw_idxs = align_cw_to_eng(cw_t_hrs, t_hrs)
    cw_T_aligned = cw_T_F[cw_idxs]  # (len(t_hrs), nD, nW)

    ti168     = int(np.abs(t_hrs - 168.0).argmin())
    cw_ti168  = int(np.abs(cw_t_hrs - 168.0).argmin())

    # ── Validity mask: di ∈ [5, 48] (per §6 of brief) ──
    VALID_DI_LO, VALID_DI_HI = 5, 49  # Python slice end-exclusive

    eng_valid = T_engine_F[:, VALID_DI_LO:VALID_DI_HI, :]  # (n_t, 44, nW)
    cw_valid  = cw_T_aligned[:, VALID_DI_LO:VALID_DI_HI, :]

    # T_max (primary gate metric)
    flat_eng = eng_valid.reshape(len(t_hrs), -1)
    flat_cw  = cw_valid.reshape(len(t_hrs), -1)

    T_max_eng = float(flat_eng.max())
    T_max_cw  = float(flat_cw.max())
    T_max_diff = T_max_eng - T_max_cw

    flat_eng_max_idx = flat_eng.argmax()
    ti_tmax_eng, di_wi_tmax_eng = divmod(int(flat_eng_max_idx), flat_eng.shape[1])
    di_tmax_eng, wi_tmax_eng = divmod(di_wi_tmax_eng, nW)
    di_tmax_eng += VALID_DI_LO

    flat_cw_max_idx = flat_cw.argmax()
    ti_tmax_cw, di_wi_tmax_cw = divmod(int(flat_cw_max_idx), flat_cw.shape[1])
    di_tmax_cw, wi_tmax_cw = divmod(di_wi_tmax_cw, nW)
    di_tmax_cw += VALID_DI_LO

    T_max_eng_t_hr = float(t_hrs[ti_tmax_eng])
    T_max_cw_t_hr  = float(cw_t_hrs[cw_idxs[ti_tmax_cw]])

    # ΔT_max (primary gate metric): max over t of [max - min over di∈[5,48], all wi]
    dT_per_t_eng = flat_eng.max(axis=1) - flat_eng.min(axis=1)  # (n_t,)
    dT_per_t_cw  = flat_cw.max(axis=1) - flat_cw.min(axis=1)

    dT_max_eng = float(dT_per_t_eng.max())
    dT_max_cw  = float(dT_per_t_cw.max())
    dT_max_diff = dT_max_eng - dT_max_cw

    # Scalar metrics at t=168
    T_core_eng_t168  = float(T_engine_F[ti168, 24, 0])
    T_core_cw_t168   = float(cw_T_F[cw_ti168, 24, 0])
    T_core_diff      = T_core_eng_t168 - T_core_cw_t168

    T_di48_eng_t168  = float(T_engine_F[ti168, 48, 0])
    T_di48_cw_t168   = float(cw_T_F[cw_ti168, 48, 0])

    alpha_core_eng_t168 = float(alpha_engine[ti168, 24, 0])

    # Diagnostic: max|R| over di∈[24,48], wi=0, t=168
    R_diag      = T_engine_F[ti168, 24:49, 0] - cw_T_F[cw_ti168, 24:49, 0]
    max_R_diag  = float(np.max(np.abs(R_diag)))
    di_at_max_R = int(np.argmax(np.abs(R_diag))) + 24

    # R_di48_sign_changes: sign changes in R(di=48, t) over all t
    R_di48_ts = T_engine_F[:, 48, 0] - cw_T_aligned[:, 48, 0]
    signs      = np.sign(R_di48_ts[R_di48_ts != 0])
    if len(signs) > 1:
        R_di48_sign_changes = int(np.sum(signs[1:] != signs[:-1]))
    else:
        R_di48_sign_changes = 0

    # ── Per-cell sanity checks (warn only) ──
    ti0 = int(np.abs(t_hrs - 0.0).argmin())
    IC_dev = float(np.max(np.abs(T_engine_F[ti0, 24:49, 0] - constr.placement_temp_F)))
    warns = []
    if IC_dev > 0.05:
        warns.append(f"IC_dev={IC_dev:.4f}°F > 0.05")
    if abs(T_di48_eng_t168 - T_so_F) > 0.5:
        warns.append(f"T_di48_t168={T_di48_eng_t168:.2f}°F vs T_soil={T_so_F} (diff={abs(T_di48_eng_t168-T_so_F):.2f})")
    if not (120.0 <= T_core_eng_t168 <= 180.0):
        warns.append(f"T_core_t168={T_core_eng_t168:.2f}°F outside [120,180]")
    if not (0.4 <= alpha_core_eng_t168 <= 0.95):
        warns.append(f"alpha_core_t168={alpha_core_eng_t168:.4f} outside [0.4,0.95]")
    for w in warns:
        print(f"    WARN: {mix_id}_{scenario_id}: {w}")

    # Save .npz
    npz_path = output_dir / f"{mix_id}_{scenario_id}.npz"
    output_dir.mkdir(parents=True, exist_ok=True)
    # Build full residual field on CW grid (use aligned CW times)
    residual_field = T_engine_F - cw_T_aligned
    np.savez_compressed(
        npz_path,
        t_hrs=t_hrs,
        T_engine_F=T_engine_F,
        T_cw_F=cw_T_aligned,
        residual_field=residual_field,
        cw_depths_m=cw_depths_m,
        cw_widths_m=cw_widths_m,
        alpha_engine=alpha_engine,
    )

    return dict(
        mix_id=mix_id, scenario_id=scenario_id,
        T_max_eng=T_max_eng, T_max_cw=T_max_cw, T_max_diff=T_max_diff,
        dT_max_eng=dT_max_eng, dT_max_cw=dT_max_cw, dT_max_diff=dT_max_diff,
        T_max_eng_t_hr=T_max_eng_t_hr, T_max_cw_t_hr=T_max_cw_t_hr,
        T_max_eng_di=di_tmax_eng, T_max_eng_wi=wi_tmax_eng,
        T_max_cw_di=di_tmax_cw, T_max_cw_wi=wi_tmax_cw,
        T_core_eng_t168=T_core_eng_t168, T_core_cw_t168=T_core_cw_t168,
        T_core_diff=T_core_diff,
        T_di48_eng_t168=T_di48_eng_t168, T_di48_cw_t168=T_di48_cw_t168,
        alpha_core_eng_t168=alpha_core_eng_t168,
        max_R_diag=max_R_diag, di_at_max_R=di_at_max_R,
        R_di48_sign_changes=R_di48_sign_changes,
        hu_factor=fac, c_val=c,
        T_pl_F=T_pl_F, T_so_F=T_so_F,
        IC_dev=IC_dev,
        warns="; ".join(warns) if warns else "",
    )


def main():
    NPZ_DIR.mkdir(parents=True, exist_ok=True)
    RES_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    cell_n  = 0
    t_sweep = time.time()

    for mix_id in MIXES:
        for scen_id in SCENARIOS:
            cell_n += 1
            t0 = time.time()
            print(f"[{cell_n:2d}/25] {mix_id}_{scen_id} ...", end=" ", flush=True)
            row = run_phase_a_cell(mix_id, scen_id, CW_DIR, NPZ_DIR)
            elapsed = time.time() - t0
            gate_T  = "PASS" if abs(row["T_max_diff"])  < 2.0 else "FAIL"
            gate_dT = "PASS" if abs(row["dT_max_diff"]) < 2.0 else "FAIL"
            print(f"done in {elapsed:.1f}s  |T_max_diff|={abs(row['T_max_diff']):.2f}°F {gate_T}  "
                  f"|dT_max_diff|={abs(row['dT_max_diff']):.2f}°F {gate_dT}")
            results.append(row)

            if cell_n % 5 == 0:
                print(f"  --- progress: {cell_n}/25 cells complete ---")

    total_elapsed = time.time() - t_sweep
    print(f"\nAll 25 cells done in {total_elapsed:.1f}s.")

    # Write aggregate CSV
    fields = [
        "mix_id", "scenario_id",
        "T_max_eng", "T_max_cw", "T_max_diff",
        "dT_max_eng", "dT_max_cw", "dT_max_diff",
        "T_max_eng_t_hr", "T_max_cw_t_hr",
        "T_max_eng_di", "T_max_eng_wi", "T_max_cw_di", "T_max_cw_wi",
        "T_core_eng_t168", "T_core_cw_t168", "T_core_diff",
        "T_di48_eng_t168", "T_di48_cw_t168",
        "alpha_core_eng_t168",
        "max_R_diag", "di_at_max_R",
        "R_di48_sign_changes",
        "hu_factor", "c_val", "T_pl_F", "T_so_F", "IC_dev", "warns",
    ]
    csv_path = RES_DIR / "phase_a_results.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in results:
            w.writerow({k: (f"{row[k]:.4f}" if isinstance(row[k], float) else row[k])
                        for k in fields})
    print(f"Written: {csv_path}")
    return results


if __name__ == "__main__":
    main()
