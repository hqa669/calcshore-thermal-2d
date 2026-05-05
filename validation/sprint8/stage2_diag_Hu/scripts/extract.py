#!/usr/bin/env python3
"""Sprint 8 Stage 2-diag-Hu §4.1–§4.3 — parse 15 CW Hu-sweep datasets.

Writes:
  data/inputs_consistency.csv   — §4.1 table
  data/Hu_residual_table.csv    — §4.3 table (headline CSV)
  data/T_core_trajectories.npz  — raw (time_hrs, T_core_F) for all runs

Run from validation/sprint8/stage2_diag_Hu/.
"""
import csv
import sys
import os
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parents[1]          # stage2_diag_Hu/
REPO = Path(__file__).resolve().parents[4]          # calcshore-thermal-2d/
sys.path.insert(0, str(REPO))

from cw_scenario_loader import parse_cw_dat, parse_cw_temp_output

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"
DATA.mkdir(exist_ok=True)

EXPECTED_T_PL   = 73.0   # °F
EXPECTED_T_SOIL = 73.0   # °F
EXPECTED_TAU    = 5.0    # hrs
EXPECTED_BETA   = 0.85

# tolerance for IC / parameter checks
TEMP_TOL  = 0.5    # °F  (slightly generous; flag but don't abort on T_pl)
PARAM_TOL = 0.01   # for tau, beta


def find_input_dat(folder: Path) -> Path:
    for name in ("input.dat", "input_test.dat"):
        p = folder / name
        if p.exists():
            return p
    raise FileNotFoundError(f"No input.dat or input_test.dat in {folder}")


def discover_folders(cw_data: Path) -> list[Path]:
    return sorted(d for d in cw_data.iterdir() if d.is_dir())


def diff_input_dats(folders: list[Path]) -> dict:
    """Diff two pairs of input files to confirm which line holds Hu and alpha_u."""
    results = {}

    au02 = [f for f in folders if "au02" in f.name]
    au08 = [f for f in folders if "au08" in f.name]

    def raw_lines(folder):
        p = find_input_dat(folder)
        with open(p, encoding="latin-1") as fh:
            return [ln.rstrip("\r\n").strip() for ln in fh]

    # Find Hu line: compare two au02 folders at very different Hu
    au02_hu1_folder = next((f for f in au02 if "Hu1" in f.name and "Hu10" not in f.name and "Hu100" not in f.name), None)
    au02_hu100k_folder = next((f for f in au02 if "Hu100000" in f.name), None)
    if au02_hu1_folder and au02_hu100k_folder:
        lines_a = raw_lines(au02_hu1_folder)
        lines_b = raw_lines(au02_hu100k_folder)
        diffs = [i+1 for i, (a, b) in enumerate(zip(lines_a, lines_b)) if a.strip() != b.strip()]
        results["Hu_diff_pair"] = (au02_hu1_folder.name, au02_hu100k_folder.name)
        results["Hu_diff_lines"] = diffs
        if len(diffs) == 1:
            results["Hu_file_line"] = diffs[0]
        else:
            results["Hu_file_line"] = f"UNEXPECTED: {diffs}"

    # Find alpha_u line: compare au02 vs au08 at same Hu
    au02_hu1k = next((f for f in au02 if "Hu1000" in f.name), None)
    au08_hu1k = next((f for f in au08 if "Hu1000" in f.name), None)
    if au02_hu1k and au08_hu1k:
        lines_a = raw_lines(au02_hu1k)
        lines_b = raw_lines(au08_hu1k)
        diffs = [i+1 for i, (a, b) in enumerate(zip(lines_a, lines_b)) if a.strip() != b.strip()]
        results["alpha_u_diff_pair"] = (au02_hu1k.name, au08_hu1k.name)
        results["alpha_u_diff_lines"] = diffs
        if len(diffs) == 1:
            results["alpha_u_file_line"] = diffs[0]
        else:
            results["alpha_u_file_line"] = f"UNEXPECTED: {diffs}"

    return results


def parse_folder(folder: Path) -> dict:
    input_path  = find_input_dat(folder)
    output_path = folder / "output.txt"

    # Parse input
    mix, geom, constr, raw = parse_cw_dat(str(input_path))

    # Parse output
    series = parse_cw_temp_output(str(output_path))
    t    = series.time_hrs          # (n_time,)
    T_F  = series.T_core_center_F   # (n_time,) — mid-depth centerline in °F

    # Restrict to t ≤ 168
    mask = t <= 168.01
    t    = t[mask]
    T_F  = T_F[mask]

    # Metrics
    T0        = float(T_F[0])
    T168      = float(np.interp(168.0, t, T_F))
    dT_168_F  = T168 - T0
    dT_max_F  = float(np.max(T_F)) - T0
    t_peak    = float(t[np.argmax(T_F)])

    row = {
        "folder":      folder.name,
        "alpha_u":     mix.alpha_u,
        "Hu_J_kg":     mix.Hu_J_kg,
        "T_pl_F":      constr.placement_temp_F,
        "T_soil_F":    constr.soil_temp_F,
        "tau_hrs":     mix.tau_hrs,
        "beta":        mix.beta,
        "T0_F":        T0,
        "T168_F":      T168,
        "dT_168_F":    dT_168_F,
        "dT_max_F":    dT_max_F,
        "t_peak_hr":   t_peak,
        "T0_C":        (T0  - 32) * 5/9,
        "T168_C":      (T168 - 32) * 5/9,
        "dT_168_C":    dT_168_F * 5/9,
        "dT_max_C":    dT_max_F  * 5/9,
    }
    return row, t, T_F


def run_assertions(rows: list[dict]) -> bool:
    ok = True
    for r in rows:
        folder = r["folder"]
        if abs(r["T_pl_F"] - EXPECTED_T_PL) > TEMP_TOL:
            print(f"  WARN: {folder} T_pl={r['T_pl_F']:.2f} (expected {EXPECTED_T_PL})", file=sys.stderr)
            ok = False
        if abs(r["T_soil_F"] - EXPECTED_T_SOIL) > TEMP_TOL:
            print(f"  WARN: {folder} T_soil={r['T_soil_F']:.2f} (expected {EXPECTED_T_SOIL})", file=sys.stderr)
            ok = False
        if abs(r["tau_hrs"] - EXPECTED_TAU) > PARAM_TOL:
            print(f"  WARN: {folder} tau={r['tau_hrs']:.4f} (expected {EXPECTED_TAU})", file=sys.stderr)
            ok = False
        if abs(r["beta"] - EXPECTED_BETA) > PARAM_TOL:
            print(f"  WARN: {folder} beta={r['beta']:.4f} (expected {EXPECTED_BETA})", file=sys.stderr)
            ok = False
        if abs(r["T0_F"] - EXPECTED_T_PL) > 1.0:
            print(f"  WARN: {folder} T_core(t=0)={r['T0_F']:.3f} °F — IC not 73 °F", file=sys.stderr)

    # Monotonicity check per alpha_u
    for au in (0.20, 0.80):
        subset = sorted([r for r in rows if abs(r["alpha_u"] - au) < 0.05], key=lambda r: r["Hu_J_kg"])
        for i in range(1, len(subset)):
            if subset[i]["dT_168_F"] < subset[i-1]["dT_168_F"] - 0.05:
                print(f"  WARN: monotonicity violation α_u≈{au}: "
                      f"Hu={subset[i-1]['Hu_J_kg']} dT={subset[i-1]['dT_168_F']:.3f}°F "
                      f"> Hu={subset[i]['Hu_J_kg']} dT={subset[i]['dT_168_F']:.3f}°F", file=sys.stderr)
                ok = False

    # Sanity: Hu≈1 at α_u=0.20 should be << 5°C core warming
    low_hu = [r for r in rows if abs(r["alpha_u"] - 0.20) < 0.05 and r["Hu_J_kg"] < 10]
    for r in low_hu:
        if r["dT_168_C"] > 5.0:
            print(f"  ABORT: {r['folder']} dT_168={r['dT_168_C']:.2f}°C at Hu≈1 — "
                  "core-extraction likely wrong. Check geometry.", file=sys.stderr)
            ok = False
    return ok


def main():
    folders = discover_folders(CW_DATA)
    if not folders:
        print("ERROR: no subfolders found in cw_data/. Copy datasets first.", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(folders)} dataset folders in cw_data/")

    # Diff analysis
    diff_info = diff_input_dats(folders)
    print(f"\ninput.dat field locations (1-based file lines):")
    print(f"  Hu line:    {diff_info.get('Hu_file_line', '?')}  "
          f"(diff between {diff_info.get('Hu_diff_pair', '?')})")
    print(f"  alpha_u line: {diff_info.get('alpha_u_file_line', '?')}  "
          f"(diff between {diff_info.get('alpha_u_diff_pair', '?')})")
    if isinstance(diff_info.get("Hu_diff_lines"), list) and len(diff_info["Hu_diff_lines"]) != 1:
        print(f"  WARNING: expected 1 differing line for Hu, got: {diff_info['Hu_diff_lines']}", file=sys.stderr)
    if isinstance(diff_info.get("alpha_u_diff_lines"), list) and len(diff_info["alpha_u_diff_lines"]) != 1:
        print(f"  WARNING: expected 1 differing line for alpha_u, got: {diff_info['alpha_u_diff_lines']}", file=sys.stderr)

    # Parse all
    rows = []
    trajs = {}
    for folder in folders:
        print(f"  Parsing {folder.name} ...", end="", flush=True)
        try:
            row, t, T_F = parse_folder(folder)
            rows.append(row)
            trajs[folder.name] = {"time_hrs": t, "T_core_F": T_F}
            print(f"  Hu={row['Hu_J_kg']:.2g} α_u={row['alpha_u']:.2f}  "
                  f"dT_168={row['dT_168_F']:.3f}°F ({row['dT_168_C']:.3f}°C)")
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            sys.exit(1)

    # Assertions
    print("\nRunning sanity assertions ...")
    ok = run_assertions(rows)
    if not ok:
        print("  WARN: some assertions failed (see above).")
    else:
        print("  All assertions passed.")

    # Write inputs_consistency.csv
    con_path = DATA / "inputs_consistency.csv"
    con_fields = ["folder", "alpha_u", "Hu_J_kg", "T_pl_F", "T_soil_F", "tau_hrs", "beta"]
    with open(con_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=con_fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in con_fields})
    print(f"\nWrote {con_path}")

    # Write Hu_residual_table.csv
    res_fields = ["folder", "alpha_u", "Hu_J_kg",
                  "T0_F", "T168_F", "dT_168_F", "dT_max_F", "t_peak_hr",
                  "T0_C", "T168_C", "dT_168_C", "dT_max_C"]
    res_path = DATA / "Hu_residual_table.csv"
    with open(res_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=res_fields)
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["alpha_u"], r["Hu_J_kg"])):
            w.writerow({k: r[k] for k in res_fields})
    print(f"Wrote {res_path}")

    # Write T_core_trajectories.npz
    npz_data = {}
    for name, d in trajs.items():
        npz_data[name + "_time_hrs"] = d["time_hrs"]
        npz_data[name + "_T_core_F"] = d["T_core_F"]
    npz_path = DATA / "T_core_trajectories.npz"
    np.savez(npz_path, **npz_data)
    print(f"Wrote {npz_path}")

    # Print summary table
    print("\n§4.3 Heat residual table (sorted by α_u, Hu):")
    hdr = f"{'folder':<36} {'α_u':>5} {'Hu':>9} {'T0°F':>7} {'T168°F':>8} {'dT168°F':>8} {'dT168°C':>8} {'dtmax°C':>8} {'t_pk hr':>8}"
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda r: (r["alpha_u"], r["Hu_J_kg"])):
        print(f"{r['folder']:<36} {r['alpha_u']:>5.2f} {r['Hu_J_kg']:>9.2f} "
              f"{r['T0_F']:>7.3f} {r['T168_F']:>8.3f} {r['dT_168_F']:>8.4f} "
              f"{r['dT_168_C']:>8.4f} {r['dT_max_C']:>8.4f} {r['t_peak_hr']:>8.1f}")

    print(f"\nDone. Files in {DATA}/")


if __name__ == "__main__":
    main()
