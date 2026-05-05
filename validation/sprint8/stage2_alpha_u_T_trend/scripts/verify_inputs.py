#!/usr/bin/env python3
"""§3.2 — Verify 32-dataset input.dat consistency for Stage 2-alpha-u-T-trend.

Checks for each folder:
  - T_pl (placement_temp_F) == T_pc encoded in folder name
  - T_soil (soil_temp_F) == T_pc
  - mix.Hu_J_kg == 1.0 J/kg
  - mix.tau_hrs == 5.0 hr
  - mix.beta == 0.85
  - geom.width_ft == 40, geom.depth_ft == 80
  - mix.alpha_u is in the expected range for the alpha-group

Writes data/inputs_consistency.csv.
Prints FAIL and exits with code 1 on any inconsistency.
"""
import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from cw_scenario_loader import parse_cw_dat

CW_DATA = HERE / "cw_data"
DATA    = HERE / "data"

ALPHA_TAG_MAP = {"02": 0.20, "04": 0.40, "06": 0.60, "08": 0.80}
ALPHA_TOLERANCE  = 0.05   # mix.alpha_u must be within ±0.05 of nominal
TEMP_TOLERANCE   = 0.5    # °F — placement vs folder-name T_pc
KINETICS_RTOL    = 0.01   # 1% relative tolerance for Hu/tau/beta
GEOM_TOLERANCE   = 0.5    # ft


def check_close(val, target, tol, name, folder, errors):
    if abs(val - target) > tol:
        errors.append(f"  {folder}: {name} = {val:.4f}, expected ≈ {target} (tol ±{tol})")


def main():
    DATA.mkdir(parents=True, exist_ok=True)

    folders = sorted(CW_DATA.glob("thermal_temp_alpha*"))
    if not folders:
        print(f"ERROR: no folders found under {CW_DATA}")
        sys.exit(1)

    print(f"§3.2 Input consistency check — {len(folders)} folders")
    print(f"{'Folder':<35} {'T_pc':>5} {'α_u_nom':>7} {'T_pl':>6} {'T_soil':>6} "
          f"{'α_u_act':>8} {'Hu':>6} {'τ':>5} {'β':>5} "
          f"{'w_ft':>5} {'d_ft':>5}  status")
    print("-" * 120)

    rows = []
    all_errors = []

    for folder in folders:
        dat = folder / "input.dat"
        if not dat.exists():
            all_errors.append(f"  {folder.name}: input.dat missing")
            continue

        # Parse folder name: thermal_temp_alpha02_70 → alpha_tag=02, T_pc=70
        m = re.match(r"thermal_temp_alpha(\d{2})_(\d+)$", folder.name)
        if not m:
            all_errors.append(f"  {folder.name}: unexpected name pattern")
            continue
        alpha_tag = m.group(1)
        T_pc      = float(m.group(2))
        alpha_nom = ALPHA_TAG_MAP.get(alpha_tag)
        if alpha_nom is None:
            all_errors.append(f"  {folder.name}: unknown alpha tag '{alpha_tag}'")
            continue

        mix, geom, constr, _ = parse_cw_dat(str(dat))

        errors = []
        check_close(constr.placement_temp_F, T_pc, TEMP_TOLERANCE,
                    "placement_temp_F", folder.name, errors)
        check_close(constr.soil_temp_F, T_pc, TEMP_TOLERANCE,
                    "soil_temp_F", folder.name, errors)
        check_close(mix.Hu_J_kg, 1.0, KINETICS_RTOL * 1.0 + 0.01,
                    "Hu_J_kg", folder.name, errors)
        check_close(mix.tau_hrs, 5.0, KINETICS_RTOL * 5.0,
                    "tau_hrs", folder.name, errors)
        check_close(mix.beta, 0.85, KINETICS_RTOL * 0.85,
                    "beta", folder.name, errors)
        check_close(geom.width_ft, 40.0, GEOM_TOLERANCE,
                    "width_ft", folder.name, errors)
        check_close(geom.depth_ft, 80.0, GEOM_TOLERANCE,
                    "depth_ft", folder.name, errors)
        check_close(mix.alpha_u, alpha_nom, ALPHA_TOLERANCE,
                    "alpha_u", folder.name, errors)

        status = "OK" if not errors else "FAIL"
        print(f"{folder.name:<35} {T_pc:>5.0f} {alpha_nom:>7.2f} "
              f"{constr.placement_temp_F:>6.1f} {constr.soil_temp_F:>6.1f} "
              f"{mix.alpha_u:>8.5f} {mix.Hu_J_kg:>6.3f} {mix.tau_hrs:>5.2f} "
              f"{mix.beta:>5.3f} {geom.width_ft:>5.1f} {geom.depth_ft:>5.1f}  [{status}]")

        if errors:
            for e in errors:
                print(e)
            all_errors.extend(errors)

        rows.append({
            "folder":           folder.name,
            "alpha_tag":        alpha_tag,
            "alpha_nom":        alpha_nom,
            "T_pc_F":           T_pc,
            "placement_temp_F": constr.placement_temp_F,
            "soil_temp_F":      constr.soil_temp_F,
            "alpha_u_actual":   mix.alpha_u,
            "Hu_J_kg":          mix.Hu_J_kg,
            "tau_hrs":          mix.tau_hrs,
            "beta":             mix.beta,
            "width_ft":         geom.width_ft,
            "depth_ft":         geom.depth_ft,
            "status":           status,
        })

    # Write CSV
    out_csv = DATA / "inputs_consistency.csv"
    fields = ["folder", "alpha_tag", "alpha_nom", "T_pc_F",
              "placement_temp_F", "soil_temp_F", "alpha_u_actual",
              "Hu_J_kg", "tau_hrs", "beta", "width_ft", "depth_ft", "status"]
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {out_csv} ({len(rows)} rows)")

    if all_errors:
        print(f"\nFAIL — {len(all_errors)} inconsistency/-ies found:")
        for e in all_errors:
            print(e)
        print("Stopping. Fix input datasets before proceeding to §3.3.")
        sys.exit(1)
    else:
        print(f"\nPASS — all {len(rows)} datasets consistent.")


if __name__ == "__main__":
    main()
