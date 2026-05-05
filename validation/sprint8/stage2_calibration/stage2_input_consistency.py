#!/usr/bin/env python3
"""Sprint 8 Stage 2 §3.1 + §3.2 — dataset copy and input.dat consistency check.

§3.1  Copies 12 folders from ~/Downloads/thermal_alpha_AFI/ into cw_data/.
§3.2  Validates each input.dat matches the Sprint 7 baseline geometry/mix.

Run:  python stage2_input_consistency.py [--no-copy]
  --no-copy  skip copy step (cw_data/ already populated)

Output: findings/input_consistency.md
"""

import argparse
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).parent
ROOT = (HERE / "../../..").resolve()
SC   = (HERE / "../../soil_calibration").resolve()

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SC))

from cw_scenario_loader import parse_cw_dat

CW_DATA  = HERE / "cw_data"
FINDINGS = HERE / "findings"
S7_RUNS  = SC / "cw_runs"
DOWNLOADS_SRC = Path.home() / "Downloads" / "thermal_alpha_AFI"

NEW_FOLDERS = [
    ("thermal_alpha02_A_73_73",  0.20,  73,  73),
    ("thermal_alpha02_F_73_45",  0.20,  73,  45),
    ("thermal_alpha02_I_100_73", 0.20, 100,  73),
    ("thermal_alpha04_A_73_73",  0.40,  73,  73),
    ("thermal_alpha04_F_73_45",  0.40,  73,  45),
    ("thermal_alpha04_I_100_73", 0.40, 100,  73),
    ("thermal_alpha06_A_73_73",  0.60,  73,  73),
    ("thermal_alpha06_F_73_45",  0.60,  73,  45),
    ("thermal_alpha06_I_100_73", 0.60, 100,  73),
    ("thermal_alpha08_A_73_73",  0.80,  73,  73),
    ("thermal_alpha08_F_73_45",  0.80,  73,  45),
    ("thermal_alpha08_I_100_73", 0.80, 100,  73),
]


# ── §3.1 copy ──────────────────────────────────────────────────────────────

def copy_datasets():
    CW_DATA.mkdir(parents=True, exist_ok=True)
    print("\n§3.1 Copying datasets from ~/Downloads/thermal_alpha_AFI/")
    print("-" * 60)
    for folder, *_ in NEW_FOLDERS:
        src = DOWNLOADS_SRC / folder
        dst = CW_DATA / folder
        if not src.exists():
            print(f"  ERROR: source not found: {src}")
            sys.exit(1)
        if dst.exists():
            print(f"  {folder:<42}  already present, skipping")
        else:
            shutil.copytree(src, dst)
            print(f"  {folder:<42}  copied")

    print("\nFile sizes (output.txt):")
    for folder, *_ in NEW_FOLDERS:
        p = CW_DATA / folder / "output.txt"
        sz = p.stat().st_size if p.exists() else 0
        flag = "" if 6_000_000 <= sz <= 9_000_000 else "  ← UNEXPECTED SIZE"
        print(f"  {folder:<42}  {sz/1e6:.2f} MB{flag}")


# ── §3.2 validate ──────────────────────────────────────────────────────────

def _extract_fields(dat_path):
    mix, geom, constr, _ = parse_cw_dat(str(dat_path))
    return {
        # geometry
        "width_ft":  geom.width_ft,
        "depth_ft":  geom.depth_ft,
        # mix
        "cement_lb_yd3":  mix.cement_type_I_II_lb_yd3,
        "water_lb_yd3":   mix.water_lb_yd3,
        "coarse_lb_yd3":  mix.coarse_agg_lb_yd3,
        "fine_lb_yd3":    mix.fine_agg_lb_yd3,
        "k_BTU":          mix.thermal_conductivity_BTU_hr_ft_F,
        # hydration kinetics
        "alpha_u": mix.alpha_u,
        "tau_hrs": mix.tau_hrs,
        "beta":    mix.beta,
        "Hu_J_kg": mix.Hu_J_kg,
        "Ea_J_mol": mix.activation_energy_J_mol,
        # construction
        "T_pl_F":   constr.placement_temp_F,
        "T_soil_F": constr.soil_temp_F,
        "form_type": constr.form_type,
    }


def validate_datasets():
    FINDINGS.mkdir(parents=True, exist_ok=True)

    s7_ref = _extract_fields(S7_RUNS / "runA_baseline" / "input.dat")

    print("\n§3.2 Input.dat consistency check")
    print("=" * 80)

    rows = []
    issues = []

    for folder, exp_au, exp_pl, exp_so in NEW_FOLDERS:
        dat = CW_DATA / folder / "input.dat"
        if not dat.exists():
            issues.append(f"MISSING: {dat}")
            continue
        f = _extract_fields(dat)

        errs = []
        # geometry must match Sprint 7
        if abs(f["width_ft"]  - s7_ref["width_ft"])  > 0.01: errs.append(f"width_ft {f['width_ft']} ≠ {s7_ref['width_ft']}")
        if abs(f["depth_ft"]  - s7_ref["depth_ft"])  > 0.01: errs.append(f"depth_ft {f['depth_ft']} ≠ {s7_ref['depth_ft']}")
        # mix design must match Sprint 7
        for key in ("cement_lb_yd3", "water_lb_yd3", "coarse_lb_yd3", "fine_lb_yd3", "k_BTU"):
            if abs(f[key] - s7_ref[key]) > 0.01:
                errs.append(f"{key} {f[key]} ≠ {s7_ref[key]}")
        # kinetics: τ=5, β=0.85, Hu=1.0 (α_u varies per dataset)
        if abs(f["tau_hrs"] - 5.0)  > 0.01: errs.append(f"tau_hrs {f['tau_hrs']} ≠ 5.0")
        if abs(f["beta"]    - 0.85) > 0.01: errs.append(f"beta {f['beta']} ≠ 0.85")
        if abs(f["Hu_J_kg"] - 1.0)  > 0.5:  errs.append(f"Hu_J_kg {f['Hu_J_kg']} ≠ 1.0")
        if abs(f["alpha_u"] - exp_au) > 0.01: errs.append(f"alpha_u {f['alpha_u']} ≠ expected {exp_au}")
        # T_pl and T_soil
        if abs(f["T_pl_F"]  - exp_pl) > 1.0: errs.append(f"T_pl_F {f['T_pl_F']} ≠ expected {exp_pl}")
        if abs(f["T_soil_F"]- exp_so) > 1.0: errs.append(f"T_soil_F {f['T_soil_F']} ≠ expected {exp_so}")
        # form type
        if f["form_type"] != s7_ref["form_type"]: errs.append(f"form_type '{f['form_type']}' ≠ '{s7_ref['form_type']}'")

        status = "OK" if not errs else "FAIL"
        if errs:
            issues.append(f"{folder}: {'; '.join(errs)}")

        print(f"  {folder:<42}  {status}  α_u={f['alpha_u']:.2f}  T_pl={f['T_pl_F']:.0f}°F  T_soil={f['T_soil_F']:.0f}°F")
        rows.append((folder, exp_au, exp_pl, exp_so, f, status, errs))

    if issues:
        print("\n*** INCONSISTENCIES FOUND — stopping ***")
        for i in issues:
            print(f"  {i}")
    else:
        print("\n✓ All 12 datasets pass consistency check.")

    _write_findings(rows, s7_ref)

    if issues:
        sys.exit(1)


def _write_findings(rows, s7_ref):
    md = ["# Sprint 8 Stage 2 §3.2 — Input.dat Consistency\n"]
    md.append("## Sprint 7 Reference (runA_baseline)\n")
    for k, v in s7_ref.items():
        md.append(f"- {k}: {v}")
    md.append("\n## New Dataset Verification\n")
    md.append("| Folder | α_u | T_pl (°F) | T_soil (°F) | τ | β | Hu | width_ft | depth_ft | k_BTU | Status |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for folder, exp_au, exp_pl, exp_so, f, status, errs in rows:
        md.append(
            f"| {folder} | {f['alpha_u']:.2f} | {f['T_pl_F']:.0f} | {f['T_soil_F']:.0f}"
            f" | {f['tau_hrs']:.1f} | {f['beta']:.2f} | {f['Hu_J_kg']:.1f}"
            f" | {f['width_ft']:.0f} | {f['depth_ft']:.0f} | {f['k_BTU']:.4f} | {status} |"
        )
        if errs:
            for e in errs:
                md.append(f"  - ERROR: {e}")
    md.append(f"\n**Verdict:** {'All 12 consistent ✓' if all(r[5]=='OK' for r in rows) else 'FAILURES FOUND — see above'}")
    p = FINDINGS / "input_consistency.md"
    p.write_text("\n".join(md))
    print(f"\nWrote: {p}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-copy", action="store_true")
    args = ap.parse_args()
    if not args.no_copy:
        copy_datasets()
    validate_datasets()


if __name__ == "__main__":
    main()
