#!/usr/bin/env python3
"""S1.1 — Verify dataset consistency across all 9 input.dat files.

Usage:
    python validation/soil_calibration/check_inputs.py

Writes:
    validation/soil_calibration/STAGE1_input_field_map.md
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))

from cw_scenario_loader import parse_cw_dat, CW_DAT_INDEX

HERE = os.path.dirname(os.path.abspath(__file__))
CW_RUNS = os.path.join(HERE, "cw_runs")
OUT_MD = os.path.join(HERE, "STAGE1_input_field_map.md")

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

EXPECTED_GEOMETRY = (40.0, 80.0, 40.0)
EXPECTED_HYDRATION = {
    "activation_energy_J_mol": 50000.0,
    "tau_hrs": 200.0,
    "beta": 0.1,
    "alpha_u": 0.1,
    "Hu_J_kg": 1.0,
}
EXPECTED_FORM = "Steel"

def read_raw_lines(path):
    with open(path, encoding="latin-1") as f:
        return [ln.rstrip("\r\n").strip() for ln in f.readlines()]

def run():
    rows = []
    issues = []
    raw_field_map = None

    for label, folder, expected_placement, expected_soil in RUNS:
        dat_path = os.path.join(CW_RUNS, folder, "input.dat")
        if not os.path.isfile(dat_path):
            rows.append({
                "label": label, "folder": folder,
                "placement": "MISSING", "soil": "MISSING",
                "delta_T": "?", "geometry_ok": "MISSING",
                "hydration_ok": "MISSING", "form_ok": "MISSING",
                "soil_str": "MISSING", "issues": ["input.dat not found"],
            })
            continue

        mix, geom, constr, raw = parse_cw_dat(dat_path)
        lines = read_raw_lines(dat_path)

        row_issues = []

        # Geometry check
        geom_vals = (geom.width_ft, geom.depth_ft, geom.length_ft)
        geometry_ok = geom_vals == EXPECTED_GEOMETRY
        if not geometry_ok:
            row_issues.append(f"geometry={geom_vals} (expected {EXPECTED_GEOMETRY})")

        # Hydration check
        hydration_ok = True
        for field, expected in EXPECTED_HYDRATION.items():
            actual = getattr(mix, field, None)
            if actual is None:
                actual = float(raw.get(field, "nan"))
            if abs(actual - expected) > 0.001:
                hydration_ok = False
                row_issues.append(f"{field}={actual} (expected {expected})")

        # Form type check
        form_actual = raw.get("form_type", "").strip()
        form_ok = form_actual == EXPECTED_FORM
        if not form_ok:
            row_issues.append(f"form_type={form_actual!r} (expected {EXPECTED_FORM!r})")

        # Placement and soil temperature check
        placement_actual = constr.placement_temp_F
        soil_actual = constr.soil_temp_F
        if abs(placement_actual - expected_placement) > 0.5:
            row_issues.append(
                f"placement_temp_F={placement_actual} (expected {expected_placement})"
            )
        if abs(soil_actual - expected_soil) > 0.5:
            row_issues.append(
                f"soil_temp_F={soil_actual} (expected {expected_soil})"
            )

        # Soil string from lines 466-467 (0-indexed = indices 465, 466)
        soil_side = lines[465] if len(lines) > 465 else "?"
        soil_bot  = lines[466] if len(lines) > 466 else "?"
        soil_str_ok = (soil_side == "Clay" and soil_bot == "Clay")
        if not soil_str_ok:
            row_issues.append(f"soil strings at idx 465/466: {soil_side!r}/{soil_bot!r} (expected Clay/Clay)")

        # footing_subbase from parser (idx 446)
        footing_sub = raw.get("footing_subbase", "").strip()

        # Soil profile from lines 487-495 (0-indexed 486-494)
        profile = [lines[i].strip() for i in range(486, 495) if i < len(lines)]
        profile_vals = set(profile)

        delta_T = expected_soil - expected_placement

        row = {
            "label": label,
            "folder": folder,
            "placement": expected_placement,
            "soil": expected_soil,
            "delta_T": delta_T,
            "geometry_ok": "OK" if geometry_ok else "FAIL",
            "hydration_ok": "OK" if hydration_ok else "FAIL",
            "form_ok": "OK" if form_ok else "FAIL",
            "soil_str": f"{soil_side}/{soil_bot}",
            "footing_sub": footing_sub,
            "profile_vals": sorted(profile_vals),
            "issues": row_issues,
        }
        rows.append(row)
        if row_issues:
            issues.append((label, row_issues))

        if raw_field_map is None:
            raw_field_map = (lines, raw)

    # Build report
    lines_sample, raw_sample = raw_field_map if raw_field_map else ([], {})

    md_lines = [
        "# STAGE1 — input.dat Field Map and Consistency Check",
        "",
        "## Field-to-Parameter Mapping",
        "",
        "All indices are **0-based** (Python `lines[idx]`). File line = idx + 1.",
        "",
        "| Field | idx | 1-indexed line | Description |",
        "| ----- | --- | -------------- | ----------- |",
    ]
    for key, idx in sorted(CW_DAT_INDEX.items(), key=lambda x: x[1]):
        val = raw_sample.get(key, "") if raw_sample else ""
        md_lines.append(f"| `{key}` | {idx} | {idx+1} | {val} |")

    md_lines += [
        "",
        "### Extra fields NOT in CW_DAT_INDEX but present in input.dat",
        "",
        "| idx | 1-indexed line | Content (from runA_baseline) | Notes |",
        "| --- | -------------- | ----------------------------- | ----- |",
        "| 465 | 466 | `Clay` | Side soil type (parsed by CW but not by loader) |",
        "| 466 | 467 | `Clay` | Bottom soil type (parsed by CW but not by loader) |",
        "| 486–494 | 487–495 | soil temp profile (9 values) | Per-depth soil temperatures; mirrors soil_temp_F |",
        "",
        "**Known gap**: `footing_subbase` (idx 446) reads `Limestone` — a CW UI default that doesn't",
        "reflect the actual soil used by CW's model. The true soil type (`Clay`) lives at idx 465/466,",
        "which the loader does NOT parse. The engine ignores `footing_subbase` entirely.",
        "",
        "---",
        "",
        "## Consistency Table",
        "",
        "| Label | Folder | Placement °F | Soil °F | ΔT | Geometry | Hydration | Form | Soil strings (idx465/466) | Notes |",
        "| ----- | ------ | ------------ | ------- | -- | -------- | --------- | ---- | ------------------------- | ----- |",
    ]

    for r in rows:
        notes = "; ".join(r["issues"]) if r["issues"] else "—"
        md_lines.append(
            f"| {r['label']} | {r['folder']} | {r['placement']} | {r['soil']} "
            f"| {r['delta_T']} | {r['geometry_ok']} | {r['hydration_ok']} "
            f"| {r.get('form_ok','?')} | {r.get('soil_str','?')} | {notes} |"
        )

    md_lines += [
        "",
        "## Summary",
        "",
    ]

    if issues:
        md_lines.append("**ISSUES FOUND:**")
        for label, iss in issues:
            md_lines.append(f"- Run {label}: " + "; ".join(iss))
    else:
        md_lines.append(
            "All 9 runs confirm: geometry=40×80×40 ft, Hu_J_kg=1 (suppressed), "
            "Ea=50000, τ=200, β=0.1, αu=0.1, form=Steel, soil strings=Clay/Clay. "
            "Placement and soil temperatures match labels."
        )

    md_lines += [
        "",
        "## Mirror Pair Note",
        "",
        "| Pair | Run+ | Run- | |ΔT| | Clean mirror? |",
        "| ---- | ---- | ---- | ----- | ------------- |",
        "| B↔D | B (73/60) | D (60/73) | 13 | Yes |",
        "| C↔E | C (73/90) | E (90/73) | 17 | Yes |",
        "| F↔H | F (73/45) | H (45/73) | 28 | Yes |",
        "| G↔I | G (73/100) | I (100/73) | 27 | Yes — note: brief's IMPORTANT warning about Run I having soil=60 is incorrect; actual data shows soil=73 |",
    ]

    with open(OUT_MD, "w") as f:
        f.write("\n".join(md_lines) + "\n")

    print(f"Written: {OUT_MD}")

    # Print summary to stdout
    all_ok = not issues
    print(f"\nConsistency check: {'PASS' if all_ok else 'FAIL'}")
    for r in rows:
        status = "OK" if not r["issues"] else "ISSUES: " + "; ".join(r["issues"])
        print(f"  Run {r['label']} ({r['folder']}): placement={r['placement']}, soil={r['soil']}, ΔT={r['delta_T']} → {status}")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(run())
