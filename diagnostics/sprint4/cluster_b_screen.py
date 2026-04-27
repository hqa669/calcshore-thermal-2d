"""Sprint 4 recon spike round 2 — hydration screen for Cluster B.

Stage 1: Compute early-window [12,36]hr hydration-fit metrics for the
Reference set (MIX-01,02,03,11,12) and Cluster B candidates (04,08,14,15).
Threshold per-metric = 1.5 × max(Reference).

Stage 2: For Cluster B candidates that pass Stage 1, diff
CWConstruction/CWGeometry/CWEnvironment fields against Reference.

Throwaway. Do not commit.
"""
import dataclasses
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CW_EXPORTS = REPO / "validation" / "cw_exports"
OUT_DIR = Path(__file__).parent
OUT_SCREEN_MD = OUT_DIR / "hydration_screen.md"
OUT_DIFF_MD   = OUT_DIR / "cluster_b_diff.md"

sys.path.insert(0, str(REPO))
from cw_scenario_loader import load_cw_scenario  # noqa: E402
from thermal_engine_2d import build_grid_half_mat, solve_hydration_2d  # noqa: E402

REFERENCE = ["MIX-01", "MIX-02", "MIX-03", "MIX-11", "MIX-12"]
CLUSTER_B = ["MIX-04", "MIX-08", "MIX-14", "MIX-15"]

EARLY_LO_HR = 12.0
EARLY_HI_HR = 36.0
THRESHOLD_MULT = 1.5


def run_engine(mix_dir: Path):
    """Replicates compare_to_cw.py:80-130 in-process. Returns scn + arrays."""
    scn = load_cw_scenario(
        str(mix_dir / "input.dat"),
        str(mix_dir / "weather.dat"),
        str(mix_dir / "output.txt"),
    )
    grid = build_grid_half_mat(scn.geometry.width_ft, scn.geometry.depth_ft)
    T0_C = (scn.construction.placement_temp_F - 32.0) * 5.0 / 9.0
    T_initial = np.full((grid.ny, grid.nx), T0_C)
    result = solve_hydration_2d(
        grid, scn.mix, T_initial,
        duration_s=168 * 3600,
        output_interval_s=1800.0,
        boundary_mode="full_2d",
        environment=scn.environment,
        construction=scn.construction,
        diagnostic_outputs=True,
    )
    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    eng_center_C = result.T_field_C[:, iy_mid, grid.nx - 1]
    eng_center_F = eng_center_C * 9.0 / 5.0 + 32.0
    return scn, result, eng_center_F


def time_to_half_peak(t_hrs, T):
    """First time the curve crosses T0 + 0.5*(T_max - T0). NaN if never."""
    T = np.asarray(T)
    T0 = float(T[0])
    Tmax = float(T.max())
    if Tmax <= T0:
        return float("nan")
    target = T0 + 0.5 * (Tmax - T0)
    above = np.where(T >= target)[0]
    if above.size == 0:
        return float("nan")
    return float(t_hrs[int(above[0])])


def peak_rise_rate(t_hrs, T, lo_hr, hi_hr):
    """Max dT/dt on [lo,hi] (°F per hr)."""
    t = np.asarray(t_hrs, dtype=float)
    T = np.asarray(T, dtype=float)
    if len(t) < 2:
        return float("nan")
    dT_dt = np.diff(T) / np.diff(t)
    t_mid = 0.5 * (t[1:] + t[:-1])
    mask = (t_mid >= lo_hr) & (t_mid <= hi_hr)
    if not mask.any():
        return float("nan")
    return float(dT_dt[mask].max())


def compute_metrics(mix_id):
    """Run engine, compute 3 hydration-fit metrics. Returns dict or {'error': ...}."""
    mix_dir = CW_EXPORTS / mix_id
    try:
        scn, result, eng_center_F = run_engine(mix_dir)
    except Exception as e:
        return {"mix": mix_id, "error": f"engine_crash: {type(e).__name__}: {e}"}

    if not np.all(np.isfinite(eng_center_F)):
        return {"mix": mix_id, "error": "engine_crash: NaN in centerline"}

    val = scn.cw_validation
    if val is None or val.T_field_F is None:
        return {"mix": mix_id, "error": "no CW T_field_F fixture"}

    cw_t_hrs = np.asarray(val.time_hrs, dtype=float)
    cw_t_s   = cw_t_hrs * 3600.0
    n_cw_t, n_cw_d, n_cw_w = val.T_field_F.shape
    cw_center_F = val.T_field_F[:, n_cw_d // 2, 0]

    # Engine interpolated to CW times
    eng_center_at_cw = np.interp(cw_t_s, result.t_s, eng_center_F)

    # Metric 1: Centerline RMS [12, 36]hr
    mask = (cw_t_hrs >= EARLY_LO_HR) & (cw_t_hrs <= EARLY_HI_HR)
    if not mask.any():
        return {"mix": mix_id, "error": "no CW samples in [12,36]hr"}
    diff = eng_center_at_cw[mask] - cw_center_F[mask]
    rms_early = float(np.sqrt(np.mean(diff ** 2)))

    # Metric 2: Centerline rise rate — peak dT/dt on [12,36], engine vs CW; report |Δ|
    eng_rise = peak_rise_rate(cw_t_hrs, eng_center_at_cw, EARLY_LO_HR, EARLY_HI_HR)
    cw_rise  = peak_rise_rate(cw_t_hrs, cw_center_F,      EARLY_LO_HR, EARLY_HI_HR)
    rise_delta = abs(eng_rise - cw_rise)

    # Metric 3: Time-to-half-peak (full series), engine vs CW; report |Δ|
    eng_t_hrs_full = result.t_s / 3600.0
    eng_thp = time_to_half_peak(eng_t_hrs_full, eng_center_F)
    cw_thp  = time_to_half_peak(cw_t_hrs, cw_center_F)
    thp_delta = abs(eng_thp - cw_thp) if (np.isfinite(eng_thp) and np.isfinite(cw_thp)) else float("nan")

    return {
        "mix":         mix_id,
        "rms_early":   rms_early,
        "eng_rise":    eng_rise,
        "cw_rise":     cw_rise,
        "rise_delta":  rise_delta,
        "eng_thp":     eng_thp,
        "cw_thp":      cw_thp,
        "thp_delta":   thp_delta,
        "scn":         scn,
    }


# ---------------------------------------------------------------------------
# Run all 9 mixes
# ---------------------------------------------------------------------------

print(f"Running engine in-process across {len(REFERENCE) + len(CLUSTER_B)} mixes...")
all_results = {}
for mix_id in REFERENCE + CLUSTER_B:
    print(f"  [{mix_id}] ", end="", flush=True)
    r = compute_metrics(mix_id)
    all_results[mix_id] = r
    if "error" in r:
        print(f"ERROR: {r['error']}", flush=True)
    else:
        print(f"rms={r['rms_early']:.2f}°F  Δrise={r['rise_delta']:.2f}°F/hr  Δt½={r['thp_delta']:.2f}hr",
              flush=True)

# ---------------------------------------------------------------------------
# Derive thresholds from Reference set
# ---------------------------------------------------------------------------

ref_ok = [all_results[m] for m in REFERENCE if "error" not in all_results[m]]
if not ref_ok:
    print("FATAL: no Reference mixes produced metrics — cannot derive thresholds.")
    sys.exit(2)

thr_rms   = max(r["rms_early"]  for r in ref_ok) * THRESHOLD_MULT
thr_rise  = max(r["rise_delta"] for r in ref_ok) * THRESHOLD_MULT
thr_thp   = max(r["thp_delta"]  for r in ref_ok) * THRESHOLD_MULT

print()
print(f"Thresholds (1.5x Reference max):")
print(f"  Centerline RMS [12,36]hr : {thr_rms:.2f}°F")
print(f"  Δ peak rise rate         : {thr_rise:.2f}°F/hr")
print(f"  Δ time-to-half-peak      : {thr_thp:.2f}hr")

# ---------------------------------------------------------------------------
# Classify Cluster B candidates
# ---------------------------------------------------------------------------

def classify(r):
    if "error" in r:
        return "engine_crash", []
    fails = []
    if r["rms_early"]  > thr_rms:  fails.append("RMS")
    if r["rise_delta"] > thr_rise: fails.append("rise")
    if r["thp_delta"]  > thr_thp:  fails.append("t_half")
    return ("hydration_pass" if not fails else "hydration_fail"), fails


classifications = {m: classify(all_results[m]) for m in CLUSTER_B}
survivors = [m for m, (cls, _) in classifications.items() if cls == "hydration_pass"]

print()
print(f"Cluster B classification:")
for m in CLUSTER_B:
    cls, fails = classifications[m]
    note = f" (failed: {','.join(fails)})" if fails else ""
    print(f"  {m}: {cls}{note}")

# ---------------------------------------------------------------------------
# Write hydration_screen.md
# ---------------------------------------------------------------------------

def fmt(x, prec=2):
    if isinstance(x, float) and not np.isfinite(x):
        return "NaN"
    return f"{x:.{prec}f}"


md = ["# Sprint 4 — Hydration-Fit Screen (Cluster B)", ""]
md.append(f"Early window: [{EARLY_LO_HR:.0f}, {EARLY_HI_HR:.0f}]hr. "
          f"Threshold = {THRESHOLD_MULT}× max(Reference) per metric.")
md.append("")
md.append("## Metrics")
md.append("")
md.append("| Mix | Set | Centerline RMS [12,36]hr (°F) | Engine rise (°F/hr) | CW rise (°F/hr) | Δ rise (°F/hr) | Engine t½ (hr) | CW t½ (hr) | Δ t½ (hr) | Result |")
md.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")

def row(mix_id, set_label):
    r = all_results[mix_id]
    if "error" in r:
        return f"| {mix_id} | {set_label} | — engine_crash — |  |  |  |  |  |  | {r['error']} |"
    if set_label == "Reference":
        result_str = "ref"
    else:
        cls, fails = classifications[mix_id]
        result_str = cls + (f" ({','.join(fails)})" if fails else "")
    return ("| {mix} | {set} | {rms} | {er} | {cr} | {dr} | {et} | {ct} | {dt} | {res} |"
            .format(mix=mix_id, set=set_label,
                    rms=fmt(r["rms_early"]), er=fmt(r["eng_rise"]), cr=fmt(r["cw_rise"]),
                    dr=fmt(r["rise_delta"]),
                    et=fmt(r["eng_thp"]), ct=fmt(r["cw_thp"]), dt=fmt(r["thp_delta"]),
                    res=result_str))

for m in REFERENCE:
    md.append(row(m, "Reference"))
for m in CLUSTER_B:
    md.append(row(m, "Cluster B"))

md.append("")
md.append("## Thresholds derived from Reference")
md.append("")
md.append(f"- Centerline RMS [12,36]hr : **{thr_rms:.2f} °F**")
md.append(f"- Δ peak rise rate         : **{thr_rise:.2f} °F/hr**")
md.append(f"- Δ time-to-half-peak      : **{thr_thp:.2f} hr**")
md.append("")
md.append("## Cluster B classification")
md.append("")
md.append("| Mix | Result | Failed metrics |")
md.append("| --- | --- | --- |")
for m in CLUSTER_B:
    cls, fails = classifications[m]
    md.append(f"| {m} | {cls} | {', '.join(fails) if fails else '—'} |")
md.append("")
if survivors:
    md.append(f"**Stage 1 survivors → Stage 2 construction-parameter diff:** "
              f"{', '.join(survivors)}")
else:
    md.append("**Stage 1: zero Cluster B mixes pass the hydration screen.** "
              "Sprint 4's evaluation set = Reference set only. "
              "All Cluster B mixes are routed to Sprint 5 (hydration scope). "
              "Stage 2 (construction-parameter diff) skipped — nothing to diff.")

OUT_SCREEN_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
print()
print(f"Wrote: {OUT_SCREEN_MD}")

# ---------------------------------------------------------------------------
# Stage 2 — only if survivors exist
# ---------------------------------------------------------------------------

if not survivors:
    print()
    print("=" * 60)
    print("STAGE 2 SKIPPED — no Cluster B mixes survived hydration screen.")
    print("=" * 60)
    print("Finding: Sprint 4's valid evaluation set = Reference set "
          f"({', '.join(REFERENCE)}).")
    print(f"Cluster B (all 4) joins Cluster A as hydration-scope (Sprint 5).")
    sys.exit(0)


# Stage 2: tabulate scenario fields
DATACLASSES = ["construction", "geometry", "environment"]


def field_value_repr(v):
    """Render a field value, skipping arrays >12 elements."""
    if isinstance(v, np.ndarray):
        if v.size > 12:
            return f"<ndarray shape={v.shape}>"
        return "[" + ",".join(f"{x:g}" for x in v.tolist()) + "]"
    if isinstance(v, list):
        if len(v) > 12:
            return f"<list len={len(v)}>"
        return "[" + ",".join(repr(x) for x in v) + "]"
    if isinstance(v, float):
        return f"{v:g}"
    return repr(v) if not isinstance(v, str) else v


def collect_fields(scn):
    out = {}
    for grp in DATACLASSES:
        obj = getattr(scn, grp)
        for f in dataclasses.fields(obj):
            v = getattr(obj, f.name)
            if isinstance(v, np.ndarray) and v.size > 12:
                continue
            out[f"{grp}.{f.name}"] = field_value_repr(v)
    return out


ref_scns = {m: all_results[m]["scn"] for m in REFERENCE}
sur_scns = {m: all_results[m]["scn"] for m in survivors}

field_table = {}  # key -> {mix: value_repr}
all_keys = set()
for m, scn in {**ref_scns, **sur_scns}.items():
    fv = collect_fields(scn)
    field_table[m] = fv
    all_keys.update(fv.keys())

ordered_keys = sorted(all_keys)


def varies_label(key):
    ref_vals = {field_table[m].get(key, "—") for m in REFERENCE}
    sur_vals = {field_table[m].get(key, "—") for m in survivors}
    ref_internal = len(ref_vals) > 1
    sur_internal = len(sur_vals) > 1
    cross = ref_vals != sur_vals if (not ref_internal and not sur_internal) else (ref_vals & sur_vals != ref_vals or ref_vals & sur_vals != sur_vals)
    if not ref_internal and not sur_internal and ref_vals == sur_vals:
        return "uniform"
    label_bits = []
    if not ref_internal and not sur_internal and ref_vals != sur_vals:
        label_bits.append("Reference≠Survivors")
    if ref_internal:
        label_bits.append("Reference-Internal")
    if sur_internal:
        label_bits.append("Survivors-Internal")
    return ", ".join(label_bits) if label_bits else "uniform"


# Markdown
md2 = ["# Sprint 4 — Cluster B Construction-Parameter Diff", ""]
md2.append(f"Reference: {', '.join(REFERENCE)}  |  Survivors: {', '.join(survivors)}")
md2.append("")
hdr = ["Field"] + REFERENCE + survivors + ["varies?"]
md2.append("| " + " | ".join(hdr) + " |")
md2.append("| " + " | ".join(["---"] * len(hdr)) + " |")
for key in ordered_keys:
    cells = [key]
    for m in REFERENCE:
        cells.append(str(field_table[m].get(key, "—")))
    for m in survivors:
        cells.append(str(field_table[m].get(key, "—")))
    cells.append(varies_label(key))
    md2.append("| " + " | ".join(cells) + " |")

OUT_DIFF_MD.write_text("\n".join(md2) + "\n", encoding="utf-8")
print(f"Wrote: {OUT_DIFF_MD}")

# Distinguishing fields summary
print()
print("=" * 60)
print("STAGE 2 SUMMARY")
print("=" * 60)
distinguishing = []
for key in ordered_keys:
    ref_vals = {field_table[m].get(key, "—") for m in REFERENCE}
    sur_vals = {field_table[m].get(key, "—") for m in survivors}
    if len(ref_vals) == 1 and (sur_vals - ref_vals):
        distinguishing.append((key, list(ref_vals)[0], sur_vals))

if not distinguishing:
    print("No fields where Reference is uniform AND Survivors differ from it.")
else:
    print(f"{len(distinguishing)} field(s) where Reference is uniform and Survivors deviate:")
    for key, refval, survals in distinguishing[:10]:
        print(f"  {key}")
        print(f"    Reference (all 5): {refval}")
        print(f"    Survivors:         {survals}")
