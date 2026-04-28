"""
batch_compare_all_mixes.py

Run engine vs CW adiabatic centerline comparison for all 15 mixes
(MIX-01 through MIX-15), using a composition-based Hu calibration
factor from the apr28 passdown:

    Hu_eff = Hu_regressed × f_type × Σᵢ kᵢ · pᵢ

Per-mix T₀ (placement temperature):
    MIX-01 .. MIX-13 → 73°F
    MIX-14           → 75°F
    MIX-15           → 45°F

Inputs (under --root, default ~/Downloads/HydrationCenter_mix01-15):
    HydrationCenter_mixNN/input.dat
    HydrationCenter_mixNN/cw_adiabatic_reference_mixNN.csv

Outputs (under --out-dir, default ./batch_results):
    engine_vs_cw_mixNN.png           (per-mix 2-panel plot)
    engine_vs_cw_mixNN.csv           (per-mix raw t/T_eng/T_cw/delta)
    summary_table.csv                (one row per mix)
    summary_table.md                 (markdown table for passdown/notes)
    batch_overview.png               (15-panel grid)

Usage:
    python batch_compare_all_mixes.py
    python batch_compare_all_mixes.py --root /path/to/HydrationCenter_mix01-15
    python batch_compare_all_mixes.py --out-dir validation/post_kinetics_correction
    python batch_compare_all_mixes.py --raw-hu      # diagnostic: turn off correction
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cw_scenario_loader import load_cw_scenario
from thermal_engine_2d import build_grid_rectangular, solve_hydration_2d


# ======================================================================
# Hu calibration (apr28 passdown)
# ======================================================================

F_TYPE = {
    "Type I":    0.9982,
    "Type I/II": 0.9897,
    "Type II":   0.9870,
    "Type V":    0.9956,
}

K_PORTLAND   = 1.0
K_FLY_ASH_F  = 0.91
K_FLY_ASH_C  = 0.883
K_SLAG       = 0.89
K_SILICA_FUME_PROVISIONAL = 0.7   # apr28: provisional until moderate-replacement runs


def compute_hu_factor(mix) -> tuple[float, str]:
    """
    Compute the apr28 calibration factor:  f_type × Σᵢ kᵢ · pᵢ

    Returns (factor, note).  `note` is a short flag string for the
    summary table (empty if everything is in-envelope).
    """
    cement_type = (mix.cement_type or "").strip()
    if cement_type not in F_TYPE:
        # Try a few common normalizations
        norm = cement_type.replace("Type", "Type ").replace("  ", " ").strip()
        if norm in F_TYPE:
            cement_type = norm
        else:
            raise KeyError(
                f"Cement type {cement_type!r} not in F_TYPE table "
                f"(known: {list(F_TYPE)})"
            )
    f_type = F_TYPE[cement_type]

    cem  = mix.cement_type_I_II_lb_yd3
    fa_f = mix.fly_ash_F_lb_yd3
    fa_c = mix.fly_ash_C_lb_yd3
    slag = mix.ggbfs_lb_yd3
    sf   = mix.silica_fume_lb_yd3
    total = cem + fa_f + fa_c + slag + sf
    if total <= 0:
        raise ValueError("total cementitious is zero; cannot compute fractions")

    p_cem  = cem  / total
    p_fa_f = fa_f / total
    p_fa_c = fa_c / total
    p_slag = slag / total
    p_sf   = sf   / total

    note_parts = []
    k_sf_used = 0.0
    if p_sf > 1e-6:
        k_sf_used = K_SILICA_FUME_PROVISIONAL
        if p_sf > 0.15:
            note_parts.append(f"SF={p_sf*100:.1f}% OUT-OF-ENVELOPE")
        else:
            note_parts.append(f"SF={p_sf*100:.1f}% k_SF=0.7 PROVISIONAL")

    k_scm_sum = (
        K_PORTLAND   * p_cem
        + K_FLY_ASH_F * p_fa_f
        + K_FLY_ASH_C * p_fa_c
        + K_SLAG      * p_slag
        + k_sf_used   * p_sf
    )
    factor = f_type * k_scm_sum
    return factor, "; ".join(note_parts)


# ======================================================================
# Per-mix table (T₀ assignment from apr28 passdown)
# ======================================================================

@dataclass
class MixSpec:
    mix_id: str          # "mix01"
    label: str           # "MIX-01"
    T0_F: float
    duration_hrs: float = 168.0


MIX_TABLE = [
    MixSpec("mix01", "MIX-01", 73.0),
    MixSpec("mix02", "MIX-02", 73.0),
    MixSpec("mix03", "MIX-03", 73.0),
    MixSpec("mix04", "MIX-04", 73.0),
    MixSpec("mix05", "MIX-05", 73.0),
    MixSpec("mix06", "MIX-06", 73.0),
    MixSpec("mix07", "MIX-07", 73.0),
    MixSpec("mix08", "MIX-08", 73.0),
    MixSpec("mix09", "MIX-09", 73.0),
    MixSpec("mix10", "MIX-10", 73.0),
    MixSpec("mix11", "MIX-11", 73.0),
    MixSpec("mix12", "MIX-12", 73.0),
    MixSpec("mix13", "MIX-13", 73.0),
    MixSpec("mix14", "MIX-14", 75.0),
    MixSpec("mix15", "MIX-15", 45.0),
]


# ======================================================================
# Helpers (lifted from plot_engine_vs_cw_centerline.py)
# ======================================================================

def f_to_c(t_f: float) -> float:
    return (t_f - 32.0) * 5.0 / 9.0


def c_to_f(t_c):
    return np.asarray(t_c) * 9.0 / 5.0 + 32.0


def load_cw_reference(csv_path: Path):
    arr = np.genfromtxt(str(csv_path), delimiter=",", skip_header=1)
    return arr[:, 0], arr[:, 1]


# ======================================================================
# One-mix engine run
# ======================================================================

def run_one_mix(input_dat: Path, T0_F: float, duration_hrs: float,
                use_calibrated_hu: bool, output_interval_s: float = 300.0):
    """
    Run engine in adiabatic mode for one mix.  Returns
    (t_hrs, T_F, summary_dict).
    """
    scenario = load_cw_scenario(str(input_dat),
                                weather_dat=None, cw_output_txt=None)
    mix = scenario.mix

    Hu_raw = mix.Hu_J_kg
    if use_calibrated_hu:
        factor, env_note = compute_hu_factor(mix)
    else:
        factor, env_note = 1.0, "raw-Hu mode"
    mix.Hu_J_kg = Hu_raw * factor

    # Tiny grid — adiabatic field stays uniform
    grid = build_grid_rectangular(Lx_m=1.0, Ly_m=1.0, nx=3, ny=3)
    T_init = np.full((grid.ny, grid.nx), f_to_c(T0_F), dtype=np.float64)

    res = solve_hydration_2d(
        grid, mix, T_init,
        duration_s=duration_hrs * 3600.0,
        output_interval_s=output_interval_s,
        boundary_mode="adiabatic",
    )

    t_hrs = res.t_s / 3600.0
    T_F = c_to_f(res.T_field_C[:, grid.ny // 2, grid.nx // 2])

    summary = {
        "cement_type": mix.cement_type,
        "cement_lb_yd3": mix.cement_type_I_II_lb_yd3,
        "fa_f_lb_yd3":   mix.fly_ash_F_lb_yd3,
        "fa_c_lb_yd3":   mix.fly_ash_C_lb_yd3,
        "slag_lb_yd3":   mix.ggbfs_lb_yd3,
        "sf_lb_yd3":     mix.silica_fume_lb_yd3,
        "total_cm_lb_yd3": mix.total_cementitious_lb_yd3,
        "Hu_raw_J_kg":   Hu_raw,
        "Hu_factor":     factor,
        "Hu_eff_J_kg":   mix.Hu_J_kg,
        "tau_hrs":       mix.tau_hrs,
        "beta":          mix.beta,
        "alpha_u":       mix.alpha_u,
        "Ea_J_mol":      mix.activation_energy_J_mol,
        "envelope_note": env_note,
    }
    return t_hrs, T_F, summary


# ======================================================================
# Per-mix plotting
# ======================================================================

def plot_one_mix(t_eng, T_eng, t_cw, T_cw, summary, mix_label, T0_F, out_png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    T_eng_on_cw = np.interp(t_cw, t_eng, T_eng)
    delta = T_eng_on_cw - T_cw
    peak_eng = float(T_eng[-1])
    peak_cw = float(T_cw[-1])
    peak_delta = peak_eng - peak_cw
    rms = float(np.sqrt(np.mean(delta ** 2)))
    max_abs = float(np.max(np.abs(delta)))

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(10, 7), sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax1.plot(t_cw, T_cw, "k-", linewidth=1.8,
             label=f"CW centerline @ T0={T0_F:.0f}°F  (peak {peak_cw:.2f}°F)")
    ax1.plot(t_eng, T_eng, "C0--", linewidth=1.5,
             label=f"Engine adiabatic, Hu×{summary['Hu_factor']:.4f}  "
                   f"(peak {peak_eng:.2f}°F)")
    ax1.set_ylabel("Temperature (°F)")
    ax1.legend(loc="lower right", fontsize=9)
    ax1.grid(True, alpha=0.3)
    title_l1 = (f"{mix_label} — T0 = {T0_F:.0f}°F adiabatic  |  "
                f"Cement: {summary['cement_type']}  "
                f"(Hu_factor = {summary['Hu_factor']:.4f})")
    title_l2 = (f"Cementitious: cem={summary['cement_lb_yd3']:.0f}, "
                f"FA-F={summary['fa_f_lb_yd3']:.0f}, "
                f"FA-C={summary['fa_c_lb_yd3']:.0f}, "
                f"slag={summary['slag_lb_yd3']:.0f}, "
                f"SF={summary['sf_lb_yd3']:.0f} lb/yd³")
    if summary["envelope_note"]:
        title_l2 += f"  [{summary['envelope_note']}]"
    ax1.set_title(title_l1 + "\n" + title_l2, fontsize=10)

    ax2.plot(t_cw, delta, "C3-", linewidth=1.0)
    ax2.axhline(0.0, color="k", linewidth=0.5)
    ax2.fill_between(t_cw, -1.0, 1.0, color="green", alpha=0.10,
                     label="±1°F gate")
    ax2.set_xlabel("Time since placement (hours)")
    ax2.set_ylabel("Engine − CW (°F)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper left", fontsize=9)

    plt.tight_layout()
    plt.savefig(str(out_png), dpi=140)
    plt.close(fig)

    return {
        "peak_engine_F": peak_eng,
        "peak_cw_F":     peak_cw,
        "peak_delta_F":  peak_delta,
        "rms_F":         rms,
        "max_abs_delta_F": max_abs,
    }


def write_one_csv(t_cw, T_cw, t_eng, T_eng, out_csv):
    T_eng_on_cw = np.interp(t_cw, t_eng, T_eng)
    delta = T_eng_on_cw - T_cw
    arr = np.column_stack([t_cw, T_eng_on_cw, T_cw, delta])
    np.savetxt(str(out_csv), arr, delimiter=",",
               header="time_hrs,T_engine_F,T_cw_F,delta_F",
               fmt=("%.4f", "%.4f", "%.4f", "%.5f"), comments="")


# ======================================================================
# 15-panel overview plot
# ======================================================================

def plot_overview(records, out_png):
    """5×3 grid of mini residual plots, one per mix."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(5, 3, figsize=(15, 18), sharex=True)
    axes = axes.flatten()
    for ax, rec in zip(axes, records):
        if rec.get("error"):
            ax.text(0.5, 0.5, f"{rec['mix_label']}\nERROR:\n{rec['error']}",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=9, color="red")
            ax.set_xticks([]); ax.set_yticks([])
            continue
        t_cw = rec["t_cw"]
        delta = rec["delta"]
        passed = abs(rec["peak_delta_F"]) <= 1.0 and rec["rms_F"] <= 1.0
        color = "C2" if passed else "C3"
        ax.plot(t_cw, delta, color=color, linewidth=1.0)
        ax.axhline(0, color="k", linewidth=0.4)
        ax.fill_between(t_cw, -1, 1, color="green", alpha=0.08)
        ax.set_title(
            f"{rec['mix_label']} @ {rec['T0_F']:.0f}°F  "
            f"f={rec['Hu_factor']:.4f}  "
            f"Δpk={rec['peak_delta_F']:+.2f}  rms={rec['rms_F']:.2f}",
            fontsize=9,
        )
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-5, 5)

    for ax in axes:
        ax.set_xlabel("Time (hr)")
        ax.set_ylabel("Engine − CW (°F)")

    fig.suptitle(
        "Engine (composition-calibrated Hu) vs CW adiabatic centerline — "
        "all 15 mixes",
        fontsize=13, y=0.995,
    )
    plt.tight_layout()
    plt.savefig(str(out_png), dpi=120)
    plt.close(fig)


# ======================================================================
# Summary table
# ======================================================================

SUMMARY_COLS = [
    "mix", "T0_F", "cement_type", "cem_pct", "fa_f_pct", "fa_c_pct",
    "slag_pct", "sf_pct", "Hu_raw_J_kg", "Hu_factor", "Hu_eff_J_kg",
    "peak_engine_F", "peak_cw_F", "peak_delta_F", "rms_F", "max_abs_F",
    "pass_pm1F", "envelope_note",
]


def write_summary_csv(records, out_csv):
    with open(out_csv, "w") as f:
        f.write(",".join(SUMMARY_COLS) + "\n")
        for r in records:
            if r.get("error"):
                f.write(f"{r['mix_label']},{r['T0_F']:.0f},,,,,,,,,,,,,,,ERROR: {r['error']}\n")
                continue
            row = [
                r["mix_label"],
                f"{r['T0_F']:.0f}",
                r["cement_type"],
                f"{r['cem_pct']:.1f}",
                f"{r['fa_f_pct']:.1f}",
                f"{r['fa_c_pct']:.1f}",
                f"{r['slag_pct']:.1f}",
                f"{r['sf_pct']:.1f}",
                f"{r['Hu_raw_J_kg']:.0f}",
                f"{r['Hu_factor']:.4f}",
                f"{r['Hu_eff_J_kg']:.0f}",
                f"{r['peak_engine_F']:.2f}",
                f"{r['peak_cw_F']:.2f}",
                f"{r['peak_delta_F']:+.2f}",
                f"{r['rms_F']:.2f}",
                f"{r['max_abs_delta_F']:.2f}",
                "Y" if r["pass_pm1F"] else "N",
                r["envelope_note"],
            ]
            f.write(",".join(row) + "\n")


def write_summary_md(records, out_md):
    lines = [
        "# Composition-calibrated Hu — engine vs CW adiabatic centerline\n",
        "Pass criterion: |peak Δ| ≤ 1.0°F and RMS ≤ 1.0°F.\n",
        "| Mix | T₀ | Cement | %Cem | %FA-F | %FA-C | %Slag | %SF | Hu_factor | Engine peak | CW peak | Δpeak | RMS | Pass | Note |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in records:
        if r.get("error"):
            lines.append(
                f"| {r['mix_label']} | {r['T0_F']:.0f}°F | — | | | | | | | | | "
                f"| | ❌ | ERROR: {r['error']} |"
            )
            continue
        passed = "✅" if r["pass_pm1F"] else "❌"
        lines.append(
            f"| {r['mix_label']} "
            f"| {r['T0_F']:.0f}°F "
            f"| {r['cement_type']} "
            f"| {r['cem_pct']:.1f} "
            f"| {r['fa_f_pct']:.1f} "
            f"| {r['fa_c_pct']:.1f} "
            f"| {r['slag_pct']:.1f} "
            f"| {r['sf_pct']:.1f} "
            f"| {r['Hu_factor']:.4f} "
            f"| {r['peak_engine_F']:.2f}°F "
            f"| {r['peak_cw_F']:.2f}°F "
            f"| {r['peak_delta_F']:+.2f}°F "
            f"| {r['rms_F']:.2f}°F "
            f"| {passed} "
            f"| {r['envelope_note']} |"
        )

    n_pass = sum(1 for r in records if not r.get("error") and r.get("pass_pm1F"))
    n_total = sum(1 for r in records if not r.get("error"))
    lines.append("")
    lines.append(f"**Overall: {n_pass}/{n_total} mixes pass ±1°F gate.**")
    Path(out_md).write_text("\n".join(lines) + "\n")


# ======================================================================
# Main batch driver
# ======================================================================

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root",
                   default=str(Path.home() / "Downloads/HydrationCenter_mix01-15"),
                   help="Root folder containing HydrationCenter_mixNN/ subfolders.")
    p.add_argument("--out-dir", default="batch_results",
                   help="Output directory.")
    p.add_argument("--raw-hu", action="store_true",
                   help="Diagnostic: disable Hu calibration "
                        "(use CW-regressed Hu as-is). Default applies the "
                        "f_type × Σ kᵢ·pᵢ correction.")
    p.add_argument("--duration-hrs", type=float, default=168.0,
                   help="Simulation duration (default 168).")
    p.add_argument("--mixes", default="all",
                   help="Comma-separated list of mix IDs (e.g. 'mix01,mix07') "
                        "or 'all' (default).")
    args = p.parse_args(argv)

    root = Path(args.root).expanduser()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.mixes == "all":
        specs = MIX_TABLE
    else:
        wanted = {s.strip() for s in args.mixes.split(",")}
        specs = [m for m in MIX_TABLE if m.mix_id in wanted]
        if not specs:
            print(f"No mixes matched {args.mixes!r}", file=sys.stderr)
            return 2

    use_calibrated = not args.raw_hu
    mode = "CALIBRATED Hu" if use_calibrated else "RAW Hu (diagnostic)"
    print(f"Batch mode: {mode}")
    print(f"Root:    {root}")
    print(f"Out:     {out_dir.resolve()}")
    print(f"Mixes:   {len(specs)}")
    print()

    records = []
    for spec in specs:
        mix_dir = root / f"HydrationCenter_{spec.mix_id}"
        input_dat = mix_dir / "input.dat"
        ref_csv = mix_dir / f"cw_adiabatic_reference_{spec.mix_id}.csv"

        rec = {
            "mix_label": spec.label,
            "mix_id":    spec.mix_id,
            "T0_F":      spec.T0_F,
        }

        if not input_dat.exists():
            rec["error"] = f"input.dat missing at {input_dat}"
            print(f"  ✗ {spec.label}: {rec['error']}")
            records.append(rec); continue
        if not ref_csv.exists():
            rec["error"] = f"reference CSV missing at {ref_csv}"
            print(f"  ✗ {spec.label}: {rec['error']}")
            records.append(rec); continue

        try:
            t_eng, T_eng, summary = run_one_mix(
                input_dat, spec.T0_F, args.duration_hrs, use_calibrated,
            )
            t_cw, T_cw = load_cw_reference(ref_csv)

            png_path = out_dir / f"engine_vs_cw_{spec.mix_id}.png"
            csv_path = out_dir / f"engine_vs_cw_{spec.mix_id}.csv"

            metrics = plot_one_mix(t_eng, T_eng, t_cw, T_cw,
                                   summary, spec.label, spec.T0_F, png_path)
            write_one_csv(t_cw, T_cw, t_eng, T_eng, csv_path)

            T_eng_on_cw = np.interp(t_cw, t_eng, T_eng)
            delta = T_eng_on_cw - T_cw
            total_cm = summary["total_cm_lb_yd3"]
            rec.update({
                "cement_type":   summary["cement_type"],
                "cem_pct":       100 * summary["cement_lb_yd3"]   / total_cm,
                "fa_f_pct":      100 * summary["fa_f_lb_yd3"]     / total_cm,
                "fa_c_pct":      100 * summary["fa_c_lb_yd3"]     / total_cm,
                "slag_pct":      100 * summary["slag_lb_yd3"]     / total_cm,
                "sf_pct":        100 * summary["sf_lb_yd3"]       / total_cm,
                "Hu_raw_J_kg":   summary["Hu_raw_J_kg"],
                "Hu_factor":     summary["Hu_factor"],
                "Hu_eff_J_kg":   summary["Hu_eff_J_kg"],
                "envelope_note": summary["envelope_note"],
                "t_cw":          t_cw,
                "delta":         delta,
                **metrics,
                "pass_pm1F": (abs(metrics["peak_delta_F"]) <= 1.0
                              and metrics["rms_F"] <= 1.0),
            })
            tag = "✓" if rec["pass_pm1F"] else "·"
            print(f"  {tag} {spec.label} @ {spec.T0_F:.0f}°F  "
                  f"f={summary['Hu_factor']:.4f}  "
                  f"Δpeak={metrics['peak_delta_F']:+.2f}°F  "
                  f"rms={metrics['rms_F']:.2f}°F"
                  f"{'  ['+summary['envelope_note']+']' if summary['envelope_note'] else ''}")
        except Exception as exc:  # noqa: BLE001
            rec["error"] = f"{type(exc).__name__}: {exc}"
            print(f"  ✗ {spec.label}: {rec['error']}")

        records.append(rec)

    # Cross-mix outputs
    print()
    print("Writing summary table…")
    write_summary_csv(records, out_dir / "summary_table.csv")
    write_summary_md(records,  out_dir / "summary_table.md")

    print("Writing 15-panel overview…")
    plot_overview(records, out_dir / "batch_overview.png")

    n_pass = sum(1 for r in records if not r.get("error") and r.get("pass_pm1F"))
    n_ok   = sum(1 for r in records if not r.get("error"))
    print()
    print("=" * 64)
    print(f"  Pass ±1°F gate (peak & RMS): {n_pass} / {n_ok}")
    print(f"  Outputs:  {out_dir.resolve()}")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
