"""Sprint 4 PR 16 — Phase 2.6: equivalent-age and hydration trajectory comparison.

Compares engine vs CW for MIX-01 (warm IC) and MIX-15 (cold IC):
  (a) Centerline temperature trajectory
  (b) Equivalent age te(t) — engine direct + CW-implied via same Arrhenius model
  (c) Hydration degree α(t)
  (d) Heat-generation proxy dα/dt

CW does not report α or te directly. CW-implied te is derived by integrating the
same Arrhenius factor af(T_CW) over the CW centerline temperature trajectory:
  te_cw(t) = te_0 + ∫₀ᵗ af(T_cw(τ)) dτ
  α_cw(t) = alpha_vec(te_cw(t), tau, beta, alpha_u)

If engine and CW use the same Arrhenius formulation and the same parameters, divergence
in the te or α trajectories would reveal WHERE the cold-IC discrepancy originates.

Comparison findings:
  1. Both mixes track → cold-IC failure downstream of equivalent age
  2. MIX-15 diverges, MIX-01 tracks → cold-IC-specific implementation difference
  3. Both diverge → baseline implementation difference, Sprint 5 priority

Output: diagnostics/sprint4/pr16_te_alpha_comparison.md

Usage (from repo root):
    python diagnostics/sprint4/pr16_te_alpha_comparison.py
"""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
CW_EXPORTS = REPO / "validation" / "cw_exports"
OUT_MD = Path(__file__).parent / "pr16_te_alpha_comparison.md"

sys.path.insert(0, str(REPO))

from cw_scenario_loader import load_cw_scenario  # noqa: E402
from thermal_engine_2d import (  # noqa: E402
    build_grid_half_mat, solve_hydration_2d,
    arrhenius_vec, hydration_alpha_vec, hydration_rate_vec,
    R_GAS, T_REF_K,
)

# Report time points (hrs) — focused on hydration window
REPORT_HRS = [0, 4, 8, 12, 16, 20, 24, 30, 36, 42, 48, 60, 72, 96, 120, 144, 168]


# ---------------------------------------------------------------------------
# Engine run
# ---------------------------------------------------------------------------

def run_engine(scn):
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
    return result, grid


# ---------------------------------------------------------------------------
# CW-implied te / α — integrate Arrhenius over CW temperature trajectory
# ---------------------------------------------------------------------------

def cw_implied_te_alpha(cw_time_hrs, cw_T_center_F, Ea, tau, beta, alpha_u, te_init=0.01):
    """Integrate the engine's Arrhenius model over the CW temperature trajectory.

    Returns (te_implied_hrs, alpha_implied) arrays at the CW timesteps.
    """
    cw_T_K = (cw_T_center_F - 32.0) * 5.0 / 9.0 + 273.15

    te_arr = np.empty(len(cw_time_hrs))
    te_arr[0] = te_init

    for i in range(1, len(cw_time_hrs)):
        dt_hrs = float(cw_time_hrs[i] - cw_time_hrs[i - 1])
        T_mid_K = float((cw_T_K[i - 1] + cw_T_K[i]) / 2.0)
        af = float(arrhenius_vec(np.array([T_mid_K]), Ea)[0])
        te_arr[i] = te_arr[i - 1] + dt_hrs * af

    alpha_arr = hydration_alpha_vec(np.maximum(te_arr, 0.01), tau, beta, alpha_u)
    return te_arr, alpha_arr


# ---------------------------------------------------------------------------
# Engine centerline extraction at report time points
# ---------------------------------------------------------------------------

def engine_centerline_at_times(result, grid, report_hrs):
    """Return engine centerline (mid-depth, x=max) T/te/α at requested hours."""
    t_hrs = result.t_s / 3600.0
    iy_mid = (grid.iy_concrete_start + grid.iy_concrete_end) // 2
    ix_cl  = grid.nx - 1   # centerline = rightmost column

    eng_T_F     = result.T_field_C[:, iy_mid, ix_cl] * 9.0 / 5.0 + 32.0
    eng_te      = result.t_e_field_hrs[:, iy_mid, ix_cl]
    eng_alpha   = result.alpha_field[:, iy_mid, ix_cl]

    out = {}
    for tgt in report_hrs:
        idx = int(np.argmin(np.abs(t_hrs - tgt)))
        out[tgt] = {
            "t_actual_hrs": float(t_hrs[idx]),
            "T_F":    float(eng_T_F[idx]),
            "te_hrs": float(eng_te[idx]),
            "alpha":  float(eng_alpha[idx]),
        }
    return out


# ---------------------------------------------------------------------------
# CW centerline extraction at report time points
# ---------------------------------------------------------------------------

def cw_centerline_at_times(scn, report_hrs, Ea, tau, beta, alpha_u):
    """Return CW centerline T, implied te, implied α at requested hours."""
    val = scn.cw_validation
    cw_t = val.time_hrs                   # shape (n_cw_t,)
    # Centerline = first width index in CW field (width stored decreasing, index 0 = max width = centerline)
    # Mid-depth = nD // 2
    n_cw_d = val.T_field_F.shape[1]
    cw_T_F = val.T_field_F[:, n_cw_d // 2, 0]   # centerline, mid-depth, °F

    te_implied, alpha_implied = cw_implied_te_alpha(cw_t, cw_T_F, Ea, tau, beta, alpha_u)

    out = {}
    for tgt in report_hrs:
        if tgt == 0:
            out[tgt] = {"T_F": float(cw_T_F[0]), "te_hrs": float(te_implied[0]),
                        "alpha": float(alpha_implied[0])}
            continue
        idx = int(np.argmin(np.abs(cw_t - tgt)))
        out[tgt] = {
            "T_F":    float(cw_T_F[idx]),
            "te_hrs": float(te_implied[idx]),
            "alpha":  float(alpha_implied[idx]),
        }
    return out


# ---------------------------------------------------------------------------
# dα/dt estimate (central difference)
# ---------------------------------------------------------------------------

def dalpha_dt(report_hrs, alpha_dict, key="alpha"):
    """Central-difference dα/dt in 1/hr at each report point."""
    out = {}
    hrs = report_hrs
    for i, t in enumerate(hrs):
        if i == 0:
            if len(hrs) > 1:
                dt = hrs[1] - hrs[0]
                da = alpha_dict[hrs[1]][key] - alpha_dict[hrs[0]][key]
                out[t] = da / dt
            else:
                out[t] = 0.0
        elif i == len(hrs) - 1:
            dt = hrs[-1] - hrs[-2]
            da = alpha_dict[hrs[-1]][key] - alpha_dict[hrs[-2]][key]
            out[t] = da / dt
        else:
            dt = hrs[i + 1] - hrs[i - 1]
            da = alpha_dict[hrs[i + 1]][key] - alpha_dict[hrs[i - 1]][key]
            out[t] = da / dt
    return out


# ---------------------------------------------------------------------------
# Markdown table builder
# ---------------------------------------------------------------------------

def build_table(report_hrs, eng, cw):
    """Build comparison table for a single mix."""
    hdr = ("| t (hr) | Eng T (°F) | CW T (°F) | ΔT (°F) "
           "| Eng te (hr) | CW te (hr) | Δte (hr) "
           "| Eng α | CW α | Δα "
           "| Eng dα/dt (×10⁻³/hr) | CW dα/dt (×10⁻³/hr) |")
    sep = "|---|---|---|---|---|---|---|---|---|---|---|---|"

    eng_dalpha = dalpha_dt(report_hrs, eng, key="alpha")
    cw_dalpha  = dalpha_dt(report_hrs, cw,  key="alpha")

    rows = [hdr, sep]
    for t in report_hrs:
        e = eng[t]
        c = cw[t]
        dT  = e["T_F"]    - c["T_F"]
        dte = e["te_hrs"] - c["te_hrs"]
        da  = e["alpha"]  - c["alpha"]
        edr = eng_dalpha[t] * 1000.0
        cdr = cw_dalpha[t]  * 1000.0
        rows.append(
            f"| {t:4d} "
            f"| {e['T_F']:7.2f} "
            f"| {c['T_F']:7.2f} "
            f"| {dT:+6.2f} "
            f"| {e['te_hrs']:8.3f} "
            f"| {c['te_hrs']:8.3f} "
            f"| {dte:+7.3f} "
            f"| {e['alpha']:.4f} "
            f"| {c['alpha']:.4f} "
            f"| {da:+.4f} "
            f"| {edr:+7.4f} "
            f"| {cdr:+7.4f} |"
        )
    return "\n".join(rows)


# ---------------------------------------------------------------------------
# Divergence detector
# ---------------------------------------------------------------------------

def detect_divergence(report_hrs, eng, cw):
    """Find first time point where ΔT exceeds 1°F and peak ΔT."""
    first_t = None
    peak_dT = 0.0
    peak_t  = None
    for t in report_hrs:
        dT = abs(eng[t]["T_F"] - cw[t]["T_F"])
        if dT > abs(peak_dT) or peak_t is None:
            peak_dT = eng[t]["T_F"] - cw[t]["T_F"]
            peak_t  = t
        if first_t is None and dT > 1.0:
            first_t = t
    return first_t, peak_dT, peak_t


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Loading scenarios ...")
    scn01 = load_cw_scenario(
        str(CW_EXPORTS / "MIX-01" / "input.dat"),
        str(CW_EXPORTS / "MIX-01" / "weather.dat"),
        str(CW_EXPORTS / "MIX-01" / "output.txt"),
    )
    scn15 = load_cw_scenario(
        str(CW_EXPORTS / "MIX-15" / "input.dat"),
        str(CW_EXPORTS / "MIX-15" / "weather.dat"),
        str(CW_EXPORTS / "MIX-15" / "output.txt"),
    )
    print(f"  MIX-01 placement_temp_F = {scn01.construction.placement_temp_F}°F")
    print(f"  MIX-15 placement_temp_F = {scn15.construction.placement_temp_F}°F")

    # Mix parameters (same for both — V2 confirmed)
    Ea      = scn01.mix.activation_energy_J_mol
    tau     = scn01.mix.tau_hrs
    beta    = scn01.mix.beta
    alpha_u = scn01.mix.alpha_u
    print(f"  Ea={Ea:.1f} J/mol  tau={tau:.3f} hr  beta={beta:.4f}  alpha_u={alpha_u:.4f}")
    print(f"  T_REF_K={T_REF_K} K = {T_REF_K - 273.15:.1f}°C")
    print()

    print("Running engine on MIX-01 ...")
    res01, grid01 = run_engine(scn01)
    print(f"  Done. Peak T = {res01.peak_T_C * 9/5 + 32:.1f}°F")

    print("Running engine on MIX-15 ...")
    res15, grid15 = run_engine(scn15)
    print(f"  Done. Peak T = {res15.peak_T_C * 9/5 + 32:.1f}°F")
    print()

    print("Extracting engine centerline trajectories ...")
    eng01 = engine_centerline_at_times(res01, grid01, REPORT_HRS)
    eng15 = engine_centerline_at_times(res15, grid15, REPORT_HRS)

    print("Computing CW-implied te/α trajectories ...")
    cw01  = cw_centerline_at_times(scn01, REPORT_HRS, Ea, tau, beta, alpha_u)
    cw15  = cw_centerline_at_times(scn15, REPORT_HRS, Ea, tau, beta, alpha_u)

    # Divergence summary
    div01_t, div01_peak, div01_peak_t = detect_divergence(REPORT_HRS, eng01, cw01)
    div15_t, div15_peak, div15_peak_t = detect_divergence(REPORT_HRS, eng15, cw15)

    print(f"  MIX-01 first |ΔT|>1°F at t={div01_t} hr, peak ΔT={div01_peak:+.2f}°F at t={div01_peak_t} hr")
    print(f"  MIX-15 first |ΔT|>1°F at t={div15_t} hr, peak ΔT={div15_peak:+.2f}°F at t={div15_peak_t} hr")
    print()

    # --- Build Markdown ---
    table01 = build_table(REPORT_HRS, eng01, cw01)
    table15 = build_table(REPORT_HRS, eng15, cw15)

    # Verdict
    both_diverge    = (div01_t is not None and div15_t is not None)
    mix15_only      = (div01_t is None     and div15_t is not None)
    neither_diverges= (div01_t is None     and div15_t is None)

    if both_diverge:
        finding = "3"
        verdict_str = ("**Finding 3: Both MIX-01 and MIX-15 show ΔT > 1°F vs CW.**\n"
                       "There is a baseline implementation difference independent of cold IC.\n"
                       "This is a Sprint 5 priority — affects the whole library.")
    elif mix15_only:
        finding = "2"
        verdict_str = ("**Finding 2: MIX-15 diverges from CW; MIX-01 tracks closely.**\n"
                       "The discrepancy is cold-IC-specific. The engine and CW agree on warm\n"
                       "placement but diverge at 45°F. Identifying where te or α first separates\n"
                       "reveals whether this is a Sprint 5 calibration or Sprint 6 physics item.")
    else:
        finding = "1"
        verdict_str = ("**Finding 1: Both MIX-01 and MIX-15 temperature trajectories track CW closely.**\n"
                       "Cold-IC failure is downstream of the temperature field itself — possibly in\n"
                       "heat-generation amplitude, specific-heat handling, or the peak extraction.\n"
                       "H4-class hypotheses exhausted; route to Sprint 6 with 'no further diagnostic' note.")

    md = [
        "# PR 16 — Phase 2.6: Equivalent-Age and Hydration Trajectory Comparison",
        "",
        "Compares engine vs CW for MIX-01 (warm IC, 60°F) and MIX-15 (cold IC, 45°F).",
        "CW-implied te/α derived by integrating engine's Arrhenius model over the CW",
        "centerline temperature trajectory (mid-depth, centerline column).",
        "",
        "If the CW and engine diverge in the same place for both mixes: baseline difference.",
        "If only MIX-15 diverges: cold-IC-specific difference.",
        "If neither diverges: cold-IC failure is downstream of equivalent age.",
        "",
        "## Mix Parameters (identical for MIX-01 and MIX-15 per V2 finding)",
        "",
        f"| Parameter | Value |",
        f"| --- | --- |",
        f"| Ea (J/mol) | {Ea:.2f} |",
        f"| tau (hr) | {tau:.4f} |",
        f"| beta | {beta:.4f} |",
        f"| alpha_u | {alpha_u:.4f} |",
        f"| T_REF_K | {T_REF_K} K ({T_REF_K - 273.15:.1f}°C = {(T_REF_K - 273.15)*9/5 + 32:.1f}°F) |",
        "",
        "## Divergence Summary",
        "",
        f"| Mix | First |ΔT| > 1°F (hr) | Peak ΔT (°F) | At t (hr) |",
        f"| --- | --- | --- | --- |",
        f"| MIX-01 (warm IC) | {'t = ' + str(div01_t) if div01_t else 'none'} | {div01_peak:+.2f} | {div01_peak_t} |",
        f"| MIX-15 (cold IC) | {'t = ' + str(div15_t) if div15_t else 'none'} | {div15_peak:+.2f} | {div15_peak_t} |",
        "",
        f"### Finding {finding}",
        "",
        verdict_str,
        "",
        "## MIX-01 — Temperature / te / α Comparison (warm IC, 60°F)",
        "",
        "CW reference: MIX-01 CW. Engine comparison: engine(MIX-01, 60°F IC).",
        "ΔT = Engine − CW (positive = engine hotter). CW-implied te/α uses engine's Arrhenius model.",
        "",
        table01,
        "",
        "## MIX-15 — Temperature / te / α Comparison (cold IC, 45°F)",
        "",
        "CW reference: MIX-15 CW. Engine comparison: engine(MIX-15, 45°F IC).",
        "ΔT = Engine − CW (positive = engine hotter). CW-implied te/α uses engine's Arrhenius model.",
        "",
        table15,
        "",
        "## Interpretation Notes",
        "",
        "### te divergence interpretation",
        "Δte > 0 at a given clock time means the engine accumulates equivalent age faster than",
        "CW-implied te (i.e., the engine is hotter than CW at that time, driving more Arrhenius",
        "factor). Δte < 0 means the engine is lagging — colder → slower hydration → lower te.",
        "",
        "### Arrhenius asymmetry at cold IC",
        f"At 45°F (280.4 K), af = {float(np.exp(-Ea/R_GAS*(1/280.37 - 1/T_REF_K))):.4f}.",
        f"At 60°F (288.7 K), af = {float(np.exp(-Ea/R_GAS*(1/288.71 - 1/T_REF_K))):.4f}.",
        f"At 73.4°F (296.15 K, T_REF), af = 1.0000.",
        "Cold concrete starts at 72% of the warm hydration rate. If CW uses a different",
        "T_REF or Ea, this ratio changes and te accumulates at a different pace.",
        "",
        "### Limitation of this diagnostic",
        "CW-implied te is computed with the ENGINE's Arrhenius model (same Ea, T_REF_K=296.15 K).",
        "If CW uses the same model, the CW-implied te is directly comparable to the engine's te.",
        "If CW uses a different T_REF or a different Arrhenius form, the comparison shows the",
        "divergence attributable to model difference (not just temperature trajectory).",
    ]

    OUT_MD.write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"Wrote: {OUT_MD}")
    print()

    # --- Stdout summary ---
    print("=" * 70)
    print("PHASE 2.6 SUMMARY — ENGINE vs CW TRAJECTORY COMPARISON")
    print("=" * 70)
    print()
    print(f"  Arrhenius at 45°F: af = {float(np.exp(-Ea/R_GAS*(1/280.37-1/T_REF_K))):.4f}")
    print(f"  Arrhenius at 60°F: af = {float(np.exp(-Ea/R_GAS*(1/288.71-1/T_REF_K))):.4f}")
    print(f"  Arrhenius at T_REF: af = 1.0000")
    print()
    print("  MIX-01  Engine vs CW temperature (mid-depth centerline):")
    for t in [0, 12, 24, 36, 48, 72, 168]:
        if t in eng01 and t in cw01:
            print(f"    t={t:4d}hr: engine={eng01[t]['T_F']:6.2f}°F  CW={cw01[t]['T_F']:6.2f}°F  "
                  f"ΔT={eng01[t]['T_F']-cw01[t]['T_F']:+.2f}°F  "
                  f"Δte={eng01[t]['te_hrs']-cw01[t]['te_hrs']:+.3f}hr  "
                  f"Δα={eng01[t]['alpha']-cw01[t]['alpha']:+.4f}")
    print()
    print("  MIX-15  Engine vs CW temperature (mid-depth centerline):")
    for t in [0, 12, 24, 36, 48, 72, 168]:
        if t in eng15 and t in cw15:
            print(f"    t={t:4d}hr: engine={eng15[t]['T_F']:6.2f}°F  CW={cw15[t]['T_F']:6.2f}°F  "
                  f"ΔT={eng15[t]['T_F']-cw15[t]['T_F']:+.2f}°F  "
                  f"Δte={eng15[t]['te_hrs']-cw15[t]['te_hrs']:+.3f}hr  "
                  f"Δα={eng15[t]['alpha']-cw15[t]['alpha']:+.4f}")
    print()
    print(f"  FINDING {finding}:")
    print(f"  {verdict_str.replace(chr(10), chr(10)+'  ')}")
    print()
    print("  WAIT — review table + interpretation, confirm decision path before Phase 3.")


if __name__ == "__main__":
    main()
