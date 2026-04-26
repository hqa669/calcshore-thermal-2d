"""
MIX-02 within-Reference field diff (pre-PR-19 mini-spike).

Loads the 5 Reference scenarios (MIX-01, 02, 03, 11, 12) and reports
which CWMixDesign / CWGeometry / CWConstruction / CWEnvironment fields
vary across them. Routes the MIX-02 outlier finding from PR 18's
recon to one of three dispositions (A/B/C) per the spike prompt.

Throwaway spike — no engine changes, no test changes.
"""
import os
import sys
from dataclasses import fields
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, _REPO)

from cw_scenario_loader import load_cw_scenario

REFERENCE_MIXES = ("MIX-01", "MIX-02", "MIX-03", "MIX-11", "MIX-12")
SCENARIO_ROOT = os.path.join(_REPO, "validation", "cw_exports")

COMPOSITION_FIELDS = {
    "cement_type_I_II_lb_yd3", "fly_ash_F_lb_yd3", "fly_ash_C_lb_yd3",
    "ggbfs_lb_yd3", "silica_fume_lb_yd3", "water_lb_yd3",
    "coarse_agg_lb_yd3", "fine_agg_lb_yd3", "air_content_pct",
    "fly_ash_CaO_pct", "cement_type", "coarse_agg_type", "fine_agg_type",
}
KINETICS_FIELDS = {
    "activation_energy_J_mol", "tau_hrs", "beta", "alpha_u", "Hu_J_kg",
}
HEAT_THERMAL_PROP_FIELDS = {
    "thermal_conductivity_BTU_hr_ft_F", "aggregate_Cp_BTU_lb_F",
}
STRESS_ONLY_FIELDS = {
    "CTE_microstrain_F",   # used by stress/cracking analysis only — not by the heat equation
}


def load_all() -> dict:
    out = {}
    for mix in REFERENCE_MIXES:
        d = os.path.join(SCENARIO_ROOT, mix)
        out[mix] = load_cw_scenario(
            os.path.join(d, "input.dat"),
            os.path.join(d, "weather.dat"),
            os.path.join(d, "output.txt"),
        )
    return out


def scalar_diff(scenarios: dict, attr_name: str) -> dict:
    """Return {field_name: {mix: value}} for non-array fields that vary."""
    obj0 = getattr(scenarios[REFERENCE_MIXES[0]], attr_name)
    out = {}
    for f in fields(obj0):
        per_mix = {}
        skip = False
        for mix in REFERENCE_MIXES:
            v = getattr(getattr(scenarios[mix], attr_name), f.name)
            if isinstance(v, (np.ndarray, list)):
                skip = True
                break
            per_mix[mix] = v
        if skip:
            continue
        vals = list(per_mix.values())
        if not all(v == vals[0] for v in vals):
            out[f.name] = per_mix
    return out


def array_diff(scenarios: dict) -> dict:
    """For CWEnvironment ndarray fields: max abs diff vs MIX-01."""
    env0 = scenarios[REFERENCE_MIXES[0]].environment
    out = {}
    for f in fields(env0):
        v0 = getattr(env0, f.name)
        if not isinstance(v0, np.ndarray):
            continue
        per_mix = {}
        differs = False
        for mix in REFERENCE_MIXES:
            v = getattr(scenarios[mix].environment, f.name)
            if v.shape != v0.shape:
                per_mix[mix] = f"SHAPE: {v.shape}"
                differs = True
            else:
                d = float(np.max(np.abs(v - v0))) if v0.size > 0 else 0.0
                per_mix[mix] = d
                if d > 0.0:
                    differs = True
        if differs:
            out[f.name] = per_mix
    return out


def list_diff(scenarios: dict) -> dict:
    """For CWEnvironment List[float] fields."""
    env0 = scenarios[REFERENCE_MIXES[0]].environment
    out = {}
    for f in fields(env0):
        v0 = getattr(env0, f.name)
        if not isinstance(v0, list):
            continue
        per_mix = {}
        differs = False
        for mix in REFERENCE_MIXES:
            v = getattr(scenarios[mix].environment, f.name)
            if v != v0:
                differs = True
                per_mix[mix] = "differs"
            else:
                per_mix[mix] = "match"
        if differs:
            out[f.name] = per_mix
    return out


def fmt_scalar_table(varying: dict, attr_name: str) -> str:
    if not varying:
        return f"_No scalar fields vary in `{attr_name}` across the 5 Reference mixes._\n"
    header = "| Field | " + " | ".join(REFERENCE_MIXES) + " |"
    sep = "|---" * (len(REFERENCE_MIXES) + 1) + "|"
    lines = [header, sep]
    for fname, per_mix in varying.items():
        row_vals = " | ".join(str(per_mix[m]) for m in REFERENCE_MIXES)
        lines.append(f"| `{fname}` | {row_vals} |")
    return "\n".join(lines) + "\n"


def fmt_composition_table(scenarios: dict) -> str:
    header = "| Mix | cement | fly_ash_F | fly_ash_C | ggbfs | silica_fume | water | total | wcm |"
    sep = "|---|---|---|---|---|---|---|---|---|"
    lines = [header, sep]
    for mix in REFERENCE_MIXES:
        m = scenarios[mix].mix
        lines.append(
            f"| {mix} | {m.cement_type_I_II_lb_yd3:.1f} | {m.fly_ash_F_lb_yd3:.1f}"
            f" | {m.fly_ash_C_lb_yd3:.1f} | {m.ggbfs_lb_yd3:.1f}"
            f" | {m.silica_fume_lb_yd3:.1f} | {m.water_lb_yd3:.1f}"
            f" | {m.total_cementitious_lb_yd3:.1f} | {m.wcm:.4f} |"
        )
    return "\n".join(lines) + "\n"


def fmt_kinetics_table(scenarios: dict) -> str:
    header = "| Mix | Ea_J_mol | tau_hrs | beta | alpha_u | Hu_J_kg |"
    sep = "|---|---|---|---|---|---|"
    lines = [header, sep]
    for mix in REFERENCE_MIXES:
        m = scenarios[mix].mix
        lines.append(
            f"| {mix} | {m.activation_energy_J_mol:.1f} | {m.tau_hrs:.4f}"
            f" | {m.beta:.5f} | {m.alpha_u:.5f} | {m.Hu_J_kg:.1f} |"
        )
    return "\n".join(lines) + "\n"


def fmt_thermal_props_table(scenarios: dict) -> str:
    header = "| Mix | k_BTU_hr_ft_F | agg_Cp_BTU_lb_F | CTE_microstrain_F |"
    sep = "|---|---|---|---|"
    lines = [header, sep]
    for mix in REFERENCE_MIXES:
        m = scenarios[mix].mix
        lines.append(
            f"| {mix} | {m.thermal_conductivity_BTU_hr_ft_F:.4f}"
            f" | {m.aggregate_Cp_BTU_lb_F:.4f}"
            f" | {m.CTE_microstrain_F:.4f} |"
        )
    return "\n".join(lines) + "\n"


def classify(mix_diff, geom_diff, constr_diff, env_scalar_diff,
             env_array_diff, env_list_diff) -> tuple[str, str]:
    """Return (disposition_letter, evidence_string)."""
    # Partition mix_diff by category
    composition_vary = {f: v for f, v in mix_diff.items() if f in COMPOSITION_FIELDS}
    kinetics_vary   = {f: v for f, v in mix_diff.items() if f in KINETICS_FIELDS}
    heat_thermal_vary = {f: v for f, v in mix_diff.items() if f in HEAT_THERMAL_PROP_FIELDS}
    stress_only_vary = {f: v for f, v in mix_diff.items() if f in STRESS_ONLY_FIELDS}
    uncategorized_vary = {f: v for f, v in mix_diff.items()
                          if f not in COMPOSITION_FIELDS
                          and f not in KINETICS_FIELDS
                          and f not in HEAT_THERMAL_PROP_FIELDS
                          and f not in STRESS_ONLY_FIELDS}

    # Disposition A: any non-composition, non-kinetics, non-stress-only divergence
    # that the engine's heat equation actually consumes.
    a_triggers = (geom_diff or constr_diff or env_scalar_diff
                  or env_array_diff or env_list_diff
                  or heat_thermal_vary or uncategorized_vary)

    if a_triggers:
        bits = []
        if geom_diff: bits.append(f"geometry: {sorted(geom_diff)}")
        if constr_diff: bits.append(f"construction: {sorted(constr_diff)}")
        if env_scalar_diff: bits.append(f"environment scalars: {sorted(env_scalar_diff)}")
        if env_array_diff: bits.append(f"environment arrays: {sorted(env_array_diff)}")
        if env_list_diff: bits.append(f"environment lists: {sorted(env_list_diff)}")
        if heat_thermal_vary: bits.append(f"heat-equation thermal properties: {sorted(heat_thermal_vary)}")
        if uncategorized_vary: bits.append(f"uncategorized mix fields: {sorted(uncategorized_vary)}")
        return "A", "; ".join(bits)

    if kinetics_vary:
        evidence = f"hydration regression params vary across mixes: {sorted(kinetics_vary)}"
        if stress_only_vary:
            evidence += (f"; stress-only fields also vary ({sorted(stress_only_vary)})"
                         " but do not affect the heat equation")
        return "C", evidence

    if composition_vary:
        evidence = (f"composition varies as expected ({sorted(composition_vary)});"
                    " no upstream divergence in geometry, construction, environment,"
                    " heat-equation thermal properties, or hydration regression parameters")
        if stress_only_vary:
            evidence += (f"; stress-only fields vary ({sorted(stress_only_vary)})"
                         " but do not affect the heat equation")
        return "B", evidence

    return "B", ("no varying fields found in ANY dataclass — even composition is"
                 " identical across the 5 Reference mixes (would be the strongest"
                 " possible Disposition B)")


def detect_kinetics_subclusters(scenarios: dict) -> dict:
    """Group mixes by (Ea, tau, beta, Hu) tuple. alpha_u excluded — varies
    within the cluster on the Reference set even where Ea/tau/beta/Hu match.
    Returns {(Ea, tau, beta, Hu) tuple: [mix names]}.
    """
    groups: dict = {}
    for mix in REFERENCE_MIXES:
        m = scenarios[mix].mix
        key = (
            round(m.activation_energy_J_mol, 4),
            round(m.tau_hrs, 6),
            round(m.beta, 6),
            round(m.Hu_J_kg, 4),
        )
        groups.setdefault(key, []).append(mix)
    return groups


def build_routing(disposition: str, evidence: str, scenarios: dict = None) -> str:
    if disposition == "A":
        return (
            f"**Disposition A — hidden field divergence.** {evidence}. "
            "PR 19 should run ablation on the 4 homogeneous Reference mixes "
            "(MIX-01, 03, 11, 12); MIX-02 routes to PR 20+ with the divergent "
            "field characterized as the candidate driver of the outlier behavior."
        )
    if disposition == "C":
        cluster_text = ""
        if scenarios:
            groups = detect_kinetics_subclusters(scenarios)
            # Largest group = the kinetics cluster; remainder are outliers
            by_size = sorted(groups.items(), key=lambda x: len(x[1]), reverse=True)
            largest_key, largest_mixes = by_size[0]
            outlier_mixes = [m for m in REFERENCE_MIXES if m not in largest_mixes]
            cluster_repr_mix = largest_mixes[0]
            crm = scenarios[cluster_repr_mix].mix
            cluster_text = (
                f" Kinetics sub-cluster analysis: cluster {{{', '.join(largest_mixes)}}}"
                f" shares Ea={crm.activation_energy_J_mol:.1f},"
                f" tau={crm.tau_hrs:.4f} hr, beta={crm.beta:.5f},"
                f" Hu={crm.Hu_J_kg:.1f} J/kg (alpha_u varies within cluster —"
                f" see kinetics table). Outlier mixes with divergent"
                f" (Ea, tau, beta, Hu): {{{', '.join(outlier_mixes)}}}."
                f" Note: MIX-03 also has divergent kinetics but was NOT flagged as a"
                f" diagnostic outlier in PR 18 (D1/D3/D4 within the cluster signature);"
                f" PR 19 should target cluster {{{', '.join(largest_mixes)}}}"
                f" for the BC ablation and report MIX-02 and MIX-03 separately."
            )
        return (
            f"**Disposition C — hydration regression divergence.** {evidence}."
            f"{cluster_text} MIX-02's outlier behavior in PR 18 (D1 negative mean,"
            " D3 inverse amplitude, D4 bottom-concentrated profile) is"
            " hydration-kinetics-driven, not boundary-physics-driven."
            " MIX-02 routes to the hydration sprint series."
            " R9 implication: differing CW-computed kinetics on mixes with identical"
            " parsed composition (MIX-01/02/03 share composition exactly) indicates"
            " CW's regression consumes inputs not captured in our parsed CWMixDesign"
            " fields — likely cement Bogue compounds or admixture chemistry."
        )
    return (
        f"**Disposition B — composition-driven outlier with no upstream divergence.** "
        f"{evidence}. MIX-02's outlier behavior in PR 18 is rooted purely in "
        "composition-driven thermal physics (likely propagating through hydration "
        "heat curve shape from differing fly_ash_F/fly_ash_C/ggbfs ratios at the "
        "same total SCM%). PR 19 should run ablation on all 5 Reference mixes and "
        "report the cluster-vs-MIX-02 split honestly; the surface BC ablation will "
        "land on the cluster cleanly while MIX-02 is reported as a "
        "composition-driven offset that BC adjustments cannot resolve."
    )


def main():
    print("Loading 5 Reference scenarios...")
    scenarios = load_all()
    print(f"Loaded: {list(scenarios.keys())}")

    mix_diff      = scalar_diff(scenarios, "mix")
    geom_diff     = scalar_diff(scenarios, "geometry")
    constr_diff   = scalar_diff(scenarios, "construction")
    env_scalar_d  = scalar_diff(scenarios, "environment")
    env_array_d   = array_diff(scenarios)
    env_list_d    = list_diff(scenarios)

    disposition, evidence = classify(
        mix_diff, geom_diff, constr_diff, env_scalar_d, env_array_d, env_list_d,
    )

    sections = [
        "# MIX-02 within-Reference field diff",
        "",
        "Pre-PR-19 mini-spike. Compares CWMixDesign / CWGeometry / CWConstruction"
        " / CWEnvironment fields across the 5 Reference mixes (MIX-01, 02, 03,"
        " 11, 12) to characterize MIX-02's outlier behavior in PR 18's recon.",
        "",
        "## Disposition: " + disposition,
        "",
        build_routing(disposition, evidence, scenarios),
        "",
        "## Construction fields that vary",
        "",
        fmt_scalar_table(constr_diff, "construction"),
        "## Geometry fields that vary",
        "",
        fmt_scalar_table(geom_diff, "geometry"),
        "## Environment scalar fields that vary",
        "",
        fmt_scalar_table(env_scalar_d, "environment"),
        "## Environment ndarray fields that vary (max abs diff vs MIX-01)",
        "",
        fmt_scalar_table(env_array_d, "environment.ndarray") if env_array_d
            else "_All hourly weather arrays bit-identical across the 5 Reference mixes._\n",
        "## Environment list fields that vary",
        "",
        fmt_scalar_table(env_list_d, "environment.list") if env_list_d
            else "_All daily-summary lists bit-identical across the 5 Reference mixes._\n",
        "## Mix design — composition (lb/yd³)",
        "",
        fmt_composition_table(scenarios),
        "## Mix design — hydration regression parameters",
        "",
        fmt_kinetics_table(scenarios),
        "## Mix design — derived thermal properties",
        "",
        fmt_thermal_props_table(scenarios),
    ]

    md = "\n".join(sections)
    out_path = os.path.join(_HERE, "mix02_recon.md")
    with open(out_path, "w") as f:
        f.write(md)

    print(f"\nDisposition: {disposition}")
    print(f"Evidence: {evidence}")
    print(f"\nWritten: {out_path}")


if __name__ == "__main__":
    main()
