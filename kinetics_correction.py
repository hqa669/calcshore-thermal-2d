"""kinetics_correction.py

Composition-based Hu calibration for the CalcShore engine.

Implements the apr28 formula:

    Hu_eff = Hu_regressed × f_type × Σᵢ kᵢ · pᵢ

where pᵢ is the mass fraction of constituent i in total cementitious
content, f_type is a cement-type-specific factor, and kᵢ is an
SCM-class-specific factor.

Validated against CW adiabatic centerline curves on 13 of 15 mixes in
the engine v3 library at ±1°F. See validation/kinetics_correction/.

Module is dependency-free (pure-Python arithmetic on a CWMixDesign-shaped
object) so it can be imported from the loader without pulling numpy.
"""

from __future__ import annotations

import warnings


# ----------------------------------------------------------------------
# Calibration tables (apr28 cross-validation, T₀ = 73°F adiabatic)
# ----------------------------------------------------------------------

# Cement-type lookup. Values isolated via cement=1 lb/yd³ trace runs in
# CW so the SCM-weighted k-sum is identically 1.0 and f_type is the only
# remaining factor.
F_TYPE = {
    "Type I":    0.9982,
    "Type I/II": 0.9897,
    "Type II":   0.9870,
    "Type V":    0.9956,
    # "Type III": NOT CHARACTERIZED — uncommon for mass concrete.
    # If a Type III mix appears, run a CW reference and add the entry.
}

# SCM coefficients (moderate-replacement DOT envelope: 15–40% per SCM).
K_PORTLAND  = 1.0
K_FLY_ASH_F = 0.91
K_FLY_ASH_C = 0.883   # SINGLE high-replacement point; needs moderate-replacement validation
K_SLAG      = 0.89

# Provisional k_SF for SF replacement < 15%. Single MIX-13 data point
# suggests k_SF may actually amplify (>1.0) at 5% replacement; current
# value is the apr28 placeholder pending binary CW runs at 5/8/12% SF.
K_SILICA_FUME_PROVISIONAL = 0.7


class KineticsCorrectionWarning(UserWarning):
    """Emitted when an out-of-envelope SCM is encountered."""


# ----------------------------------------------------------------------
# Public API
# ----------------------------------------------------------------------

def compute_hu_factor(mix) -> tuple[float, str]:
    """Compute (factor, env_note) for a CWMixDesign-shaped object.

    `factor` is the multiplicative correction Hu_eff / Hu_regressed.
    `env_note` is empty if all constituents are in-envelope, otherwise
    a short flag string for diagnostics.

    The argument is duck-typed: any object with attributes
    ``cement_type``, ``cement_type_I_II_lb_yd3``, ``fly_ash_F_lb_yd3``,
    ``fly_ash_C_lb_yd3``, ``ggbfs_lb_yd3``, ``silica_fume_lb_yd3`` works.

    Raises
    ------
    KeyError
        Cement type not in F_TYPE table.
    ValueError
        Total cementitious is zero.
    NotImplementedError
        Cement = 0 with slag > 0 (CW α_u sentinel regime; engine doesn't
        currently support this).
    """
    cement_type = _normalize_cement_type(mix.cement_type)
    if cement_type not in F_TYPE:
        raise KeyError(
            f"Cement type {cement_type!r} not in F_TYPE table "
            f"(known: {sorted(F_TYPE)}). Run a CW adiabatic reference "
            f"with this cement type and add the entry."
        )
    f_type = F_TYPE[cement_type]

    cem  = mix.cement_type_I_II_lb_yd3
    fa_f = mix.fly_ash_F_lb_yd3
    fa_c = mix.fly_ash_C_lb_yd3
    slag = mix.ggbfs_lb_yd3
    sf   = mix.silica_fume_lb_yd3

    total = cem + fa_f + fa_c + slag + sf
    if total <= 0.0:
        raise ValueError("Total cementitious is zero; cannot compute Hu factor.")

    # CW α_u sentinel regime: cement=0 with slag>0 triggers a discontinuous
    # step in CW (α_u drops 0.99 → 0.28). The engine does not currently
    # support this regime; bail out rather than silently produce wrong Hu.
    if cem == 0.0 and slag > 0.0:
        raise NotImplementedError(
            "Pure-slag mix (cement=0) is in CW's α_u sentinel regime "
            "and is not supported by engine v3 kinetics. Customer "
            "should add even a small Portland fraction (≥5%)."
        )

    p_cem  = cem  / total
    p_fa_f = fa_f / total
    p_fa_c = fa_c / total
    p_slag = slag / total
    p_sf   = sf   / total

    notes: list[str] = []
    k_sf_used, sf_note = _silica_fume_handler(p_sf)
    if sf_note:
        notes.append(sf_note)

    k_weighted_sum = (
        K_PORTLAND  * p_cem
        + K_FLY_ASH_F * p_fa_f
        + K_FLY_ASH_C * p_fa_c
        + K_SLAG      * p_slag
        + k_sf_used   * p_sf
    )
    factor = f_type * k_weighted_sum
    return factor, "; ".join(notes)


# ----------------------------------------------------------------------
# Internal helpers
# ----------------------------------------------------------------------

def _normalize_cement_type(raw) -> str:
    """Tolerate small variations: 'I/II', 'Type I/II', '  Type I/II  '.

    The CW input.dat parser sometimes returns the bare designator
    (e.g. 'I/II') rather than the prefixed form ('Type I/II').
    """
    if not raw:
        return ""
    s = str(raw).strip()
    if s.startswith("Type "):
        return s
    return f"Type {s}"


def _silica_fume_handler(p_sf: float) -> tuple[float, str]:
    """Returns (k_sf, note) for the SF mass fraction.

    - p_sf ≈ 0      : k=0, no note (SF contributes nothing to the sum).
    - 0 < p_sf ≤ 15%: k=K_SILICA_FUME_PROVISIONAL, PROVISIONAL warning + note.
    - p_sf > 15%    : k=K_SILICA_FUME_PROVISIONAL, OUT-OF-ENVELOPE warning + note.
                      (CW returns sentinel Hu near zero in this regime,
                      so the engine result is provisional/likely wrong.)
    """
    if p_sf < 1e-6:
        return 0.0, ""
    if p_sf > 0.15:
        warnings.warn(
            f"Silica fume {p_sf*100:.1f}% exceeds the 15% envelope. "
            "CW returns sentinel Hu near zero in this regime. The engine "
            "result is provisional and likely wrong.",
            KineticsCorrectionWarning,
            stacklevel=3,
        )
        return K_SILICA_FUME_PROVISIONAL, f"SF={p_sf*100:.1f}% OUT-OF-ENVELOPE"
    warnings.warn(
        f"Silica fume {p_sf*100:.1f}%: k_SF={K_SILICA_FUME_PROVISIONAL} is "
        f"a provisional placeholder pending binary-mix CW runs. MIX-13 "
        f"empirically suggests k_SF may be ≈2.0 at 5% replacement.",
        KineticsCorrectionWarning,
        stacklevel=3,
    )
    return (
        K_SILICA_FUME_PROVISIONAL,
        f"SF={p_sf*100:.1f}% k_SF={K_SILICA_FUME_PROVISIONAL} PROVISIONAL",
    )
