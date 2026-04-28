"""tests/test_kinetics_correction.py

Pin the apr28 calibration:
- 4 cross-validation points from the apr28 passdown (composition → factor)
- 14-mix engine v3 library factors (validated by batch_compare_all_mixes.py
  against CW adiabatic centerline curves)
- Envelope guards (SF, pure-slag, unknown cement type, zero cementitious)
- Cement-type string normalization

If any coefficient in kinetics_correction.py changes, these break.
"""

import pytest

from kinetics_correction import (
    F_TYPE,
    K_FLY_ASH_C,
    K_FLY_ASH_F,
    K_PORTLAND,
    K_SLAG,
    KineticsCorrectionWarning,
    compute_hu_factor,
)


class FakeMix:
    """Minimal CWMixDesign-shaped object for testing the correction in isolation."""

    def __init__(self, *, cement_type, cem=0.0, fa_f=0.0, fa_c=0.0, slag=0.0, sf=0.0):
        self.cement_type             = cement_type
        self.cement_type_I_II_lb_yd3 = cem
        self.fly_ash_F_lb_yd3        = fa_f
        self.fly_ash_C_lb_yd3        = fa_c
        self.ggbfs_lb_yd3            = slag
        self.silica_fume_lb_yd3      = sf


# ----------------------------------------------------------------------
# apr28 cross-validation table
# (kinetics_passdown_apr28.md, "Cross-validation across 4 representative mixes")
# ----------------------------------------------------------------------

@pytest.mark.parametrize("label, mix, expected", [
    ("MIX-00A 100% I/II",
     FakeMix(cement_type="Type I/II", cem=575),
     0.9897),
    ("MIX-04-apr28 78% I/II + 22% FA-F",
     FakeMix(cement_type="Type I/II", cem=448.5, fa_f=126.5),
     0.9704),
    ("MIX-01-apr28 61% I/II + 22% FA-F + 17% slag",
     FakeMix(cement_type="Type I/II", cem=350.75, fa_f=126.5, slag=97.75),
     0.9514),
    ("MIX-07-apr28 9% I/II + 22% FA-F + 70% slag",
     FakeMix(cement_type="Type I/II", cem=51.75, fa_f=126.5, slag=402.5),
     0.8948),
])
def test_apr28_cross_validation(label, mix, expected):
    factor, _ = compute_hu_factor(mix)
    assert abs(factor - expected) < 5e-4, f"{label}: {factor} vs {expected}"


# ----------------------------------------------------------------------
# 14-mix engine v3 library factors (validated by batch_compare_all_mixes.py)
# ----------------------------------------------------------------------

@pytest.mark.parametrize("mix_id, mix, expected", [
    # Reference cluster — 39.1% SCM, varying cement types
    ("MIX-01", FakeMix(cement_type="Type I/II", cem=350, fa_f=125, slag=100), 0.9514),
    ("MIX-02", FakeMix(cement_type="Type V",    cem=350, fa_f=125, slag=100), 0.9571),
    ("MIX-03", FakeMix(cement_type="Type I",    cem=350, fa_f=125, slag=100), 0.9596),
    # B1 — low-SCM
    ("MIX-04", FakeMix(cement_type="Type I",    cem=450, fa_f=125),           0.9787),
    ("MIX-08", FakeMix(cement_type="Type I/II", cem=475, slag=100),           0.9708),
    # Cluster A — high-SCM slag-heavy
    ("MIX-05", FakeMix(cement_type="Type I/II", cem=275, fa_f=125, slag=175), 0.9372),
    ("MIX-06", FakeMix(cement_type="Type I/II", cem=165, fa_f=125, slag=285), 0.9164),
    ("MIX-07", FakeMix(cement_type="Type I/II", cem=50,  fa_f=125, slag=400), 0.8946),
    # FA-C variant
    ("MIX-09", FakeMix(cement_type="Type I/II", cem=350, fa_c=125, slag=100), 0.9456),
    # High-FA variant
    ("MIX-10", FakeMix(cement_type="Type I/II", cem=275, fa_f=200, slag=100), 0.9398),
    # Reference cluster duplicates (different w/cm)
    ("MIX-11", FakeMix(cement_type="Type I/II", cem=350, fa_f=125, slag=100), 0.9514),
    ("MIX-12", FakeMix(cement_type="Type I/II", cem=350, fa_f=125, slag=100), 0.9514),
    ("MIX-14", FakeMix(cement_type="Type I/II", cem=350, fa_f=125, slag=100), 0.9514),
    ("MIX-15", FakeMix(cement_type="Type I/II", cem=350, fa_f=125, slag=100), 0.9514),
])
def test_v3_library_hu_factors(mix_id, mix, expected):
    factor, _ = compute_hu_factor(mix)
    assert abs(factor - expected) < 5e-4, f"{mix_id}: {factor} vs {expected}"


# ----------------------------------------------------------------------
# Envelope guards
# ----------------------------------------------------------------------

def test_silica_fume_moderate_emits_warning():
    """At ≤15% SF replacement: warning + 'PROVISIONAL' flag in the note."""
    mix = FakeMix(cement_type="Type I/II", cem=525, sf=50)  # 8.7% SF
    with pytest.warns(KineticsCorrectionWarning, match="provisional placeholder"):
        factor, note = compute_hu_factor(mix)
    assert "PROVISIONAL" in note


def test_silica_fume_high_emits_oob_warning():
    """Above 15% SF replacement: warning + 'OUT-OF-ENVELOPE' flag in the note."""
    mix = FakeMix(cement_type="Type I/II", cem=400, sf=175)  # 30% SF
    with pytest.warns(KineticsCorrectionWarning, match="exceeds the 15% envelope"):
        factor, note = compute_hu_factor(mix)
    assert "OUT-OF-ENVELOPE" in note


def test_pure_slag_raises():
    mix = FakeMix(cement_type="Type I/II", cem=0, slag=575)
    with pytest.raises(NotImplementedError, match="sentinel"):
        compute_hu_factor(mix)


def test_unknown_cement_type_raises():
    mix = FakeMix(cement_type="Type III", cem=575)
    with pytest.raises(KeyError, match="Type III"):
        compute_hu_factor(mix)


def test_zero_total_raises():
    mix = FakeMix(cement_type="Type I/II")
    with pytest.raises(ValueError, match="zero"):
        compute_hu_factor(mix)


def test_normalization_tolerates_bare_designator():
    """The loader sometimes returns 'I/II' instead of 'Type I/II'."""
    mix = FakeMix(cement_type="I/II", cem=575)
    factor, _ = compute_hu_factor(mix)
    assert abs(factor - 0.9897) < 5e-4


# ----------------------------------------------------------------------
# Sanity: lookup tables haven't been edited
# ----------------------------------------------------------------------

def test_f_type_table_pinned():
    assert F_TYPE == {
        "Type I":    0.9982,
        "Type I/II": 0.9897,
        "Type II":   0.9870,
        "Type V":    0.9956,
    }


def test_k_coefficients_pinned():
    assert K_PORTLAND  == 1.0
    assert K_FLY_ASH_F == 0.91
    assert K_FLY_ASH_C == 0.883
    assert K_SLAG      == 0.89
