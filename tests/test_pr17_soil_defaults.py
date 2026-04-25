"""
PR 17 contract test: soil_lag_hrs and soil_damping default to no-op values
per R4 disposition.
"""
from cw_scenario_loader import CWConstruction


def test_soil_lag_hrs_default_is_zero():
    """Per R4 disposition: default lag is 0.0 (no-op, reduces to T_amb)."""
    c = CWConstruction()
    assert c.soil_lag_hrs == 0.0


def test_soil_damping_default_is_one():
    """Per R4 disposition: default damping is 1.0 (no-op, reduces to T_amb)."""
    c = CWConstruction()
    assert c.soil_damping == 1.0


def test_soil_default_pair_is_no_op():
    """Per ADR-08 (PR 17 update): default pair (lag=0.0, damping=1.0) reduces
    to T_amb-pinned ground behavior."""
    c = CWConstruction()
    # The combination lag=0 + damping=1 means T_ground(t) = T_amb_daily_mean +
    # (T_amb(t-0) - T_amb_daily_mean) * 1.0 = T_amb(t).
    # This is the pre-Sprint-3 behavior.
    assert c.soil_lag_hrs == 0.0 and c.soil_damping == 1.0
