"""
Tests for cw_scenario_loader: dew_point_C, sky_temp_C, and CWEnvironment
computed-field population.
"""

import os

import numpy as np
import pytest

from cw_scenario_loader import dew_point_C, sky_temp_C

MIX01_DIR = os.path.join(
    os.path.dirname(__file__), "..", "validation", "cw_exports", "MIX-01"
)
MIX01_INPUT = os.path.join(MIX01_DIR, "input.dat")
MIX01_WEATHER = os.path.join(MIX01_DIR, "weather.dat")
MIX01_OUTPUT = os.path.join(MIX01_DIR, "output.txt")


def test_dew_point_magnus_reference_values():
    # Reference points from published meteorology tables (Alduchov & Eskridge 1996)
    T = np.array([25.0, 20.0, 0.0, 30.0])
    RH = np.array([50.0, 80.0, 100.0, 30.0])
    dp = dew_point_C(T, RH)
    assert abs(dp[0] - 13.86) < 0.1,  f"(25°C, 50%) → expected ~13.86°C, got {dp[0]:.3f}"
    assert abs(dp[1] - 16.44) < 0.1,  f"(20°C, 80%) → expected ~16.44°C, got {dp[1]:.3f}"
    assert abs(dp[2] - 0.0)   < 0.05, f"(0°C, 100%) → expected ~0.00°C, got {dp[2]:.3f}"
    assert abs(dp[3] - 10.5)  < 0.2,  f"(30°C, 30%) → expected ~10.5°C, got {dp[3]:.3f}"


def test_sky_temp_berdahl_martin_clear_sky():
    T_air = np.array([25.0])
    T_dp = np.array([15.0])
    cloud = np.array([0.0])
    T_sky = sky_temp_C(T_air, T_dp, cloud)
    assert T_sky.shape == (1,)
    delta = float(T_air[0]) - float(T_sky[0])
    assert delta > 5.0, f"Clear sky should be >5°C colder than air; got ΔT={delta:.2f}°C"


def test_sky_temp_berdahl_martin_overcast():
    T_air = np.array([25.0])
    T_dp = np.array([15.0])
    cloud = np.array([8.0])
    T_sky = sky_temp_C(T_air, T_dp, cloud)
    delta = abs(float(T_sky[0]) - float(T_air[0]))
    assert delta < 2.0, f"Overcast sky should be within 2°C of air; got |ΔT|={delta:.2f}°C"


@pytest.fixture
def mix01_scenario():
    pytest.importorskip("cw_scenario_loader", reason="cw_scenario_loader not available")
    if not os.path.exists(MIX01_INPUT):
        pytest.skip(f"MIX-01 fixture not found at {MIX01_DIR}")
    from cw_scenario_loader import load_cw_scenario
    return load_cw_scenario(MIX01_INPUT, MIX01_WEATHER, MIX01_OUTPUT)


def test_cw_environment_has_hourly_sky_and_dewpoint(mix01_scenario):
    env = mix01_scenario.environment
    n = 24 * 7  # 7-day analysis (input.dat line 8 = "7 days")
    assert len(env.T_dp_C) == n,  f"T_dp_C length {len(env.T_dp_C)} != {n}"
    assert len(env.T_sky_C) == n, f"T_sky_C length {len(env.T_sky_C)} != {n}"
    assert len(env.T_dp_C) == len(env.T_air_F), "T_dp_C not aligned with T_air_F"
    assert 5.0 < env.T_dp_C.mean() < 30.0, (
        f"T_dp_C mean {env.T_dp_C.mean():.1f}°C not in expected Austin-July range 5-30°C"
    )
    assert -10.0 < env.T_sky_C.mean() < 30.0, (
        f"T_sky_C mean {env.T_sky_C.mean():.1f}°C not in expected range -10 to 30°C"
    )
    T_air_C_mean = (env.T_air_F.mean() - 32.0) * 5.0 / 9.0
    assert env.T_sky_C.mean() < T_air_C_mean, (
        f"T_sky_C mean {env.T_sky_C.mean():.1f}°C should be below T_air_C mean {T_air_C_mean:.1f}°C"
    )
