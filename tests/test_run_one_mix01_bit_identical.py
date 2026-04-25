"""
PR 13 contract test: run_one() must be bit-identical to Sprint 3 close on MIX-01.
Any numeric deviation is a regression in the refactor, not a feature.
"""
import json
import os
import pytest

pytest.importorskip("cw_scenario_loader")
pytest.importorskip("thermal_engine_2d")

from compare_to_cw import run_one  # noqa: E402

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "mix01_sprint3_baseline.json")
SCENARIO = "validation/cw_exports/MIX-01"


def _load_baseline():
    with open(FIXTURE) as f:
        return json.load(f)


def test_mix01_bit_identical():
    if not os.path.isdir(SCENARIO):
        pytest.skip("MIX-01 scenario directory not found")

    expected = _load_baseline()
    result = run_one(SCENARIO, png_path=None)

    assert not result["skipped"], "run_one returned skipped for MIX-01"

    # Exact float equality — refactor must not change a single bit
    for key in ("peak_max_F", "peak_max_hr", "peak_grad_F", "peak_grad_hr"):
        assert result["engine"][key] == expected["engine"][key], (
            f"engine.{key}: got {result['engine'][key]!r}, expected {expected['engine'][key]!r}"
        )

    for key in ("peak_max_F", "peak_max_hr", "peak_grad_F", "peak_grad_hr"):
        assert result["cw"][key] == expected["cw"][key], (
            f"cw.{key}: got {result['cw'][key]!r}, expected {expected['cw'][key]!r}"
        )

    for key in ("peak_max_F", "peak_grad_F"):
        assert result["deltas"][key] == expected["deltas"][key], (
            f"deltas.{key}: got {result['deltas'][key]!r}, expected {expected['deltas'][key]!r}"
        )

    for key in ("field_F", "center_F", "corner_F"):
        assert result["rms"][key] == expected["rms"][key], (
            f"rms.{key}: got {result['rms'][key]!r}, expected {expected['rms'][key]!r}"
        )

    assert result["s0_overall"] == expected["s0_overall"]
    assert result["s0_pass_count"] == expected["s0_pass_count"]
