"""
PR 17 regression test: run_one() on MIX-01 must match post-PR-17 baseline.

PR 17 changed soil defaults from (lag=5.0, damping=0.7) to (lag=0.0,
damping=1.0), shifting MIX-01 numbers. This file replaces
test_run_one_mix01_bit_identical.py per PR 11 + PR 17 precedent: when a
calibration default changes by design, the bit-identical contract is moot
and the baseline resets. The original sprint-3-complete fixture is preserved
at tests/fixtures/mix01_sprint3_baseline.json for historical reference.
"""
import json
import os
import pytest

pytest.importorskip("cw_scenario_loader")
pytest.importorskip("thermal_engine_2d")

from compare_to_cw import run_one  # noqa: E402

FIXTURE = os.path.join(os.path.dirname(__file__), "fixtures", "mix01_pr17_baseline.json")
SCENARIO = "validation/cw_exports/MIX-01"


def _load_baseline():
    with open(FIXTURE) as f:
        return json.load(f)


def test_mix01_bit_identical_post_pr17():
    if not os.path.isdir(SCENARIO):
        pytest.skip("MIX-01 scenario directory not found")

    expected = _load_baseline()
    result = run_one(SCENARIO, png_path=None)

    assert not result["skipped"], "run_one returned skipped for MIX-01"

    # Exact float equality — any numeric deviation is a regression
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
