"""
PR 13: run_one() must return a skipped sentinel dict for MIX-13 (no output.txt),
without running the engine.
"""
import time
import os
import pytest

pytest.importorskip("cw_scenario_loader")

from compare_to_cw import run_one  # noqa: E402

SCENARIO = "validation/cw_exports/MIX-13"


def test_mix13_skipped_sentinel():
    if not os.path.isdir(SCENARIO):
        pytest.skip("MIX-13 scenario directory not found")

    t0 = time.perf_counter()
    result = run_one(SCENARIO, png_path=None)
    elapsed = time.perf_counter() - t0

    assert result["skipped"] is True
    assert result["skip_reason"] == "no_cw_output"
    # No engine run — should complete in milliseconds, not seconds
    assert elapsed < 0.5, f"sentinel path took {elapsed:.2f}s — engine may have run"
    # Engine key must be absent (or None) — refactor must not sneak in a run
    assert "engine" not in result or result.get("engine") is None
