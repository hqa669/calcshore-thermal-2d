"""
PR 14 contract test: validation/sprint4_baseline.md must exist and
contain the expected partitioning.
"""
import os
import pytest

BASELINE = "validation/sprint4_baseline.md"


def test_baseline_file_exists():
    assert os.path.isfile(BASELINE), f"{BASELINE} not committed"


def test_baseline_has_15_rows():
    with open(BASELINE) as f:
        content = f.read()
    # Expect at least 15 mix references (MIX-01..15)
    for n in range(1, 16):
        assert f"MIX-{n:02d}" in content, f"MIX-{n:02d} missing from baseline"


def test_reference_set_passes_5_of_5():
    with open(BASELINE) as f:
        content = f.read()
    reference_mixes = ["MIX-01", "MIX-02", "MIX-03", "MIX-11", "MIX-12"]
    for mix in reference_mixes:
        # Find the row for this mix and confirm "5/5" appears in it
        for line in content.splitlines():
            if line.startswith(f"| {mix}") or f"| {mix} " in line:
                assert "5/5" in line, f"{mix} not 5/5 in baseline: {line}"
                break
        else:
            pytest.fail(f"{mix} row not found in baseline")


def test_mix13_skipped_row_present():
    with open(BASELINE) as f:
        content = f.read()
    assert "MIX-13" in content
    # Skipped row should have "skipped" or "no_cw_output" marker
    for line in content.splitlines():
        if "MIX-13" in line:
            assert "skipped" in line.lower() or "no_cw_output" in line, \
                f"MIX-13 row missing skipped marker: {line}"
            break
