"""
PR 13: structural tests for GROUPS dict in run_all.py — no engine runs.
"""
from run_all import GROUPS


def test_group_keys():
    assert set(GROUPS) == {"reference", "b1", "b2", "cluster_a", "all", "evaluation_set"}


def test_group_sizes():
    assert len(GROUPS["reference"]) == 5
    assert len(GROUPS["b1"]) == 2
    assert len(GROUPS["b2"]) == 1
    assert len(GROUPS["cluster_a"]) == 5
    assert len(GROUPS["evaluation_set"]) == 8
    assert len(GROUPS["all"]) == 14


def test_evaluation_set_is_reference_plus_b1_plus_b2():
    assert set(GROUPS["evaluation_set"]) == (
        set(GROUPS["reference"]) | set(GROUPS["b1"]) | set(GROUPS["b2"])
    )


def test_mix13_absent_from_all_named_groups():
    for grp_name, members in GROUPS.items():
        assert "MIX-13" not in members, f"MIX-13 unexpectedly found in group '{grp_name}'"


def test_all_excludes_mix13():
    assert "MIX-13" not in GROUPS["all"]


def test_all_contains_14_mixes():
    assert len(GROUPS["all"]) == 14
