"""
PR 13: print_gate_table() must handle a skipped sentinel gracefully —
one-line output mentioning 'skipped' and skip_reason, no engine or gate tokens.
"""
import io
import contextlib

from compare_to_cw import print_gate_table


def test_print_gate_table_skipped_one_line():
    sentinel = {
        "scenario_dir": "validation/cw_exports/MIX-13",
        "skipped": True,
        "skip_reason": "no_cw_output",
        "scenario_name": "MIX-13",
    }
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print_gate_table(sentinel)
    out = buf.getvalue()

    assert "skipped" in out
    assert "no_cw_output" in out
    # Skipped path must not contain engine or gate tokens
    assert "Engine" not in out
    assert "PASS" not in out
    assert "FAIL" not in out
