"""
PR 13 smoke test: run_all with --group reference must complete, write valid markdown,
and show 5/5 S0 for all Reference mixes.
"""
import os
import subprocess
import sys
import tempfile

import pytest

pytest.importorskip("cw_scenario_loader")
pytest.importorskip("thermal_engine_2d")


def test_run_all_reference_smoke():
    if not os.path.isdir("validation/cw_exports/MIX-01"):
        pytest.skip("CW export scenarios not found")

    with tempfile.NamedTemporaryFile(suffix=".md", delete=False) as tf:
        md_path = tf.name

    try:
        proc = subprocess.run(
            [sys.executable, "run_all.py",
             "--group", "reference",
             "--quiet",
             "--output-md", md_path],
            capture_output=True, text=True,
        )
        assert proc.returncode == 0, (
            f"run_all.py exited {proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )

        assert os.path.isfile(md_path), "Markdown file was not written"
        md = open(md_path).read()

        # 1 header row + separator + 5 data rows
        lines = [l for l in md.splitlines() if l.startswith("|")]
        assert len(lines) == 7, f"Expected 7 table lines (header+sep+5 data), got {len(lines)}"

        # Every Reference mix must pass all S0 gates
        data_lines = lines[2:]  # skip header and separator
        for line in data_lines:
            assert "5/5" in line, f"Reference mix does not show 5/5 S0: {line}"
    finally:
        if os.path.isfile(md_path):
            os.unlink(md_path)
