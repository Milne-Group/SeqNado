"""Smoke tests for the seqnado benchmark command."""

import subprocess


def test_seqnado_benchmark_help_exits_zero():
    res = subprocess.run(["seqnado", "benchmark", "--help"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr
    assert "benchmark" in res.stdout.lower() or "usage" in res.stdout.lower()

