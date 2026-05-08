"""Tests for pipeline CLI helpers."""

from seqnado.cli.commands.pipeline import _ensure_default_snakemake_flag


def test_ensure_default_snakemake_flag_adds_missing_flag() -> None:
    options = ["--dry-run", "--unlock"]

    result = _ensure_default_snakemake_flag(options, "--benchmark-extended")

    assert result == ["--benchmark-extended", "--dry-run", "--unlock"]


def test_ensure_default_snakemake_flag_preserves_existing_flag() -> None:
    options = ["--benchmark-extended", "--dry-run"]

    result = _ensure_default_snakemake_flag(options, "--benchmark-extended")

    assert result == ["--benchmark-extended", "--dry-run"]
