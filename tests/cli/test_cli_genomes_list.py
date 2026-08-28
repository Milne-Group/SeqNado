"""Tests for `seqnado genomes list` CLI command."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path


def _write_genome_config(tmp_root: Path, genome_data: dict) -> Path:
    cfg_dir = tmp_root / ".config" / "seqnado"
    cfg_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = cfg_dir / "genome_config.json"
    cfg_path.write_text(json.dumps(genome_data, indent=2))
    return cfg_path


def _make_bt2_index(tmp_root: Path, name: str) -> str:
    bt2_dir = tmp_root / "bt2"
    bt2_dir.mkdir(parents=True, exist_ok=True)
    prefix = bt2_dir / name
    for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
        (bt2_dir / f"{name}{suffix}").touch()
    return str(prefix)


def _make_star_index(tmp_root: Path, name: str) -> str:
    star_dir = tmp_root / f"STAR_{name}"
    star_dir.mkdir(parents=True, exist_ok=True)
    (star_dir / "SA").write_text("mock")
    return str(star_dir)


def _run_genomes_list(assay: str, tmp_path: Path, monkeypatch) -> subprocess.CompletedProcess:
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
    return subprocess.run(
        ["seqnado", "genomes", "list", assay],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )


def test_genomes_list_valid_config(tmp_path: Path, monkeypatch):
    bt2_index = _make_bt2_index(tmp_path, "hg38")
    _write_genome_config(tmp_path, {"hg38": {"bt2_index": bt2_index}})

    result = _run_genomes_list("atac", tmp_path, monkeypatch)

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert "hg38" in result.stdout
    assert "Bowtie2" in result.stdout


def test_genomes_list_rna_shows_star_index_not_bowtie(tmp_path: Path, monkeypatch):
    """Regression test: `genomes list rna` must show the STAR index, not silently
    fall back to Bowtie2 due to the assay str/enum mismatch."""
    star_index = _make_star_index(tmp_path, "hg38")
    _write_genome_config(tmp_path, {"hg38": {"star_index": star_index}})

    result = _run_genomes_list("rna", tmp_path, monkeypatch)

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert "hg38" in result.stdout
    assert "STAR" in result.stdout
    assert star_index in result.stdout
    assert "Bowtie2" not in result.stdout


def test_genomes_list_atac_flags_missing_bowtie_index_as_partial(tmp_path: Path, monkeypatch):
    """A genome with only star_index set has no Bowtie2 index for ATAC -- this
    must be flagged, not silently reported as if the genome were fully configured."""
    star_index = _make_star_index(tmp_path, "hg38")
    _write_genome_config(tmp_path, {"hg38": {"star_index": star_index}})

    result = _run_genomes_list("atac", tmp_path, monkeypatch)

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert "hg38" in result.stdout
    assert "No Bowtie2 index configured" in result.stdout


def test_genomes_list_mixed_valid_and_invalid_entries(tmp_path: Path, monkeypatch):
    bt2_index = _make_bt2_index(tmp_path, "hg38")
    _write_genome_config(
        tmp_path,
        {
            "hg38": {"bt2_index": bt2_index},
            "broken": {"bt2_index": bt2_index, "fasta": "/does/not/exist.fa"},
        },
    )

    result = _run_genomes_list("atac", tmp_path, monkeypatch)

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert "hg38" in result.stdout
    assert "broken" in result.stdout
    assert "failed to load" in result.stdout
    assert "Genome file not found" in result.stdout
    assert "genome_config.json" in result.stdout
    assert "seqnado genomes edit" in result.stdout


def test_genomes_list_all_entries_invalid(tmp_path: Path, monkeypatch):
    _write_genome_config(
        tmp_path,
        {"broken": {"fasta": "/does/not/exist.fa"}},
    )

    result = _run_genomes_list("atac", tmp_path, monkeypatch)

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    assert "broken" in result.stdout
    assert "failed to load" in result.stdout
    assert "No genome config found" not in result.stdout


def test_genomes_list_unknown_assay_errors(tmp_path: Path, monkeypatch):
    _write_genome_config(tmp_path, {})

    result = _run_genomes_list("not-a-real-assay", tmp_path, monkeypatch)

    assert result.returncode == 2, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
