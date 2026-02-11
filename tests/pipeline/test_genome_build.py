"""Tests for seqnado genomes build workflow."""

from __future__ import annotations

import gzip
import shutil
import subprocess
from pathlib import Path

import pytest
from helpers import GenomeResources

GENOME_NAME = "chr21"


@pytest.mark.pipeline
@pytest.mark.snakemake
class TestGenomeBuildDryRun:
    """Verify the genome build Snakefile produces a valid DAG via dry-run."""

    def test_single_genome_dry_run(self, tmp_path: Path) -> None:
        """Dry-run for a single genome produces expected rules."""
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "--name",
                "hg38",
                "--outdir",
                "output",
                "--preset",
                "le",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=tmp_path,
        )

        assert result.returncode == 0, (
            f"Dry-run failed (rc={result.returncode}):\n{result.stderr}"
        )

        # Verify key rules appear in the dry-run plan
        expected_rules = [
            "download_fasta",
            "download_gtf",
            "download_chrom_sizes",
            "download_blacklist",
            "samtools_faidx",
            "build_bowtie2",
            "build_star",
        ]
        output = result.stdout + result.stderr
        for rule in expected_rules:
            assert rule in output, f"Expected rule '{rule}' not found in dry-run output"

    def test_multi_genome_dry_run(self, tmp_path: Path) -> None:
        """Dry-run for comma-separated genomes plans jobs for each."""
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "--name",
                "hg38,mm39",
                "--outdir",
                "output",
                "--preset",
                "le",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=tmp_path,
        )

        assert result.returncode == 0, (
            f"Dry-run failed (rc={result.returncode}):\n{result.stderr}"
        )

        output = result.stdout + result.stderr
        # Both genomes should appear in the plan
        assert "hg38" in output, "hg38 not found in dry-run output"
        assert "mm39" in output, "mm39 not found in dry-run output"

    def test_spikein_dry_run(self, tmp_path: Path) -> None:
        """Dry-run with --spikein plans composite genome rules."""
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "--name",
                "hg38",
                "--spikein",
                "dm6",
                "--outdir",
                "output",
                "--preset",
                "le",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=tmp_path,
        )

        assert result.returncode == 0, (
            f"Dry-run failed (rc={result.returncode}):\n{result.stderr}"
        )

        output = result.stdout + result.stderr
        assert "hg38_dm6" in output, "Composite genome name not in dry-run output"


def _stage_genome(resources: GenomeResources, out_dir: Path, genome_name: str) -> None:
    """Pre-stage downloaded files so Snakemake only runs indexing rules.

    Copies FASTA, GTF, chrom.sizes, and blacklist into the directory structure
    that the Snakefile_genome download rules would produce.
    """
    genome_dir = out_dir / genome_name

    (genome_dir / "sequence").mkdir(parents=True)
    (genome_dir / "genes").mkdir(parents=True)

    shutil.copy(resources.fasta, genome_dir / "sequence" / f"{genome_name}.fa")
    shutil.copy(
        resources.chromosome_sizes,
        genome_dir / "sequence" / f"{genome_name}.chrom.sizes",
    )
    shutil.copy(
        resources.gtf, genome_dir / "genes" / f"{genome_name}.ncbiRefSeq.gtf"
    )

    # Blacklist must be gzipped to match download rule output
    blacklist_gz = genome_dir / f"{genome_name}-blacklist.bed.gz"
    with open(resources.blacklist, "rb") as f_in, gzip.open(blacklist_gz, "wb") as f_out:
        f_out.writelines(f_in)


@pytest.mark.pipeline
@pytest.mark.snakemake
@pytest.mark.requires_apptainer
@pytest.mark.slow
class TestGenomeBuild:
    """Build indices from pre-staged chr21 test data via apptainer."""

    def test_build_single_genome(self, tmp_path: Path) -> None:
        """Build bt2, STAR, and faidx indices from chr21 FASTA."""
        resources = GenomeResources.download_resources(
            tmp_path / "test_genomes", "build"
        )
        out_dir = tmp_path / "output"
        _stage_genome(resources, out_dir, GENOME_NAME)

        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "--name",
                GENOME_NAME,
                "--outdir",
                str(out_dir),
                "--preset",
                "t",
                "--cores",
                "4",
            ],
            capture_output=True,
            text=True,
            timeout=600,
            cwd=tmp_path,
        )

        assert result.returncode == 0, (
            f"Build failed (rc={result.returncode}):\n{result.stderr}"
        )

        genome_dir = out_dir / GENOME_NAME
        assert (genome_dir / "sequence" / f"{GENOME_NAME}.fa.fai").exists()
        assert (genome_dir / "bt2_index" / f"{GENOME_NAME}.1.bt2").exists()
        assert (genome_dir / "STAR_2.7.10b").is_dir()
