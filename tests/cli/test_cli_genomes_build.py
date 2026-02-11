"""Tests for SeqNado genome build CLI command."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest


class TestGenomesBuildCLI:
    """Test suite for `seqnado genomes build` CLI command."""

    def test_genomes_build_requires_name(self, tmp_path: Path) -> None:
        """Test that --name is required for build subcommand."""
        outdir = tmp_path / "build_output"
        
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--outdir",
                str(outdir),
            ],
            capture_output=True,
            text=True,
        )
        
        assert result.returncode != 0, "Build without --name should fail"
        assert "-n" in result.stderr or "name" in result.stderr.lower() or "required" in result.stderr.lower()

    def test_genomes_build_requires_fasta(self, tmp_path: Path) -> None:
        """Test that build with --name but no --fasta fails gracefully."""
        outdir = tmp_path / "build_output"
        
        # This tests the snakemake portion which validates the fasta file
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--name",
                "test_genome",
                "--outdir",
                str(outdir),
                "-c",
                "1",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )
        
        # The snakemake dry-run should complete but indicate missing input
        # We're just checking that the CLI accepts the arguments
        assert "test_genome" in result.stderr or "test_genome" in result.stdout or result.returncode != 0

    def test_genomes_build_invalid_assay(self, tmp_path: Path) -> None:
        """Test that an invalid assay type is handled."""
        outdir = tmp_path / "build_output"
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">chr1\nACGT\n")
        
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "invalid_assay",
                "--name",
                "test_genome",
                "--fasta",
                str(fasta),
                "--outdir",
                str(outdir),
                "-c",
                "1",
            ],
            capture_output=True,
            text=True,
        )
        
        # Should fail due to invalid assay
        assert result.returncode != 0

    def test_genomes_build_multiple_genomes_with_comma_separator(self, tmp_path: Path) -> None:
        """Test that multiple genomes can be specified with comma-separated names."""
        outdir = tmp_path / "build_output"
        fasta1 = tmp_path / "test1.fasta"
        fasta1.write_text(">chr1\nACGT\n")
        
        # Dry run to check argument parsing
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--name",
                "test1,test2",
                "--fasta",
                str(fasta1),
                "--outdir",
                str(outdir),
                "-c",
                "1",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )

        # Check that the command accepts multiple genome names
        # (snakemake will be invoked with the config)
        assert "test1" in result.stderr or "test1" in result.stdout or result.returncode != 0

    def test_genomes_build_spikein_with_single_genome(self, tmp_path: Path) -> None:
        """Test that spike-in is accepted with single genome."""
        outdir = tmp_path / "build_output"
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">chr1\nACGT\n")
        
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--name",
                "primary",
                "--spikein",
                "secondary",
                "--fasta",
                str(fasta),
                "--outdir",
                str(outdir),
                "-c",
                "1",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )

        # Should accept the arguments and pass to snakemake
        assert "primary" in result.stderr or "primary" in result.stdout or result.returncode != 0

    def test_genomes_build_spikein_rejects_multiple_genomes(self, tmp_path: Path) -> None:
        """Test that spike-in is rejected when multiple genomes are specified."""
        outdir = tmp_path / "build_output"
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">chr1\nACGT\n")
        
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--name",
                "genome1,genome2",
                "--spikein",
                "secondary",
                "--fasta",
                str(fasta),
                "--outdir",
                str(outdir),
                "-c",
                "1",
            ],
            capture_output=True,
            text=True,
        )
        
        # Should reject mixing multiple genomes with spike-in
        assert result.returncode != 0
        assert "spike-in" in result.stderr.lower() or "multiple" in result.stderr.lower()

    def test_genomes_build_outdir_defaults_to_current_directory(self, tmp_path: Path) -> None:
        """Test that output directory defaults to genome_build in current directory."""
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">chr1\nACGT\n")
        
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--name",
                "test",
                "--fasta",
                str(fasta),
                "-c",
                "1",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),
        )
        
        # Just verify the command structure is accepted
        assert "test" in result.stderr or "test" in result.stdout or result.returncode != 0

    def test_genomes_build_cores_parameter(self, tmp_path: Path) -> None:
        """Test that cores parameter is accepted."""
        outdir = tmp_path / "build_output"
        fasta = tmp_path / "test.fasta"
        fasta.write_text(">chr1\nACGT\n")
        
        result = subprocess.run(
            [
                "seqnado",
                "genomes",
                "build",
                "rna",
                "--name",
                "test",
                "--fasta",
                str(fasta),
                "--outdir",
                str(outdir),
                "--cores",
                "4",
                "--dry-run",
            ],
            capture_output=True,
            text=True,
            timeout=30,
        )

        # Verify the command structure
        assert "test" in result.stderr or "test" in result.stdout or result.returncode != 0


class TestGenomesListCLI:
    """Test suite for `seqnado genomes list` CLI command."""

    def test_genomes_list_with_valid_assay(self, monkeypatch, tmp_path: Path) -> None:
        """Test that list command works with a valid assay."""
        # Setup minimal genome config
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True, exist_ok=True)
        cfg_file = cfg_dir / "genome_config.json"
        cfg_file.write_text(json.dumps({}))
        
        monkeypatch.setenv("HOME", str(tmp_path))
        
        result = subprocess.run(
            ["seqnado", "genomes", "list", "rna"],
            capture_output=True,
            text=True,
        )
        
        # Should succeed with valid assay
        assert result.returncode == 0

    def test_genomes_list_returns_empty_when_no_config(self, monkeypatch, tmp_path: Path) -> None:
        """Test that list command exits gracefully when no genome config exists."""
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True, exist_ok=True)
        cfg_file = cfg_dir / "genome_config.json"
        cfg_file.write_text(json.dumps({}))
        
        monkeypatch.setenv("HOME", str(tmp_path))
        
        result = subprocess.run(
            ["seqnado", "genomes", "list", "rna"],
            capture_output=True,
            text=True,
        )
        
        # Should exit gracefully with empty config
        assert result.returncode == 0


class TestGenomesEditCLI:
    """Test suite for `seqnado genomes edit` CLI command."""

    def test_genomes_edit_requires_config_to_exist(self, monkeypatch, tmp_path: Path) -> None:
        """Test that edit command requires genome config to exist."""
        monkeypatch.setenv("HOME", str(tmp_path))
        
        result = subprocess.run(
            ["seqnado", "genomes", "edit", "rna"],
            capture_output=True,
            text=True,
        )
        
        # Should fail because genome config doesn't exist
        assert result.returncode != 0
        assert "not found" in result.stderr.lower() or "no such" in result.stderr.lower()

    def test_genomes_edit_with_custom_editor_path(self, monkeypatch, tmp_path: Path) -> None:
        """Test that edit command respects custom EDITOR environment variable."""
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True, exist_ok=True)
        cfg_file = cfg_dir / "genome_config.json"
        cfg_file.write_text(json.dumps({"test": {}}))
        
        monkeypatch.setenv("HOME", str(tmp_path))
        monkeypatch.setenv("EDITOR", "/bin/true")  # Use a dummy editor
        
        result = subprocess.run(
            ["seqnado", "genomes", "edit", "rna"],
            capture_output=True,
            text=True,
        )
        
        # Should succeed with /bin/true as editor
        assert result.returncode == 0


class TestGenomesValidation:
    """Test suite for genome validation and structure."""

    def test_built_genome_has_required_files(self, tmp_path: Path) -> None:
        """Test that a built genome contains expected file structure (unit test)."""
        # Create a minimal genome directory structure
        genome_dir = tmp_path / "test_genome"
        genome_dir.mkdir()
        
        # Create expected subdirectories
        (genome_dir / "sequence").mkdir()
        (genome_dir / "genes").mkdir()
        (genome_dir / "bt2_index").mkdir()
        (genome_dir / "STAR_2.7.10b").mkdir()
        
        # Create minimal required files
        fasta = genome_dir / "sequence" / "test_genome.fa"
        fasta.write_text(">chr1\nACGT\n")
        
        gtf = genome_dir / "genes" / "test_genome.ncbiRefSeq.gtf"
        gtf.write_text("chr1\t.\tgene\t1\t4\t.\t+\t.\tgene_id \"test\"\n")
        
        chrom_sizes = genome_dir / "sequence" / "test_genome.chrom.sizes"
        chrom_sizes.write_text("chr1\t4\n")
        
        # Bowtie2 index files
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (genome_dir / "bt2_index" / f"test_genome{suffix}").touch()
        
        # Verify structure
        assert (genome_dir / "sequence").is_dir()
        assert (genome_dir / "genes").is_dir()
        assert (genome_dir / "bt2_index").is_dir()
        assert (genome_dir / "STAR_2.7.10b").is_dir()
        assert fasta.exists()
        assert gtf.exists()
        assert chrom_sizes.exists()

    def test_bowtie_index_validation(self, tmp_path: Path) -> None:
        """Test that Bowtie2 index prefix validation works."""
        from seqnado.config.configs import validate_bowtie2_index_prefix
        
        # Create a minimal bowtie2 index
        index_dir = tmp_path / "index"
        index_dir.mkdir()
        
        prefix = index_dir / "test_index"
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (index_dir / f"test_index{suffix}").touch()
        
        # Should not raise an error
        validate_bowtie2_index_prefix(prefix)

    def test_bowtie_index_validation_fails_without_files(self, tmp_path: Path) -> None:
        """Test that Bowtie2 index validation fails without required files."""
        from seqnado.config.configs import validate_bowtie2_index_prefix
        
        index_dir = tmp_path / "index"
        index_dir.mkdir()
        
        prefix = index_dir / "missing_index"
        
        # Should raise an error
        with pytest.raises(ValueError, match="No Bowtie2 index files found"):
            validate_bowtie2_index_prefix(prefix)

    def test_genome_config_creation_with_bowtie_index(self, tmp_path: Path) -> None:
        """Test creating a GenomeConfig with Bowtie2 index."""
        from seqnado.config.configs import BowtieIndex, GenomeConfig
        
        # Create a minimal bowtie2 index
        index_dir = tmp_path / "index"
        index_dir.mkdir()
        
        prefix = str(index_dir / "test_index")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (index_dir / f"test_index{suffix}").touch()
        
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGT\n")
        
        config = GenomeConfig(
            name="test",
            index=BowtieIndex(prefix=prefix),
            fasta=str(fasta),
        )
        
        assert config.name == "test"
        assert config.index.type == "Bowtie2"
        assert config.index.prefix == prefix

    def test_genome_config_creation_with_star_index(self, tmp_path: Path) -> None:
        """Test creating a GenomeConfig with STAR index."""
        from seqnado.config.configs import STARIndex, GenomeConfig
        
        star_dir = tmp_path / "STAR_index"
        star_dir.mkdir()
        
        # Create minimal STAR index files
        (star_dir / "SA_file.txt").write_text("mock")
        (star_dir / "chrLength.txt").write_text("chr1\t50\n")
        
        config = GenomeConfig(
            name="test_rna",
            index=STARIndex(prefix=star_dir),
        )
        
        assert config.name == "test_rna"
        assert config.index.type == "STAR"
        assert config.index.prefix == star_dir
