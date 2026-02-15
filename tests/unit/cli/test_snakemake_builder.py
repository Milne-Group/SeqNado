"""Unit tests for Snakemake command builder."""

import pytest
from pathlib import Path
from seqnado.cli.snakemake_builder import SnakemakeCommandBuilder


class TestSnakemakeCommandBuilder:
    """Tests for SnakemakeCommandBuilder class."""

    def test_builder_basic_command(self, tmp_path):
        """Test building a basic Snakemake command."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile, cores=4)
        cmd = builder.build()

        assert "snakemake" in cmd
        assert "--snakefile" in cmd
        assert str(snakefile) in cmd
        assert "-c" in cmd
        assert "4" in cmd

    def test_builder_add_configfile(self, tmp_path):
        """Test adding configfile to command."""
        snakefile = tmp_path / "Snakefile"
        config = tmp_path / "config.yaml"
        snakefile.touch()
        config.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_configfile(config).build()

        assert "--configfile" in cmd
        assert str(config) in cmd

    def test_builder_add_config(self, tmp_path):
        """Test adding config key-value pairs."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_config(genome="hg38", output_dir="/tmp").build()

        assert "--config" in cmd
        assert "genome=hg38" in cmd
        assert "output_dir=/tmp" in cmd

    def test_builder_add_dry_run(self, tmp_path):
        """Test adding --dry-run flag."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_dry_run().build()

        assert "--dry-run" in cmd

    def test_builder_add_unlock(self, tmp_path):
        """Test adding --unlock flag."""
        snakefile = tmp_path / "Snakefile"
        snakefile.touch()

        builder = SnakemakeCommandBuilder(snakefile)
        cmd = builder.add_unlock().build()

        assert "--unlock" in cmd

    def test_builder_chaining(self, tmp_path):
        """Test fluent API chaining."""
        snakefile = tmp_path / "Snakefile"
        config = tmp_path / "config.yaml"
        snakefile.touch()
        config.touch()

        cmd = (
            SnakemakeCommandBuilder(snakefile, cores=8)
            .add_configfile(config)
            .add_dry_run()
            .add_config(test="value")
            .build()
        )

        assert len(cmd) > 0
        assert "snakemake" in cmd
        assert "--dry-run" in cmd


# TODO: Add tests for profile resolution and pass-through args
