"""Unit tests for genome build and configuration."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from seqnado import Assay
from seqnado.config.configs import BowtieIndex, GenomeConfig, STARIndex
from seqnado.config.user_input import load_genome_configs


class TestGenomeConfigCreation:
    """Test suite for creating GenomeConfig objects."""

    def test_genome_config_with_valid_bowtie_index(self, tmp_path: Path) -> None:
        """Test creating GenomeConfig with valid Bowtie2 index."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        # Create minimal Bowtie2 index files
        prefix = str(bt2_dir / "test")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"test{suffix}").touch()
        
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGT\n")
        
        config = GenomeConfig(
            name="test_genome",
            index=BowtieIndex(prefix=prefix),
            fasta=str(fasta),
        )
        
        assert config.name == "test_genome"
        assert config.index.type == "Bowtie2"
        assert config.index.prefix == prefix
        assert config.fasta is not None

    def test_genome_config_with_optional_fields(self, tmp_path: Path) -> None:
        """Test creating GenomeConfig with optional fields."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "test")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"test{suffix}").touch()
        
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGT\n")
        
        gtf = tmp_path / "test.gtf"
        gtf.write_text("chr1\t.\tgene\t1\t4\t.\t+\t.\tgene_id \"test\"\n")
        
        chrom_sizes = tmp_path / "test.chrom.sizes"
        chrom_sizes.write_text("chr1\t4\n")
        
        config = GenomeConfig(
            name="test_genome",
            index=BowtieIndex(prefix=prefix),
            fasta=str(fasta),
            gtf=str(gtf),
            chromosome_sizes=str(chrom_sizes),
            organism="Homo sapiens",
            version="hg38",
            bin_size=5000,
        )
        
        assert config.organism == "Homo sapiens"
        assert config.version == "hg38"
        assert config.bin_size == 5000
        assert config.gtf is not None

    def test_genome_config_with_star_index(self, tmp_path: Path) -> None:
        """Test creating GenomeConfig with STAR index."""
        star_dir = tmp_path / "STAR"
        star_dir.mkdir()
        
        # Create minimal STAR index files
        (star_dir / "SA").write_text("mock")
        (star_dir / "chrLength.txt").write_text("chr1\t50\n")
        
        config = GenomeConfig(
            name="test_rna",
            index=STARIndex(prefix=star_dir),
        )
        
        assert config.name == "test_rna"
        assert config.index.type == "STAR"
        assert config.index.prefix == star_dir

    def test_genome_config_none_string_conversion(self, tmp_path: Path) -> None:
        """Test that 'none' strings are converted to None."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "test")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"test{suffix}").touch()
        
        config = GenomeConfig(
            name="test",
            index=BowtieIndex(prefix=prefix),
            gtf="none",
        )
        
        assert config.gtf is None

    def test_genome_config_default_bin_size(self, tmp_path: Path) -> None:
        """Test that bin_size defaults to 1000."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "test")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"test{suffix}").touch()
        
        config = GenomeConfig(
            name="test",
            index=BowtieIndex(prefix=prefix),
        )
        
        assert config.bin_size == 1000


class TestBowtieIndexValidation:
    """Test suite for Bowtie2 index validation."""

    def test_valid_bowtie_index_passes(self, tmp_path: Path) -> None:
        """Test that valid Bowtie2 index passes validation."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "genome")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"genome{suffix}").touch()
        
        index = BowtieIndex(prefix=prefix)
        assert index.prefix == prefix
        assert len(index.files) > 0

    def test_bowtie_index_with_missing_files_raises_error(self, tmp_path: Path) -> None:
        """Test that Bowtie2 index with missing files raises error."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "missing")
        
        with pytest.raises(ValueError, match="No Bowtie2 index files found"):
            BowtieIndex(prefix=prefix)

    def test_bowtie_index_files_property(self, tmp_path: Path) -> None:
        """Test that files property returns all index files."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "genome")
        suffixes = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]
        for suffix in suffixes:
            (bt2_dir / f"genome{suffix}").touch()
        
        index = BowtieIndex(prefix=prefix)
        files = index.files
        
        assert len(files) == 6
        assert all(f.exists() for f in files)

    def test_bowtie_index_accepts_bt2l_format(self, tmp_path: Path) -> None:
        """Test that Bowtie2 index accepts .bt2l format (large genome)."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "large_genome")
        for suffix in [".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"]:
            (bt2_dir / f"large_genome{suffix}").touch()
        
        index = BowtieIndex(prefix=prefix)
        assert len(index.files) == 6

    def test_bowtie_index_with_nonexistent_directory(self, tmp_path: Path) -> None:
        """Test that Bowtie2 index with nonexistent directory raises error."""
        prefix = str(tmp_path / "nonexistent" / "genome")
        
        with pytest.raises(ValueError, match="does not exist"):
            BowtieIndex(prefix=prefix)

    def test_bowtie_index_none_prefix_is_valid(self) -> None:
        """Test that None prefix is valid for Bowtie2 index."""
        index = BowtieIndex(prefix=None)
        assert index.prefix is None
        assert index.files == []


class TestSTARIndexValidation:
    """Test suite for STAR index validation."""

    def test_valid_star_index_passes(self, tmp_path: Path) -> None:
        """Test that valid STAR index passes validation."""
        star_dir = tmp_path / "STAR"
        star_dir.mkdir()
        
        # Create minimal STAR files
        (star_dir / "SA").write_text("mock")
        (star_dir / "chrLength.txt").write_text("chr1\t50\n")
        
        index = STARIndex(prefix=star_dir)
        assert index.prefix == star_dir
        assert len(index.files) > 0

    def test_star_index_with_nonexistent_directory(self) -> None:
        """Test that STAR index with nonexistent directory raises error."""
        from pathlib import Path
        
        nonexistent = Path("/nonexistent/path/to/star")
        
        with pytest.raises(ValueError, match="does not exist or is not a directory"):
            STARIndex(prefix=nonexistent)

    def test_star_index_files_property(self, tmp_path: Path) -> None:
        """Test that files property returns all index files."""
        star_dir = tmp_path / "STAR"
        star_dir.mkdir()
        
        (star_dir / "SA").write_text("mock")
        (star_dir / "chrLength.txt").write_text("chr1\t50\n")
        (star_dir / "genomeParameters.txt").write_text("genome params")
        
        index = STARIndex(prefix=star_dir)
        files = index.files
        
        assert len(files) >= 3

    def test_star_index_none_prefix_is_valid(self) -> None:
        """Test that None prefix is valid for STAR index."""
        index = STARIndex(prefix=None)
        assert index.prefix is None
        assert index.files == []


class TestLoadGenomeConfigs:
    """Test suite for loading genome configurations."""

    def test_load_genome_configs_for_rna(self, tmp_path: Path, monkeypatch) -> None:
        """Test loading genome configs for RNA assay (STAR index)."""
        # Setup config
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True)
        cfg_file = cfg_dir / "genome_config.json"
        
        # Create STAR index
        star_dir = tmp_path / "STAR"
        star_dir.mkdir()
        (star_dir / "SA").write_text("mock")
        
        genome_data = {
            "hg38": {
                "star_index": str(star_dir),
            }
        }
        
        cfg_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
        
        result = load_genome_configs(Assay.RNA)
        
        assert "hg38" in result
        assert isinstance(result["hg38"], GenomeConfig)
        assert isinstance(result["hg38"].index, STARIndex)

    def test_load_genome_configs_for_chip(self, tmp_path: Path, monkeypatch) -> None:
        """Test loading genome configs for ChIP assay (Bowtie index)."""
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True)
        cfg_file = cfg_dir / "genome_config.json"
        
        # Create Bowtie2 index
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        prefix = str(bt2_dir / "hg38")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"hg38{suffix}").touch()
        
        genome_data = {
            "hg38": {
                "bt2_index": prefix,
            }
        }
        
        cfg_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
        
        result = load_genome_configs(Assay.CHIP)
        
        assert "hg38" in result
        assert isinstance(result["hg38"].index, BowtieIndex)

    def test_load_genome_configs_multiple_genomes(self, tmp_path: Path, monkeypatch) -> None:
        """Test loading multiple genome configurations."""
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True)
        cfg_file = cfg_dir / "genome_config.json"
        
        # Create indices
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        genome_data = {}
        for genome_name in ["hg38", "mm10", "dm6"]:
            bt2_prefix = str(bt2_dir / genome_name)
            for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
                (bt2_dir / f"{genome_name}{suffix}").touch()
            genome_data[genome_name] = {"bt2_index": bt2_prefix}
        
        cfg_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
        
        result = load_genome_configs(Assay.CHIP)
        
        assert len(result) == 3
        assert all(name in result for name in ["hg38", "mm10", "dm6"])

    def test_load_genome_configs_sets_genome_name(self, tmp_path: Path, monkeypatch) -> None:
        """Test that genome name is correctly set from config key."""
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True)
        cfg_file = cfg_dir / "genome_config.json"
        
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        prefix = str(bt2_dir / "custom")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"custom{suffix}").touch()
        
        genome_data = {
            "my_custom_genome": {"bt2_index": prefix}
        }
        
        cfg_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
        
        result = load_genome_configs(Assay.CHIP)
        
        assert result["my_custom_genome"].name == "my_custom_genome"

    def test_load_genome_configs_with_missing_config_file(self, tmp_path: Path, monkeypatch) -> None:
        """Test error handling when config file doesn't exist."""
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
        
        with pytest.raises(SystemExit):
            load_genome_configs(Assay.CHIP)

    def test_load_genome_configs_with_invalid_assay(self, tmp_path: Path, monkeypatch) -> None:
        """Test that invalid assay is handled gracefully."""
        cfg_dir = tmp_path / ".config" / "seqnado"
        cfg_dir.mkdir(parents=True)
        cfg_file = cfg_dir / "genome_config.json"
        cfg_file.write_text(json.dumps({}))
        
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))
        
        # Should handle gracefully - either return empty dict or raise error
        try:
            result = load_genome_configs("INVALID_ASSAY")  # type: ignore
            assert isinstance(result, dict)
        except (ValueError, AttributeError, SystemExit):
            # Expected for invalid assay
            pass


class TestGenomeConfigSerializationAndDeserialization:
    """Test suite for genome config JSON serialization."""

    def test_genome_config_json_round_trip(self, tmp_path: Path) -> None:
        """Test that genome config can be serialized and deserialized."""
        bt2_dir = tmp_path / "bt2"
        bt2_dir.mkdir()
        
        prefix = str(bt2_dir / "test")
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"test{suffix}").touch()
        
        fasta = tmp_path / "test.fa"
        fasta.write_text(">chr1\nACGT\n")
        
        # Create original config
        original = GenomeConfig(
            name="test",
            index=BowtieIndex(prefix=prefix),
            fasta=str(fasta),
            organism="Homo sapiens",
            version="test_v1",
        )
        
        # Convert to dict and back
        config_dict = original.model_dump()
        reconstructed = GenomeConfig(**config_dict)
        
        assert reconstructed.name == original.name
        assert reconstructed.organism == original.organism
        assert reconstructed.version == original.version
