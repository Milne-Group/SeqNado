"""Unit tests for seqnado.config.multiomics module."""

import pytest
from pathlib import Path

from seqnado.config.multiomics import MultiomicsConfig, none_str_to_none


class TestNoneStrToNone:
    """Tests for the none_str_to_none validator function."""

    def test_converts_string_none_to_none(self):
        """Test that string 'none' is converted to None."""
        assert none_str_to_none("none") is None

    def test_converts_uppercase_none_to_none(self):
        """Test that uppercase 'NONE' is converted to None."""
        assert none_str_to_none("NONE") is None

    def test_converts_mixed_case_none_to_none(self):
        """Test that mixed case 'NoNe' is converted to None."""
        assert none_str_to_none("NoNe") is None

    def test_converts_none_with_whitespace_to_none(self):
        """Test that 'none' with whitespace is converted to None."""
        assert none_str_to_none("  none  ") is None
        assert none_str_to_none("\tnone\n") is None

    def test_returns_non_none_strings_unchanged(self):
        """Test that non-'none' strings are returned unchanged."""
        assert none_str_to_none("something") == "something"
        assert none_str_to_none("") == ""
        assert none_str_to_none("not_none") == "not_none"

    def test_returns_non_string_values_unchanged(self):
        """Test that non-string values are returned unchanged."""
        assert none_str_to_none(None) is None
        assert none_str_to_none(42) == 42
        assert none_str_to_none([1, 2, 3]) == [1, 2, 3]
        assert none_str_to_none({"key": "value"}) == {"key": "value"}


class TestMultiomicsConfig:
    """Tests for the MultiomicsConfig model."""

    def test_default_values(self):
        """Test MultiomicsConfig with default values."""
        # create_dataset=False to avoid requiring fasta_index
        config = MultiomicsConfig(create_dataset=False)

        assert config.output_dir == "seqnado_output/"
        assert config.binsize == 1000
        assert config.assays == []
        assert config.assay_configs == {}
        assert config.create_heatmaps is True
        assert config.create_dataset is False
        assert config.create_summary is True
        assert config.regions_bed is None

    def test_custom_output_dir(self):
        """Test MultiomicsConfig with custom output directory."""
        config = MultiomicsConfig(output_dir="custom_output/", create_dataset=False)
        assert config.output_dir == "custom_output/"

    def test_custom_binsize(self):
        """Test MultiomicsConfig with custom binsize."""
        config = MultiomicsConfig(binsize=500, create_dataset=False)
        assert config.binsize == 500

    def test_custom_assays(self):
        """Test MultiomicsConfig with custom assays list."""
        config = MultiomicsConfig(assays=["atac", "chip", "rna"], create_dataset=False)
        assert config.assays == ["atac", "chip", "rna"]

    def test_custom_flags(self):
        """Test MultiomicsConfig with custom analysis flags."""
        config = MultiomicsConfig(
            create_heatmaps=False,
            create_dataset=False,
            create_summary=False,
        )
        assert config.create_heatmaps is False
        assert config.create_dataset is False
        assert config.create_summary is False

    def test_regions_bed_path(self):
        """Test MultiomicsConfig with regions BED file."""
        bed_path = Path("/path/to/regions.bed")
        config = MultiomicsConfig(regions_bed=bed_path)
        assert config.regions_bed == bed_path

    def test_assay_configs_dict(self):
        """Test MultiomicsConfig with assay configs."""
        assay_configs = {
            "atac": {"assay": "atac", "project": {"name": "test"}},
            "chip": {"assay": "chip", "project": {"name": "test"}},
        }
        config = MultiomicsConfig(assay_configs=assay_configs, create_dataset=False)
        assert config.assay_configs == assay_configs

    def test_model_serialization(self):
        """Test that model can be serialized to dict."""
        config = MultiomicsConfig(
            output_dir="test_output/",
            binsize=200,
            assays=["atac", "rna"],
            create_dataset=False,
        )

        data = config.model_dump()

        assert "output_dir" in data
        assert data["output_dir"] == "test_output/"
        assert data["binsize"] == 200
        assert data["assays"] == ["atac", "rna"]

    def test_model_from_dict(self):
        """Test that model can be created from dict."""
        data = {
            "output_dir": "custom/",
            "binsize": 200,
            "assays": ["chip"],
            "create_heatmaps": False,
            "create_dataset": False,
        }
        config = MultiomicsConfig(**data)

        assert config.output_dir == "custom/"
        assert config.binsize == 200
        assert config.assays == ["chip"]
        assert config.create_heatmaps is False
