"""Unit tests for seqnado.config.multiomics module."""

from seqnado.config.multiomics import MultiomicsConfig


class TestMultiomicsConfig:
    """Tests for the MultiomicsConfig model."""

    def test_default_values(self):
        """Test MultiomicsConfig with default values."""
        config = MultiomicsConfig()

        assert config.assays == []
        assert config.assay_configs == {}
        assert config.create_heatmaps is True
        assert config.create_dataset is True
        assert config.create_summary is True

    def test_custom_assays(self):
        """Test MultiomicsConfig with custom assays list."""
        config = MultiomicsConfig(assays=["atac", "chip", "rna"])
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

    def test_assay_configs_dict(self):
        """Test MultiomicsConfig with assay configs."""
        assay_configs = {
            "atac": {"assay": "atac", "project": {"name": "test"}},
            "chip": {"assay": "chip", "project": {"name": "test"}},
        }
        config = MultiomicsConfig(assay_configs=assay_configs)
        assert config.assay_configs == assay_configs

    def test_model_serialization(self):
        """Test that model can be serialized to dict."""
        config = MultiomicsConfig(assays=["atac", "rna"])

        data = config.model_dump()

        assert data["assays"] == ["atac", "rna"]

    def test_model_from_dict(self):
        """Test that model can be created from dict."""
        data = {
            "assays": ["chip"],
            "create_heatmaps": False,
            "create_dataset": False,
        }
        config = MultiomicsConfig(**data)

        assert config.assays == ["chip"]
        assert config.create_heatmaps is False
