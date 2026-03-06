"""Tests for configuration validation in seqnado.config.configs module."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from seqnado import Assay, GenomicCoordinate
from seqnado.config.configs import (
    UCSCHubConfig,
    PCRDuplicatesConfig,
    PCRDuplicateHandling,
    PCRDuplicateTool,
    ProjectConfig,
    QCConfig,
    BigwigConfig,
    PlottingConfig,
    PileupMethod,
)


class TestUCSCHubConfig:
    """Tests for UCSCHubConfig validation and class methods."""

    def test_default_ucsc_hub_config(self):
        """Test default UCSC hub config creation."""
        config = UCSCHubConfig()

        assert config.directory == "seqnado_output/hub/"
        assert config.name == "seqnado_hub"
        assert config.genome == "hg38"
        assert config.email == "test@example.com"
        assert config.color_by == ['samplename']
        assert config.subgroup_by == ['method', 'norm']

    def test_ucsc_hub_config_for_rna_assay(self):
        """Test UCSC hub config for RNA assay."""
        config = UCSCHubConfig.for_assay(Assay.RNA)

        assert config.subgroup_by == ["method", "norm", "strand"]
        assert config.overlay_by == ["samplename", "method", "norm"]

    def test_ucsc_hub_config_for_mcc_assay(self):
        """Test UCSC hub config for MCC assay."""
        config = UCSCHubConfig.for_assay(Assay.MCC)

        assert config.color_by == ['viewpoint', 'samplename']
        assert config.subgroup_by == ['viewpoint']

    def test_ucsc_hub_config_for_other_assays(self):
        """Test UCSC hub config for other assays (ATAC, ChIP, etc.)."""
        # CHIP and CAT have antibody in color_by and subgroup_by, others don't
        for assay in [Assay.ATAC, Assay.SNP]:
            config = UCSCHubConfig.for_assay(assay)
            assert config.color_by == ['samplename']
            assert config.subgroup_by == ['method', 'norm']
        
        # CHIP and CAT include antibody in color_by and subgroup_by
        for assay in [Assay.CHIP, Assay.CAT]:
            config = UCSCHubConfig.for_assay(assay)
            assert config.color_by == ['antibody', 'samplename']
            assert config.subgroup_by == ['method', 'norm', 'antibody']

    def test_ucsc_hub_config_invalid_name(self):
        """Test that invalid hub name raises error."""
        with pytest.raises(ValidationError, match="alphanumeric"):
            UCSCHubConfig(name="invalid name with spaces!")

    def test_ucsc_hub_config_empty_name(self):
        """Test that empty hub name raises error."""
        with pytest.raises(ValidationError, match="must not be empty"):
            UCSCHubConfig(name="")

    def test_ucsc_hub_config_overlay_by_none_string(self):
        """Test that 'None' string is coerced to None for overlay_by."""
        config = UCSCHubConfig(overlay_by="None")
        assert config.overlay_by is None

    def test_ucsc_hub_config_overlay_by_comma_separated(self):
        """Test that comma-separated string is parsed into list."""
        config = UCSCHubConfig(overlay_by="samplename,method,norm")
        assert config.overlay_by == ["samplename", "method", "norm"]

    def test_ucsc_hub_config_overlay_by_single_value(self):
        """Test that single string becomes single-element list."""
        config = UCSCHubConfig(overlay_by="samplename")
        assert config.overlay_by == ["samplename"]

    def test_ucsc_hub_config_overlay_by_list(self):
        """Test that list is preserved as list."""
        config = UCSCHubConfig(overlay_by=["samplename", "method"])
        assert config.overlay_by == ["samplename", "method"]

    def test_ucsc_hub_config_supergroup_by_none(self):
        """Test that None is preserved for supergroup_by."""
        config = UCSCHubConfig(supergroup_by=None)
        assert config.supergroup_by is None

    def test_ucsc_hub_config_supergroup_by_string(self):
        """Test that string is coerced to list for supergroup_by."""
        config = UCSCHubConfig(supergroup_by="treatment")
        assert config.supergroup_by == ["treatment"]

    def test_ucsc_hub_config_directory_validation_absolute_path(self, tmp_path):
        """Test directory validation for absolute paths."""
        # Parent directory must exist for absolute paths
        hub_dir = tmp_path / "hub"
        config = UCSCHubConfig(directory=str(hub_dir))
        assert config.directory == str(hub_dir)

    def test_ucsc_hub_config_directory_validation_relative_path(self):
        """Test directory validation skips validation for relative paths."""
        # Relative paths are allowed (used during config generation)
        config = UCSCHubConfig(directory="seqnado_output/hub/")
        assert config.directory == "seqnado_output/hub/"

    def test_ucsc_hub_config_genomic_coordinate_default(self):
        """Test default genomic position."""
        config = UCSCHubConfig()
        assert isinstance(config.default_position, GenomicCoordinate)
        assert config.default_position.chromosome == "chr1"
        assert config.default_position.start == 1
        assert config.default_position.end == 1000000


class TestPCRDuplicatesConfig:
    """Tests for PCRDuplicatesConfig."""

    def test_default_pcr_duplicates_config(self):
        """Test default PCR duplicates config."""
        config = PCRDuplicatesConfig()
        assert config.strategy == PCRDuplicateHandling.NONE
        assert config.tool is None

    def test_pcr_duplicates_config_with_tool(self):
        """Test PCR duplicates config with tool specified."""
        config = PCRDuplicatesConfig(
            strategy=PCRDuplicateHandling.REMOVE,
            tool=PCRDuplicateTool.PICARD
        )
        assert config.strategy == PCRDuplicateHandling.REMOVE
        assert config.tool == PCRDuplicateTool.PICARD

    def test_pcr_duplicates_config_mark_strategy(self):
        """Test PCR duplicates config with MARK strategy."""
        config = PCRDuplicatesConfig(
            strategy=PCRDuplicateHandling.MARK,
            tool=PCRDuplicateTool.SAMTOOLS
        )
        assert config.strategy == PCRDuplicateHandling.MARK
        assert config.tool == PCRDuplicateTool.SAMTOOLS


class TestProjectConfig:
    """Tests for ProjectConfig."""

    def test_project_config_minimal(self):
        """Test minimal project config."""
        config = ProjectConfig(name="test_project")
        assert config.name == "test_project"
        assert config.description is None
        # Date should be set to today by default
        from datetime import date
        assert config.date == date.today()

    def test_project_config_with_description(self):
        """Test project config with description."""
        config = ProjectConfig(
            name="my_project",
            description="A test project"
        )
        assert config.name == "my_project"
        assert config.description == "A test project"

    def test_project_config_with_custom_date(self):
        """Test project config with custom date."""
        from datetime import date
        custom_date = date(2024, 1, 1)
        config = ProjectConfig(name="test", date=custom_date)
        assert config.date == custom_date


class TestQCConfig:
    """Tests for QCConfig."""

    def test_default_qc_config(self):
        """Test default QC config."""
        config = QCConfig()
        assert config.run_fastq_screen is True
        assert config.calculate_library_complexity is False
        assert config.calculate_fraction_of_reads_in_peaks is False

    def test_qc_config_all_enabled(self):
        """Test QC config with all features enabled."""
        config = QCConfig(
            run_fastq_screen=True,
            calculate_library_complexity=True,
            calculate_fraction_of_reads_in_peaks=True
        )
        assert config.run_fastq_screen is True
        assert config.calculate_library_complexity is True
        assert config.calculate_fraction_of_reads_in_peaks is True


class TestBigwigConfig:
    """Tests for BigwigConfig."""

    def test_default_bigwig_config(self):
        """Test default bigwig config."""
        config = BigwigConfig()
        assert config.pileup_method is None
        assert config.binsize is None

    def test_bigwig_config_with_methods(self):
        """Test bigwig config with pileup methods."""
        config = BigwigConfig(
            pileup_method=[PileupMethod.DEEPTOOLS, PileupMethod.HOMER],
            binsize=50
        )
        assert PileupMethod.DEEPTOOLS in config.pileup_method
        assert PileupMethod.HOMER in config.pileup_method
        assert config.binsize == 50


class TestPlottingConfig:
    """Tests for PlottingConfig."""

    def test_default_plotting_config(self):
        """Test default plotting config."""
        config = PlottingConfig()
        assert config.coordinates is None
        assert config.file_format == "pdf"

    def test_plotting_config_png_format(self):
        """Test plotting config with PNG format."""
        config = PlottingConfig(file_format="png")
        assert config.file_format == "png"

    def test_plotting_config_with_coordinates(self):
        """Test plotting config with coordinates."""
        config = PlottingConfig(coordinates="chr1:1000-2000")
        assert config.coordinates == "chr1:1000-2000"


