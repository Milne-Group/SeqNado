"""Tests for seqnado.config.third_party_tools module."""

from seqnado import Assay
from seqnado.config.third_party_tools import Bamnado, ThirdPartyToolsConfig


class TestBamnado:
    """Tests for the Bamnado tool config class."""

    def test_default_bam_normalize_config(self):
        """bam_normalize has a sensible default bin size for background-window scaling."""
        config = Bamnado()
        assert "--bin-size 10000" in str(config.bam_normalize.command_line_arguments)

    def test_bam_coverage_and_bam_normalize_binsize_are_independent(self):
        """bam_coverage (track resolution) and bam_normalize (background window) bin
        sizes are different concepts and must not be conflated."""
        config = Bamnado()
        assert "--bin-size 10" in str(config.bam_coverage.command_line_arguments)
        assert "--bin-size 10000" in str(config.bam_normalize.command_line_arguments)


class TestThirdPartyToolsConfigBinsizeInheritance:
    """Tests for BigwigConfig.binsize propagating into tool CLI options via for_assay()."""

    def test_binsize_propagates_to_deeptools_and_bamnado_bam_coverage(self):
        config = ThirdPartyToolsConfig.for_assay(Assay.CHIP, binsize=50)
        assert "--binSize 50" in str(config.deeptools.bam_coverage.command_line_arguments)
        assert "--bin-size 50" in str(config.bamnado.bam_coverage.command_line_arguments)

    def test_binsize_does_not_propagate_to_bam_normalize(self):
        """bam_normalize's bin size is a separate concept (background window size)
        and must not be overwritten by the bigwig track-resolution binsize."""
        config = ThirdPartyToolsConfig.for_assay(Assay.CHIP, binsize=50)
        assert "--bin-size 10000" in str(config.bamnado.bam_normalize.command_line_arguments)

    def test_no_binsize_keeps_defaults(self):
        config = ThirdPartyToolsConfig.for_assay(Assay.CHIP)
        assert "--binSize 10" in str(config.deeptools.bam_coverage.command_line_arguments)
        assert "--bin-size 10" in str(config.bamnado.bam_coverage.command_line_arguments)
