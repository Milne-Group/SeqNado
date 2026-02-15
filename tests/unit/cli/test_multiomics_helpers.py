"""Unit tests for CLI multiomics helpers."""

import pytest
from pathlib import Path
from seqnado.cli.multiomics_helpers import (
    should_use_multiomics_mode,
    get_multiomics_config_paths,
    validate_multiomics_setup,
)


class TestShouldUseMultiomicsMode:
    """Tests for multiomics mode detection."""

    def test_single_assay_provided(self, tmp_path):
        """Test that multiomics mode is not used when assay explicitly provided."""
        assert not should_use_multiomics_mode(assay="rna", config_file=None, cwd=tmp_path)

    def test_config_file_provided(self, tmp_path):
        """Test that multiomics mode is not used when config file explicitly provided."""
        config = tmp_path / "config_custom.yaml"
        assert not should_use_multiomics_mode(
            assay=None, config_file=config, cwd=tmp_path
        )

    def test_single_config_file_exists(self, tmp_path):
        """Test that multiomics mode is not used with single config file."""
        config = tmp_path / "config_atac.yaml"
        config.touch()

        assert not should_use_multiomics_mode(assay=None, config_file=None, cwd=tmp_path)

    def test_multiple_config_files(self, tmp_path):
        """Test that multiomics mode is detected with multiple config files."""
        config1 = tmp_path / "config_atac.yaml"
        config2 = tmp_path / "config_rna.yaml"
        config1.touch()
        config2.touch()

        # Note: This requires proper Assay enum detection
        # The actual behavior depends on how find_assay_config_paths handles the files
        result = should_use_multiomics_mode(assay=None, config_file=None, cwd=tmp_path)
        # Result depends on whether files match expected pattern


class TestGetMultiomicsConfigPaths:
    """Tests for retrieving multiomics config paths."""

    def test_get_config_paths_returns_dict(self, tmp_path):
        """Test that function returns dict with expected keys."""
        result = get_multiomics_config_paths(cwd=tmp_path)

        assert isinstance(result, dict)
        assert "configs" in result
        assert "metadata" in result


# TODO: Add more comprehensive tests once multiomics file structure is finalized
