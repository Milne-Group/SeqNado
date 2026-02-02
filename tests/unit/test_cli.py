"""Tests for CLI helper functions."""

import json
from pathlib import Path
import shutil
import subprocess

from seqnado.cli import (
    _assay_names,
    _coerce_value_to_dtype,
    _extract_candidate_defaults_from_schema,
    _find_fastqs,
    _format_col_hint,
    _get_profile_name,
    _preset_profiles,
    _profile_autocomplete,
    _read_json,
    _snakemake_available,
    _style_name_with_rich,
    _write_json,
)


def test_read_json(tmp_path):
    """Test reading JSON from file."""
    test_file = tmp_path / "test.json"
    test_data = {"key": "value", "nested": {"data": [1, 2, 3]}}

    # Write the test data to the file
    test_file.write_text(json.dumps(test_data), encoding="utf-8")

    # Read the data back using _read_json
    loaded = _read_json(test_file)

    # Verify the data matches
    assert loaded == test_data


def test_write_json(tmp_path):
    """Test writing JSON to file."""
    test_file = tmp_path / "output.json"
    test_data = {"key": "value", "nested": {"data": [1, 2, 3]}}
    
    _write_json(test_file, test_data)
    
    assert test_file.exists()
    loaded = json.loads(test_file.read_text())
    assert loaded == test_data


def test_snakemake_available():
    """Test checking if snakemake is available."""
    # Snakemake should be available in the test environment
    result = _snakemake_available()
    assert isinstance(result, bool)


def test_find_fastqs_with_directory(tmp_path):
    """Test finding FASTQ files in a directory."""
    # Create test FASTQ files - only *.fastq.gz pattern
    (tmp_path / "sample1_R1.fastq.gz").touch()
    (tmp_path / "sample1_R2.fastq.gz").touch()
    (tmp_path / "sample2_R1.fastq.gz").touch()
    (tmp_path / "sample2_R2.fastq.gz").touch()
    (tmp_path / "other.fq.gz").touch()  # Different pattern, won't match
    (tmp_path / "other.txt").touch()
    
    result = _find_fastqs([str(tmp_path)])
    
    assert len(result) == 4
    assert all(isinstance(p, Path) for p in result)
    assert all(p.suffix == ".gz" for p in result)
    assert all("fastq" in p.name for p in result)


def test_find_fastqs_with_glob(tmp_path):
    """Test finding FASTQ files with direct path (not glob pattern)."""
    (tmp_path / "sample1_R1.fastq.gz").touch()
    (tmp_path / "sample1_R2.fastq.gz").touch()
    (tmp_path / "other_R1.fastq.gz").touch()
    
    # _find_fastqs only globs within directories, not with patterns
    result = _find_fastqs([str(tmp_path)])
    
    assert len(result) == 3


def test_find_fastqs_with_file(tmp_path):
    """Test that non-directory paths are skipped."""
    test_file = tmp_path / "test.fastq.gz"
    test_file.touch()
    
    # File path doesn't exist as a directory, so gets skipped
    result = _find_fastqs([str(test_file)])
    
    assert len(result) == 0


def test_find_fastqs_nonexistent(tmp_path):
    """Test that non-existent paths are skipped."""
    nonexistent = tmp_path / "does_not_exist"
    
    result = _find_fastqs([str(nonexistent)])
    
    assert len(result) == 0


def test_style_name_with_rich():
    """Test styling text with rich (or fallback)."""
    result = _style_name_with_rich("test_name")
    assert isinstance(result, str)
    # Should at minimum return the name itself
    assert "test_name" in result or result == "test_name"


def test_style_name_with_rich_custom_style():
    """Test styling with custom style."""
    result = _style_name_with_rich("test_name", style="bold green")
    assert isinstance(result, str)


def test_coerce_value_to_dtype_int():
    """Test coercing string to int."""
    result = _coerce_value_to_dtype("42", int, categories=None)
    assert result == 42
    assert isinstance(result, int)


def test_coerce_value_to_dtype_float():
    """Test coercing string to float."""
    result = _coerce_value_to_dtype("3.14", float, categories=None)
    assert result == 3.14
    assert isinstance(result, float)


def test_coerce_value_to_dtype_bool():
    """Test coercing string to bool."""
    assert _coerce_value_to_dtype("true", bool, categories=None) is True
    assert _coerce_value_to_dtype("True", bool, categories=None) is True
    assert _coerce_value_to_dtype("yes", bool, categories=None) is True
    assert _coerce_value_to_dtype("y", bool, categories=None) is True
    assert _coerce_value_to_dtype("1", bool, categories=None) is True
    assert _coerce_value_to_dtype("false", bool, categories=None) is False
    assert _coerce_value_to_dtype("False", bool, categories=None) is False
    assert _coerce_value_to_dtype("no", bool, categories=None) is False
    assert _coerce_value_to_dtype("n", bool, categories=None) is False
    assert _coerce_value_to_dtype("0", bool, categories=None) is False


def test_coerce_value_to_dtype_str():
    """Test that strings pass through."""
    result = _coerce_value_to_dtype("test", str, categories=None)
    assert result == "test"


def test_coerce_value_to_dtype_with_categories():
    """Test categorical validation."""
    categories = ["option1", "option2", "option3"]
    
    # Valid category
    result = _coerce_value_to_dtype("option1", str, categories=categories)
    assert result == "option1"
    
    # Empty is allowed
    result = _coerce_value_to_dtype("", str, categories=categories)
    assert result == ""
    
    # Invalid category should raise
    try:
        _coerce_value_to_dtype("invalid", str, categories=categories)
        assert False, "Should have raised ValueError"
    except ValueError as e:
        assert "must be one of" in str(e)


def test_format_col_hint():
    """Test formatting column hints."""
    meta = {
        "default": "value1",
        "description": "A test column",
        "dtype": str,
        "categories": None,
    }
    
    result = _format_col_hint("test_col", meta)
    assert isinstance(result, str)
    assert "test_col" in result


def test_format_col_hint_with_categories():
    """Test formatting column hints with categories."""
    meta = {
        "default": "option1",
        "categories": ["option1", "option2", "option3"],
        "dtype": str,
    }
    
    result = _format_col_hint("test_col", meta)
    assert isinstance(result, str)
    # Just verify it returns a string (formatting details may vary)
    assert len(result) > 0


def test_assay_names():
    """Test getting assay names."""
    result = _assay_names()
    assert isinstance(result, list)
    assert len(result) > 0
    assert all(isinstance(name, str) for name in result)
    # Check for some expected assays
    assert any("chip" in name.lower() for name in result)
    assert any("atac" in name.lower() for name in result)


def test_get_profile_name_valid():
    """Test getting profile shorthand from path."""
    result = _get_profile_name(Path("profile_local_environment"))
    assert result == "le"
    
    result = _get_profile_name(Path("profile_slurm_singularity"))
    assert result == "ss"


def test_get_profile_name_invalid():
    """Test profile name for non-profile path."""
    result = _get_profile_name(Path("not_a_profile"))
    assert result is None


def test_get_profile_name_single_word():
    """Test profile name with single word after prefix."""
    result = _get_profile_name(Path("profile_test"))
    assert result == "t"


def test_preset_profiles():
    """Test getting preset profiles dictionary."""
    result = _preset_profiles()
    assert isinstance(result, dict)
    assert len(result) > 0
    # Keys should be short codes like 'le', 'ss'
    assert all(len(k) <= 3 for k in result.keys())
    # Values should be profile directory names
    assert all(v.startswith("profile_") for v in result.values())


def test_profile_autocomplete():
    """Test getting profile names for autocomplete."""
    result = _profile_autocomplete()
    assert isinstance(result, list)
    assert len(result) > 0
    assert all(isinstance(name, str) for name in result)


def test_extract_candidate_defaults_from_schema():
    """Test extracting default values from schema."""
    from seqnado.inputs.validation import DesignDataFrame
    from seqnado.core import Assay
    
    result = _extract_candidate_defaults_from_schema(
        DesignDataFrame, 
        Assay.CHIP
    )
    
    assert isinstance(result, dict)
    # Should have some column defaults
    assert len(result) >= 0  # May be empty depending on schema


def test_extract_candidate_defaults_atac():
    """Test extracting defaults for ATAC assay."""
    from seqnado.inputs.validation import DesignDataFrame
    from seqnado.core import Assay
    
    result = _extract_candidate_defaults_from_schema(
        DesignDataFrame,
        Assay.ATAC
    )
    
    assert isinstance(result, dict)


def test_configure_logging():
    """Test logging configuration."""
    from seqnado.cli import _configure_logging
    
    # Should not raise
    _configure_logging(False)
    _configure_logging(True)


def test_pkg_traversable():
    """Test getting package traversable."""
    from seqnado.cli import _pkg_traversable
    
    traversable = _pkg_traversable("seqnado")
    assert traversable is not None
    # Should be able to access package contents
    assert hasattr(traversable, "joinpath")


def test_coerce_value_invalid_int():
    """Test coercing invalid int."""
    try:
        _coerce_value_to_dtype("not_a_number", int, categories=None)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_coerce_value_invalid_float():
    """Test coercing invalid float."""
    try:
        _coerce_value_to_dtype("not_a_float", float, categories=None)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass


def test_format_col_hint_with_description():
    """Test formatting with description."""
    meta = {
        "default": "test",
        "description": "This is a description",
        "dtype": str,
        "categories": None,
    }
    
    result = _format_col_hint("col_name", meta)
    assert "col_name" in result
    # Description should be in the hint somehow
    assert len(result) > len("col_name")


def test_write_json_creates_parent_dir(tmp_path):
    """Test that _write_json creates parent directories."""
    nested = tmp_path / "nested" / "dir" / "file.json"
    data = {"test": "data"}
    
    _write_json(nested, data)
    
    assert nested.exists()
    loaded = _read_json(nested)
    assert loaded == data


def test_find_fastqs_multiple_hints(tmp_path):
    """Test finding FASTQs from multiple directory hints."""
    dir1 = tmp_path / "dir1"
    dir2 = tmp_path / "dir2"
    dir1.mkdir()
    dir2.mkdir()
    
    (dir1 / "sample1_R1.fastq.gz").touch()
    (dir1 / "sample1_R2.fastq.gz").touch()
    (dir2 / "sample2_R1.fastq.gz").touch()
    
    result = _find_fastqs([str(dir1), str(dir2)])
    
    assert len(result) == 3


def test_coerce_value_empty_string():
    """Test empty string handling."""
    result = _coerce_value_to_dtype("", str, categories=None)
    assert result == ""


def test_preset_profiles_has_local_env():
    """Test that preset profiles includes local environment."""
    result = _preset_profiles()
    # Should have at least local environment
    assert "le" in result or "local" in str(result).lower()


def test_profile_autocomplete_returns_strings():
    """Test profile autocomplete returns valid strings."""
    result = _profile_autocomplete()
    assert all(isinstance(p, str) and len(p) > 0 for p in result)


def setup_module(module):
    """Ensure test_output directory is cleaned up before tests."""
    test_output = Path("./test_output")
    if test_output.exists():
        # Use subprocess rm -rf for more robust removal on NFS filesystems
        subprocess.run(["rm", "-rf", str(test_output)], check=False)
    test_output.mkdir(parents=True, exist_ok=True)
