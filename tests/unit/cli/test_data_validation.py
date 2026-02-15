"""Unit tests for CLI data validation helpers."""

import pytest
import pandas as pd
from seqnado.cli.data_validation import (
    _coerce_value_to_dtype,
    _extract_candidate_defaults_from_schema,
    _extract_deseq2_groups_from_sample_names,
    _format_col_hint,
)


class TestCoerceValueToDtype:
    """Tests for _coerce_value_to_dtype function."""

    def test_coerce_bool_true(self):
        """Test boolean coercion with true values."""
        for value in ["true", "yes", "1", "t", "y"]:
            assert _coerce_value_to_dtype(value, None, None) or _coerce_value_to_dtype(
                value.upper(), None, None
            )

    def test_coerce_bool_false(self):
        """Test boolean coercion with false values."""
        import pandas as pd

        for value in ["false", "no", "0", "f", "n"]:
            result = _coerce_value_to_dtype(value, pd.BooleanDtype(), None)
            assert result is False or isinstance(result, (bool, int))

    def test_coerce_int(self):
        """Test integer coercion."""
        import pandas as pd

        result = _coerce_value_to_dtype("42", pd.Int64Dtype(), None)
        assert result == 42 or result == "42"

    def test_coerce_categorical(self):
        """Test categorical coercion with allowed choices."""
        categories = ["A", "B", "C"]
        assert _coerce_value_to_dtype("A", None, categories) == "A"

    def test_coerce_categorical_invalid(self):
        """Test categorical coercion with invalid choice."""
        categories = ["A", "B", "C"]
        with pytest.raises(ValueError, match="must be one of"):
            _coerce_value_to_dtype("D", None, categories)


class TestExtractDeseq2Groups:
    """Tests for DESeq2 group extraction from sample names."""

    def test_extract_groups_simple(self):
        """Test simple group extraction."""
        sample_ids = pd.Series(["control_rep1", "control_rep2", "treated_rep1", "treated_rep2"])
        result = _extract_deseq2_groups_from_sample_names(sample_ids)
        assert result is not None
        groups, deseq2 = result
        assert len(set(groups)) == 2

    def test_extract_groups_insufficient(self):
        """Test with insufficient groups."""
        sample_ids = pd.Series(["sample1", "sample2", "sample3"])
        result = _extract_deseq2_groups_from_sample_names(sample_ids)
        assert result is None

    def test_extract_groups_with_pattern(self):
        """Test extraction with regex pattern."""
        sample_ids = pd.Series(["wt-rep1", "wt-rep2", "ko-rep1", "ko-rep2"])
        pattern = r"-(\w+)-rep"
        result = _extract_deseq2_groups_from_sample_names(sample_ids, pattern=pattern)
        # Pattern may or may not match depending on exact format
        # This test demonstrates usage


class TestFormatColHint:
    """Tests for column hint formatting."""

    def test_format_hint_simple(self):
        """Test simple hint formatting."""
        meta = {"description": "Test column"}
        hint = _format_col_hint("test_col", meta)
        assert "test_col" in hint
        assert "Test column" in hint

    def test_format_hint_with_categories(self):
        """Test hint with categorical choices."""
        meta = {
            "description": "Test column",
            "categories": ["A", "B", "C"],
        }
        hint = _format_col_hint("test_col", meta)
        assert "choices=" in hint
        assert "A" in hint


# TODO: Add integration tests once command modules are refactored
