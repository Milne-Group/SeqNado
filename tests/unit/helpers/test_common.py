"""Tests for seqnado.workflow.helpers.common module."""

import pytest

from seqnado import DataScalingTechnique
from seqnado.workflow.helpers.common import (
    define_memory_requested,
    define_time_requested,
    get_scale_method,
)


class TestDefineMemoryRequested:
    """Tests for define_memory_requested function."""

    def test_default_memory(self):
        """Test default memory allocation (1 attempt, 1G initial)."""
        result = define_memory_requested()
        assert result == "1G"

    def test_memory_scaling_with_attempts(self):
        """Test memory doubles with each attempt."""
        assert define_memory_requested(attempts=1, initial_value=2) == "2G"
        assert define_memory_requested(attempts=2, initial_value=2) == "4G"
        assert define_memory_requested(attempts=3, initial_value=2) == "8G"
        assert define_memory_requested(attempts=4, initial_value=2) == "16G"

    def test_memory_with_scale_factor(self):
        """Test memory with scale factor."""
        result = define_memory_requested(attempts=1, initial_value=4, scale=2.5)
        assert result == "10G"

    def test_memory_avoids_decimals(self):
        """Test that memory values are always integers."""
        result = define_memory_requested(attempts=1, initial_value=3, scale=1.5)
        assert result == "4G"  # int(3 * 1.5) = 4
        assert "." not in result


class TestDefineTimeRequested:
    """Tests for define_time_requested function."""

    def test_default_time(self):
        """Test default time allocation (1 attempt, 1 hour)."""
        result = define_time_requested()
        assert result == "1.0h"

    def test_time_scaling_with_attempts(self):
        """Test time doubles with each attempt."""
        assert define_time_requested(attempts=1, initial_value=2) == "2.0h"
        assert define_time_requested(attempts=2, initial_value=2) == "4.0h"
        assert define_time_requested(attempts=3, initial_value=2) == "8.0h"

    def test_time_with_scale_factor(self):
        """Test time with scale factor."""
        result = define_time_requested(attempts=1, initial_value=4, scale=2.5)
        assert result == "10.0h"

    def test_time_allows_decimals(self):
        """Test that time can have decimal values."""
        result = define_time_requested(attempts=1, initial_value=3, scale=1.5)
        assert result == "4.5h"


class TestGetScaleMethod:
    """Tests for get_scale_method function."""

    def test_scale_method_default(self):
        """Test default scale method (scaled-per-sample only)."""
        config = {}
        result = get_scale_method(config)
        assert result == ["scaled-per-sample"]  # Returned as lowercase string from .value

    def test_scale_method_with_spikein(self):
        """Test scale method with spikein."""
        config = {"spikein": True}
        result = get_scale_method(config)
        assert "scaled-per-sample" in result
        assert "spikein" in result
        assert len(result) == 2

    def test_scale_method_with_scale(self):
        """Test scale method with scaled-per-group scaling."""
        config = {"scale": True}
        result = get_scale_method(config)
        assert "scaled-per-sample" in result
        assert "scaled-per-group" in result
        assert len(result) == 2

    def test_scale_method_spikein_takes_precedence(self):
        """Test that spikein takes precedence over scale."""
        config = {"spikein": True, "scale": True}
        result = get_scale_method(config)
        assert "spikein" in result
        assert "scaled-per-group" not in result
