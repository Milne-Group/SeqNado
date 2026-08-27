"""Tests for seqnado.workflow.helpers.peaks module."""

import re
from unittest.mock import MagicMock, Mock

import pytest

from seqnado import FileType, DataScalingTechnique
from seqnado.inputs import FastqCollection, FastqCollectionForIP
from seqnado.config.third_party_tools import CommandLineArguments
from seqnado.workflow.helpers.peaks import (
    get_control_file,
    correct_macs_options,
    get_lanceotron_call_peaks_threshold,
)


class TestGetControlFile:
    """Tests for get_control_file function."""

    def test_fastq_collection_returns_no_control(self):
        """Test that FastqCollection returns NO_CONTROL_PRESENT_FOR_FILE."""
        wildcards = Mock(sample_id="sample1")
        input_files = Mock(spec=FastqCollection)
        output_dir = "/test/output"

        result = get_control_file(wildcards, FileType.BAM, input_files, output_dir)
        assert result == "NO_CONTROL_PRESENT_FOR_FILE"

    def test_fastq_collection_for_ip_no_control(self):
        """Test FastqCollectionForIP with no control returns NO_CONTROL_PRESENT_FOR_FILE."""
        wildcards = Mock(sample_id="sample1")
        input_files = Mock(spec=FastqCollectionForIP)
        input_files.get_control_performed.return_value = None
        output_dir = "/test/output"

        result = get_control_file(wildcards, FileType.BAM, input_files, output_dir)
        assert result == "NO_CONTROL_PRESENT_FOR_FILE"

    def test_get_control_bam_file(self):
        """Test getting control BAM file path."""
        wildcards = Mock(sample_id="treatment1")
        input_files = Mock(spec=FastqCollectionForIP)
        input_files.get_control_performed.return_value = "control1"
        output_dir = "/test/output"

        result = get_control_file(wildcards, FileType.BAM, input_files, output_dir)
        assert result == "/test/output/aligned/control1.bam"

    def test_get_control_bigwig_file(self):
        """Test getting control BigWig file path."""
        wildcards = Mock(sample_id="treatment1")
        input_files = Mock(spec=FastqCollectionForIP)
        input_files.get_control_performed.return_value = "control1"
        output_dir = "/test/output"

        result = get_control_file(wildcards, FileType.BIGWIG, input_files, output_dir)
        expected = f"/test/output/bigwigs/deeptools/{DataScalingTechnique.PER_SAMPLE.value}/control1.bigWig"
        assert result == expected

    def test_get_control_tag_directory(self):
        """Test getting control tag directory path."""
        wildcards = Mock(sample_id="treatment1")
        input_files = Mock(spec=FastqCollectionForIP)
        input_files.get_control_performed.return_value = "control1"
        output_dir = "/test/output"

        result = get_control_file(wildcards, FileType.TAG_DIRECTORY, input_files, output_dir)
        assert result == "/test/output/tag_dirs/control1"

    def test_unsupported_file_type_raises_error(self):
        """Test that unsupported file type raises ValueError."""
        wildcards = Mock(sample_id="treatment1")
        input_files = Mock(spec=FastqCollectionForIP)
        input_files.get_control_performed.return_value = "control1"
        output_dir = "/test/output"

        # Use a mock FileType that's not handled
        unsupported_type = Mock()
        unsupported_type.__str__ = lambda x: "UNSUPPORTED"

        with pytest.raises(ValueError, match="Unsupported file type"):
            get_control_file(wildcards, unsupported_type, input_files, output_dir)


class TestCorrectMacsOptions:
    """Tests for correct_macs_options function."""

    def test_single_sample_paired_end(self):
        """Test single sample with paired-end data keeps -f flag."""
        wildcards = Mock(sample_id="sample1", spec=["sample_id"])
        options = CommandLineArguments(value="-f BAMPE --other-option")
        input_files = Mock()
        input_files.is_paired_end.return_value = True
        sample_groupings = Mock()

        result = correct_macs_options(wildcards, options, input_files, sample_groupings)
        assert "-f" in str(result) or "BAMPE" in str(result)

    def test_single_sample_single_end(self):
        """Test single sample with single-end data removes -f flag."""
        wildcards = Mock(sample_id="sample1", spec=["sample_id"])
        options = CommandLineArguments(value="-f BAMPE --other-option")
        input_files = Mock()
        input_files.is_paired_end.return_value = False
        sample_groupings = Mock()

        result = correct_macs_options(wildcards, options, input_files, sample_groupings)
        result_str = str(result)
        assert "-f" not in result_str or "BAMPE" not in result_str

    def test_grouped_samples_all_paired(self):
        """Test grouped samples where all are paired-end keeps -f flag."""
        wildcards = Mock(group="group1", sample_id="sample1")
        options = CommandLineArguments(value="-f BAMPE --other-option")

        group_mock = Mock()
        group_mock.samples = ["sample1", "sample2"]

        grouping_mock = Mock()
        grouping_mock.get_group.return_value = group_mock

        sample_groupings = Mock()
        sample_groupings.get_grouping.return_value = grouping_mock

        input_files = Mock()
        input_files.is_paired_end.return_value = True

        result = correct_macs_options(wildcards, options, input_files, sample_groupings)
        assert "-f" in str(result) or "BAMPE" in str(result)

    def test_grouped_samples_mixed_pairing(self):
        """Test grouped samples with mixed pairing removes -f flag."""
        wildcards = Mock(group="group1", sample_id="sample1")
        options = CommandLineArguments(value="-f BAMPE --other-option")

        group_mock = Mock()
        group_mock.samples = ["sample1", "sample2"]

        grouping_mock = Mock()
        grouping_mock.get_group.return_value = group_mock

        sample_groupings = Mock()
        sample_groupings.get_grouping.return_value = grouping_mock

        input_files = Mock()
        # First sample paired, second single
        input_files.is_paired_end.side_effect = [True, False]

        result = correct_macs_options(wildcards, options, input_files, sample_groupings)
        result_str = str(result)
        assert "-f" not in result_str or "BAMPE" not in result_str


class TestGetLanceotronCallPeaksThreshold:
    """Tests for get_lanceotron_call_peaks_threshold function."""

    def test_extract_threshold_with_c_flag(self):
        """Test extracting threshold with -c flag."""
        wildcards = Mock()
        config = Mock()
        config.third_party_tools.lanceotron.call_peaks.command_line_arguments = "-c 0.8 --other-option"

        result = get_lanceotron_call_peaks_threshold(wildcards, config)
        assert result == "0.8"

    def test_extract_threshold_with_integer(self):
        """Test extracting integer threshold."""
        wildcards = Mock()
        config = Mock()
        config.third_party_tools.lanceotron.call_peaks.command_line_arguments = "-c 1 --other-option"

        result = get_lanceotron_call_peaks_threshold(wildcards, config)
        assert result == "1"

    def test_no_threshold_returns_default(self):
        """Test that missing threshold returns default 0.5."""
        wildcards = Mock()
        config = Mock()
        config.third_party_tools.lanceotron.call_peaks.command_line_arguments = "--other-option"

        result = get_lanceotron_call_peaks_threshold(wildcards, config)
        assert result == "0.5"

    def test_threshold_without_space(self):
        """Test extracting threshold when directly attached to -c flag."""
        wildcards = Mock()
        config = Mock()
        config.third_party_tools.lanceotron.call_peaks.command_line_arguments = "-c0.75"

        result = get_lanceotron_call_peaks_threshold(wildcards, config)
        assert result == "0.75"
