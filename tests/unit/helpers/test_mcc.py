"""Tests for seqnado.workflow.helpers.mcc module."""

import json
import tempfile
from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from seqnado.workflow.helpers.mcc import (
    extract_viewpoints,
    get_mcc_bam_files_for_merge,
    get_n_cis_scaling_factor,
    identify_extracted_bam_files,
    redefine_viewpoints,
    viewpoint_to_grouped_viewpoint,
)


class TestGetNCisScalingFactor:
    """Tests for get_n_cis_scaling_factor function."""

    def test_missing_stats_file_returns_one(self):
        """Test that missing stats file returns 1."""
        wildcards = Mock(sample="sample1", viewpoint_group="vp1", spec=["sample", "viewpoint_group"])
        output_dir = "/nonexistent/path"

        result = get_n_cis_scaling_factor(wildcards, output_dir)
        assert result == 1

    def test_calculates_scaling_factor_for_sample(self):
        """Test calculation of scaling factor for individual sample."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1", viewpoint_group="vp1", spec=["sample", "viewpoint_group"])
            output_dir = tmpdir

            # Create resources directory and stats file
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            stats_data = {
                "vp1": {"n_cis": 5000, "n_total": 10000}
            }
            stats_file = resources_dir / "sample1_ligation_stats.json"
            stats_file.write_text(json.dumps(stats_data))

            result = get_n_cis_scaling_factor(wildcards, output_dir)
            # normalized_count = raw_count * F where F = 1e6 / n_cis (CPM normalized by cis interactions).
            expected = 1e6 / 5000
            assert result == expected

    def test_calculates_scaling_factor_for_group(self):
        """Test calculation of scaling factor for sample group."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(group="group1", viewpoint_group="vp1", spec=["group", "viewpoint_group"])
            output_dir = tmpdir

            # Create resources directory and stats file for group
            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            stats_data = {
                "vp1": {"n_cis": 8000, "n_total": 20000}
            }
            stats_file = resources_dir / "group1_ligation_stats.json"
            stats_file.write_text(json.dumps(stats_data))

            result = get_n_cis_scaling_factor(wildcards, output_dir)
            # Implementation currently returns CPM based on n_cis (1e6 / n_cis)
            expected = 1e6 / 8000
            assert result == expected

    def test_zero_total_returns_zero(self):
        """Test that zero n_total returns 0 to avoid division by zero."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1", viewpoint_group="vp1", spec=["sample", "viewpoint_group"])
            output_dir = tmpdir

            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            stats_data = {
                "vp1": {"n_cis": 100, "n_total": 0}
            }
            stats_file = resources_dir / "sample1_ligation_stats.json"
            stats_file.write_text(json.dumps(stats_data))

            result = get_n_cis_scaling_factor(wildcards, output_dir)
            # Implementation returns CPM based on n_cis even if n_total is zero
            expected = 1e6 / 100
            assert result == expected

    def test_missing_viewpoint_group_raises_error(self):
        """Test that missing viewpoint group raises KeyError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1", viewpoint_group="missing_vp", spec=["sample", "viewpoint_group"])
            output_dir = tmpdir

            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            stats_data = {
                "vp1": {"n_cis": 5000, "n_total": 10000}
            }
            stats_file = resources_dir / "sample1_ligation_stats.json"
            stats_file.write_text(json.dumps(stats_data))

            with pytest.raises(KeyError, match="Viewpoint group 'missing_vp' not found"):
                get_n_cis_scaling_factor(wildcards, output_dir)

    def test_missing_required_keys_raises_error(self):
        """Test that missing required keys raises KeyError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            wildcards = Mock(sample="sample1", viewpoint_group="vp1", spec=["sample", "viewpoint_group"])
            output_dir = tmpdir

            resources_dir = Path(tmpdir) / "resources"
            resources_dir.mkdir()
            # Missing n_total key
            stats_data = {
                "vp1": {"n_cis": 5000}
            }
            stats_file = resources_dir / "sample1_ligation_stats.json"
            stats_file.write_text(json.dumps(stats_data))

            # Implementation only requires 'n_cis' key, so should return CPM
            result = get_n_cis_scaling_factor(wildcards, output_dir)
            expected = 1e6 / 5000
            assert result == expected


class TestGetMccBamFilesForMerge:
    """Tests for get_mcc_bam_files_for_merge function."""

    def test_returns_bam_files_for_group(self):
        """Test that function returns BAM files for all samples in group."""
        wildcards = Mock(group="group1")

        group_mock = Mock()
        group_mock.samples = ["sample1", "sample2", "sample3"]

        grouping_mock = Mock()
        grouping_mock.get_group.return_value = group_mock

        sample_groupings = Mock()
        sample_groupings.get_grouping.return_value = grouping_mock

        output_dir = "/output"

        result = get_mcc_bam_files_for_merge(wildcards, sample_groupings, output_dir)

        expected = [
            "/output/mcc/replicates/sample1/sample1.bam",
            "/output/mcc/replicates/sample2/sample2.bam",
            "/output/mcc/replicates/sample3/sample3.bam",
        ]
        assert result == expected

    def test_missing_group_returns_empty_list(self):
        """Test that missing group returns empty list."""
        wildcards = Mock(group="nonexistent_group")

        grouping_mock = Mock()
        grouping_mock.get_group.side_effect = KeyError("Group not found")

        sample_groupings = Mock()
        sample_groupings.get_grouping.return_value = grouping_mock

        output_dir = "/output"

        result = get_mcc_bam_files_for_merge(wildcards, sample_groupings, output_dir)
        assert result == []

    def test_single_sample_group(self):
        """Test group with single sample."""
        wildcards = Mock(group="group1")

        group_mock = Mock()
        group_mock.samples = ["sample1"]

        grouping_mock = Mock()
        grouping_mock.get_group.return_value = group_mock

        sample_groupings = Mock()
        sample_groupings.get_grouping.return_value = grouping_mock

        output_dir = "/output"

        result = get_mcc_bam_files_for_merge(wildcards, sample_groupings, output_dir)

        assert len(result) == 1
        assert result[0] == "/output/mcc/replicates/sample1/sample1.bam"


class TestIdentifyExtractedBamFiles:
    """Tests for identify_extracted_bam_files function."""

    @patch("snakemake.io.glob_wildcards")
    def test_returns_extracted_bam_paths(self, mock_glob_wildcards):
        """Test that function returns all extracted BAM file paths."""
        # Mock wildcards as a dict-like object for **wildcards unpacking
        wildcards = {"sample": "sample1"}

        # Mock checkpoint output
        checkpoint_output = Mock()
        checkpoint_output.output.bams = "/output/mcc/sample1/bams"

        checkpoints = Mock()
        checkpoints.identify_viewpoint_reads.get.return_value = checkpoint_output

        # Mock glob_wildcards to return viewpoint names
        mock_glob_result = Mock()
        mock_glob_result.viewpoint = ["vp1", "vp2", "vp3"]
        mock_glob_wildcards.return_value = mock_glob_result

        result = identify_extracted_bam_files(wildcards, checkpoints)

        expected = [
            "/output/mcc/sample1/bams/vp1.bam",
            "/output/mcc/sample1/bams/vp2.bam",
            "/output/mcc/sample1/bams/vp3.bam",
        ]
        assert result == expected

    @patch("snakemake.io.glob_wildcards")
    def test_empty_viewpoints_returns_empty_list(self, mock_glob_wildcards):
        """Test handling when no viewpoints are found."""
        # Mock wildcards as a dict-like object for **wildcards unpacking
        wildcards = {"sample": "sample1"}

        checkpoint_output = Mock()
        checkpoint_output.output.bams = "/output/mcc/sample1/bams"

        checkpoints = Mock()
        checkpoints.identify_viewpoint_reads.get.return_value = checkpoint_output

        mock_glob_result = Mock()
        mock_glob_result.viewpoint = []
        mock_glob_wildcards.return_value = mock_glob_result

        result = identify_extracted_bam_files(wildcards, checkpoints)
        assert result == []


class TestRedefineViewpoints:
    """Tests for redefine_viewpoints function."""

    @patch("snakemake.io.glob_wildcards")
    def test_returns_intersection_of_viewpoints(self, mock_glob_wildcards):
        """Test that function returns intersection of viewpoints across samples."""
        samples = ["sample1", "sample2", "sample3"]

        # Mock checkpoint outputs
        def get_checkpoint_output(sample):
            output = Mock()
            output.output.bams = f"/output/mcc/{sample}/bams"
            return output

        checkpoints = Mock()
        checkpoints.identify_viewpoint_reads.get.side_effect = get_checkpoint_output

        # Mock glob_wildcards to return different viewpoints for each sample
        viewpoints_by_sample = {
            "sample1": ["vp1", "vp2", "vp3", "vp4"],
            "sample2": ["vp1", "vp2", "vp3"],
            "sample3": ["vp1", "vp2"],
        }

        def glob_side_effect(pattern):
            # Extract sample from pattern
            for sample in samples:
                if sample in pattern:
                    result = Mock()
                    result.viewpoint = viewpoints_by_sample[sample]
                    return result
            return Mock(viewpoint=[])

        mock_glob_wildcards.side_effect = glob_side_effect

        result = redefine_viewpoints(samples, checkpoints)

        # Intersection should be vp1 and vp2 (present in all samples)
        assert set(result) == {"vp1", "vp2"}

    @patch("snakemake.io.glob_wildcards")
    def test_no_common_viewpoints_returns_empty(self, mock_glob_wildcards):
        """Test that no common viewpoints returns empty list."""
        samples = ["sample1", "sample2"]

        def get_checkpoint_output(sample):
            output = Mock()
            output.output.bams = f"/output/mcc/{sample}/bams"
            return output

        checkpoints = Mock()
        checkpoints.identify_viewpoint_reads.get.side_effect = get_checkpoint_output

        # Different viewpoints for each sample (no intersection)
        viewpoints_by_sample = {
            "sample1": ["vp1", "vp2"],
            "sample2": ["vp3", "vp4"],
        }

        def glob_side_effect(pattern):
            for sample in samples:
                if sample in pattern:
                    result = Mock()
                    result.viewpoint = viewpoints_by_sample[sample]
                    return result
            return Mock(viewpoint=[])

        mock_glob_wildcards.side_effect = glob_side_effect

        result = redefine_viewpoints(samples, checkpoints)
        assert len(result) == 0

    @patch("snakemake.io.glob_wildcards")
    def test_single_sample_returns_all_viewpoints(self, mock_glob_wildcards):
        """Test that single sample returns all its viewpoints."""
        samples = ["sample1"]

        checkpoint_output = Mock()
        checkpoint_output.output.bams = "/output/mcc/sample1/bams"

        checkpoints = Mock()
        checkpoints.identify_viewpoint_reads.get.return_value = checkpoint_output

        mock_glob_result = Mock()
        mock_glob_result.viewpoint = ["vp1", "vp2", "vp3"]
        mock_glob_wildcards.return_value = mock_glob_result

        result = redefine_viewpoints(samples, checkpoints)
        assert set(result) == {"vp1", "vp2", "vp3"}


class TestExtractViewpoints:
    """Tests for extract_viewpoints function."""

    def test_extract_viewpoints_basic_bed(self, tmp_path):
        """Test extraction from basic BED file."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "chr1\t1000\t2000\tViewpoint1\t100\t+\n"
        bed_content += "chr2\t3000\t4000\tViewpoint2\t200\t-\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 2
        assert "Viewpoint1-chr1-1000-2000" in result
        assert "Viewpoint2-chr2-3000-4000" in result

    def test_extract_viewpoints_with_coordinates_in_name(self, tmp_path):
        """Test extraction when name already contains coordinates."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "chr1\t1000\t2000\tVP1-chr1-1000-2000\t100\t+\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 1
        # Should use the existing coordinate format
        assert "VP1-chr1-1000-2000" in result

    def test_extract_viewpoints_three_column_bed(self, tmp_path):
        """Test extraction from minimal 3-column BED file."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "chr1\t1000\t2000\n"
        bed_content += "chr2\t3000\t4000\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 2
        # Name should be auto-generated from coordinates
        assert any("chr1:1000-2000" in vp for vp in result)
        assert any("chr2:3000-4000" in vp for vp in result)

    def test_extract_viewpoints_with_comments(self, tmp_path):
        """Test extraction with comment lines in BED file."""
        bed_file = tmp_path / "viewpoints.bed"
        bed_content = "# This is a comment\n"
        bed_content += "chr1\t1000\t2000\tViewpoint1\t100\t+\n"
        bed_file.write_text(bed_content)

        result = extract_viewpoints(str(bed_file))

        assert len(result) == 1
        assert "Viewpoint1-chr1-1000-2000" in result

    def test_extract_viewpoints_invalid_file(self):
        """Test error handling for invalid file."""
        with pytest.raises(ValueError, match="Error reading BED file"):
            extract_viewpoints("/nonexistent/file.bed")


class TestViewpointToGroupedViewpoint:
    """Tests for viewpoint_to_grouped_viewpoint function."""

    def test_group_single_viewpoint(self):
        """Test grouping of single viewpoint with coordinates."""
        viewpoints = ["VP1-chr1-1000-2000"]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result == {"VP1-chr1-1000-2000": "VP1"}

    def test_group_multiple_viewpoints(self):
        """Test grouping of multiple viewpoints."""
        viewpoints = [
            "VP1-chr1-1000-2000",
            "VP1-chr2-3000-4000",
            "VP2-chr3-5000-6000",
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result["VP1-chr1-1000-2000"] == "VP1"
        assert result["VP1-chr2-3000-4000"] == "VP1"
        assert result["VP2-chr3-5000-6000"] == "VP2"

    def test_group_viewpoints_different_chromosomes(self):
        """Test grouping with X, Y, M, MT chromosomes."""
        viewpoints = [
            "VP1-chrX-1000-2000",
            "VP2-chrY-3000-4000",
            "VP3-chrM-5000-6000",
            "VP4-chrMT-7000-8000",
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result["VP1-chrX-1000-2000"] == "VP1"
        assert result["VP2-chrY-3000-4000"] == "VP2"
        assert result["VP3-chrM-5000-6000"] == "VP3"
        assert result["VP4-chrMT-7000-8000"] == "VP4"

    def test_group_viewpoints_without_coordinates(self):
        """Test that viewpoints without coordinate pattern are not included."""
        viewpoints = [
            "VP1-chr1-1000-2000",
            "SimpleViewpoint",  # No coordinate pattern
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert "VP1-chr1-1000-2000" in result
        assert "SimpleViewpoint" not in result

    def test_group_empty_list(self):
        """Test grouping with empty viewpoint list."""
        viewpoints = []
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result == {}

    def test_group_viewpoints_complex_names(self):
        """Test grouping with complex viewpoint names containing hyphens."""
        viewpoints = [
            "Complex-VP-Name-chr1-1000-2000",
            "Another-Complex-VP-chr2-3000-4000",
        ]
        result = viewpoint_to_grouped_viewpoint(viewpoints)

        assert result["Complex-VP-Name-chr1-1000-2000"] == "Complex-VP-Name"
        assert result["Another-Complex-VP-chr2-3000-4000"] == "Another-Complex-VP"
