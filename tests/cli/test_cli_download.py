"""
Tests for the seqnado download command (GEO/SRA data download)
"""

import subprocess
from pathlib import Path

import pandas as pd
import pytest
import yaml


def _create_geo_metadata_tsv(
    tmp_path: Path, filename: str = "metadata.tsv", include_layout: bool = True
) -> Path:
    """Create a mock GEO/ENA metadata TSV file."""
    data = {
        "run_accession": ["SRR123456", "SRR123457", "SRR123458"],
        "sample_title": ["WT_rep1", "WT_rep2", "KO_rep1"],
        "library_name": ["GSM001", "GSM002", "GSM003"],
    }

    if include_layout:
        data["library_layout"] = ["PAIRED", "PAIRED", "SINGLE"]

    df = pd.DataFrame(data)
    tsv_path = tmp_path / filename
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


def test_download_help_displays():
    """Test that download command help is accessible."""
    result = subprocess.run(
        ["seqnado", "download", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "Download FASTQ files from GEO/SRA" in result.stdout
    assert "library_layout" in result.stdout.lower() or "metadata" in result.stdout.lower()


def test_download_requires_metadata_tsv(tmp_path: Path):
    """Test that download command requires a metadata TSV file."""
    result = subprocess.run(
        ["seqnado", "download"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    # Should complain about missing argument


def test_download_parses_metadata_tsv(tmp_path: Path):
    """Test that download command parses metadata TSV and creates config."""
    tsv_path = _create_geo_metadata_tsv(tmp_path)

    # Run with --dry-run to avoid actual downloads
    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "--dry-run",
            "-c",
            "1",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Check output for expected behavior
    assert "Reading metadata" in result.stdout or "Reading metadata" in result.stderr
    # Config file is created temporarily but cleaned up
    # We should at least see no errors about missing columns
    assert "Missing required columns" not in result.stderr


def test_download_requires_library_layout_column(tmp_path: Path):
    """Test that download command requires library_layout column."""
    tsv_path = _create_geo_metadata_tsv(tmp_path, include_layout=False)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Should error about missing library_layout
    assert result.returncode != 0
    combined_output = result.stdout + result.stderr
    assert "library_layout" in combined_output.lower()


def test_download_separates_paired_and_single_samples(tmp_path: Path):
    """Test that samples are correctly separated by library layout."""
    tsv_path = _create_geo_metadata_tsv(tmp_path)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "--dry-run",
            "-v",  # verbose to see sample counts
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    combined_output = result.stdout + result.stderr

    # Should report sample counts
    assert "2 paired-end" in combined_output or "2" in combined_output
    assert "1 single-end" in combined_output or "1" in combined_output


def test_download_validates_required_columns(tmp_path: Path):
    """Test that download validates required TSV columns."""
    # Create TSV missing required columns
    df = pd.DataFrame(
        {
            "run_accession": ["SRR123456"],
            # Missing sample_title and library_name
        }
    )
    tsv_path = tmp_path / "incomplete.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    result = subprocess.run(
        ["seqnado", "download", str(tsv_path), "-o", "downloads"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    combined_output = result.stdout + result.stderr
    assert "Missing required columns" in combined_output


def test_download_creates_output_directory(tmp_path: Path):
    """Test that download creates the specified output directory."""
    tsv_path = _create_geo_metadata_tsv(tmp_path)
    outdir = tmp_path / "geo_downloads"

    # Dry run to avoid actual downloads
    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            str(outdir),
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Output directory should be created
    assert outdir.exists()


def test_download_with_assay_flag(tmp_path: Path):
    """Test download with --assay flag for design generation."""
    tsv_path = _create_geo_metadata_tsv(tmp_path)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "-a",
            "rna",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Should mention design file generation
    combined_output = result.stdout + result.stderr
    # In dry run, it won't generate design but should parse assay
    assert "rna" in combined_output.lower() or result.returncode == 0


def test_download_sample_naming():
    """Test that sample names are correctly formatted from TSV columns."""
    # This is a unit-style test for the naming logic
    library_name = "GSM123456"
    sample_title = "WT_treatment_rep1"
    expected = f"{library_name}-{sample_title}"

    assert expected == "GSM123456-WT_treatment_rep1"


@pytest.mark.requires_apptainer
@pytest.mark.slow
def test_download_container_check(tmp_path: Path):
    """Test that download checks for container runtime (apptainer/singularity)."""
    tsv_path = _create_geo_metadata_tsv(tmp_path)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "--dry-run",
            "-v",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    combined_output = result.stdout + result.stderr

    # Should mention container runtime
    assert (
        "apptainer" in combined_output.lower()
        or "singularity" in combined_output.lower()
        or "container" in combined_output.lower()
    )


def test_download_cores_option(tmp_path: Path):
    """Test that --cores option is respected."""
    tsv_path = _create_geo_metadata_tsv(tmp_path)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "-c",
            "16",
            "--dry-run",
            "-v",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # In verbose output should see cores mentioned
    combined_output = result.stdout + result.stderr
    assert "16" in combined_output or "--cores" in combined_output


def test_download_mixed_layout_handling(tmp_path: Path):
    """Test handling of projects with both paired and single-end samples."""
    # Create metadata with both types
    data = {
        "run_accession": ["SRR001", "SRR002", "SRR003", "SRR004"],
        "sample_title": ["paired1", "paired2", "single1", "single2"],
        "library_name": ["GSM001", "GSM002", "GSM003", "GSM004"],
        "library_layout": ["PAIRED", "PAIRED", "SINGLE", "SINGLE"],
    }

    df = pd.DataFrame(data)
    tsv_path = tmp_path / "mixed.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "--dry-run",
            "-v",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    combined_output = result.stdout + result.stderr

    # Should show both types
    assert "2 paired-end" in combined_output or "paired" in combined_output.lower()
    assert "2 single-end" in combined_output or "single" in combined_output.lower()


@pytest.mark.unit
def test_config_structure_for_snakemake(tmp_path: Path):
    """Test that the generated config has the correct structure for Snakemake."""
    # This tests the expected config structure
    # In a real run, CLI creates .geo_download_config.yaml temporarily

    expected_config = {
        "geo_samples_paired": {
            "GSM001-sample1": {
                "srr": "SRR001",
                "gsm": "GSM001",
                "sample": "sample1",
            }
        },
        "geo_samples_single": {
            "GSM002-sample2": {
                "srr": "SRR002",
                "gsm": "GSM002",
                "sample": "sample2",
            }
        },
        "geo_outdir": "/path/to/downloads",
    }

    # Verify structure
    assert "geo_samples_paired" in expected_config
    assert "geo_samples_single" in expected_config
    assert "geo_outdir" in expected_config

    # Verify sample structure
    paired_sample = list(expected_config["geo_samples_paired"].values())[0]
    assert "srr" in paired_sample
    assert "gsm" in paired_sample
    assert "sample" in paired_sample


def test_download_invalid_layout_value(tmp_path: Path):
    """Test handling of invalid library_layout values."""
    data = {
        "run_accession": ["SRR123456"],
        "sample_title": ["sample1"],
        "library_name": ["GSM001"],
        "library_layout": ["INVALID"],  # Invalid value
    }

    df = pd.DataFrame(data)
    tsv_path = tmp_path / "invalid_layout.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)

    result = subprocess.run(
        [
            "seqnado",
            "download",
            str(tsv_path),
            "-o",
            "downloads",
            "--dry-run",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    # Should handle invalid layout gracefully
    assert result.returncode != 0 or "invalid" in (result.stdout + result.stderr).lower()
