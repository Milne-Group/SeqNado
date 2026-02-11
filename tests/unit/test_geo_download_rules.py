"""
Unit tests for GEO download Snakemake rules
"""

from pathlib import Path

import pytest
import yaml


@pytest.fixture
def mock_geo_config():
    """Provide a mock configuration for GEO download rules."""
    return {
        "geo_samples_paired": {
            "GSM001-WT_rep1": {
                "srr": "SRR123456",
                "gsm": "GSM001",
                "sample": "WT_rep1",
            },
            "GSM002-WT_rep2": {
                "srr": "SRR123457",
                "gsm": "GSM002",
                "sample": "WT_rep2",
            },
        },
        "geo_samples_single": {
            "GSM003-KO_rep1": {
                "srr": "SRR123458",
                "gsm": "GSM003",
                "sample": "KO_rep1",
            },
        },
        "geo_outdir": "geo_data",
    }


@pytest.fixture
def download_rules_file():
    """Return the path to the download.smk rules file."""
    return Path(__file__).parent.parent.parent / "seqnado" / "workflow" / "rules" / "geo" / "download.smk"


def test_download_rules_file_exists(download_rules_file):
    """Test that the download.smk rules file exists."""
    assert download_rules_file.exists(), f"Rules file not found: {download_rules_file}"


def test_download_rules_has_paired_rule(download_rules_file):
    """Test that download.smk contains the geo_download_paired rule."""
    content = download_rules_file.read_text()
    assert "rule geo_download_paired:" in content
    assert "geo_samples_paired" in content


def test_download_rules_has_single_rule(download_rules_file):
    """Test that download.smk contains the geo_download_single rule."""
    content = download_rules_file.read_text()
    assert "rule geo_download_single:" in content
    assert "geo_samples_single" in content


def test_download_rules_has_all_rule(download_rules_file):
    """Test that download.smk contains the geo_download_all rule."""
    content = download_rules_file.read_text()
    assert "rule geo_download_all:" in content


def test_paired_rule_output_structure(download_rules_file):
    """Test that paired-end rule has correct output structure."""
    content = download_rules_file.read_text()

    # Paired rule should output both R1 and R2
    assert "_R1.fastq.gz" in content
    assert "_R2.fastq.gz" in content


def test_single_rule_output_structure(download_rules_file):
    """Test that single-end rule has correct output structure."""
    content = download_rules_file.read_text()

    # Single rule should NOT create empty R2 files
    # Check that single rule exists and has single output
    lines = content.split("\n")
    in_single_rule = False
    single_rule_section = []

    for line in lines:
        if "rule geo_download_single:" in line:
            in_single_rule = True
        elif in_single_rule:
            if line.startswith("rule "):
                break
            single_rule_section.append(line)

    single_rule_text = "\n".join(single_rule_section)

    # Should have output but NOT create empty R2
    assert "output:" in single_rule_text
    assert "touch" not in single_rule_text or "_R2" not in single_rule_text


def test_rules_use_container(download_rules_file):
    """Test that rules specify an Apptainer/Singularity container."""
    content = download_rules_file.read_text()
    assert "container:" in content
    assert "sra-tools" in content


def test_rules_have_retry_logic(download_rules_file):
    """Test that rules include retry logic for downloads."""
    content = download_rules_file.read_text()
    assert "max_retries" in content
    assert "prefetch" in content
    assert "attempt" in content


def test_rules_use_pigz(download_rules_file):
    """Test that rules use pigz for parallel compression."""
    content = download_rules_file.read_text()
    assert "pigz" in content


def test_rules_cleanup_sra_files(download_rules_file):
    """Test that rules clean up SRA files after extraction."""
    content = download_rules_file.read_text()
    assert "rm -rf" in content and ".sra" in content


def test_config_structure_paired(mock_geo_config):
    """Test that paired samples config has correct structure."""
    paired = mock_geo_config["geo_samples_paired"]
    assert len(paired) == 2

    for sample_name, sample_info in paired.items():
        assert "srr" in sample_info
        assert "gsm" in sample_info
        assert "sample" in sample_info
        assert sample_info["srr"].startswith("SRR")
        assert sample_info["gsm"].startswith("GSM")


def test_config_structure_single(mock_geo_config):
    """Test that single samples config has correct structure."""
    single = mock_geo_config["geo_samples_single"]
    assert len(single) == 1

    for sample_name, sample_info in single.items():
        assert "srr" in sample_info
        assert "gsm" in sample_info
        assert "sample" in sample_info


def test_config_has_outdir(mock_geo_config):
    """Test that config includes output directory."""
    assert "geo_outdir" in mock_geo_config
    assert isinstance(mock_geo_config["geo_outdir"], str)


def test_sample_name_format():
    """Test that sample names follow expected format."""
    gsm = "GSM123456"
    sample = "WT_treatment_rep1"
    expected_name = f"{gsm}-{sample}"

    assert expected_name == "GSM123456-WT_treatment_rep1"
    assert "-" in expected_name
    assert expected_name.startswith("GSM")


def test_config_yaml_serializable(mock_geo_config, tmp_path):
    """Test that config can be serialized to YAML (for Snakemake)."""
    config_file = tmp_path / "test_config.yaml"

    with open(config_file, "w") as f:
        yaml.dump(mock_geo_config, f)

    assert config_file.exists()

    # Read it back and verify structure
    with open(config_file) as f:
        loaded = yaml.safe_load(f)

    assert loaded == mock_geo_config


def test_paired_rule_resources(download_rules_file):
    """Test that rules specify appropriate resources."""
    content = download_rules_file.read_text()

    # Should specify threads, memory, and runtime
    assert "threads:" in content
    assert "mem_mb" in content or "mem" in content
    assert "runtime" in content or "time" in content


def test_rules_have_logging(download_rules_file):
    """Test that rules specify log files."""
    content = download_rules_file.read_text()
    assert "log:" in content
    assert "logs/geo_download" in content


@pytest.mark.unit
def test_no_empty_placeholder_files(download_rules_file):
    """Test that rules do NOT create empty R2 placeholder files."""
    content = download_rules_file.read_text()

    # Should NOT have any logic that creates empty R2 files
    # This was the bad pattern we removed
    assert "touch {wildcards.sample_name}_R2.fastq.gz" not in content
    assert "touch {params.outdir}/{wildcards.sample_name}_R2.fastq.gz" not in content


def test_fasterq_dump_options(download_rules_file):
    """Test that fasterq-dump is called with appropriate options."""
    content = download_rules_file.read_text()

    assert "fasterq-dump" in content

    # Paired should use --split-files
    lines = content.split("\n")
    paired_section = []
    capture = False

    for line in lines:
        if "rule geo_download_paired:" in line:
            capture = True
        elif capture and line.startswith("rule "):
            break
        elif capture:
            paired_section.append(line)

    paired_text = "\n".join(paired_section)
    assert "--split-files" in paired_text
    assert "fasterq-dump" in paired_text


def test_prefetch_max_size(download_rules_file):
    """Test that prefetch specifies maximum download size."""
    content = download_rules_file.read_text()
    assert "prefetch" in content
    assert "--max-size" in content
    # Should allow large files (50G or more)
    assert "50G" in content or "100G" in content
