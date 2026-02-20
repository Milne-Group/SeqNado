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


@pytest.fixture
def rules_content(download_rules_file):
    """Return the content of download.smk."""
    return download_rules_file.read_text()


def _get_rule_section(content, rule_name):
    """Extract a rule's section from the smk file content."""
    lines = content.split("\n")
    section = []
    capture = False
    for line in lines:
        if f"rule {rule_name}:" in line:
            capture = True
        elif capture and (line.startswith("rule ") or line.startswith("ruleorder:")):
            break
        elif capture:
            section.append(line)
    return "\n".join(section)


def test_download_rules_file_exists(download_rules_file):
    """Test that the download.smk rules file exists."""
    assert download_rules_file.exists(), f"Rules file not found: {download_rules_file}"


# --- Rule existence tests ---


def test_has_prefetch_rule(rules_content):
    """Test that download.smk contains the geo_prefetch rule."""
    assert "rule geo_prefetch:" in rules_content


def test_has_paired_rule(rules_content):
    """Test that download.smk contains the geo_fastq_dump_paired rule."""
    assert "rule geo_fastq_dump_paired:" in rules_content
    assert "geo_samples_paired" in rules_content


def test_has_single_rule(rules_content):
    """Test that download.smk contains the geo_fastq_dump_single rule."""
    assert "rule geo_fastq_dump_single:" in rules_content
    assert "geo_samples_single" in rules_content


def test_has_compress_rule(rules_content):
    """Test that download.smk contains the geo_compress_fastq_files rule."""
    assert "rule geo_compress_fastq_files:" in rules_content


def test_has_download_all_rule(rules_content):
    """Test that download.smk contains the geo_download_all rule."""
    assert "rule geo_download_all:" in rules_content


# --- Prefetch rule tests ---


def test_prefetch_uses_temp_output(rules_content):
    """Test that prefetch output is marked as temp for auto-cleanup."""
    section = _get_rule_section(rules_content, "geo_prefetch")
    assert "temp(" in section
    assert "sra_cache" in section


def test_prefetch_max_size(rules_content):
    """Test that prefetch specifies maximum download size."""
    section = _get_rule_section(rules_content, "geo_prefetch")
    assert "prefetch" in section
    assert "--max-size" in section
    assert "50G" in section or "100G" in section


def test_prefetch_uses_container(rules_content):
    """Test that prefetch uses sra-tools container."""
    section = _get_rule_section(rules_content, "geo_prefetch")
    assert "container:" in section
    assert "sra-tools" in section


# --- Paired rule tests ---


def test_paired_rule_outputs_r1_r2(rules_content):
    """Test that paired rule outputs R1 and R2 fastq files."""
    section = _get_rule_section(rules_content, "geo_fastq_dump_paired")
    assert "_R1.fastq" in section
    assert "_R2.fastq" in section


def test_paired_rule_uses_split_3(rules_content):
    """Test that paired rule uses fasterq-dump --split-3."""
    section = _get_rule_section(rules_content, "geo_fastq_dump_paired")
    assert "fasterq-dump" in section
    assert "--split-3" in section


def test_paired_rule_has_wildcard_constraints(rules_content):
    """Test that paired rule constrains sample_name wildcard."""
    section = _get_rule_section(rules_content, "geo_fastq_dump_paired")
    assert "wildcard_constraints:" in section
    assert "geo_samples_paired" in section


def test_paired_rule_uses_container(rules_content):
    """Test that paired rule uses sra-tools container."""
    section = _get_rule_section(rules_content, "geo_fastq_dump_paired")
    assert "container:" in section
    assert "sra-tools" in section


# --- Single rule tests ---


def test_single_rule_output_structure(rules_content):
    """Test that single-end rule has single output without R1/R2."""
    section = _get_rule_section(rules_content, "geo_fastq_dump_single")
    assert "output:" in section
    assert "_R1" not in section
    assert "_R2" not in section


def test_single_rule_has_wildcard_constraints(rules_content):
    """Test that single rule constrains sample_name wildcard."""
    section = _get_rule_section(rules_content, "geo_fastq_dump_single")
    assert "wildcard_constraints:" in section
    assert "geo_samples_single" in section


# --- Compress rule tests ---
def test_compress_uses_pigz(rules_content):
    """Test that compress rule uses pigz for parallel compression."""
    section = _get_rule_section(rules_content, "geo_compress_fastq_files")
    assert "pigz" in section


def test_compress_uses_threads(rules_content):
    """Test that compress rule uses multiple threads."""
    section = _get_rule_section(rules_content, "geo_compress_fastq_files")
    assert "threads:" in section


def test_compress_input_output_pattern(rules_content):
    """Test that compress rule converts .fastq to .fastq.gz."""
    section = _get_rule_section(rules_content, "geo_compress_fastq_files")
    assert ".fastq.gz" in section
    assert "{filename}" in section


# --- Download all rule tests ---


def test_download_all_requests_gz_files(rules_content):
    """Test that download_all expects compressed outputs."""
    section = _get_rule_section(rules_content, "geo_download_all")
    assert "_R1.fastq.gz" in section
    assert "_R2.fastq.gz" in section
    assert ".fastq.gz" in section

# --- General rule tests ---
def test_ruleorder_exists(rules_content):
    """Test that ruleorder is defined for paired vs single."""
    assert "ruleorder: geo_fastq_dump_paired > geo_fastq_dump_single" in rules_content


@pytest.mark.unit
def test_no_empty_placeholder_files(rules_content):
    """Test that rules do NOT create empty R2 placeholder files."""
    assert "touch {wildcards.sample_name}_R2.fastq.gz" not in rules_content
    assert "touch {params.outdir}/{wildcards.sample_name}_R2.fastq.gz" not in rules_content


# --- Config structure tests ---
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

    with open(config_file) as f:
        loaded = yaml.safe_load(f)

    assert loaded == mock_geo_config


def test_rules_have_resources(rules_content):
    """Test that rules specify appropriate resources."""
    assert "mem" in rules_content
    assert "runtime" in rules_content
