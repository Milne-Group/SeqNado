import subprocess
from pathlib import Path
import pandas as pd


def _write_fastq(tmp: Path, name: str) -> Path:
    p = tmp / name
    p.touch()
    return p


def test_cli_design_generates_csv(tmp_path: Path, monkeypatch):
    # create a couple fastqs
    r1 = _write_fastq(tmp_path, "chipA_H3K27ac_L001_R1_001.fastq.gz")
    r2 = _write_fastq(tmp_path, "chipA_H3K27ac_L001_R2_001.fastq.gz")

    # Invoke CLI via subprocess with explicit file paths
    result = subprocess.run(
        [
            "seqnado",
            "design",
            "chip",
            "-o",
            "metadata.csv",
            "--no-interactive",
            "--accept-all-defaults",
            str(r1),
            str(r2),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    out = tmp_path / "metadata.csv"
    assert out.exists(), "metadata.csv should be created"
    text = out.read_text()
    # basic sanity checks on columns
    assert "sample" in text.lower()  # Could be sample_id or sample_name
    assert "r1" in text.lower()


def test_cli_design_includes_default_columns(tmp_path: Path):
    """Test that design command with --accept-all-defaults includes scaling_group column."""
    # Create dummy fastq files
    r1 = _write_fastq(tmp_path, "test_sample_R1.fastq.gz")
    r2 = _write_fastq(tmp_path, "test_sample_R2.fastq.gz")

    # Run seqnado design with both flags
    result = subprocess.run(
        [
            "seqnado",
            "design",
            "rna",
            "--no-interactive",
            "--accept-all-defaults",
            "-o",
            "metadata.csv",
            str(r1),
            str(r2),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"

    # Read the generated CSV
    metadata_file = tmp_path / "metadata.csv"
    assert metadata_file.exists(), "metadata.csv should be created"

    df = pd.read_csv(metadata_file)

    # Check that scaling_group column is present with default value
    assert "scaling_group" in df.columns, "scaling_group column should be present"
    assert all(df["scaling_group"] == "default"), "scaling_group should have default value 'default'"


def test_cli_design_consensus_by_ip(tmp_path: Path):
    """--consensus-by ip should group ChIP samples by antibody into consensus_group."""
    _write_fastq(tmp_path, "sample1_H3K27ac_R1.fastq.gz")
    _write_fastq(tmp_path, "sample1_H3K27ac_R2.fastq.gz")
    _write_fastq(tmp_path, "sample2_H3K27ac_R1.fastq.gz")
    _write_fastq(tmp_path, "sample2_H3K27ac_R2.fastq.gz")
    _write_fastq(tmp_path, "sample1_Menin_R1.fastq.gz")
    _write_fastq(tmp_path, "sample1_Menin_R2.fastq.gz")

    result = subprocess.run(
        [
            "seqnado",
            "design",
            "chip",
            "--no-interactive",
            "--accept-all-defaults",
            "--consensus-by",
            "ip",
            "-o",
            "metadata.csv",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    df = pd.read_csv(tmp_path / "metadata.csv")

    assert "consensus_group" in df.columns
    # every row's consensus_group should match its ip (antibody)
    assert (df["consensus_group"] == df["ip"]).all()
    assert set(df["consensus_group"]) == {"H3K27ac", "Menin"}


def test_cli_design_condition_by_regex(tmp_path: Path):
    """--condition-by should populate 'condition' from a regex extracted from sample names."""
    _write_fastq(tmp_path, "sample-control-rep1_R1.fastq.gz")
    _write_fastq(tmp_path, "sample-control-rep1_R2.fastq.gz")
    _write_fastq(tmp_path, "sample-treated-rep1_R1.fastq.gz")
    _write_fastq(tmp_path, "sample-treated-rep1_R2.fastq.gz")

    result = subprocess.run(
        [
            "seqnado",
            "design",
            "atac",
            "--no-interactive",
            "--accept-all-defaults",
            "--condition-by",
            r"-(control|treated)-",
            "-o",
            "metadata.csv",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, f"stderr:\n{result.stderr}\nstdout:\n{result.stdout}"
    df = pd.read_csv(tmp_path / "metadata.csv")

    assert "condition" in df.columns
    assert set(df["condition"]) == {"control", "treated"}


def test_cli_design_extracts_named_filename_metadata(tmp_path: Path):
    """Named captures become columns and can drive existing grouping options."""
    for name in [
        "SEM-RES5-VTP-Rep1_MLL_R1.fastq.gz",
        "SEM-RES5-VTP-Rep1_MLL_R2.fastq.gz",
        "SEM-STOCK-VTP-Rep2_MLL_R1.fastq.gz",
        "SEM-STOCK-VTP-Rep2_MLL_R2.fastq.gz",
    ]:
        _write_fastq(tmp_path, name)

    pattern = (
        r"^(?P<cell_line>SEM-(?:RES5|STOCK))-(?P<condition>VTP)-"
        r"(?P<replicate>Rep\d+)_(?P<target>[^_]+)_R[12]\.fastq\.gz$"
    )
    result = subprocess.run(
        [
            "seqnado",
            "design",
            "chip",
            "--no-interactive",
            "--accept-all-defaults",
            "--metadata-regex",
            pattern,
            "--consensus-by",
            "target",
            "-o",
            "metadata.csv",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    df = pd.read_csv(tmp_path / "metadata.csv")
    assert set(df["cell_line"]) == {"SEM-RES5", "SEM-STOCK"}
    assert set(df["replicate"]) == {"Rep1", "Rep2"}
    assert set(df["condition"]) == {"VTP"}
    assert set(df["target"]) == {"MLL"}
    assert (df["consensus_group"] == df["target"]).all()


def test_cli_design_accepts_pcre_named_capture_syntax(tmp_path: Path):
    _write_fastq(tmp_path, "SEM-RES5-VTP-Rep1_MLL_R1.fastq.gz")
    _write_fastq(tmp_path, "SEM-RES5-VTP-Rep1_MLL_R2.fastq.gz")
    result = subprocess.run(
        [
            "seqnado",
            "design",
            "chip",
            "--no-interactive",
            "--accept-all-defaults",
            "--metadata-regex",
            r"^(?<cell_line>SEM)-(?<state>RES5)-(?<treatment>VTP)-Rep\d+_(?<target>[^_]+)_R[12]\.fastq\.gz$",
            "-o",
            "metadata.csv",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    df = pd.read_csv(tmp_path / "metadata.csv")
    assert df.loc[0, "cell_line"] == "SEM"
    assert df.loc[0, "state"] == "RES5"
    assert df.loc[0, "treatment"] == "VTP"
    assert df.loc[0, "target"] == "MLL"


def test_cli_design_extracts_explicit_field_and_keeps_unmatched_blank(tmp_path: Path):
    _write_fastq(tmp_path, "SEM-RES5-VTP-Rep1_R1.fastq.gz")
    _write_fastq(tmp_path, "SEM-RES5-VTP-Rep1_R2.fastq.gz")
    _write_fastq(tmp_path, "other-sample_R1.fastq.gz")
    _write_fastq(tmp_path, "other-sample_R2.fastq.gz")

    result = subprocess.run(
        [
            "seqnado",
            "design",
            "atac",
            "--no-interactive",
            "--accept-all-defaults",
            "--metadata-field",
            r"cell_line=^SEM-(RES5|STOCK)-",
            "-o",
            "metadata.csv",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    df = pd.read_csv(tmp_path / "metadata.csv")
    assert set(df.loc[df["cell_line"].notna(), "cell_line"]) == {"RES5"}
    assert df["cell_line"].isna().sum() == 1


@pytest.mark.parametrize(
    ("options", "error"),
    [
        (["--metadata-regex", "["], "Invalid --metadata-regex"),
        (["--metadata-regex", "(not_named)"], "named capture group"),
        (["--metadata-field", "batch=(one)(two)"], "exactly one capture group"),
        (
            ["--metadata-regex", "(?P<batch>one)", "--metadata-field", "batch=(two)"],
            "Duplicate metadata destination",
        ),
        (["--metadata-regex", "(?P<r1>sample)"], "protected FASTQ path"),
    ],
)
def test_cli_design_rejects_invalid_metadata_extractors(tmp_path: Path, options, error):
    _write_fastq(tmp_path, "sample_R1.fastq.gz")
    _write_fastq(tmp_path, "sample_R2.fastq.gz")
    result = subprocess.run(
        [
            "seqnado",
            "design",
            "atac",
            "--no-interactive",
            "--accept-all-defaults",
            *options,
            "-o",
            "metadata.csv",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 2
    assert error in result.stderr
