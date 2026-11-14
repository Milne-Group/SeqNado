"""Tests for seqnado.workflow.helpers.spikein contig classification."""

import shutil
import subprocess
from pathlib import Path

import pytest

from seqnado import SpikeInMethod
from seqnado.config.configs import SpikeInConfig
from seqnado.workflow.helpers.spikein import (
    DEFAULT_REFERENCE_CONTIG_REGEX,
    classify_contigs,
    ere_escape,
    reference_contig_regex,
    spikein_contig_regex,
    spikein_regexes_from_config,
)

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    "contigs, ref_prefix, exo_prefix, expected",
    [
        pytest.param(
            ["chr1", "chrX", "chr1_KI270706v1_random", "chrM", "dm6_chr2L", "dm6_chrX"],
            None,
            ["dm6"],
            (
                ["chr1", "chrX", "chr1_KI270706v1_random", "chrM"],
                ["dm6_chr2L", "dm6_chrX"],
                [],
            ),
            id="ucsc-hg38-with-dm6",
        ),
        pytest.param(
            ["1", "22", "X", "Y", "MT", "GL000220.1", "dm6_2L", "dm6_X"],
            None,
            ["dm6"],
            (["1", "22", "X", "Y", "MT"], ["dm6_2L", "dm6_X"], ["GL000220.1"]),
            id="ensembl-primary-with-dm6",
        ),
        pytest.param(
            ["NC_000001.11", "NC_000023.11", "spike_lambda"],
            ["NC_"],
            ["spike"],
            (["NC_000001.11", "NC_000023.11"], ["spike_lambda"], []),
            id="custom-endogenous-prefix",
        ),
        pytest.param(
            ["chr21", "AmpR", "Cas9_3p", "Cas9_5p"],
            None,
            ["AmpR", "Cas9_3p", "Cas9_5p"],
            (["chr21"], ["AmpR", "Cas9_3p", "Cas9_5p"], []),
            id="explicit-contig-list",
        ),
        pytest.param(
            ["chr21", "Lambda", "2kb-unmod", "Lambda_NEB"],
            None,
            ["Lambda", "2kb-unmod"],
            (["chr21"], ["Lambda", "2kb-unmod", "Lambda_NEB"], []),
            id="exact-and-underscore-prefixed",
        ),
    ],
)
def test_classify_contigs(contigs, ref_prefix, exo_prefix, expected):
    ref_re = reference_contig_regex(ref_prefix)
    exo_re = spikein_contig_regex(exo_prefix)
    assert classify_contigs(contigs, ref_re, exo_re) == expected


@pytest.mark.parametrize(
    "exo_prefix, contig",
    [
        ("dm6", "dm6x_chr1"),  # no delimiter
        ("dm6", "dm61"),
        ("ERCC.v2", "ERCCxv2"),  # regex metachar is literal
        ("a+b", "aab"),
    ],
)
def test_spikein_regex_does_not_overmatch(exo_prefix, contig):
    _, exo, _ = classify_contigs([contig], reference_contig_regex(), spikein_contig_regex([exo_prefix]))
    assert exo == []


def test_spikein_regex_takes_precedence_over_reference():
    # A spike-in prefixed with "chr" must still be classified as spike-in.
    ref, exo, _ = classify_contigs(
        ["chr1", "chrSpike_1"], reference_contig_regex(), spikein_contig_regex(["chrSpike"])
    )
    assert ref == ["chr1"]
    assert exo == ["chrSpike_1"]


def test_reference_regex_default():
    assert reference_contig_regex(None) == DEFAULT_REFERENCE_CONTIG_REGEX
    assert reference_contig_regex([]) == DEFAULT_REFERENCE_CONTIG_REGEX


def test_spikein_regex_requires_prefix():
    with pytest.raises(ValueError):
        spikein_contig_regex([])


def test_ere_escape():
    assert ere_escape("a.b(c)") == r"a\.b\(c\)"
    assert ere_escape("250bp-v1") == "250bp-v1"


@pytest.mark.parametrize(
    "kwargs, exo_contig",
    [
        ({"exogenous_genome": "dm6"}, "dm6_chr2L"),
        ({"exogenous_genome": "plasmid", "exogenous_contig_prefix": ["AmpR"]}, "AmpR"),
    ],
)
def test_regexes_from_config(kwargs, exo_contig):
    cfg = SpikeInConfig(method=[SpikeInMethod.ORLANDO], **kwargs)
    ref_re, exo_re = spikein_regexes_from_config(cfg)
    assert classify_contigs(["chr1", exo_contig], ref_re, exo_re) == (["chr1"], [exo_contig], [])


def test_regexes_from_config_requires_exogenous():
    cfg = SpikeInConfig(method=[SpikeInMethod.ORLANDO])
    with pytest.raises(ValueError, match="exogenous"):
        spikein_regexes_from_config(cfg)


SCRIPT = (
    Path(__file__).resolve().parents[3]
    / "seqnado" / "workflow" / "scripts" / "split_bam_by_contig.sh"
)

requires_samtools = pytest.mark.skipif(
    shutil.which("samtools") is None, reason="samtools not available"
)

HEADER_CONTIGS = [
    ("1", 1000), ("X", 1000), ("MT", 500), ("GL000220.1", 500),
    ("dm6_2L", 1000), ("dm6_chrX", 1000), ("ERCC.v2", 100), ("ERCCxv2", 100),
]
# (contig, number of reads)
READS = [("1", 2), ("X", 1), ("MT", 1), ("GL000220.1", 1), ("dm6_2L", 1),
         ("dm6_chrX", 2), ("ERCC.v2", 1), ("ERCCxv2", 1)]


@pytest.fixture
def synthetic_bam(tmp_path: Path) -> Path:
    lines = ["@HD\tVN:1.6\tSO:coordinate"]
    lines += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in HEADER_CONTIGS]
    n = 0
    for contig, count in READS:
        for j in range(count):
            n += 1
            lines.append(f"r{n}\t0\t{contig}\t{10 + j}\t60\t4M\t*\t0\t0\tACGT\tIIII")
    sam = tmp_path / "in.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp_path / "in.bam"
    subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    return bam


def run_split(bam: Path, tmp_path: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["bash", str(SCRIPT), "--bam", str(bam), "--log", str(tmp_path / "split.log"), *args],
        capture_output=True,
        text=True,
    )


def contigs_in(bam: Path) -> list[str]:
    out = subprocess.run(
        ["samtools", "view", str(bam)], capture_output=True, text=True, check=True
    ).stdout
    return [line.split("\t")[2] for line in out.splitlines()]


@requires_samtools
@pytest.mark.parametrize(
    "ref_prefix, exo_prefix, expected_ref, expected_exo",
    [
        pytest.param(None, ["dm6"], ["1", "1", "X", "MT"], ["dm6_2L", "dm6_chrX", "dm6_chrX"], id="default-ref-dm6"),
        pytest.param(["X"], ["dm6"], ["X"], ["dm6_2L", "dm6_chrX", "dm6_chrX"], id="custom-ref-prefix"),
        pytest.param(None, ["ERCC.v2"], ["1", "1", "X", "MT"], ["ERCC.v2"], id="metachar-literal"),
    ],
)
def test_script_splits_bam(synthetic_bam, tmp_path, ref_prefix, exo_prefix, expected_ref, expected_exo):
    ref_out, exo_out = tmp_path / "ref.bam", tmp_path / "exo.bam"
    result = run_split(
        synthetic_bam, tmp_path,
        "--ref-regex", reference_contig_regex(ref_prefix),
        "--exo-regex", spikein_contig_regex(exo_prefix),
        "--ref-out", str(ref_out), "--exo-out", str(exo_out),
    )
    assert result.returncode == 0, result.stderr + (tmp_path / "split.log").read_text()
    assert contigs_in(ref_out) == expected_ref
    assert contigs_in(exo_out) == expected_exo
    log = (tmp_path / "split.log").read_text()
    assert "unassigned" in log and "done" in log


@requires_samtools
def test_script_fails_when_no_spikein_contigs(synthetic_bam, tmp_path):
    result = run_split(
        synthetic_bam, tmp_path,
        "--ref-regex", reference_contig_regex(),
        "--exo-regex", spikein_contig_regex(["plasmid"]),
        "--ref-out", str(tmp_path / "ref.bam"), "--exo-out", str(tmp_path / "exo.bam"),
    )
    assert result.returncode != 0
    assert "exogenous_contig_prefix" in result.stderr
    assert "dm6_2L" in result.stderr  # header contigs listed to help the user


@requires_samtools
def test_script_fails_when_no_reference_contigs(synthetic_bam, tmp_path):
    result = run_split(
        synthetic_bam, tmp_path,
        "--ref-regex", reference_contig_regex(["NC_"]),
        "--exo-regex", spikein_contig_regex(["dm6"]),
        "--ref-out", str(tmp_path / "ref.bam"), "--exo-out", str(tmp_path / "exo.bam"),
    )
    assert result.returncode != 0
    assert "endogenous_contig_prefix" in result.stderr


@requires_samtools
def test_script_warn_mode_writes_empty_indexed_bam(synthetic_bam, tmp_path):
    out = tmp_path / "lambda.bam"
    result = run_split(
        synthetic_bam, tmp_path,
        "--exo-regex", spikein_contig_regex(["Lambda"]),
        "--exo-out", str(out), "--on-empty-exo", "warn", "--index",
    )
    assert result.returncode == 0, result.stderr
    assert contigs_in(out) == []
    assert Path(str(out) + ".bai").exists()
    assert "WARNING" in (tmp_path / "split.log").read_text()


@requires_samtools
def test_script_methylation_reference_excludes_spikeins(synthetic_bam, tmp_path):
    """Reference-only output with spike-in exclusion, as used by methylation_bam_splits."""
    out = tmp_path / "ref.bam"
    result = run_split(
        synthetic_bam, tmp_path,
        "--ref-regex", reference_contig_regex(),
        "--exo-regex", spikein_contig_regex(["dm6", "ERCC.v2"]),
        "--ref-out", str(out),
    )
    assert result.returncode == 0, result.stderr
    assert contigs_in(out) == ["1", "1", "X", "MT"]


def test_script_rejects_unknown_argument(tmp_path):
    result = subprocess.run(["bash", str(SCRIPT), "--nope"], capture_output=True, text=True)
    assert result.returncode == 2
