"""Unit tests for seqnado.analysis.SeqNadoProject."""

import pytest
import pandas as pd
from pathlib import Path

from seqnado.analysis import SeqNadoProject
from seqnado.analysis.project import _parse_bigwig_path, _parse_peak_path


# ---------------------------------------------------------------------------
# Path-parsing helpers
# ---------------------------------------------------------------------------

BIGWIG_DIR = Path("seqnado_output/bigwigs")


@pytest.mark.parametrize("rel, expected", [
    (
        "deeptools/unscaled/sample1.bigWig",
        {"sample": "sample1", "method": "deeptools", "scale": "unscaled",
         "spikein_method": None, "merged": False, "strand": None},
    ),
    (
        "deeptools/csaw/sample1.bigWig",
        {"sample": "sample1", "method": "deeptools", "scale": "csaw",
         "spikein_method": None, "merged": False, "strand": None},
    ),
    (
        "deeptools/merged/unscaled/group1.bigWig",
        {"sample": "group1", "method": "deeptools", "scale": "unscaled",
         "spikein_method": None, "merged": True, "strand": None},
    ),
    (
        "deeptools/spikein/orlando/sample1.bigWig",
        {"sample": "sample1", "method": "deeptools", "scale": "spikein",
         "spikein_method": "orlando", "merged": False, "strand": None},
    ),
    (
        "deeptools/merged/spikein/orlando/group1.bigWig",
        {"sample": "group1", "method": "deeptools", "scale": "spikein",
         "spikein_method": "orlando", "merged": True, "strand": None},
    ),
    (
        "deeptools/unscaled/sample1_plus.bigWig",
        {"sample": "sample1", "method": "deeptools", "scale": "unscaled",
         "spikein_method": None, "merged": False, "strand": "plus"},
    ),
    (
        "deeptools/unscaled/sample1_minus.bigWig",
        {"sample": "sample1", "method": "deeptools", "scale": "unscaled",
         "spikein_method": None, "merged": False, "strand": "minus"},
    ),
    (
        "deeptools/aggregated/treated.bigWig",
        {"sample": "treated", "method": "deeptools", "scale": "aggregated",
         "spikein_method": None, "merged": False, "strand": None},
    ),
    (
        "taps/sample1_hg38.bigWig",
        {"sample": "sample1_hg38", "method": "taps", "scale": "unscaled",
         "spikein_method": None, "merged": False, "strand": None},
    ),
])
@pytest.mark.unit
def test_parse_bigwig_path(rel, expected):
    p = BIGWIG_DIR / rel
    result = _parse_bigwig_path(p, BIGWIG_DIR)
    assert result is not None
    for key, val in expected.items():
        assert result[key] == val, f"{key}: got {result[key]!r}, want {val!r}"


@pytest.mark.unit
def test_parse_peak_path():
    peaks_dir = Path("seqnado_output/peaks")
    p = peaks_dir / "macs2" / "sample1.bed"
    rec = _parse_peak_path(p, peaks_dir)
    assert rec["method"] == "macs2"
    assert rec["sample"] == "sample1"
    assert rec["merged"] is False

    p_merged = peaks_dir / "macs2" / "merged" / "group1.bed"
    rec_m = _parse_peak_path(p_merged, peaks_dir)
    assert rec_m["merged"] is True
    assert rec_m["sample"] == "group1"


# ---------------------------------------------------------------------------
# SeqNadoProject with fake output directory
# ---------------------------------------------------------------------------

@pytest.fixture()
def fake_project(tmp_path):
    """Create a minimal fake SeqNado output directory."""
    out = tmp_path / "seqnado_output"

    # Aligned BAMs
    bam_dir = out / "aligned"
    bam_dir.mkdir(parents=True)
    for name in ["sampleA_H3K27ac", "sampleA_Input", "sampleB_H3K27ac"]:
        (bam_dir / f"{name}.bam").touch()
        (bam_dir / f"{name}.bam.bai").touch()

    # BigWigs
    bw_base = out / "bigwigs" / "deeptools"
    (bw_base / "unscaled").mkdir(parents=True)
    (bw_base / "unscaled" / "sampleA_H3K27ac.bigWig").touch()
    (bw_base / "unscaled" / "sampleB_H3K27ac.bigWig").touch()
    (bw_base / "csaw").mkdir(parents=True)
    (bw_base / "csaw" / "sampleA_H3K27ac.bigWig").touch()
    (bw_base / "merged" / "unscaled").mkdir(parents=True)
    (bw_base / "merged" / "unscaled" / "H3K27ac_merged.bigWig").touch()

    # Peaks
    pk_base = out / "peaks" / "macs2"
    pk_base.mkdir(parents=True)
    (pk_base / "sampleA_H3K27ac.bed").touch()
    (pk_base / "sampleB_H3K27ac.bed").touch()
    (pk_base / "merged").mkdir()
    (pk_base / "merged" / "consensus.bed").touch()

    # Counts
    counts_dir = out / "readcounts" / "feature_counts"
    counts_dir.mkdir(parents=True)
    counts_path = counts_dir / "read_counts.tsv"
    df = pd.DataFrame({"sampleA_H3K27ac": [10, 20], "sampleB_H3K27ac": [15, 25]},
                      index=["geneA", "geneB"])
    df.to_csv(counts_path, sep="\t")

    # Design file
    design = pd.DataFrame({
        "sample_id": ["sampleA", "sampleB"],
        "condition": ["treated", "control"],
        "ip": ["H3K27ac", "H3K27ac"],
        "scaling_group": ["grp1", "grp1"],
        "r1": ["sampleA_R1.fq.gz", "sampleB_R1.fq.gz"],
        "r2": ["sampleA_R2.fq.gz", "sampleB_R2.fq.gz"],
    })
    design.to_csv(tmp_path / "design.tsv", sep="\t", index=False)

    return SeqNadoProject(out, design=tmp_path / "design.tsv")


@pytest.mark.unit
def test_samples(fake_project):
    # Inferred from BAM files
    assert "sampleA_H3K27ac" in fake_project.samples


@pytest.mark.unit
def test_conditions(fake_project):
    assert set(fake_project.conditions) == {"treated", "control"}


@pytest.mark.unit
def test_antibodies(fake_project):
    assert "H3K27ac" in fake_project.antibodies


@pytest.mark.unit
def test_bigwigs_all(fake_project):
    bws = fake_project.bigwigs()
    assert len(bws) == 4  # 2 unscaled + 1 csaw + 1 merged


@pytest.mark.unit
def test_bigwigs_filter_scale(fake_project):
    bws = fake_project.bigwigs(scale="unscaled", merged=False)
    assert len(bws) == 2
    assert all("unscaled" in str(p) for p in bws)


@pytest.mark.unit
def test_bigwigs_filter_merged(fake_project):
    bws = fake_project.bigwigs(merged=True)
    assert len(bws) == 1
    assert "merged" in str(bws[0])


@pytest.mark.unit
def test_bigwigs_filter_condition(fake_project):
    bws = fake_project.bigwigs(condition="treated", scale="unscaled", merged=False)
    assert all("sampleA" in str(p) for p in bws)
    assert not any("sampleB" in str(p) for p in bws)


@pytest.mark.unit
def test_peaks_all(fake_project):
    pks = fake_project.peaks()
    assert len(pks) == 3


@pytest.mark.unit
def test_peaks_merged(fake_project):
    pks = fake_project.peaks(merged=True)
    assert len(pks) == 1
    assert pks[0].name == "consensus.bed"


@pytest.mark.unit
def test_bams(fake_project):
    bams = fake_project.bams()
    assert len(bams) == 3
    assert all(p.suffix == ".bam" for p in bams)


@pytest.mark.unit
def test_bams_filter_antibody(fake_project):
    bams = fake_project.bams(antibody="H3K27ac")
    assert all("H3K27ac" in p.name for p in bams)
    assert not any("Input" in p.name for p in bams)


@pytest.mark.unit
def test_counts_path(fake_project):
    p = fake_project.counts()
    assert p is not None
    assert p.exists()


@pytest.mark.unit
def test_load_counts(fake_project):
    df = fake_project.load_counts()
    assert "sampleA_H3K27ac" in df.columns
    assert "geneA" in df.index


@pytest.mark.unit
def test_filter_view(fake_project):
    view = fake_project.filter(condition="treated")
    assert "sampleA_H3K27ac" in view.samples
    assert "sampleB_H3K27ac" not in view.samples

    bws = view.bigwigs(scale="unscaled", merged=False)
    assert all("sampleA" in str(p) for p in bws)


@pytest.mark.unit
def test_filter_chained(fake_project):
    view = fake_project.filter(condition="treated").filter(antibody="H3K27ac")
    assert len(view.samples) >= 1


@pytest.mark.unit
def test_summary(fake_project):
    df = fake_project.summary()
    assert "File type" in df.columns
    assert "Present" in df.columns


@pytest.mark.unit
def test_repr(fake_project):
    r = repr(fake_project)
    assert "SeqNadoProject" in r
