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
    bed_content = "chr1\t1000\t2000\tpeak1\t100\t.\n"
    (pk_base / "sampleA_H3K27ac.bed").write_text(bed_content)
    (pk_base / "sampleB_H3K27ac.bed").write_text(bed_content)
    (pk_base / "merged").mkdir()
    (pk_base / "merged" / "consensus.bed").write_text(bed_content)

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


# ---------------------------------------------------------------------------
# Enum params
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_bigwigs_enum_method(fake_project):
    from seqnado.core import PileupMethod, DataScalingTechnique
    bws_str  = fake_project.bigwigs(method="deeptools", scale="unscaled", merged=False)
    bws_enum = fake_project.bigwigs(method=PileupMethod.DEEPTOOLS,
                                    scale=DataScalingTechnique.UNSCALED, merged=False)
    assert bws_str == bws_enum


@pytest.mark.unit
def test_peaks_enum_method(fake_project):
    from seqnado.core import PeakCallingMethod
    pks_str  = fake_project.peaks(method="macs2")
    pks_enum = fake_project.peaks(method=PeakCallingMethod.MACS2)
    assert pks_str == pks_enum


# ---------------------------------------------------------------------------
# Discovery helpers
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_pileup_methods(fake_project):
    assert "deeptools" in fake_project.pileup_methods


@pytest.mark.unit
def test_scales(fake_project):
    assert "unscaled" in fake_project.scales
    assert "csaw" in fake_project.scales


@pytest.mark.unit
def test_peak_methods(fake_project):
    assert "macs2" in fake_project.peak_methods


@pytest.mark.unit
def test_consensus_groups(fake_project):
    # merged bigWig and peak created in fixture
    groups = fake_project.consensus_groups
    assert len(groups) >= 1


@pytest.mark.unit
def test_samples_for_condition(fake_project):
    treated = fake_project.samples_for(condition="treated")
    assert all("sampleA" in s for s in treated)
    assert not any("sampleB" in s for s in treated)


@pytest.mark.unit
def test_samples_for_antibody(fake_project):
    ab = fake_project.samples_for(antibody="H3K27ac")
    assert len(ab) >= 1
    assert all("H3K27ac" in s for s in ab)


@pytest.mark.unit
def test_metadata_for(fake_project):
    meta = fake_project.metadata_for("sampleA_H3K27ac")
    assert meta.get("condition") == "treated"
    assert meta.get("antibody") == "H3K27ac"


@pytest.mark.unit
def test_metadata_for_missing(fake_project):
    assert fake_project.metadata_for("nonexistent") == {}


@pytest.mark.unit
def test_what_exists(fake_project):
    df = fake_project.what_exists()
    assert "sample" in df.columns
    assert "bam" in df.columns
    assert "bigwig" in df.columns
    assert "peaks" in df.columns
    # at least one sample with BAM present
    assert df["bam"].any()


@pytest.mark.unit
def test_condition_pairs(fake_project):
    pairs = fake_project.condition_pairs()
    assert ("treated", "control") in pairs
    assert ("control", "treated") in pairs


# ---------------------------------------------------------------------------
# New loaders (smoke test with fake data)
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_load_peaks(fake_project):
    df = fake_project.load_peaks()
    assert not df.empty
    assert "chrom" in df.columns
    assert "sample" in df.columns
    assert "peak_method" in df.columns


@pytest.mark.unit
def test_load_peaks_filtered(fake_project):
    df = fake_project.load_peaks(condition="treated")
    assert all("sampleA" in r for r in df["sample"])


@pytest.mark.unit
def test_load_alignment_stats(fake_project, tmp_path):
    # Create fake alignment_stats.tsv
    p = fake_project.output_dir / "qc" / "alignment_stats.tsv"
    p.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"sample": ["sampleA_H3K27ac"], "Retained": [1000]}).to_csv(p, sep="\t", index=False)
    df = fake_project.load_alignment_stats()
    assert "sample" in df.columns


@pytest.mark.unit
def test_logs(fake_project):
    log_dir = fake_project.output_dir / "logs" / "bowtie2"
    log_dir.mkdir(parents=True, exist_ok=True)
    (log_dir / "sampleA_H3K27ac.log").touch()
    logs = fake_project.logs(rule="bowtie2")
    assert len(logs) == 1
    assert fake_project.log_rules == ["bowtie2"]


@pytest.mark.unit
def test_load_benchmarks(fake_project):
    bench_dir = fake_project.output_dir / ".benchmark" / "bowtie2"
    bench_dir.mkdir(parents=True, exist_ok=True)
    bdf = pd.DataFrame({"s": [1.2], "h:m:s": ["0:00:01"], "max_rss": [100.0]})
    bdf.to_csv(bench_dir / "sampleA_H3K27ac.tsv", sep="\t", index=False)
    df = fake_project.load_benchmarks(rule="bowtie2")
    assert "rule" in df.columns
    assert "sample" in df.columns
    assert df["rule"].iloc[0] == "bowtie2"
    assert fake_project.benchmark_rules == ["bowtie2"]


@pytest.mark.unit
def test_protocol(fake_project):
    p = fake_project.output_dir / "protocol.txt"
    p.write_text("Methods: ...")
    assert fake_project.protocol() == p


@pytest.mark.unit
def test_bai_for(fake_project):
    bam = fake_project.bams()[0]
    bai = Path(str(bam) + ".bai")
    assert fake_project.bai_for(bam) == bai


@pytest.mark.unit
def test_sample_files(fake_project):
    result = fake_project.sample_files("sampleA_H3K27ac")
    assert "bam" in result
    assert "bigwigs" in result
    assert "peaks" in result


# ---------------------------------------------------------------------------
# New gaps: reload, select_files, enrich, load_deseq2, load_mageck,
#           _FilteredProject discovery props, summary(detail=True)
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_reload_clears_cache(fake_project):
    _ = fake_project.samples          # populate cache
    _ = fake_project.pileup_methods
    fake_project.reload()
    # After reload the cached values should be gone from __dict__
    assert "samples" not in fake_project.__dict__
    assert "_bw_index" not in fake_project.__dict__
    # Re-access still works
    assert "sampleA_H3K27ac" in fake_project.samples


@pytest.mark.unit
def test_select_files_non_recursive(fake_project):
    # recursive=False → glob only at root, no BAMs there
    files = fake_project.select_files("*.bam", recursive=False)
    assert len(files) == 0


@pytest.mark.unit
def test_select_files_recursive(fake_project):
    # recursive=True (default) → rglob finds all 3 BAMs
    files = fake_project.select_files("*.bam")
    assert len(files) == 3


@pytest.mark.unit
def test_select_files_subdir_pattern(fake_project):
    # explicit subdir pattern without **
    files = fake_project.select_files("aligned/*.bam", recursive=False)
    assert len(files) == 3


@pytest.mark.unit
def test_enrich(fake_project):
    bws = fake_project.bigwigs(scale="unscaled", merged=False)
    df = fake_project.enrich(bws)
    assert "path" in df.columns
    assert "sample" in df.columns
    assert "condition" in df.columns
    assert len(df) == 2
    row_a = df[df["sample"] == "sampleA_H3K27ac"].iloc[0]
    assert row_a["condition"] == "treated"


@pytest.mark.unit
def test_enrich_strand_stripped(fake_project):
    bw_dir = fake_project.output_dir / "bigwigs" / "deeptools" / "unscaled"
    strand_bw = bw_dir / "sampleA_H3K27ac_plus.bigWig"
    strand_bw.touch()
    fake_project.reload()
    df = fake_project.enrich([strand_bw])
    assert df.iloc[0]["sample"] == "sampleA_H3K27ac"
    assert df.iloc[0]["condition"] == "treated"


@pytest.mark.unit
def test_load_deseq2(fake_project):
    de_dir = fake_project.output_dir / "deseq2_results"
    de_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": ["geneA", "geneB"], "log2FoldChange": [1.0, -1.0]}).to_csv(
        de_dir / "DEseq2_treated_vs_control.csv", index=False
    )
    df = fake_project.load_deseq2()
    assert "contrast" in df.columns
    assert df["contrast"].iloc[0] == "treated_vs_control"


@pytest.mark.unit
def test_load_deseq2_single_contrast(fake_project):
    de_dir = fake_project.output_dir / "deseq2_results"
    de_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": ["geneA"], "log2FoldChange": [2.0]}).to_csv(
        de_dir / "DEseq2_treated_vs_control.csv", index=False
    )
    df = fake_project.load_deseq2(contrast="treated_vs_control")
    assert len(df) == 1
    assert df["contrast"].iloc[0] == "treated_vs_control"


@pytest.mark.unit
def test_load_deseq2_not_found(fake_project):
    with pytest.raises(FileNotFoundError):
        fake_project.load_deseq2(contrast="nonexistent")


@pytest.mark.unit
def test_load_mageck(fake_project):
    mageck_dir = fake_project.output_dir / "readcounts" / "mageck"
    mageck_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"id": ["gene1"], "num": [5], "neg|score": [0.01]}).to_csv(
        mageck_dir / "mageck_test.gene_summary.txt", sep="\t", index=False
    )
    df = fake_project.load_mageck()
    assert "id" in df.columns


@pytest.mark.unit
def test_load_mageck_not_found(fake_project):
    with pytest.raises(FileNotFoundError):
        fake_project.load_mageck(analysis="mle", type="sgrna")


@pytest.mark.unit
def test_filtered_project_conditions(fake_project):
    view = fake_project.filter(condition="treated")
    assert view.conditions == ["treated"]
    assert "control" not in view.conditions


@pytest.mark.unit
def test_filtered_project_antibodies(fake_project):
    view = fake_project.filter(condition="treated")
    assert "H3K27ac" in view.antibodies


@pytest.mark.unit
def test_filtered_project_pileup_methods(fake_project):
    view = fake_project.filter(condition="treated")
    assert "deeptools" in view.pileup_methods


@pytest.mark.unit
def test_filtered_project_scales(fake_project):
    view = fake_project.filter(condition="treated")
    assert "unscaled" in view.scales


@pytest.mark.unit
def test_filtered_project_peak_methods(fake_project):
    view = fake_project.filter(condition="treated")
    assert "macs2" in view.peak_methods


@pytest.mark.unit
def test_summary_detail(fake_project):
    df = fake_project.summary(detail=True)
    assert "Detail" in df.columns
    assert "Count" in df.columns
    bw_rows = df[df["File type"] == "BigWig"]
    assert len(bw_rows) > 0
    assert all(bw_rows["Present"] == "✓")


@pytest.mark.unit
def test_summary_simple(fake_project):
    df = fake_project.summary()
    assert "File type" in df.columns
    assert "Present" in df.columns
    assert "Detail" not in df.columns
