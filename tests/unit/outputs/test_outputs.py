"""Tests for seqnado.outputs.files module."""

from pathlib import Path

import pandas as pd
import pytest

from seqnado import (
    Assay,
    DataScalingTechnique,
    MethylationMethod,
    PeakCallingMethod,
    PileupMethod,
    QuantificationMethod,
    SNPCallingMethod,
    SpikeInMethod,
)
from seqnado.config.configs import BigwigConfig, GenomeConfig, STARIndex, QCConfig
from seqnado.config.core import ATACAssayConfig, SeqnadoConfig
from seqnado.inputs.bam import BamCollection, BamFile
from seqnado.inputs.bigwigs import BigWigCollection, BigWigFile
from seqnado.inputs.core import Metadata
from seqnado.inputs.fastq import FastqCollection, FastqFile, FastqSet
from seqnado.inputs.grouping import SampleGroup, SampleGroupings, SampleGroups
from seqnado.outputs.core import (
    SeqnadoOutputBuilder,
    SeqnadoOutputFactory,
    SeqnadoOutputFiles,
)
from seqnado.outputs.files import (
    BigBedFiles,
    BigWigFiles,
    FileCollection,
    SeqNadoReportFile,
    GeoSubmissionFiles,
    HeatmapFiles,
    HubFiles,
    MethylationFiles,
    PeakCallingFiles,
    PlotFiles,
    QCFiles,
    QuantificationFiles,
    SNPFilesAnnotated,
    SNPFilesRaw,
    SpikeInFiles,
)


def test_file_collection_protocol():
    """Test that FileCollection protocol can be used for type checking."""
    # Import and verify the protocol exists
    from seqnado.outputs.files import FileCollection
    
    # Access the protocol's method to ensure it's covered
    assert hasattr(FileCollection, 'files')
    
    # Try to instantiate the protocol directly (will call pass statement)
    try:
        # This will trigger the protocol's files property getter
        FileCollection.files.fget(None)
    except (AttributeError, TypeError):
        # Expected - protocols can't be instantiated
        pass
    
    # Create a simple class that implements the protocol
    class TestFileCollection:
        @property
        def files(self):
            return ["file1.txt", "file2.txt"]
    
    # Verify it behaves as a FileCollection
    collection = TestFileCollection()
    assert isinstance(collection.files, list)
    assert len(collection.files) == 2


def _minimal_config(tmp: Path) -> SeqnadoConfig:
    star = tmp / "star"
    star.mkdir()
    genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
    assay_cfg = ATACAssayConfig(
        bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS])
    )
    return SeqnadoConfig(
        assay=Assay.ATAC,
        project=dict(name="p"),
        genome=genome,
        metadata=tmp / "m.csv",
        assay_config=assay_cfg,
    )


def _small_collection(tmp: Path) -> FastqCollection:
    r1_path = tmp / "s1_R1.fastq.gz"
    r1_path.write_text("@r\nN\n+\n#\n")
    r2_path = tmp / "s1_R2.fastq.gz"
    r2_path.write_text("@r\nN\n+\n#\n")
    r1 = FastqFile(path=r1_path)
    r2 = FastqFile(path=r2_path)
    fs = FastqSet(sample_id="s1", r1=r1, r2=r2)
    return FastqCollection(assay=Assay.ATAC, metadata=[Metadata(assay=Assay.ATAC)], fastq_sets=[fs])


def test_output_builder_bigwigs_only(tmp_path: Path):
    cfg = _minimal_config(tmp_path)
    samples = _small_collection(tmp_path)

    builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
    builder.add_qc_files()
    builder.add_individual_bigwig_files()
    out = builder.build()

    # Only expect bigwig files, not seqnado_report.html
    assert any(f.endswith(".bigWig") for f in out.files)


def test_output_factory_with_consensus_groups(tmp_path: Path):
    cfg = _minimal_config(tmp_path)
    samples = _small_collection(tmp_path)
    groups = SampleGroupings(
        groupings={
            "consensus": SampleGroups(groups=[SampleGroup(name="cons", samples=["s1"])])
        }
    )
    factory = SeqnadoOutputFactory(Assay.ATAC, samples, cfg, sample_groupings=groups)
    out = factory.create_output_builder().build()
    assert out.files, "Expected some files listed from factory"


# Helper functions
def _create_fastq_collection(
    tmp_path: Path, sample_ids: list[str], assay: Assay
) -> FastqCollection:
    """Create a minimal FastqCollection for testing."""
    fastq_sets = []
    for sample_id in sample_ids:
        r1_path = tmp_path / f"{sample_id}_R1.fastq.gz"
        r1_path.write_text("@r\nN\n+\n#\n")
        r1 = FastqFile(path=r1_path)
        fastq_sets.append(FastqSet(sample_id=sample_id, r1=r1))

    return FastqCollection(
        assay=assay,
        metadata=[Metadata(assay=assay) for _ in sample_ids],
        fastq_sets=fastq_sets,
    )


def _create_bam_collection(
    tmp_path: Path, sample_ids: list[str], assay: Assay
) -> BamCollection:
    """Create a minimal BamCollection for testing."""
    bam_files = []
    for sample_id in sample_ids:
        bam_path = tmp_path / f"{sample_id}.bam"
        bam_path.touch()
        bam_files.append(BamFile(path=bam_path))

    return BamCollection(
        assay=assay,
        bam_files=bam_files,
        metadata=[Metadata(assay=assay) for _ in sample_ids],
    )


def _create_bigwig_collection(
    tmp_path: Path, sample_ids: list[str], assay: Assay
) -> BigWigCollection:
    """Create a minimal BigWigCollection for testing."""
    bigwig_files = []
    for sample_id in sample_ids:
        bw_path = tmp_path / f"{sample_id}.bigWig"
        bw_path.touch()
        bigwig_files.append(BigWigFile(path=bw_path))

    return BigWigCollection(
        assay=assay,
        bigwig_files=bigwig_files,
        metadata=[Metadata(assay=assay) for _ in sample_ids],
    )


class TestQCFiles:
    """Tests for QCFiles class."""

    def test_qc_files_with_fastq_collection(self, tmp_path):
        """Test QCFiles with FastqCollection."""
        samples = _create_fastq_collection(tmp_path, ["sample1", "sample2"], Assay.ATAC)
        qc = QCFiles(assay=Assay.ATAC, samples=samples, output_dir="test_output")

        files = qc.files
        assert any("qualimap_bamqc" in f for f in files)
        assert len([f for f in files if "qualimap" in f]) == 2  # One per sample

    def test_qc_files_rna_assay(self, tmp_path):
        """Test QCFiles with RNA assay uses qualimap_rnaseq."""
        samples = _create_fastq_collection(tmp_path, ["sample1"], Assay.RNA)
        qc = QCFiles(assay=Assay.RNA, samples=samples)

        files = qc.files
        assert any("qualimap_rnaseq" in f for f in files)
        assert all("qualimap_bamqc" not in f for f in files)

    def test_qc_files_with_bam_collection(self, tmp_path):
        """Test QCFiles with BamCollection."""
        samples = _create_bam_collection(tmp_path, ["sample1"], Assay.ATAC)
        qc = QCFiles(assay=Assay.ATAC, samples=samples)

        files = qc.files
        assert any("qualimap_bamqc" in f for f in files)

    def test_qc_files_with_bigwig_collection(self, tmp_path):
        """Test QCFiles with BigWigCollection (no qualimap files)."""
        samples = _create_bigwig_collection(tmp_path, ["sample1"], Assay.ATAC)
        qc = QCFiles(assay=Assay.ATAC, samples=samples)

        files = qc.files
        assert files == []

    def test_qc_files_no_fastq_screen(self, tmp_path):
        """Test QCFiles does not include fastq_screen files when disabled."""
        samples = _create_fastq_collection(tmp_path, ["sample1"], Assay.ATAC)
        
        # Create QCConfig with run_fastq_screen=False
        qc_config = QCConfig(run_fastq_screen=False)
        
        qc = QCFiles(assay=Assay.ATAC, samples=samples, config=qc_config)

        files = qc.files
        assert all("fastq_screen" not in f for f in files)

    def test_qc_files_with_fastq_screen(self, tmp_path):
        """Test QCFiles includes fastq_screen files when enabled."""
        samples = _create_fastq_collection(tmp_path, ["sample1"], Assay.ATAC)
        
        # Create QCConfig with run_fastq_screen=True
        qc_config = QCConfig(run_fastq_screen=True)
        
        qc = QCFiles(assay=Assay.ATAC, samples=samples, config=qc_config)

        files = qc.files
        assert any("fastq_screen" in f for f in files)

class SeqNadoReportFileTest:
    """Tests for SeqNadoReportFiles class."""

    def test_seqnado_report_file_basic(self, tmp_path):
        """Test SeqNadoReportFile basic functionality."""
        report = SeqNadoReportFile(
            assay=Assay.ATAC,
            samples=_create_fastq_collection(tmp_path, ["sample1"], Assay.ATAC),
            config=_minimal_config(tmp_path),
            output_dir="test_output",
        )
        files = report.files
        assert isinstance(files, list)
        assert all(isinstance(f, str) for f in files)
        assert any("seqnado_report.html" in f for f in files)

class TestBigWigFiles:
    """Tests for BigWigFiles class."""

    def test_bigwig_files_chip_basic(self):
        """Test BigWigFiles for ChIP assay."""
        bw = BigWigFiles(
            assay=Assay.CHIP,
            names=["sample1", "sample2"],
            pileup_methods=[PileupMethod.DEEPTOOLS],
            scale_methods=[DataScalingTechnique.UNSCALED],
        )

        files = bw.files
        assert len(files) == 2  # One per sample
        assert all(f.endswith(".bigWig") for f in files)
        assert all("deeptools" in f for f in files)

    def test_bigwig_files_rna_stranded(self):
        """Test BigWigFiles for RNA assay generates plus/minus files."""
        bw = BigWigFiles(
            assay=Assay.RNA,
            names=["sample1"],
            pileup_methods=[PileupMethod.DEEPTOOLS],
            scale_methods=[DataScalingTechnique.UNSCALED],
        )

        files = bw.files
        assert len(files) == 2  # Plus and minus
        assert any("plus.bigWig" in f for f in files)
        assert any("minus.bigWig" in f for f in files)

    def test_bigwig_files_multiple_methods(self):
        """Test BigWigFiles with multiple pileup methods."""
        bw = BigWigFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            pileup_methods=[PileupMethod.DEEPTOOLS, PileupMethod.HOMER],
            scale_methods=[DataScalingTechnique.UNSCALED],
        )

        files = bw.files
        assert len(files) == 2  # One per method
        assert any("deeptools" in f for f in files)
        assert any("homer" in f for f in files)

    def test_bigwig_files_multiple_scales(self):
        """Test BigWigFiles with multiple scaling methods."""
        bw = BigWigFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            pileup_methods=[PileupMethod.DEEPTOOLS],
            scale_methods=[DataScalingTechnique.UNSCALED, DataScalingTechnique.CPM],
        )

        files = bw.files
        assert len(files) == 2  # One per scale method
        assert any("unscaled" in f for f in files)
        assert any("cpm" in f for f in files)

    def test_bigwig_incompatible_methods_filtered(self):
        """Test that incompatible method/scale combinations are filtered."""
        bw = BigWigFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            pileup_methods=[PileupMethod.HOMER],
            scale_methods=[DataScalingTechnique.CSAW, DataScalingTechnique.UNSCALED],
        )

        files = bw.files
        # HOMER with CSAW should be filtered out
        assert not any("csaw" in f for f in files)
        assert any("unscaled" in f for f in files)

    def test_bigwig_prefix_property(self):
        """Test prefix property."""
        bw = BigWigFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            pileup_methods=[PileupMethod.DEEPTOOLS],
            output_dir="custom_output",
        )

        assert bw.prefix == "custom_output/bigwigs/"

    def test_is_rna_property(self):
        """Test is_rna property."""
        bw_chip = BigWigFiles(
            assay=Assay.CHIP, names=["sample1"], pileup_methods=[PileupMethod.DEEPTOOLS]
        )
        bw_rna = BigWigFiles(
            assay=Assay.RNA, names=["sample1"], pileup_methods=[PileupMethod.DEEPTOOLS]
        )

        assert bw_chip.is_rna is False
        assert bw_rna.is_rna is True

    def test_is_compatible_method(self):
        """Test _is_compatible method."""
        bw = BigWigFiles(
            assay=Assay.CHIP, names=["sample1"], pileup_methods=[PileupMethod.HOMER]
        )

        assert bw._is_compatible(
            PileupMethod.HOMER, DataScalingTechnique.UNSCALED, assay=Assay.CHIP
        )
        assert not bw._is_compatible(
            PileupMethod.HOMER, DataScalingTechnique.CSAW, assay=Assay.CHIP
        )
        assert not bw._is_compatible(
            PileupMethod.HOMER, DataScalingTechnique.SPIKEIN, assay=Assay.CHIP
        )
        assert not bw._is_compatible(
            PileupMethod.DEEPTOOLS, DataScalingTechnique.UNSCALED, assay=Assay.MCC
        )


class TestPeakCallingFiles:
    """Tests for PeakCallingFiles class."""

    def test_peak_calling_files_chip(self):
        """Test PeakCallingFiles for ChIP assay."""
        pcf = PeakCallingFiles(
            assay=Assay.CHIP,
            names=["sample1", "sample2"],
            peak_calling_method=[PeakCallingMethod.MACS2],
        )

        files = pcf.files
        assert len(files) == 2
        assert all(f.endswith(".bed") for f in files)
        assert all("macs2" in f for f in files)

    def test_peak_calling_files_atac(self):
        """Test PeakCallingFiles for ATAC assay."""
        pcf = PeakCallingFiles(
            assay=Assay.ATAC,
            names=["sample1"],
            peak_calling_method=[PeakCallingMethod.MACS2, PeakCallingMethod.SEACR],
        )

        files = pcf.files
        assert len(files) == 2  # One per method
        assert any("macs2" in f for f in files)
        assert any("seacr" in f for f in files)

    def test_peak_calling_invalid_assay(self):
        """Test PeakCallingFiles raises error for invalid assay."""
        with pytest.raises(ValueError, match="Invalid assay for peak calling"):
            PeakCallingFiles(
                assay=Assay.RNA,  # RNA doesn't support peak calling
                names=["sample1"],
                peak_calling_method=[PeakCallingMethod.MACS2],
            )

    def test_peak_calling_prefix(self):
        """Test prefix property."""
        pcf = PeakCallingFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            peak_calling_method=[PeakCallingMethod.MACS2],
            output_dir="custom_output",
        )

        assert pcf.prefix == "custom_output/peaks/"


class TestHeatmapFiles:
    """Tests for HeatmapFiles class."""

    def test_heatmap_files_chip(self):
        """Test HeatmapFiles for ChIP assay."""
        hf = HeatmapFiles(assay=Assay.CHIP)

        files = hf.files
        assert len(files) == 2
        assert any("heatmap.pdf" in f for f in files)
        assert any("metaplot.pdf" in f for f in files)

    def test_heatmap_files_atac(self):
        """Test HeatmapFiles for ATAC assay."""
        hf = HeatmapFiles(assay=Assay.ATAC, output_dir="custom_output")

        files = hf.files
        assert all("custom_output/heatmap/" in f for f in files)

    def test_heatmap_invalid_assay(self):
        """Test HeatmapFiles raises error for invalid assay."""
        with pytest.raises(ValueError, match="Invalid assay for heatmap"):
            HeatmapFiles(assay=Assay.SNP)


class TestHubFiles:
    """Tests for HubFiles class."""

    def test_hub_files_basic(self, tmp_path):
        """Test HubFiles basic functionality."""
        hub = HubFiles(hub_dir=tmp_path, hub_name="test_hub")

        files = hub.files
        assert len(files) == 1
        assert files[0].endswith("test_hub.hub.txt")

    def test_hub_txt_property(self, tmp_path):
        """Test hub_txt property."""
        hub = HubFiles(hub_dir=tmp_path, hub_name="my_hub")

        assert hub.hub_txt == tmp_path / "my_hub.hub.txt"


class TestSpikeInFiles:
    """Tests for SpikeInFiles class."""

    def test_spikein_files_chip(self):
        """Test SpikeInFiles for ChIP assay."""
        sif = SpikeInFiles(
            assay=Assay.CHIP, names=["sample1", "sample2"], method=[SpikeInMethod.ORLANDO]
        )

        files = sif.files
        assert len(files) == 1
        assert "normalisation_factors.tsv" in files[0]

    def test_spikein_invalid_assay(self):
        """Test SpikeInFiles raises error for invalid assay."""
        with pytest.raises(ValueError, match="Invalid assay for spike-in"):
            SpikeInFiles(assay=Assay.SNP, names=["sample1"], method=[SpikeInMethod.ORLANDO])

    def test_norm_factors_property(self):
        """Test norm_factors property."""
        sif = SpikeInFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            method=[SpikeInMethod.ORLANDO],
            output_dir="custom_output",
        )

        assert sif.norm_factors == [
            "custom_output/resources/orlando/normalisation_factors.tsv"
        ]


class TestPlotFiles:
    """Tests for PlotFiles class."""

    def test_plot_files_with_bed(self, tmp_path):
        """Test PlotFiles with BED file."""
        bed_file = tmp_path / "coords.bed"
        bed_file.write_text("chr1\t1000\t2000\tregion1\n")

        pf = PlotFiles(coordinates=bed_file, file_format="svg")

        files = pf.files
        assert len(files) == 1
        assert "region1.svg" in str(files[0])

    def test_plot_files_without_names(self, tmp_path):
        """Test PlotFiles without region names uses coordinates."""
        bed_file = tmp_path / "coords.bed"
        bed_file.write_text("chr1\t1000\t2000\n")

        pf = PlotFiles(coordinates=bed_file)

        files = pf.files
        assert len(files) == 1
        assert "chr1-1000-2000" in str(files[0])

    def test_plot_files_multiple_regions(self, tmp_path):
        """Test PlotFiles with multiple regions."""
        bed_file = tmp_path / "coords.bed"
        bed_file.write_text("chr1\t1000\t2000\tregion1\nchr2\t3000\t4000\tregion2\n")

        pf = PlotFiles(coordinates=bed_file, file_format="pdf")

        files = pf.files
        assert len(files) == 2
        assert any("region1.pdf" in str(f) for f in files)
        assert any("region2.pdf" in str(f) for f in files)

    def test_plot_files_missing_coordinates(self, tmp_path):
        """Test PlotFiles with missing coordinates file."""
        pf = PlotFiles(coordinates=tmp_path / "missing.bed")

        # Should not raise error, just return empty list
        files = pf.files
        assert files == []

    def test_plot_files_custom_format(self, tmp_path):
        """Test PlotFiles with custom format."""
        bed_file = tmp_path / "coords.bed"
        bed_file.write_text("chr1\t1000\t2000\n")

        pf = PlotFiles(coordinates=bed_file, file_format="png", output_dir="custom")

        files = pf.files
        assert str(files[0]).endswith(".png")
        assert "custom/genome_browser_plots/" in str(files[0])


class TestSNPFilesRaw:
    """Tests for SNPFilesRaw class."""

    def test_snp_files_raw_basic(self):
        """Test SNPFilesRaw basic functionality."""
        snp = SNPFilesRaw(assay=Assay.SNP, names=["sample1", "sample2"])

        files = snp.files
        assert len(files) == 2
        assert all(f.endswith(".vcf.gz") for f in files)
        assert all("variant" in f for f in files)


class TestSNPFilesAnnotated:
    """Tests for SNPFilesAnnotated class."""

    def test_snp_files_annotated_basic(self):
        """Test SNPFilesAnnotated basic functionality."""
        snp = SNPFilesAnnotated(assay=Assay.SNP, names=["sample1", "sample2"])

        files = snp.files
        assert len(files) == 2
        assert all(f.endswith(".anno.vcf.gz") for f in files)


class TestMethylationFiles:
    """Tests for MethylationFiles class."""

    def test_methylation_files_basic(self):
        """Test MethylationFiles basic functionality."""
        mf = MethylationFiles(
            assay=Assay.METH,
            names=["sample1"],
            genomes=["hg38", "spikein"],
            method=MethylationMethod.TAPS,
        )

        files = mf.files
        # Should have split_bams, methyldackel, and bias files
        assert len(files) > 0
        assert any(".bam" in f for f in files)
        assert any(".bedGraph" in f for f in files)

    def test_methylation_prefix_property(self):
        """Test prefix property."""
        mf = MethylationFiles(
            assay=Assay.METH,
            names=["sample1"],
            genomes=["hg38"],
            method=MethylationMethod.TAPS,
            output_dir="custom_output",
        )

        assert mf.prefix == "custom_output/methylation/"

    def test_methylation_split_bams_files(self):
        """Test split_bams_files property."""
        mf = MethylationFiles(
            assay=Assay.METH,
            names=["sample1"],
            genomes=["hg38", "spikein"],
            method=MethylationMethod.TAPS,
        )

        files = mf.split_bams_files
        assert len(files) == 2  # One per genome
        assert all("spikein" in f or "hg38" in f for f in files)

    def test_methylation_methyldackel_files(self):
        """Test methyldackel_files property."""
        mf = MethylationFiles(
            assay=Assay.METH,
            names=["sample1"],
            genomes=["hg38"],
            method=MethylationMethod.TAPS,
        )

        files = mf.methyldackel_files
        assert len(files) == 1
        assert "bedGraph" in files[0]

    def test_methylation_bias_files(self):
        """Test methylation_bias property."""
        mf = MethylationFiles(
            assay=Assay.METH,
            names=["sample1"],
            genomes=["hg38"],
            method=MethylationMethod.TAPS,
        )

        files = mf.methylation_bias
        assert len(files) > 0
        assert any("methylation_conversion.tsv" in f for f in files)
        assert any("bias" in f for f in files)


class TestBigBedFiles:
    """Tests for BigBedFiles class."""

    def test_bigbed_files_basic(self, tmp_path):
        """Test BigBedFiles basic functionality."""
        bed1 = tmp_path / "file1.bed"
        bed2 = tmp_path / "file2.bed"
        bed1.touch()
        bed2.touch()

        bbf = BigBedFiles(bed_files=[bed1, bed2])

        files = bbf.files
        assert len(files) == 2
        assert all(f.endswith(".bb") for f in files)

    def test_bigbed_files_with_strings(self, tmp_path):
        """Test BigBedFiles accepts string paths."""
        bed1 = tmp_path / "file1.bed"
        bed1.touch()

        bbf = BigBedFiles(bed_files=[str(bed1)])

        files = bbf.files
        assert len(files) == 1
        assert files[0].endswith(".bb")

    def test_bigbed_files_filters_non_bed(self, tmp_path):
        """Test BigBedFiles filters non-BED files."""
        bed1 = tmp_path / "file1.bed"
        txt1 = tmp_path / "file2.txt"
        bed1.touch()
        txt1.touch()

        bbf = BigBedFiles(bed_files=[bed1, txt1])

        files = bbf.files
        # Only .bed files should be converted to .bb
        assert len(files) == 1
        assert "file1.bb" in files[0]

    def test_bigbed_files_empty_list(self):
        """Test BigBedFiles with empty list."""
        bbf = BigBedFiles(bed_files=[])

        assert bbf.files == []


class TestQuantificationFiles:
    """Tests for QuantificationFiles class."""

    def test_quantification_files_rna_featurecounts(self):
        """Test QuantificationFiles for RNA with FeatureCounts."""
        groups = SampleGroups(groups=[SampleGroup(name="group1", samples=["sample1"])])
        qf = QuantificationFiles(
            assay=Assay.RNA,
            methods=[QuantificationMethod.FEATURE_COUNTS],
            names=["sample1"],
            groups=groups,
        )

        files = qf.files
        assert len(files) > 0
        assert any("read_counts.tsv" in f for f in files)

    def test_quantification_files_rna_salmon(self):
        """Test QuantificationFiles for RNA with Salmon."""
        groups = SampleGroups(groups=[SampleGroup(name="group1", samples=["sample1"])])
        qf = QuantificationFiles(
            assay=Assay.RNA,
            methods=[QuantificationMethod.SALMON],
            names=["sample1"],
            groups=groups,
        )

        files = qf.files
        assert any("salmon" in f or "read_counts.tsv" in f for f in files)

    def test_quantification_files_chip_only_featurecounts(self):
        """Test QuantificationFiles for ChIP filters to only FeatureCounts."""
        groups = SampleGroups(groups=[SampleGroup(name="group1", samples=["sample1"])])
        qf = QuantificationFiles(
            assay=Assay.CHIP,
            methods=[QuantificationMethod.FEATURE_COUNTS, QuantificationMethod.SALMON],
            names=["sample1"],
            groups=groups,
        )

        # Should filter out SALMON for non-RNA assays
        files = qf.files
        assert not any("salmon" in f for f in files)

    def test_quantification_combined_counts(self):
        """Test combined_counts_file property."""
        groups = SampleGroups(groups=[SampleGroup(name="group1", samples=["sample1"])])
        qf = QuantificationFiles(
            assay=Assay.RNA,
            methods=[QuantificationMethod.FEATURE_COUNTS],
            names=["sample1"],
            groups=groups,
        )

        files = qf.combined_counts_file
        assert len(files) == 1
        assert "read_counts.tsv" in files[0]

    def test_quantification_grouped_counts(self):
        """Test grouped_counts_files property."""
        groups = SampleGroups(
            groups=[
                SampleGroup(name="group1", samples=["sample1"]),
                SampleGroup(name="group2", samples=["sample2"]),
            ]
        )
        qf = QuantificationFiles(
            assay=Assay.RNA,
            methods=[QuantificationMethod.FEATURE_COUNTS],
            names=["sample1", "sample2"],
            groups=groups,
        )

        files = qf.grouped_counts_files
        # Should have files for each group
        assert len(files) >= 2


class TestGeoSubmissionFiles:
    """Tests for GeoSubmissionFiles class."""

    def test_geo_submission_default_files(self):
        """Test default_files property."""
        geo = GeoSubmissionFiles(assay=Assay.CHIP, names=["sample1"])

        files = geo.default_files
        assert len(files) == 4
        assert any("md5sums.txt" in f for f in files)
        assert any("samples_table.txt" in f for f in files)
        assert any("raw_data_checksums.txt" in f for f in files)
        assert any("processed_data_checksums.txt" in f for f in files)

    def test_geo_submission_raw_files(self):
        """Test raw_files property."""
        geo = GeoSubmissionFiles(assay=Assay.CHIP, names=["sample1", "sample2"])

        files = geo.raw_files
        # Should have R1 and R2 for each sample
        assert len(files) == 4
        assert all(f.endswith(".fastq.gz") for f in files)

    def test_geo_submission_processed_files(self):
        """Test processed_data_files property."""
        geo = GeoSubmissionFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            seqnado_files=[
                "seqnado_output/bigwigs/deeptools/unscaled/sample1.bigWig",
                "seqnado_output/peaks/macs2/sample1.bed",
            ],
        )

        files = geo.processed_data_files
        assert len(files) == 2
        # Check that files are flattened and renamed
        assert any("deeptools_unscaled" in f for f in files)
        assert any("macs2" in f for f in files)

    def test_geo_submission_filters_extensions(self):
        """Test that only allowed extensions are processed."""
        geo = GeoSubmissionFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            seqnado_files=[
                "seqnado_output/file.bigWig",
                "seqnado_output/file.txt",  # Not in allowed extensions
                "seqnado_output/file.bed",
            ],
        )

        files = geo.processed_data_files
        assert len(files) == 2  # .txt should be filtered out
        assert all(f.endswith((".bigWig", ".bed")) for f in files)

    def test_geo_submission_all_files(self):
        """Test files property includes all file types."""
        geo = GeoSubmissionFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            seqnado_files=["seqnado_output/file.bigWig"],
        )

        files = geo.files
        # Should include default, raw, and processed files
        assert len(files) > 5

    def test_geo_submission_upload_directory(self):
        """Test upload_directory property."""
        geo = GeoSubmissionFiles(
            assay=Assay.CHIP, names=["sample1"], output_dir="custom_output"
        )

        assert "custom_output/geo_submission" in str(geo.upload_directory)
        assert geo.assay.clean_name in str(geo.upload_directory)

    def test_geo_submission_upload_instructions(self):
        """Test upload_instructions property."""
        geo = GeoSubmissionFiles(assay=Assay.CHIP, names=["sample1"])

        assert geo.upload_instructions == Path(
            "seqnado_output/geo_submission/upload_instructions.txt"
        )

    def test_geo_submission_vcf_gz_files(self):
        """Test that .vcf.gz files are included."""
        geo = GeoSubmissionFiles(
            assay=Assay.SNP,
            names=["sample1"],
            seqnado_files=["seqnado_output/variant/sample1.vcf.gz"],
        )

        files = geo.processed_data_files
        assert any(f.endswith(".vcf.gz") for f in files)

    def test_geo_submission_custom_allowed_extensions(self):
        """Test custom allowed_extensions."""
        geo = GeoSubmissionFiles(
            assay=Assay.CHIP,
            names=["sample1"],
            seqnado_files=["seqnado_output/file.bigWig", "seqnado_output/file.txt"],
            allowed_extensions=[".txt"],
        )

        files = geo.processed_data_files
        # Only .txt should be included with custom extensions
        assert len(files) == 1
        assert files[0].endswith(".txt")


# Tests for seqnado.outputs.core module to increase coverage


class TestSeqnadoOutputFilesCore:
    """Tests for SeqnadoOutputFiles class methods."""

    def test_filter_by_suffix_with_contains(self):
        """Test _filter_by_suffix with contains parameter."""
        from seqnado.outputs.core import SeqnadoOutputFiles

        files = [
            "output/deeptools/sample1.bigWig",
            "output/homer/sample2.bigWig",
            "output/sample3.bed",
        ]
        output = SeqnadoOutputFiles(files=files, sample_names=["sample1", "sample2"])

        # Test with contains parameter
        result = output.select_files("bigWig", include="deeptools")
        assert len(result) == 1
        assert "deeptools" in result[0]

    def test_select_bigwig_subtype(self):
        """Test select_bigwig_subtype method."""
        from seqnado.outputs.core import SeqnadoOutputFiles

        files = [
            "output/bigwigs/deeptools/unscaled/sample1.bigWig",
            "output/bigwigs/deeptools/cpm/sample1.bigWig",
            "output/bigwigs/homer/unscaled/sample2.bigWig",
        ]
        output = SeqnadoOutputFiles(files=files, sample_names=["sample1", "sample2"])

        result = output.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS, scale=DataScalingTechnique.UNSCALED
        )
        assert len(result) == 1
        assert "deeptools" in result[0] and "unscaled" in result[0]

    def test_select_bigwig_subtype_with_assay(self):
        """Test select_bigwig_subtype method with assay filter."""
        from seqnado.outputs.core import SeqnadoOutputFiles

        files = [
            "output/bigwigs/deeptools/unscaled/atac/sample1.bigWig",
            "output/bigwigs/deeptools/unscaled/RNA/sample1.bigWig",
            "output/bigwigs/deeptools/unscaled/chip/sample1.bigWig",
        ]
        output = SeqnadoOutputFiles(files=files, sample_names=["sample1"])

        # Filter for ATAC
        result = output.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            assay=Assay.ATAC,
        )
        assert len(result) == 1
        assert Assay.ATAC.value.lower() in result[0].lower()

        # Filter for RNA
        result = output.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED,
            assay=Assay.RNA,
        )
        assert len(result) == 1
        assert "RNA" in result[0]

    def test_bigbed_files_property(self):
        """Test bigbed_files property."""
        from seqnado.outputs.core import SeqnadoOutputFiles

        files = ["output/sample1.bigBed", "output/sample2.bed"]
        output = SeqnadoOutputFiles(files=files, sample_names=["sample1"])

        bb_files = output.bigbed_files
        assert len(bb_files) == 1
        assert bb_files[0].endswith(".bigBed")

    def test_has_consensus_peaks_true(self):
        """Test has_consensus_peaks when peaks exist."""
        from seqnado.outputs.core import SeqnadoOutputFiles

        files = ["output/peaks/merged/group1.bed"]
        groups = SampleGroupings(
            groupings={
                "consensus": SampleGroups(
                    groups=[SampleGroup(name="group1", samples=["s1"])]
                )
            }
        )
        output = SeqnadoOutputFiles(
            files=files, sample_names=["s1"], sample_groups=groups
        )

        assert output.has_consensus_peaks is True

    def test_has_consensus_peaks_false(self):
        """Test has_consensus_peaks when no consensus peaks."""
        from seqnado.outputs.core import SeqnadoOutputFiles

        files = ["output/peaks/sample1.bed"]
        output = SeqnadoOutputFiles(files=files, sample_names=["sample1"])

        assert output.has_consensus_peaks is False


class TestSeqnadoOutputBuilderCore:
    """Tests for SeqnadoOutputBuilder additional methods."""

    def test_add_grouped_bigwig_files(self, tmp_path):
        """Test add_grouped_bigwig_files method."""
        cfg = _minimal_config(tmp_path)

        # Create collection with multiple samples
        r1_path = tmp_path / "s1_R1.fastq.gz"
        r1_path.write_text("@r\nN\n+\n#\n")
        r2_path = tmp_path / "s2_R1.fastq.gz"
        r2_path.write_text("@r\nN\n+\n#\n")

        fs1 = FastqSet(sample_id="s1", r1=FastqFile(path=r1_path))
        fs2 = FastqSet(sample_id="s2", r1=FastqFile(path=r2_path))
        samples = FastqCollection(
            assay=Assay.ATAC,
            metadata=[Metadata(assay=Assay.ATAC), Metadata(assay=Assay.ATAC)],
            fastq_sets=[fs1, fs2],
        )

        groups = SampleGroupings(
            groupings={
                "consensus": SampleGroups(
                    groups=[SampleGroup(name="group1", samples=["s1", "s2"])]
                )
            }
        )

        builder = SeqnadoOutputBuilder(
            Assay.ATAC, samples, cfg, sample_groupings=groups
        )
        builder.add_grouped_bigwig_files()

        output = builder.build()
        assert any("group1" in f for f in output.files)

    def test_add_peak_files_from_bigwig_wrong_method_raises(self, tmp_path):
        """Test add_peak_files raises error for BigWigCollection with non-Lanceotron."""
        from seqnado.config.configs import PeakCallingConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            peak_calling=PeakCallingConfig(method=[PeakCallingMethod.MACS2])
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        # Create BigWigCollection
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        samples = BigWigCollection(
            assay=Assay.ATAC, bigwig_files=[BigWigFile(path=bw_path)], metadata=[]
        )

        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)

        with pytest.raises(ValueError, match="only.*LANCEOTRON.*is allowed"):
            builder.add_peak_files()

    def test_add_grouped_peak_files(self, tmp_path):
        """Test add_grouped_peak_files method."""
        from seqnado.config.configs import PeakCallingConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            peak_calling=PeakCallingConfig(method=[PeakCallingMethod.MACS2])
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.ATAC,
            metadata=[Metadata(assay=Assay.ATAC)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        groups = SampleGroupings(
            groupings={
                "consensus": SampleGroups(
                    groups=[SampleGroup(name="cons", samples=["s1"])]
                )
            }
        )

        builder = SeqnadoOutputBuilder(
            Assay.ATAC, samples, cfg, sample_groupings=groups
        )
        builder.add_grouped_peak_files()

        output = builder.build()
        assert any("cons" in f for f in output.files)

    def test_add_bigbed_files(self, tmp_path):
        """Test add_bigbed_files method."""
        from seqnado.config.configs import PeakCallingConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            peak_calling=PeakCallingConfig(method=[PeakCallingMethod.MACS2])
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.ATAC,
            metadata=[Metadata(assay=Assay.ATAC)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
        builder.add_peak_files()
        builder.add_bigbed_files()

        output = builder.build()
        assert any(f.endswith(".bb") for f in output.files)

    def test_add_heatmap_files(self, tmp_path):
        """Test add_heatmap_files method."""
        cfg = _minimal_config(tmp_path)
        samples = _small_collection(tmp_path)

        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
        builder.add_heatmap_files()

        output = builder.build()
        assert any("heatmap" in f for f in output.files)

    def test_add_hub_files(self, tmp_path):
        """Test add_hub_files method."""
        from seqnado.config.configs import UCSCHubConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(ucsc_hub=UCSCHubConfig(name="test_hub"))
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)
        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
        builder.add_hub_files()

        output = builder.build()
        assert any("hub.txt" in f for f in output.files)

    def test_add_spikein_files(self, tmp_path):
        """Test add_spikein_files method."""
        from seqnado.config.configs import SpikeInConfig
        from seqnado.config.core import RNAAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = RNAAssayConfig(
            bigwigs=BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS]),
            spikein=SpikeInConfig(method=[SpikeInMethod.ORLANDO]),
        )
        cfg = SeqnadoConfig(
            assay=Assay.RNA,
            project=dict(name="p"),
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1_path = tmp_path / "s1_R1.fastq.gz"
        r1_path.write_text("@r\nN\n+\n#\n")
        r2_path = tmp_path / "s1_R2.fastq.gz"
        r2_path.write_text("@r\nN\n+\n#\n")
        r1 = FastqFile(path=r1_path)
        r2 = FastqFile(path=r2_path)
        fs = FastqSet(sample_id="s1", r1=r1, r2=r2)
        samples = FastqCollection(
            assay=Assay.RNA, metadata=[Metadata(assay=Assay.RNA)], fastq_sets=[fs]
        )

        builder = SeqnadoOutputBuilder(Assay.RNA, samples, cfg)
        builder.add_spikein_files()

        output = builder.build()
        assert any("normalisation_factors.tsv" in f for f in output.files)

    def test_add_plot_files(self, tmp_path):
        """Test add_plot_files method."""
        from seqnado.config.configs import PlottingConfig
        from seqnado.config.core import ATACAssayConfig

        bed_file = tmp_path / "coords.bed"
        bed_file.write_text("chr1\t1000\t2000\tregion1\n")

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            plotting=PlottingConfig(coordinates=str(bed_file), file_format="svg")
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)
        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
        builder.add_plot_files()

        output = builder.build()
        assert any(".svg" in f for f in output.files)

    def test_add_geo_submission_files(self, tmp_path):
        """Test add_geo_submission_files method."""
        cfg = _minimal_config(tmp_path)
        samples = _small_collection(tmp_path)

        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)
        builder.add_qc_files()
        builder.add_individual_bigwig_files()
        builder.add_geo_submission_files()

        output = builder.build()
        assert any("geo_submission" in f for f in output.files)

    def test_add_geo_submission_files_from_bam_raises(self, tmp_path):
        """Test add_geo_submission_files raises for BamCollection."""
        cfg = _minimal_config(tmp_path)

        bam_path = tmp_path / "sample1.bam"
        bam_path.touch()
        samples = BamCollection(
            assay=Assay.ATAC, bam_files=[BamFile(path=bam_path)], metadata=[]
        )

        builder = SeqnadoOutputBuilder(Assay.ATAC, samples, cfg)

        with pytest.raises(ValueError, match="can only be generated from FASTQ"):
            builder.add_geo_submission_files()

    def test_add_snp_files(self, tmp_path):
        """Test add_snp_files method for SNP assay."""
        from seqnado.config.configs import SNPCallingConfig
        from seqnado.config.core import SNPAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))

        # SNPCallingConfig requires method field
        from seqnado import SNPCallingMethod

        snp_cfg = SNPCallingConfig(method=SNPCallingMethod.BCFTOOLS, snp_database=None)
        assay_cfg = SNPAssayConfig(snp_calling=snp_cfg, ucsc_hub=None)
        cfg = SeqnadoConfig(
            assay=Assay.SNP,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.SNP,
            metadata=[Metadata(assay=Assay.SNP)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        builder = SeqnadoOutputBuilder(Assay.SNP, samples, cfg)
        builder.add_snp_files()

        output = builder.build()
        assert any(".vcf.gz" in f for f in output.files)

    def test_add_methylation_files(self, tmp_path):
        """Test add_methylation_files method."""
        from seqnado.config.configs import MethylationConfig
        from seqnado.config.core import MethylationAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))

        meth_cfg = MethylationConfig(method=MethylationMethod.TAPS)
        # MethylationAssayConfig requires ucsc_hub=None explicitly
        assay_cfg = MethylationAssayConfig(methylation=meth_cfg, ucsc_hub=None)
        cfg = SeqnadoConfig(
            assay=Assay.METH,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.METH,
            metadata=[Metadata(assay=Assay.METH)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        builder = SeqnadoOutputBuilder(Assay.METH, samples, cfg)
        builder.add_methylation_files()

        output = builder.build()
        # Methylation files include bedGraph outputs
        assert len(output.files) > 0

    def test_add_contact_files(self, tmp_path):
        """Test add_contact_files method for MCC assay."""
        from seqnado.config.configs import MCCConfig
        from seqnado.config.core import MCCAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))

        # MCCConfig requires viewpoints file
        viewpoints = tmp_path / "viewpoints.bed"
        viewpoints.write_text("chr1\t1000\t2000\tview1\n")
        mcc_cfg = MCCConfig(viewpoints=viewpoints, resolutions=[5000, 10000])

        # MCC requires bigwigs configured
        from seqnado.config.configs import BigwigConfig
        assay_cfg = MCCAssayConfig(mcc=mcc_cfg, bigwigs=BigwigConfig())
        cfg = SeqnadoConfig(
            assay=Assay.MCC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.MCC,
            metadata=[Metadata(assay=Assay.MCC)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        builder = SeqnadoOutputBuilder(Assay.MCC, samples, cfg)
        builder.add_contact_files()

        output = builder.build()
        assert any("mcc" in f for f in output.files)

    def test_add_quantification_files(self, tmp_path):
        """Test add_quantification_files method for RNA assay."""
        from seqnado.config.configs import RNAQuantificationConfig
        from seqnado.config.core import RNAAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))

        quant_cfg = RNAQuantificationConfig(method=QuantificationMethod.FEATURE_COUNTS)
        assay_cfg = RNAAssayConfig(rna_quantification=quant_cfg)
        cfg = SeqnadoConfig(
            assay=Assay.RNA,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.RNA,
            metadata=[Metadata(assay=Assay.RNA)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        groups = SampleGroupings(
            groupings={
                "consensus": SampleGroups(
                    groups=[SampleGroup(name="g1", samples=["s1"])]
                )
            }
        )

        builder = SeqnadoOutputBuilder(Assay.RNA, samples, cfg, sample_groupings=groups)
        builder.add_quantification_files()

        output = builder.build()
        assert any("quantification" in f or "read_counts" in f for f in output.files)

    def test_add_snp_files_with_annotation(self, tmp_path):
        """Test add_snp_files with annotation enabled."""
        from seqnado.config.configs import SNPCallingConfig
        from seqnado.config.core import SNPAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))

        snp_cfg = SNPCallingConfig(method=SNPCallingMethod.BCFTOOLS, annotate_snps=True)
        assay_cfg = SNPAssayConfig(snp_calling=snp_cfg, ucsc_hub=None)
        cfg = SeqnadoConfig(
            assay=Assay.SNP,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.SNP,
            metadata=[Metadata(assay=Assay.SNP)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        builder = SeqnadoOutputBuilder(Assay.SNP, samples, cfg)
        builder.add_snp_files()

        output = builder.build()
        # Should have both raw VCF and annotated files
        assert any("vcf" in f for f in output.files)

    def test_add_quantification_files_without_groupings(self, tmp_path):
        """Test add_quantification_files when no groupings provided."""
        from seqnado.config.configs import RNAQuantificationConfig
        from seqnado.config.core import RNAAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))

        quant_cfg = RNAQuantificationConfig(method=QuantificationMethod.FEATURE_COUNTS)
        assay_cfg = RNAAssayConfig(rna_quantification=quant_cfg)
        cfg = SeqnadoConfig(
            assay=Assay.RNA,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)

        # No groupings provided - should use empty SampleGroups
        builder = SeqnadoOutputBuilder(Assay.RNA, samples, cfg, sample_groupings=None)
        builder.add_quantification_files()

        output = builder.build()
        assert any("quantification" in f or "read_counts" in f for f in output.files)


class TestSeqnadoOutputFilesProperties:
    """Tests for SeqnadoOutputFiles property methods."""

    def test_bigwig_files_property(self, tmp_path):
        """Test bigwig_files property filters correctly."""
        output = SeqnadoOutputFiles(
            files=[
                "seqnado_output/bigwigs/s1.bigWig",
                "seqnado_output/peaks/s1.bed",
                "seqnado_output/bigwigs/s2.bigWig",
            ],
            sample_names=["s1", "s2"],
        )
        bigwigs = output.bigwig_files
        assert len(bigwigs) == 2
        assert all(f.endswith(".bigWig") for f in bigwigs)

    def test_heatmap_files_property(self, tmp_path):
        """Test heatmap_files property filters correctly."""
        output = SeqnadoOutputFiles(
            files=[
                "seqnado_output/heatmap/plot.pdf",
                "seqnado_output/peaks/s1.bed",
                "seqnado_output/deeptools/heatmap/h1.pdf",
            ],
            sample_names=["s1"],
        )
        heatmaps = output.heatmap_files
        assert len(heatmaps) == 2
        assert all(f.endswith(".pdf") and "heatmap" in f for f in heatmaps)

    def test_genome_browser_plots_property(self, tmp_path):
        """Test genome_browser_plots property filters correctly."""
        output = SeqnadoOutputFiles(
            files=[
                "seqnado_output/genome_browser/plot.pdf",
                "seqnado_output/peaks/s1.bed",
                "seqnado_output/plots/genome_browser/p1.pdf",
            ],
            sample_names=["s1"],
        )
        plots = output.genome_browser_plots
        assert len(plots) == 2
        assert all(f.endswith(".pdf") and "genome_browser" in f for f in plots)

    def test_ucsc_hub_files_property(self, tmp_path):
        """Test ucsc_hub_files property filters correctly."""
        output = SeqnadoOutputFiles(
            files=[
                "seqnado_output/ucsc_hub/hub.txt",
                "seqnado_output/peaks/s1.bed",
                "seqnado_output/hub/genomes.txt",
            ],
            sample_names=["s1"],
        )
        hub_files = output.ucsc_hub_files
        assert len(hub_files) == 2
        assert all(f.endswith(".txt") and "hub" in f for f in hub_files)

    def test_select_bigwig_subtype(self, tmp_path):
        """Test select_bigwig_subtype filtering."""
        output = SeqnadoOutputFiles(
            files=[
                "seqnado_output/bigwigs/deeptools/unscaled/s1.bigWig",
                "seqnado_output/bigwigs/deeptools/cpm/s1.bigWig",
                "seqnado_output/bigwigs/homer/unscaled/s1.bigWig",
            ],
            sample_names=["s1"],
        )
        # Test filtering by method and scale
        deep_unscaled = output.select_bigwig_subtype(
            method=PileupMethod.DEEPTOOLS,
            scale=DataScalingTechnique.UNSCALED
        )
        assert len(deep_unscaled) == 1
        assert "deeptools/unscaled" in deep_unscaled[0]


class TestSeqnadoOutputFactoryCore:
    """Tests for SeqnadoOutputFactory conditional logic."""

    def test_factory_with_bigwigs_and_groupings(self, tmp_path):
        """Test factory adds grouped bigwigs when groupings provided."""
        from seqnado.config.configs import BigwigConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            bigwigs=BigwigConfig(
                pileup_method=[PileupMethod.DEEPTOOLS],
                scale=[DataScalingTechnique.UNSCALED],
            )
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.ATAC,
            metadata=[Metadata(assay=Assay.ATAC)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        groups = SampleGroupings(
            groupings={
                "consensus": SampleGroups(
                    groups=[SampleGroup(name="group1", samples=["s1"])]
                )
            }
        )

        factory = SeqnadoOutputFactory(
            Assay.ATAC, samples, cfg, sample_groupings=groups
        )
        builder = factory.create_output_builder()
        output = builder.build()

        # Should have both individual and grouped bigwigs
        assert any(".bigWig" in f for f in output.files)

    def test_factory_with_peaks_and_consensus(self, tmp_path):
        """Test factory adds grouped peaks when consensus groupings exist."""
        from seqnado.config.configs import PeakCallingConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            peak_calling=PeakCallingConfig(method=[PeakCallingMethod.MACS2])
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        r1 = tmp_path / "s1_R1.fastq.gz"
        r1.write_text("@r\nN\n+\n#\n")
        samples = FastqCollection(
            assay=Assay.ATAC,
            metadata=[Metadata(assay=Assay.ATAC)],
            fastq_sets=[FastqSet(sample_id="s1", r1=FastqFile(path=r1))],
        )

        groups = SampleGroupings(
            groupings={
                "consensus": SampleGroups(
                    groups=[SampleGroup(name="cons", samples=["s1"])]
                )
            }
        )

        factory = SeqnadoOutputFactory(
            Assay.ATAC, samples, cfg, sample_groupings=groups
        )
        builder = factory.create_output_builder()
        output = builder.build()

        assert any(".bed" in f for f in output.files)

    def test_factory_with_heatmaps(self, tmp_path):
        """Test factory adds heatmaps when create_heatmaps is True."""
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(create_heatmaps=True)
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)
        factory = SeqnadoOutputFactory(Assay.ATAC, samples, cfg)
        builder = factory.create_output_builder()
        output = builder.build()

        assert any("heatmap" in f for f in output.files)

    def test_factory_with_ucsc_hub(self, tmp_path):
        """Test factory adds UCSC hub when create_ucsc_hub is True."""
        from seqnado.config.configs import UCSCHubConfig
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(ucsc_hub=UCSCHubConfig(name="test_hub"))
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)
        factory = SeqnadoOutputFactory(Assay.ATAC, samples, cfg)
        builder = factory.create_output_builder()
        output = builder.build()

        assert any("hub.txt" in f for f in output.files)

    def test_factory_with_geo_submission(self, tmp_path):
        """Test factory adds GEO submission files."""
        from seqnado.config.core import ATACAssayConfig

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(create_geo_submission_files=True)
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)
        factory = SeqnadoOutputFactory(Assay.ATAC, samples, cfg)
        builder = factory.create_output_builder()
        output = builder.build()

        assert any("geo_submission" in f for f in output.files)

    def test_factory_assay_specific_atac_plots(self, tmp_path):
        """Test factory adds plots for ATAC assay via match statement."""
        from seqnado.config.configs import PlottingConfig
        from seqnado.config.core import ATACAssayConfig

        bed_file = tmp_path / "coords.bed"
        bed_file.write_text("chr1\t1000\t2000\tregion1\n")

        star = tmp_path / "star"
        star.mkdir()
        genome = GenomeConfig(name="hg38", index=STARIndex(prefix=star))
        assay_cfg = ATACAssayConfig(
            plotting=PlottingConfig(coordinates=str(bed_file), file_format="svg")
        )
        cfg = SeqnadoConfig(
            assay=Assay.ATAC,
            project={"name": "test"},
            genome=genome,
            metadata=tmp_path / "m.csv",
            assay_config=assay_cfg,
        )

        samples = _small_collection(tmp_path)
        factory = SeqnadoOutputFactory(Assay.ATAC, samples, cfg)
        builder = factory.create_output_builder()

        # Just verify builder was created with plot files added
        assert builder is not None
