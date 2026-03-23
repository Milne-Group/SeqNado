"""Tests for SeqnadoConfig class in seqnado.config.core module."""

from pathlib import Path
from datetime import date

import pytest

from seqnado import Assay, PCRDuplicateHandling, PCRDuplicateTool
from seqnado.config.configs import (
    PCRDuplicatesConfig,
    ProjectConfig,
    QCConfig,
)
from seqnado.config.core import (
    ATACAssayConfig,
    ChIPAssayConfig,
    CATAssayConfig,
    RNAAssayConfig,
    SNPAssayConfig,
    MethylationAssayConfig,
    MCCAssayConfig,
    CRISPRAssayConfig,
    SeqnadoConfig,
)
from seqnado.config.configs import GenomeConfig, STARIndex, BowtieIndex


@pytest.fixture
def tmp_metadata(tmp_path: Path) -> Path:
    """Create a temporary metadata file."""
    metadata = tmp_path / "metadata.csv"
    metadata.write_text("samplename,fastq1,fastq2\nsample1,/path/to/read1.fq,/path/to/read2.fq\n")
    return metadata


@pytest.fixture
def tmp_star_index(tmp_path: Path) -> Path:
    """Create a temporary STAR index directory."""
    star_dir = tmp_path / "star_index"
    star_dir.mkdir()
    (star_dir / "SA").write_text("mock")
    return star_dir


@pytest.fixture
def tmp_bowtie_index(tmp_path: Path) -> Path:
    """Create a temporary Bowtie2 index."""
    idx_prefix = tmp_path / "bowtie_index" / "genome"
    idx_prefix.parent.mkdir(parents=True)
    for suffix in [".1.bt2", ".2.bt2", ".rev.1.bt2"]:
        (Path(str(idx_prefix) + suffix)).write_text("")
    return idx_prefix


class TestPCRDuplicateDefaults:
    """Tests for PCR duplicate handling defaults based on assay type."""

    def test_atac_defaults_to_remove(self, tmp_metadata, tmp_bowtie_index):
        """Test ATAC assay defaults to REMOVE duplicates."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE

    def test_chip_defaults_to_remove(self, tmp_metadata, tmp_bowtie_index):
        """Test ChIP assay defaults to REMOVE duplicates."""
        config = SeqnadoConfig(
            assay=Assay.CHIP,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE

    def test_cat_defaults_to_remove(self, tmp_metadata, tmp_bowtie_index):
        """Test CAT assay defaults to REMOVE duplicates."""
        config = SeqnadoConfig(
            assay=Assay.CAT,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE

    def test_snp_defaults_to_remove(self, tmp_metadata, tmp_bowtie_index):
        """Test SNP assay defaults to REMOVE duplicates."""
        config = SeqnadoConfig(
            assay=Assay.SNP,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE

    def test_meth_defaults_to_remove(self, tmp_metadata, tmp_bowtie_index):
        """Test METH assay defaults to REMOVE duplicates."""
        config = SeqnadoConfig(
            assay=Assay.METH,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE

    def test_rna_defaults_to_keep(self, tmp_metadata, tmp_star_index):
        """Test RNA assay defaults to KEEP duplicates."""
        config = SeqnadoConfig(
            assay=Assay.RNA,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=STARIndex(prefix=tmp_star_index)),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.NONE

    def test_mcc_defaults_to_keep(self, tmp_metadata, tmp_bowtie_index):
        """Test MCC assay defaults to KEEP duplicates."""
        config = SeqnadoConfig(
            assay=Assay.MCC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.NONE

    def test_crispr_defaults_to_keep(self, tmp_metadata, tmp_bowtie_index):
        """Test CRISPR assay defaults to KEEP duplicates."""
        config = SeqnadoConfig(
            assay=Assay.CRISPR,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.NONE

    def test_explicit_pcr_duplicates_override_default(self, tmp_metadata, tmp_bowtie_index):
        """Test that explicitly setting pcr_duplicates overrides the default."""
        # ATAC normally defaults to REMOVE, but we explicitly set it to NONE
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            pcr_duplicates=PCRDuplicatesConfig(strategy=PCRDuplicateHandling.NONE),
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.NONE

    def test_explicit_pcr_duplicates_with_tool(self, tmp_metadata, tmp_star_index):
        """Test that explicitly setting pcr_duplicates with tool works."""
        # RNA normally defaults to KEEP, but we explicitly set it to REMOVE with samtools
        config = SeqnadoConfig(
            assay=Assay.RNA,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=STARIndex(prefix=tmp_star_index)),
            metadata=tmp_metadata,
            pcr_duplicates=PCRDuplicatesConfig(
                strategy=PCRDuplicateHandling.REMOVE,
                tool=PCRDuplicateTool.SAMTOOLS
            ),
        )

        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE
        assert config.pcr_duplicates.tool == PCRDuplicateTool.SAMTOOLS


class TestSeqnadoConfigValidation:
    """Tests for SeqnadoConfig validation and default values."""

    def test_minimal_config_creation(self, tmp_metadata, tmp_bowtie_index):
        """Test creating a minimal SeqnadoConfig."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test_project"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.assay == Assay.ATAC
        assert config.project.name == "test_project"
        assert config.genome.name == "hg38"
        assert config.metadata == tmp_metadata
        assert isinstance(config.qc, QCConfig)
        assert isinstance(config.pcr_duplicates, PCRDuplicatesConfig)

    def test_default_qc_config(self, tmp_metadata, tmp_bowtie_index):
        """Test default QC config is created."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
        )

        assert config.qc.run_fastq_screen is True
        assert config.qc.calculate_library_complexity is False

    def test_custom_qc_config(self, tmp_metadata, tmp_bowtie_index):
        """Test custom QC config can be provided."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            qc=QCConfig(run_fastq_screen=False, calculate_library_complexity=True),
        )

        assert config.qc.run_fastq_screen is False
        assert config.qc.calculate_library_complexity is True


class TestAssayConfigMatching:
    """Tests for assay_config validation to match assay type."""

    def test_atac_assay_config(self, tmp_metadata, tmp_bowtie_index):
        """Test ATAC assay with ATACAssayConfig."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=ATACAssayConfig(tn5_shift=True, create_heatmaps=False),
        )

        assert isinstance(config.assay_config, ATACAssayConfig)
        assert config.assay_config.tn5_shift is True

    def test_rna_assay_config(self, tmp_metadata, tmp_star_index):
        """Test RNA assay with RNAAssayConfig."""
        config = SeqnadoConfig(
            assay=Assay.RNA,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=STARIndex(prefix=tmp_star_index)),
            metadata=tmp_metadata,
            assay_config=RNAAssayConfig(create_heatmaps=True),
        )

        assert isinstance(config.assay_config, RNAAssayConfig)
        assert config.assay_config.create_heatmaps is True

    def test_chip_assay_config(self, tmp_metadata, tmp_bowtie_index):
        """Test ChIP assay with ChIPAssayConfig."""
        config = SeqnadoConfig(
            assay=Assay.CHIP,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=ChIPAssayConfig(create_heatmaps=True),
        )

        assert isinstance(config.assay_config, ChIPAssayConfig)

    def test_assay_config_from_dict(self, tmp_metadata, tmp_bowtie_index):
        """Test assay_config can be created from dict."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config={"tn5_shift": True, "create_heatmaps": False},
        )

        assert isinstance(config.assay_config, ATACAssayConfig)
        assert config.assay_config.tn5_shift is True


class TestSeqnadoConfigFromYAML:
    """Tests for loading SeqnadoConfig from YAML."""

    def test_from_yaml_basic(self, tmp_path, tmp_metadata, tmp_bowtie_index):
        """Test loading basic config from YAML."""
        yaml_content = f"""
assay: ATAC
project:
  name: test_project
  date: 2024-01-01
genome:
  name: hg38
  index:
    type: Bowtie2
    prefix: {tmp_bowtie_index}
metadata: {tmp_metadata}
"""
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text(yaml_content)

        config = SeqnadoConfig.from_yaml(yaml_file)

        assert config.assay == Assay.ATAC
        assert config.project.name == "test_project"
        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE  # default for ATAC

    def test_from_yaml_with_explicit_pcr_duplicates(self, tmp_path, tmp_metadata, tmp_bowtie_index):
        """Test loading config from YAML with explicit pcr_duplicates."""
        yaml_content = f"""
assay: RNA
project:
  name: test_project
genome:
  name: hg38
  index:
    type: Bowtie2
    prefix: {tmp_bowtie_index}
metadata: {tmp_metadata}
pcr_duplicates:
  strategy: remove
  tool: samtools
"""
        yaml_file = tmp_path / "config.yaml"
        yaml_file.write_text(yaml_content)

        config = SeqnadoConfig.from_yaml(yaml_file)

        assert config.assay == Assay.RNA
        assert config.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE
        assert config.pcr_duplicates.tool == PCRDuplicateTool.SAMTOOLS


class TestSpikeinConfig:
    """Tests for spikein field in BaseAssayConfig."""

    def test_spikein_none_default_for_atac(self, tmp_metadata, tmp_bowtie_index):
        """Test spikein defaults to None for ATAC assay."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=ATACAssayConfig(),
        )

        assert config.assay_config.spikein is None

    def test_spikein_none_default_for_snp(self, tmp_metadata, tmp_bowtie_index):
        """Test spikein defaults to None for SNP assay."""
        config = SeqnadoConfig(
            assay=Assay.SNP,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=SNPAssayConfig(),
        )

        assert config.assay_config.spikein is None

    def test_spikein_none_default_for_methylation(self, tmp_metadata, tmp_bowtie_index):
        """Test spikein defaults to None for Methylation assay."""
        config = SeqnadoConfig(
            assay=Assay.METH,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=MethylationAssayConfig(ucsc_hub=None),
        )

        assert config.assay_config.spikein is None

    def test_spikein_none_default_for_crispr(self, tmp_metadata, tmp_bowtie_index):
        """Test spikein defaults to None for CRISPR assay."""
        config = SeqnadoConfig(
            assay=Assay.CRISPR,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=CRISPRAssayConfig(),
        )

        assert config.assay_config.spikein is None

    def test_has_spikein_false_when_spikein_none(self, tmp_metadata, tmp_bowtie_index):
        """Test has_spikein computed field is False when spikein is None."""
        config = SeqnadoConfig(
            assay=Assay.CHIP,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=ChIPAssayConfig(spikein=None),
        )

        assert config.assay_config.spikein is None
        assert config.assay_config.has_spikein is False

    def test_has_spikein_false_for_atac(self, tmp_metadata, tmp_bowtie_index):
        """Test has_spikein is False for ATAC with no spike-in config."""
        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=ProjectConfig(name="test"),
            genome=GenomeConfig(name="hg38", index=BowtieIndex(prefix=str(tmp_bowtie_index))),
            metadata=tmp_metadata,
            assay_config=ATACAssayConfig(),
        )

        assert config.assay_config.has_spikein is False
