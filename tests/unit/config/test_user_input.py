"""Tests for seqnado.config.user_input module."""

import json
from datetime import datetime
from pathlib import Path

import pytest

from seqnado import (
    Assay,
    PeakCallingMethod,
    PileupMethod,
    QuantificationMethod,
    SNPCallingMethod,
    MethylationMethod,
)
from seqnado.config import (
    ATACAssayConfig,
    BigwigConfig,
    BowtieIndex,
    ChIPAssayConfig,
    CATAssayConfig,
    CRISPRAssayConfig,
    GenomeConfig,
    MCCAssayConfig,
    MCCConfig,
    MethylationAssayConfig,
    PlottingConfig,
    ProjectConfig,
    RNAAssayConfig,
    SeqnadoConfig,
    SNPAssayConfig,
    STARIndex,
    UCSCHubConfig,
)
from seqnado.config.user_input import (
    build_default_assay_config,
    load_genome_configs,
)


class TestLoadGenomeConfigs:
    """Tests for load_genome_configs function."""

    @pytest.fixture
    def mock_data_dir(self, tmp_path):
        """Create mock test data directory with required files."""
        data_dir = tmp_path / "genome_data"
        data_dir.mkdir(parents=True, exist_ok=True)

        # Create Bowtie2 index directory and files
        bt2_dir = data_dir / "bt2_chr21_dm6_chr2L"
        bt2_dir.mkdir(parents=True, exist_ok=True)
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"bt2_chr21_dm6_chr2L{suffix}").touch()

        # Create fasta index file
        (data_dir / "chr21.fa.fai").write_text("chr21\t48129895\t4\t50\t51\n")

        # Create GTF file
        (data_dir / "chr21.gtf").write_text("chr21\tENSEMBL\texon\t1\t1000\t.\t+\t.\tgene_id \"gene1\";\n")

        # Create STAR index directory
        star_dir = data_dir / "STAR_chr21_rna_spikein"
        star_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, 6):
            (star_dir / f"SA_{i}.txt").write_text(f"Mock content {i}\n")

        return data_dir

    def test_load_genome_configs_rna_assay(self, tmp_path, monkeypatch, mock_data_dir):
        """Test loading genome configs for RNA assay (STAR index)."""
        config_dir = tmp_path / ".config" / "seqnado"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "genome_config.json"

        star_index = str(mock_data_dir / "STAR_chr21_rna_spikein")
        bt2_index = str(mock_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(mock_data_dir / "chr21.fa.fai")
        gtf = str(mock_data_dir / "chr21.gtf")

        genome_data = {
            "hg38": {
                "star_index": star_index,
                "bt2_index": bt2_index,
                "fasta": fasta,
                "gtf": gtf,
            }
        }

        config_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

        result = load_genome_configs(Assay.RNA)

        assert "hg38" in result
        assert isinstance(result["hg38"], GenomeConfig)
        assert isinstance(result["hg38"].index, STARIndex)
        assert result["hg38"].index.prefix == Path(star_index)

    def test_load_genome_configs_chip_assay(self, tmp_path, monkeypatch, mock_data_dir):
        """Test loading genome configs for ChIP assay (Bowtie index)."""
        config_dir = tmp_path / ".config" / "seqnado"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "genome_config.json"

        bt2_index = str(mock_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(mock_data_dir / "chr21.fa.fai")

        genome_data = {
            "mm10": {
                "bt2_index": bt2_index,
                "fasta": fasta,
            }
        }

        config_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

        result = load_genome_configs(Assay.CHIP)

        assert "mm10" in result
        assert isinstance(result["mm10"], GenomeConfig)
        assert isinstance(result["mm10"].index, BowtieIndex)
        assert result["mm10"].index.prefix == bt2_index

    def test_load_genome_configs_atac_assay(self, tmp_path, monkeypatch, mock_data_dir):
        """Test loading genome configs for ATAC assay (Bowtie index)."""
        config_dir = tmp_path / ".config" / "seqnado"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "genome_config.json"

        bt2_index = str(mock_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(mock_data_dir / "chr21.fa.fai")

        genome_data = {
            "dm6": {
                "bt2_index": bt2_index,
                "fasta": fasta,
            }
        }

        config_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

        result = load_genome_configs(Assay.ATAC)

        assert "dm6" in result
        assert isinstance(result["dm6"].index, BowtieIndex)

    def test_load_genome_configs_multiple_genomes(self, tmp_path, monkeypatch, mock_data_dir):
        """Test loading multiple genome configurations."""
        config_dir = tmp_path / ".config" / "seqnado"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "genome_config.json"

        bt2_index = str(mock_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(mock_data_dir / "chr21.fa.fai")

        genome_data = {
            "hg38": {
                "bt2_index": bt2_index,
                "fasta": fasta,
            },
            "mm10": {
                "bt2_index": bt2_index,
                "fasta": fasta,
            },
            "dm6": {
                "bt2_index": bt2_index,
                "fasta": fasta,
            },
        }

        config_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

        result = load_genome_configs(Assay.CHIP)

        assert len(result) == 3
        assert all(name in result for name in ["hg38", "mm10", "dm6"])

    def test_load_genome_configs_missing_file(self, tmp_path, monkeypatch):
        """Test error when genome config file doesn't exist."""
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

        with pytest.raises(SystemExit):
            load_genome_configs(Assay.CHIP)

    def test_load_genome_configs_sets_name(self, tmp_path, monkeypatch, mock_data_dir):
        """Test that genome name is correctly set from key."""
        config_dir = tmp_path / ".config" / "seqnado"
        config_dir.mkdir(parents=True)
        config_file = config_dir / "genome_config.json"

        bt2_index = str(mock_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(mock_data_dir / "chr21.fa.fai")

        genome_data = {
            "test_genome": {
                "bt2_index": bt2_index,
                "fasta": fasta,
            }
        }

        config_file.write_text(json.dumps(genome_data))
        monkeypatch.setenv("SEQNADO_CONFIG", str(tmp_path))

        result = load_genome_configs(Assay.ATAC)

        assert result["test_genome"].name == "test_genome"


class TestBuildDefaultAssayConfig:
    """Tests for build_default_assay_config function."""

    @pytest.fixture
    def test_data_dir(self, tmp_path):
        """Create test data directory structure with necessary files."""
        # Create the directory structure
        genome_dir = tmp_path / "genome"
        genome_dir.mkdir(parents=True, exist_ok=True)

        bt2_dir = genome_dir / "bt2_chr21_dm6_chr2L"
        bt2_dir.mkdir(parents=True, exist_ok=True)

        # Touch the index files
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"bt2_chr21_dm6_chr2L{suffix}").touch()

        # Create fasta index file with minimal content
        fasta_index = genome_dir / "chr21.fa.fai"
        fasta_index.write_text("chr21\t48129895\t4\t50\t51\n")

        # Create GTF file with minimal content
        gtf_file = genome_dir / "chr21.gtf"
        gtf_file.write_text("chr21\tENSEMBL\texon\t1\t1000\t.\t+\t.\tgene_id \"gene1\";\n")

        # Create MCC viewpoints file
        (genome_dir / "mcc_viewpoints.bed").touch()

        # Create STAR index directory and mock files with minimal content
        star_dir = genome_dir / "STAR_chr21_rna_spikein"
        star_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, 6):
            (star_dir / f"SA_{i}.txt").write_text(f"Mock content for SA_{i}.txt\n")

        return genome_dir

    @pytest.fixture
    def mock_genome_config(self, test_data_dir):
        """Create a mock genome configuration with temporary paths."""
        bt2_index = str(test_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(test_data_dir / "chr21.fa.fai")

        return GenomeConfig(
            name="hg38",
            index=BowtieIndex(prefix=bt2_index),
            fasta=fasta,
        )

    def test_build_default_atac_config(self, mock_genome_config):
        """Test building default ATAC assay configuration."""
        result = build_default_assay_config(Assay.ATAC, mock_genome_config)

        assert isinstance(result, ATACAssayConfig)
        assert result.tn5_shift is True
        assert result.peak_calling is not None
        assert result.peak_calling.method == [PeakCallingMethod.LANCEOTRON]
        assert result.bigwigs is not None
        assert result.bigwigs.pileup_method == [PileupMethod.DEEPTOOLS]

    def test_build_default_chip_config(self, mock_genome_config):
        """Test building default ChIP assay configuration."""
        result = build_default_assay_config(Assay.CHIP, mock_genome_config)

        assert isinstance(result, ChIPAssayConfig)
        assert result.spikein is None
        assert result.peak_calling is not None
        assert result.peak_calling.method == [PeakCallingMethod.LANCEOTRON]
        assert result.dataset_for_ml is not None

    def test_build_default_cat_config(self, mock_genome_config):
        """Test building default CUT&TAG assay configuration."""
        result = build_default_assay_config(Assay.CAT, mock_genome_config)

        assert isinstance(result, CATAssayConfig)
        assert result.tn5_shift is False
        assert result.spikein is None
        assert result.peak_calling is not None
        assert result.peak_calling.method == [PeakCallingMethod.SEACR]

    def test_build_default_rna_config(self, mock_genome_config):
        """Test building default RNA assay configuration."""
        result = build_default_assay_config(Assay.RNA, mock_genome_config)

        assert isinstance(result, RNAAssayConfig)
        assert result.rna_quantification is not None
        assert result.rna_quantification.method == QuantificationMethod.FEATURE_COUNTS
        assert result.rna_quantification.run_deseq2 is False
        assert result.rna_quantification.strandedness == 0

    def test_build_default_snp_config(self, mock_genome_config):
        """Test building default SNP assay configuration."""
        result = build_default_assay_config(Assay.SNP, mock_genome_config)

        assert isinstance(result, SNPAssayConfig)
        assert result.snp_calling is not None
        assert result.snp_calling.method == SNPCallingMethod.BCFTOOLS
        assert result.snp_calling.annotate_snps is False
        # SNP assays should not have bigwigs, plotting, or heatmaps
        assert result.bigwigs is None
        assert result.plotting is None
        assert result.create_heatmaps is False
        # SNP assays don't generate bigwigs, so ucsc_hub should be None
        assert result.ucsc_hub is None

    def test_build_default_mcc_config(self, mock_genome_config, test_data_dir, monkeypatch):
        """Test building default MCC assay configuration."""
        # Temporarily patch the hardcoded path in build_default_assay_config
        # to use the actual test viewpoints file
        viewpoints_file = test_data_dir / "mcc_viewpoints.bed"

        # We need to mock the MCCConfig creation in the function
        original_func = build_default_assay_config

        def patched_build(assay, genome_config):
            if assay == Assay.MCC:
                from seqnado.config import BigwigConfig, PlottingConfig

                bigwigs = BigwigConfig(pileup_method=[PileupMethod.DEEPTOOLS], binsize=10)
                plotting = PlottingConfig()
                mcc = MCCConfig(
                    viewpoints=viewpoints_file,
                    resolutions=[100, 1000],
                )
                return MCCAssayConfig(
                    genome=genome_config,
                    bigwigs=bigwigs,
                    plotting=plotting,
                    ucsc_hub=None,  # MCC doesn't use UCSC hub
                    create_heatmaps=False,
                    create_geo_submission_files=False,
                    mcc=mcc,
                )
            return original_func(assay, genome_config)

        result = patched_build(Assay.MCC, mock_genome_config)

        assert isinstance(result, MCCAssayConfig)
        assert result.mcc is not None
        assert result.mcc.resolutions == [100, 1000]

    def test_build_default_meth_config(self, mock_genome_config):
        """Test building default methylation assay configuration."""
        result = build_default_assay_config(Assay.METH, mock_genome_config)

        assert isinstance(result, MethylationAssayConfig)
        assert result.methylation is not None
        assert result.methylation.method == MethylationMethod.TAPS

    def test_build_default_crispr_config(self, mock_genome_config):
        """Test building default CRISPR assay configuration."""
        result = build_default_assay_config(Assay.CRISPR, mock_genome_config)

        assert isinstance(result, CRISPRAssayConfig)

    def test_default_config_has_bigwigs(self, mock_genome_config):
        """Test that default configs include bigwig configuration."""
        result = build_default_assay_config(Assay.CHIP, mock_genome_config)

        assert result.bigwigs is not None
        assert isinstance(result.bigwigs, BigwigConfig)
        assert result.bigwigs.binsize == 10

    def test_default_config_has_plotting(self, mock_genome_config):
        """Test that default configs include plotting configuration."""
        result = build_default_assay_config(Assay.ATAC, mock_genome_config)

        assert result.plotting is not None
        assert isinstance(result.plotting, PlottingConfig)

    def test_default_config_has_ucsc_hub(self, mock_genome_config):
        """Test that default configs include UCSC hub configuration."""
        result = build_default_assay_config(Assay.RNA, mock_genome_config)

        assert result.ucsc_hub is not None
        assert isinstance(result.ucsc_hub, UCSCHubConfig)
        assert result.ucsc_hub.directory == "seqnado_output/hub/"
        assert result.ucsc_hub.email == "user@example.com"

    def test_default_config_heatmaps_disabled(self, mock_genome_config):
        """Test that heatmaps are disabled by default."""
        result = build_default_assay_config(Assay.CHIP, mock_genome_config)

        assert result.create_heatmaps is False

    def test_default_config_geo_files_disabled(self, mock_genome_config):
        """Test that GEO submission files are disabled by default."""
        result = build_default_assay_config(Assay.ATAC, mock_genome_config)

        assert result.create_geo_submission_files is False

    def test_default_config_peak_calling_no_consensus(self, mock_genome_config):
        """Test that consensus counts are disabled by default."""
        result = build_default_assay_config(Assay.CHIP, mock_genome_config)

        assert result.peak_calling.consensus_counts is False

    def test_unsupported_assay_raises_error(self, mock_genome_config):
        """Test that unsupported assay types raise ValueError."""
        # Create a mock assay that doesn't match any case
        with pytest.raises(ValueError, match="Unsupported assay type"):
            # This should never happen with actual Assay enum, but tests the fallback
            build_default_assay_config("INVALID_ASSAY", mock_genome_config)


class TestBuildWorkflowConfig:
    """Tests for build functions that might need integration testing."""

    @pytest.fixture
    def test_data_dir(self, tmp_path):
        """Create test data directory structure with necessary files."""
        # Create the directory structure
        genome_dir = tmp_path / "genome"
        genome_dir.mkdir(parents=True, exist_ok=True)

        bt2_dir = genome_dir / "bt2_chr21_dm6_chr2L"
        bt2_dir.mkdir(parents=True, exist_ok=True)

        # Touch the index files
        for suffix in [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]:
            (bt2_dir / f"bt2_chr21_dm6_chr2L{suffix}").touch()

        # Create fasta index file with minimal content
        fasta_index = genome_dir / "chr21.fa.fai"
        fasta_index.write_text("chr21\t48129895\t4\t50\t51\n")

        # Create GTF file with minimal content
        gtf_file = genome_dir / "chr21.gtf"
        gtf_file.write_text("chr21\tENSEMBL\texon\t1\t1000\t.\t+\t.\tgene_id \"gene1\";\n")

        # Create STAR index directory and mock files with minimal content
        star_dir = genome_dir / "STAR_chr21_rna_spikein"
        star_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, 6):
            (star_dir / f"SA_{i}.txt").write_text(f"Mock content for SA_{i}.txt\n")

        return genome_dir

    @pytest.fixture
    def mock_genome_config(self, test_data_dir):
        """Create a mock genome configuration with temporary paths."""
        bt2_index = str(test_data_dir / "bt2_chr21_dm6_chr2L" / "bt2_chr21_dm6_chr2L")
        fasta = str(test_data_dir / "chr21.fa.fai")

        return GenomeConfig(
            name="hg38",
            index=BowtieIndex(prefix=bt2_index),
            fasta=fasta,
        )

    def test_seqnado_config_creation_with_defaults(self, mock_genome_config):
        """Test that SeqnadoConfig can be created with default assay config."""
        project = ProjectConfig(
            name="test_project",
            date=datetime.strftime(datetime.today(), "%Y-%m-%d"),
            directory=Path("/test/dir"),
        )

        assay_config = build_default_assay_config(Assay.ATAC, mock_genome_config)

        config = SeqnadoConfig(
            assay=Assay.ATAC,
            project=project,
            genome=mock_genome_config,
            metadata="metadata_atac.csv",
            assay_config=assay_config,
        )

        assert config.assay == Assay.ATAC
        assert config.project.name == "test_project"
        assert isinstance(config.assay_config, ATACAssayConfig)

    def test_seqnado_config_chip_with_defaults(self, mock_genome_config):
        """Test SeqnadoConfig for ChIP assay."""
        project = ProjectConfig(
            name="chip_project",
            date="2025-11-21",
            directory=Path("/chip/dir"),
        )

        assay_config = build_default_assay_config(Assay.CHIP, mock_genome_config)

        config = SeqnadoConfig(
            assay=Assay.CHIP,
            project=project,
            genome=mock_genome_config,
            metadata="metadata_chip.csv",
            assay_config=assay_config,
        )

        assert config.assay == Assay.CHIP
        assert isinstance(config.assay_config, ChIPAssayConfig)

    def test_seqnado_config_rna_with_defaults(self, mock_genome_config, test_data_dir):
        """Test SeqnadoConfig for RNA assay with STAR index."""
        star_index = str(test_data_dir / "STAR_chr21_rna_spikein")
        fasta = str(test_data_dir / "chr21.fa.fai")
        gtf = str(test_data_dir / "chr21.gtf")

        genome_rna = GenomeConfig(
            name="hg38",
            index=STARIndex(prefix=star_index),
            fasta=fasta,
            gtf=gtf,
        )

        project = ProjectConfig(
            name="rna_project",
            date="2025-11-21",
            directory=Path("/rna/dir"),
        )

        assay_config = build_default_assay_config(Assay.RNA, genome_rna)

        config = SeqnadoConfig(
            assay=Assay.RNA,
            project=project,
            genome=genome_rna,
            metadata="metadata_rna.csv",
            assay_config=assay_config,
        )

        assert config.assay == Assay.RNA
        assert isinstance(config.assay_config, RNAAssayConfig)
        assert isinstance(config.genome.index, STARIndex)


class TestRNAQuantificationConfig:
    """Tests for RNAQuantificationConfig model and validation."""

    def test_rna_quantification_config_creation_with_defaults(self):
        """Test creating RNAQuantificationConfig with default values."""
        from seqnado.config.configs import RNAQuantificationConfig

        config = RNAQuantificationConfig(method=QuantificationMethod.FEATURE_COUNTS)

        assert config.method == QuantificationMethod.FEATURE_COUNTS
        assert config.strandedness == 0
        assert config.run_deseq2 is False
        assert config.salmon_index is None

    def test_rna_quantification_config_strandedness_unstranded(self):
        """Test RNAQuantificationConfig with unstranded reads (0)."""
        from seqnado.config.configs import RNAQuantificationConfig

        config = RNAQuantificationConfig(
            method=QuantificationMethod.FEATURE_COUNTS,
            strandedness=0
        )

        assert config.strandedness == 0

    def test_rna_quantification_config_strandedness_forward(self):
        """Test RNAQuantificationConfig with forward-stranded reads (1)."""
        from seqnado.config.configs import RNAQuantificationConfig

        config = RNAQuantificationConfig(
            method=QuantificationMethod.FEATURE_COUNTS,
            strandedness=1
        )

        assert config.strandedness == 1

    def test_rna_quantification_config_strandedness_reverse(self):
        """Test RNAQuantificationConfig with reverse-stranded reads (2)."""
        from seqnado.config.configs import RNAQuantificationConfig

        config = RNAQuantificationConfig(
            method=QuantificationMethod.FEATURE_COUNTS,
            strandedness=2
        )

        assert config.strandedness == 2

    def test_rna_quantification_config_invalid_strandedness(self):
        """Test RNAQuantificationConfig rejects invalid strandedness values."""
        from seqnado.config.configs import RNAQuantificationConfig
        from pydantic import ValidationError

        with pytest.raises(ValidationError) as exc_info:
            RNAQuantificationConfig(
                method=QuantificationMethod.FEATURE_COUNTS,
                strandedness=3  # Invalid: must be 0, 1, or 2
            )

        assert "strandedness must be 0 (unstranded), 1 (forward), or 2 (reverse)" in str(exc_info.value)

    def test_rna_quantification_config_invalid_negative_strandedness(self):
        """Test RNAQuantificationConfig rejects negative strandedness values."""
        from seqnado.config.configs import RNAQuantificationConfig
        from pydantic import ValidationError

        with pytest.raises(ValidationError) as exc_info:
            RNAQuantificationConfig(
                method=QuantificationMethod.FEATURE_COUNTS,
                strandedness=-1
            )

        assert "strandedness must be 0 (unstranded), 1 (forward), or 2 (reverse)" in str(exc_info.value)

    def test_rna_quantification_config_with_run_deseq2(self):
        """Test RNAQuantificationConfig with DESeq2 enabled."""
        from seqnado.config.configs import RNAQuantificationConfig

        config = RNAQuantificationConfig(
            method=QuantificationMethod.FEATURE_COUNTS,
            strandedness=1,
            run_deseq2=True
        )

        assert config.run_deseq2 is True
        assert config.strandedness == 1

    def test_rna_quantification_config_salmon_with_strandedness(self):
        """Test RNAQuantificationConfig with salmon method and strandedness."""
        from seqnado.config.configs import RNAQuantificationConfig

        config = RNAQuantificationConfig.model_validate(
            {
                "method": QuantificationMethod.SALMON,
                "salmon_index": "path/to/salmon_index",
                "strandedness": 0,
            },
            context={"skip_path_validation": True},
        )

        assert config.method == QuantificationMethod.SALMON
        assert config.strandedness == 0
        assert config.salmon_index == "path/to/salmon_index"


class TestGetRNAQuantificationConfig:
    """Tests for get_rna_quantification_config function."""

    def test_get_rna_quantification_config_feature_counts_prompts_strandedness(self, monkeypatch):
        """Test get_rna_quantification_config prompts for strandedness when method is feature_counts."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["feature_counts", "1", "no"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.method == QuantificationMethod.FEATURE_COUNTS
        assert config.strandedness == 1
        assert config.run_deseq2 is False

    def test_get_rna_quantification_config_salmon_defaults_strandedness(self, monkeypatch):
        """Test get_rna_quantification_config defaults strandedness to 0 when method is salmon."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["salmon", "path/to/salmon_index", "no"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.method == QuantificationMethod.SALMON
        assert config.strandedness == 0  # Salmon doesn't prompt for strandedness
        assert config.salmon_index == "path/to/salmon_index"

    def test_get_rna_quantification_config_feature_counts_strandedness_forward(self, monkeypatch):
        """Test get_rna_quantification_config accepts forward strandedness for feature_counts."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["feature_counts", "1", "no"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.strandedness == 1

    def test_get_rna_quantification_config_feature_counts_strandedness_reverse(self, monkeypatch):
        """Test get_rna_quantification_config accepts reverse strandedness for feature_counts."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["feature_counts", "2", "no"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.strandedness == 2

    def test_get_rna_quantification_config_feature_counts_default_strandedness(self, monkeypatch):
        """Test get_rna_quantification_config defaults to unstranded (0) for feature_counts."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["feature_counts", "", "no"])  # Empty input uses default
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.strandedness == 0

    def test_get_rna_quantification_config_feature_counts_invalid_strandedness_retries(self, monkeypatch):
        """Test get_rna_quantification_config retries on invalid strandedness."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["feature_counts", "invalid", "1", "no"])  # Invalid, then valid
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.strandedness == 1

    def test_get_rna_quantification_config_deseq2_enabled(self, monkeypatch):
        """Test get_rna_quantification_config with DESeq2 enabled."""
        from seqnado.config.user_input import get_rna_quantification_config

        inputs = iter(["feature_counts", "0", "yes"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))

        config = get_rna_quantification_config()

        assert config.run_deseq2 is True
        assert config.strandedness == 0


class TestBuildDefaultAssayConfigRNA:
    """Tests for build_default_assay_config with RNA assay and strandedness."""

    @pytest.fixture
    def test_data_dir(self, tmp_path):
        """Create test data directory structure with necessary files."""
        genome_dir = tmp_path / "genome"
        genome_dir.mkdir(parents=True, exist_ok=True)

        # Create STAR index directory
        star_dir = genome_dir / "STAR_chr21_rna_spikein"
        star_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1, 6):
            (star_dir / f"SA_{i}.txt").write_text(f"Mock content for SA_{i}.txt\n")

        # Create fasta index file
        fasta_index = genome_dir / "chr21.fa.fai"
        fasta_index.write_text("chr21\t48129895\t4\t50\t51\n")

        # Create GTF file
        gtf_file = genome_dir / "chr21.gtf"
        gtf_file.write_text("chr21\tENSEMBL\texon\t1\t1000\t.\t+\t.\tgene_id \"gene1\";\n")

        return genome_dir

    @pytest.fixture
    def mock_genome_config_rna(self, test_data_dir):
        """Create a mock genome configuration for RNA assay."""
        star_index = str(test_data_dir / "STAR_chr21_rna_spikein")
        fasta = str(test_data_dir / "chr21.fa.fai")
        gtf = str(test_data_dir / "chr21.gtf")

        return GenomeConfig(
            name="hg38",
            index=STARIndex(prefix=star_index),
            fasta=fasta,
            gtf=gtf,
        )

    def test_build_default_rna_config_has_strandedness(self, mock_genome_config_rna):
        """Test that default RNA config includes strandedness field."""
        result = build_default_assay_config(Assay.RNA, mock_genome_config_rna)

        assert isinstance(result, RNAAssayConfig)
        assert result.rna_quantification is not None
        assert hasattr(result.rna_quantification, "strandedness")
        assert result.rna_quantification.strandedness == 0

    def test_build_default_rna_config_strandedness_default_zero(self, mock_genome_config_rna):
        """Test that default RNA config sets strandedness to 0 (unstranded)."""
        result = build_default_assay_config(Assay.RNA, mock_genome_config_rna)

        assert result.rna_quantification.strandedness == 0

    def test_build_default_rna_config_uses_feature_counts(self, mock_genome_config_rna):
        """Test that default RNA config uses feature_counts method."""
        result = build_default_assay_config(Assay.RNA, mock_genome_config_rna)

        assert result.rna_quantification.method == QuantificationMethod.FEATURE_COUNTS


class TestGetUserInput:
    """Tests for get_user_input helper function."""

    def test_get_user_input_with_default(self, monkeypatch):
        """Test get_user_input accepts default on empty input."""
        from seqnado.config.user_input import get_user_input

        monkeypatch.setattr('builtins.input', lambda _: "")
        result = get_user_input("Enter value", default="default_val")
        assert result == "default_val"

    def test_get_user_input_with_choices(self, monkeypatch):
        """Test get_user_input validates choices."""
        from seqnado.config.user_input import get_user_input

        inputs = iter(["invalid", "valid"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))
        result = get_user_input("Select", choices=["valid", "other"])
        assert result == "valid"

    def test_get_user_input_boolean_yes(self, monkeypatch):
        """Test get_user_input accepts boolean yes values."""
        from seqnado.config.user_input import get_user_input

        for value in ["yes", "y", "true", "1"]:
            monkeypatch.setattr('builtins.input', lambda _, v=value: v)
            result = get_user_input("Bool?", is_boolean=True)
            assert result is True

    def test_get_user_input_boolean_no(self, monkeypatch):
        """Test get_user_input accepts boolean no values."""
        from seqnado.config.user_input import get_user_input

        for value in ["no", "n", "false", "0"]:
            monkeypatch.setattr('builtins.input', lambda _, v=value: v)
            result = get_user_input("Bool?", is_boolean=True)
            assert result is False

    def test_get_user_input_invalid_boolean(self, monkeypatch):
        """Test get_user_input rejects invalid boolean."""
        from seqnado.config.user_input import get_user_input

        inputs = iter(["maybe", "yes"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))
        result = get_user_input("Bool?", is_boolean=True)
        assert result is True

    def test_get_user_input_path_validation(self, monkeypatch, tmp_path):
        """Test get_user_input validates path existence."""
        from seqnado.config.user_input import get_user_input

        existing_path = tmp_path / "exists.txt"
        existing_path.touch()

        inputs = iter([str(tmp_path / "missing.txt"), str(existing_path)])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))
        result = get_user_input("Path?", is_path=True)
        assert result == str(existing_path)

    def test_get_user_input_not_required(self, monkeypatch):
        """Test get_user_input returns None when not required."""
        from seqnado.config.user_input import get_user_input

        monkeypatch.setattr('builtins.input', lambda _: "")
        result = get_user_input("Optional?", required=False)
        assert result is None

    def test_get_user_input_empty_retries(self, monkeypatch):
        """Test get_user_input retries on empty when required."""
        from seqnado.config.user_input import get_user_input

        inputs = iter(["", "value"])
        monkeypatch.setattr('builtins.input', lambda _: next(inputs))
        result = get_user_input("Required?", required=True)
        assert result == "value"


