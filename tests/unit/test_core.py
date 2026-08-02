"""Tests for seqnado.core module."""

import pytest
from pydantic import ValidationError

from seqnado.config import GenomicCoordinate
from seqnado.core import (
    Assay,
    FileType,
    PileupMethod,
    DataScalingTechnique,
    PeakCallingMethod,
    PCRDuplicateHandling,
    PCRDuplicateTool,
    SpikeInMethod,
    SNPCallingMethod,
    QuantificationMethod,
    MethylationMethod,
    Molecule,
    Organism,
    LibraryType,
    NONE_VALUES,
)


class TestNoneValues:
    """Tests for NONE_VALUES constant."""

    def test_none_values_contains_expected_strings(self):
        """Test that NONE_VALUES contains all expected null-like values."""
        assert None in NONE_VALUES
        assert "None" in NONE_VALUES
        assert "none" in NONE_VALUES
        assert "null" in NONE_VALUES
        assert "Null" in NONE_VALUES
        assert "NULL" in NONE_VALUES
        assert "." in NONE_VALUES
        assert "" in NONE_VALUES
        assert "NA" in NONE_VALUES


class TestAssayEnum:
    """Tests for Assay enum."""

    def test_all_assays_returns_list(self):
        """Test all_assays class method returns all assays."""
        assays = Assay.all_assays()
        assert isinstance(assays, list)
        assert Assay.RNA in assays
        assert Assay.ATAC in assays

    def test_clean_name_for_all_assays(self):
        """Test clean_name property for all assay types."""
        assert Assay.RNA.clean_name == "rna"
        assert Assay.ATAC.clean_name == "atac"
        assert Assay.SNP.clean_name == "snp"
        assert Assay.CHIP.clean_name == "chip"
        assert Assay.CAT.clean_name == "cat"
        assert Assay.METH.clean_name == "meth"
        assert Assay.MCC.clean_name == "mcc"
        assert Assay.CRISPR.clean_name == "crispr"
        assert Assay.MULTIOMICS.clean_name == "multiomics"

    def test_clean_name_unknown_assay_raises_error(self):
        """Test clean_name raises error for unmapped assay."""
        # Create a mock assay instance that won't be in the short_names dict
        from unittest.mock import MagicMock
        mock_assay = MagicMock(spec=Assay)
        mock_assay.__class__ = Assay
        
        # Get the clean_name property function
        clean_name_prop = Assay.clean_name.fget
        
        # Call it with the mock - should raise ValueError
        with pytest.raises(ValueError, match="Unknown assay type"):
            clean_name_prop(mock_assay)

    def test_from_clean_name(self):
        """Test from_clean_name class method."""
        assert Assay.from_clean_name("rna") == Assay.RNA
        assert Assay.from_clean_name("atac") == Assay.ATAC
        assert Assay.from_clean_name("chip") == Assay.CHIP
        assert Assay.from_clean_name("multiomics") == Assay.MULTIOMICS

    def test_from_clean_name_invalid_raises_error(self):
        """Test from_clean_name raises error for invalid name."""
        with pytest.raises(ValueError, match="Unknown clean name"):
            Assay.from_clean_name("invalid_assay")

    def test_all_assay_clean_names(self):
        """Test all_assay_clean_names returns all clean names."""
        clean_names = Assay.all_assay_clean_names()
        assert "rna" in clean_names
        assert "atac" in clean_names
        assert "chip" in clean_names
        assert len(clean_names) == len(Assay)

    def test_non_ip_assays(self):
        """Test non_ip_assays returns assays without IP requirement."""
        non_ip = Assay.non_ip_assays()
        assert Assay.RNA in non_ip
        assert Assay.ATAC in non_ip
        assert Assay.SNP in non_ip
        assert Assay.CHIP not in non_ip
        assert Assay.CAT not in non_ip

    def test_ip_assays(self):
        """Test ip_assays returns assays requiring IP."""
        ip_assays = Assay.ip_assays()
        assert Assay.CHIP in ip_assays
        assert Assay.CAT in ip_assays
        assert len(ip_assays) == 2

    def test_non_multiomics_assays(self):
        """Test non_multiomics_assays excludes MULTIOMICS."""
        non_multi = Assay.non_multiomics_assays()
        assert Assay.MULTIOMICS not in non_multi
        assert Assay.RNA in non_multi
        assert Assay.ATAC in non_multi
        assert Assay.CHIP in non_multi
        assert len(non_multi) == len(Assay) - 1


class TestFileTypeEnum:
    """Tests for FileType enum."""

    def test_file_type_values(self):
        """Test FileType enum has expected values."""
        assert FileType.BAM.value == "bam"
        assert FileType.BED.value == "bed"
        assert FileType.BIGWIG.value == "bigwig"
        assert FileType.FASTQ.value == "fastq"
        assert FileType.TAG_DIRECTORY.value == "homer_tag_directory"


class TestPileupMethodEnum:
    """Tests for PileupMethod enum."""

    def test_pileup_method_values(self):
        """Test PileupMethod enum has expected values."""
        assert PileupMethod.DEEPTOOLS.value == "deeptools"
        assert PileupMethod.HOMER.value == "homer"
        assert PileupMethod.BAMNADO.value == "bamnado"


class TestDataScalingTechniqueEnum:
    """Tests for DataScalingTechnique enum."""

    def test_data_scaling_values(self):
        """Test DataScalingTechnique enum has expected values."""
        assert DataScalingTechnique.UNSCALED.value == "unscaled"
        assert DataScalingTechnique.CSAW.value == "csaw"
        assert DataScalingTechnique.CPM.value == "cpm"
        assert DataScalingTechnique.RPKM.value == "rpkm"
        assert DataScalingTechnique.SPIKEIN.value == "spikein"
        assert DataScalingTechnique.MERGED.value == "merged"


class TestPeakCallingMethodEnum:
    """Tests for PeakCallingMethod enum."""

    def test_peak_calling_method_values(self):
        """Test PeakCallingMethod enum has expected values."""
        assert PeakCallingMethod.MACS2.value == "macs2"
        assert PeakCallingMethod.MACS3.value == "macs3"
        assert PeakCallingMethod.HOMER.value == "homer"
        assert PeakCallingMethod.LANCEOTRON.value == "lanceotron"
        assert PeakCallingMethod.SEACR.value == "seacr"
        assert PeakCallingMethod.LANCEOTRON_MCC.value == "lanceotron-mcc"


class TestPCRDuplicateHandlingEnum:
    """Tests for PCRDuplicateHandling enum."""

    def test_pcr_duplicate_handling_values(self):
        """Test PCRDuplicateHandling enum has expected values."""
        assert PCRDuplicateHandling.REMOVE.value == "remove"
        assert PCRDuplicateHandling.MARK.value == "mark"
        assert PCRDuplicateHandling.NONE.value == "keep"


class TestPCRDuplicateToolEnum:
    """Tests for PCRDuplicateTool enum."""

    def test_pcr_duplicate_tool_values(self):
        """Test PCRDuplicateTool enum has expected values."""
        assert PCRDuplicateTool.PICARD.value == "picard"
        assert PCRDuplicateTool.SAMTOOLS.value == "samtools"
        assert PCRDuplicateTool.NONE.value == "None"


class TestSpikeInMethodEnum:
    """Tests for SpikeInMethod enum."""

    def test_spike_in_method_values(self):
        """Test SpikeInMethod enum has expected values."""
        assert SpikeInMethod.ORLANDO.value == "orlando"
        assert SpikeInMethod.WITH_INPUT.value == "with_input"


class TestSNPCallingMethodEnum:
    """Tests for SNPCallingMethod enum."""

    def test_snp_calling_method_values(self):
        """Test SNPCallingMethod enum has expected values."""
        assert SNPCallingMethod.BCFTOOLS.value == "bcftools"
        assert SNPCallingMethod.DEEPVARIANT.value == "deepvariant"


class TestQuantificationMethodEnum:
    """Tests for QuantificationMethod enum."""

    def test_quantification_method_values(self):
        """Test QuantificationMethod enum has expected values."""
        assert QuantificationMethod.FEATURE_COUNTS.value == "feature_counts"
        assert QuantificationMethod.SALMON.value == "salmon"


class TestMethylationMethodEnum:
    """Tests for MethylationMethod enum."""

    def test_methylation_method_values(self):
        """Test MethylationMethod enum has expected values."""
        assert MethylationMethod.TAPS.value == "taps"
        assert MethylationMethod.BISULFITE.value == "bisulfite"


class TestMoleculeEnum:
    """Tests for Molecule enum."""

    def test_molecule_values(self):
        """Test Molecule enum has expected values."""
        assert Molecule.rna_total.value == "total RNA"
        assert Molecule.rna_polya.value == "polyA RNA"
        assert Molecule.dna_genomic.value == "genomic DNA"
        assert Molecule.protein.value == "protein"


class TestOrganismEnum:
    """Tests for Organism enum."""

    def test_organism_values(self):
        """Test Organism enum has expected values."""
        assert Organism.HUMAN.value == "Homo sapiens"
        assert Organism.MOUSE.value == "Mus musculus"
        assert Organism.RAT.value == "Rattus norvegicus"
        assert Organism.DROSOPHILA.value == "Drosophila melanogaster"
        assert Organism.UNKNOWN.value == "Unknown"


class TestLibraryTypeEnum:
    """Tests for LibraryType enum."""

    def test_library_type_values(self):
        """Test LibraryType enum has expected values."""
        assert LibraryType.SINGLE.value == "single-end"
        assert LibraryType.PAIRED.value == "paired-end"


class TestGenomicCoordinate:
    """Tests for GenomicCoordinate model."""

    def test_valid_genomic_coordinate(self):
        """Test creating a valid genomic coordinate."""
        coord = GenomicCoordinate(chromosome="chr1", start=1000, end=2000)
        assert coord.chromosome == "chr1"
        assert coord.start == 1000
        assert coord.end == 2000

    def test_negative_start_raises_error(self):
        """Test that negative start coordinate raises error."""
        with pytest.raises(ValidationError):
            GenomicCoordinate(chromosome="chr1", start=-100, end=2000)

    def test_negative_end_raises_error(self):
        """Test that negative end coordinate raises error."""
        with pytest.raises(ValidationError, match="non-negative"):
            GenomicCoordinate(chromosome="chr1", start=1000, end=-2000)

    def test_end_less_than_start_raises_error(self):
        """Test that end < start raises error."""
        with pytest.raises(
            ValidationError, match="End coordinate must be greater than start coordinate"
        ):
            GenomicCoordinate(chromosome="chr1", start=2000, end=1000)

    def test_from_string_valid_format(self):
        """Test creating coordinate from string representation."""
        coord = GenomicCoordinate.from_string("chr1:1000-2000")
        assert coord.chromosome == "chr1"
        assert coord.start == 1000
        assert coord.end == 2000

    def test_from_string_different_chromosome(self):
        """Test from_string with different chromosome names."""
        coord = GenomicCoordinate.from_string("chrX:5000-10000")
        assert coord.chromosome == "chrX"
        assert coord.start == 5000
        assert coord.end == 10000

    def test_from_string_validates_coordinates(self):
        """Test that from_string validates coordinates."""
        with pytest.raises(ValidationError):
            GenomicCoordinate.from_string("chr1:2000-1000")
