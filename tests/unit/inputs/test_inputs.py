"""Comprehensive tests for seqnado.inputs module.

This test file covers all input-related modules:
- core: Base classes, metadata, utility functions
- fastq: FASTQ file handling
- bam: BAM file handling
- bigwigs: BigWig file handling
- grouping: Sample grouping
- helpers: Helper functions
- interfaces: Collection protocols
- validation: Data validation schemas
"""

import csv
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from pandera.errors import SchemaError

from seqnado import Assay, Organism
from seqnado.inputs.core import (
    Metadata,
    BaseCollection,
    predict_organism,
    is_valid_path,
    clean_sample_name,
    extract_read_number,
    is_control_sample,
    ILLUMINA_FILENAME_PATTERNS,
)
from seqnado.inputs.fastq import (
    FastqFile,
    FastqSet,
    FastqCollection,
    FastqFileIP,
    FastqSetIP,
    FastqCollectionForIP,
    ExperimentIP,
)
from seqnado.inputs.bam import BamFile, BamCollection
from seqnado.inputs.bigwigs import BigWigFile, BigWigCollection
from seqnado.inputs.grouping import SampleGroup, SampleGroupings, SampleGroups
from seqnado.inputs.helpers import get_sample_collection
from seqnado.inputs.interfaces import (
    CollectionLike,
    ensure_collection,
    detect_file_type,
)
from seqnado.inputs.validation import DesignDataFrame, ViewpointsFile


# =============================================================================
# Helper Functions
# =============================================================================

def _write_fastq(tmp: Path, name: str) -> Path:
    """Create a minimal FASTQ file for testing."""
    p = tmp / name
    p.write_text("@r\nN\n+\n#\n")
    return p


# =============================================================================
# CORE MODULE TESTS
# =============================================================================

class TestMetadata:
    """Tests for Metadata class."""

    def test_metadata_basic(self):
        """Test basic Metadata creation."""
        md = Metadata(assay=Assay.CHIP)
        assert md.assay == Assay.CHIP
        assert md.scaling_group == "default"

    def test_metadata_with_all_fields(self):
        """Test Metadata with all fields set."""
        md = Metadata(
            assay=Assay.RNA,
            consensus_group="group1",
            scaling_group="batch1",
            deseq2=1
        )
        assert md.assay == Assay.RNA
        assert md.consensus_group == "group1"
        assert md.scaling_group == "batch1"
        assert md.deseq2 == 1

    def test_metadata_optional_fields(self):
        """Test that optional fields can be omitted."""
        md = Metadata(assay=Assay.CHIP)
        assert md.consensus_group is None
        assert md.deseq2 is None

    def test_metadata_pandas_na_rejected(self):
        """Test that pandas NA is rejected."""
        from pydantic import ValidationError
        with pytest.raises(ValidationError):
            Metadata(assay=Assay.CHIP, consensus_group=pd.NA)

    def test_metadata_retains_extra_fields(self):
        md = Metadata(assay=Assay.CHIP, donor="D1", treatment="VTP")
        assert md.model_dump()["donor"] == "D1"
        assert md.model_dump()["treatment"] == "VTP"


class TestUtilityFunctions:
    """Tests for utility functions in core module."""

    def test_predict_organism_human(self):
        """Test predict_organism for human genomes."""
        assert predict_organism("hg38") == Organism.HUMAN
        assert predict_organism("hg19") == Organism.HUMAN

    def test_predict_organism_mouse(self):
        """Test predict_organism for mouse genomes."""
        assert predict_organism("mm10") == Organism.MOUSE
        assert predict_organism("mm39") == Organism.MOUSE

    def test_predict_organism_unknown(self):
        """Test predict_organism for unknown genomes."""
        assert predict_organism("dm6") == Organism.UNKNOWN
        assert predict_organism("custom") == Organism.UNKNOWN

    def test_is_valid_path_exists(self, tmp_path):
        """Test is_valid_path for existing file."""
        f = tmp_path / "test.txt"
        f.touch()
        assert is_valid_path(f) is True

    def test_is_valid_path_not_exists(self, tmp_path):
        """Test is_valid_path for non-existent file."""
        f = tmp_path / "missing.txt"
        assert is_valid_path(f) is False

    def test_is_valid_path_none(self):
        """Test is_valid_path with None."""
        assert is_valid_path(None) is False

    def test_is_valid_path_invalid_strings(self):
        """Test is_valid_path with invalid string values."""
        assert is_valid_path("-") is False
        assert is_valid_path(".") is False
        assert is_valid_path("") is False
        assert is_valid_path("None") is False

    def test_is_valid_path_type_error(self):
        """Test is_valid_path with invalid types."""
        assert is_valid_path(123) is False
        assert is_valid_path([]) is False

    def test_clean_sample_name(self):
        """Test clean_sample_name with Illumina patterns."""
        assert clean_sample_name("sample_S1_L001_R1_001") == "sample"
        assert clean_sample_name("sample_S23_R2") == "sample"
        assert clean_sample_name("sample__extra") == "sample_extra"

    def test_extract_read_number(self):
        """Test extract_read_number from filenames."""
        assert extract_read_number("sample_R1_001") == 1
        assert extract_read_number("sample_R2_001") == 2
        assert extract_read_number("sample_1") == 1
        assert extract_read_number("sample_2") == 2
        assert extract_read_number("sample_single") is None

    def test_is_control_sample(self):
        """Test is_control_sample detection."""
        assert is_control_sample("Input") is True
        assert is_control_sample("mock_control") is True
        assert is_control_sample("IgG") is True
        assert is_control_sample("H3K4me3") is False


class TestBaseCollection:
    """Tests for BaseCollection class methods."""

    def test_build_metadata_callable(self):
        """Test _build_metadata with callable."""
        from seqnado.inputs.core import BaseCollection
        
        def metadata_factory(sample_name: str) -> Metadata:
            return Metadata(assay=Assay.CHIP, consensus_group=sample_name)
        
        md = BaseCollection._build_metadata(
            sample_name="sample1",
            metadata=metadata_factory,
            assay=Assay.CHIP
        )
        assert md.assay == Assay.CHIP
        assert md.consensus_group == "sample1"

    def test_build_metadata_instance(self):
        """Test _build_metadata with Metadata instance."""
        from seqnado.inputs.core import BaseCollection
        
        original_md = Metadata(assay=Assay.RNA, scaling_group="batch1")
        md = BaseCollection._build_metadata(
            sample_name="sample1",
            metadata=original_md,
            assay=Assay.RNA
        )
        assert md.assay == Assay.RNA
        assert md.scaling_group == "batch1"

    def test_build_metadata_default(self):
        """Test _build_metadata with None (default)."""
        from seqnado.inputs.core import BaseCollection
        
        md = BaseCollection._build_metadata(
            sample_name="sample1",
            metadata=None,
            assay=Assay.ATAC
        )
        assert md.assay == Assay.ATAC

    def test_build_metadata_assay_override(self):
        """Test that assay is always overridden."""
        from seqnado.inputs.core import BaseCollection
        
        original_md = Metadata(assay=Assay.RNA)
        md = BaseCollection._build_metadata(
            sample_name="sample1",
            metadata=original_md,
            assay=Assay.CHIP  # Different assay
        )
        assert md.assay == Assay.CHIP  # Should be overridden


# =============================================================================
# FASTQ MODULE TESTS
# =============================================================================

class TestFastqFile:
    """Tests for FastqFile class."""

    def test_fastqfile_basic(self, tmp_path):
        """Test basic FastqFile creation."""
        f = _write_fastq(tmp_path, "sample_R1.fastq.gz")
        fq = FastqFile(path=f)
        assert fq.path.is_absolute()
        assert fq.read_number == 1

    def test_fastqfile_parses_read_and_sample(self, tmp_path):
        """Test FastqFile parses read number and sample name."""
        f = _write_fastq(tmp_path, "chip_A_S1_L001_R1_001.fastq.gz")
        fq = FastqFile(path=f)
        assert fq.read_number == 1
        assert fq.is_paired is True
        assert clean_sample_name("chip_A_S1_L001_R1_001") == fq.sample_base

    def test_fastqfile_not_found(self, tmp_path):
        """Test FastqFile with missing file."""
        with pytest.raises(FileNotFoundError):
            FastqFile(path=tmp_path / "missing.fastq.gz")

    def test_fastqfile_stem_property(self, tmp_path):
        """Test FastqFile stem property."""
        f = _write_fastq(tmp_path, "sample.fastq.gz")
        fq = FastqFile(path=f)
        assert fq.stem == "sample"

    def test_fastqfile_comparison(self, tmp_path):
        """Test FastqFile comparison operators."""
        f1 = _write_fastq(tmp_path, "a_sample.fastq.gz")
        f2 = _write_fastq(tmp_path, "b_sample.fastq.gz")
        fq1 = FastqFile(path=f1)
        fq2 = FastqFile(path=f2)
        assert fq1 < fq2
        assert fq2 > fq1

    def test_fastqfile_equality_hash(self, tmp_path):
        """Test FastqFile equality and hash."""
        f = _write_fastq(tmp_path, "sample.fastq.gz")
        fq1 = FastqFile(path=f)
        fq2 = FastqFile(path=f)
        assert fq1 == fq2
        assert hash(fq1) == hash(fq2)

    def test_fastqfile_is_lane(self, tmp_path):
        """Test FastqFile is_lane property."""
        f_lane = _write_fastq(tmp_path, "sample_L001_R1.fastq.gz")
        f_no_lane = _write_fastq(tmp_path, "sample_R1.fastq.gz")
        
        fq_lane = FastqFile(path=f_lane)
        fq_no_lane = FastqFile(path=f_no_lane)
        
        assert fq_lane.is_lane is True
        assert fq_no_lane.is_lane is False

    def test_fastqfile_use_resolved_name(self, tmp_path):
        """Test FastqFile with use_resolved_name."""
        f = _write_fastq(tmp_path, "sample.fastq.gz")
        fq = FastqFile(path=f, use_resolved_name=True)
        assert fq.path.is_absolute()
        assert fq.path.exists()


class TestFastqFileIP:
    """Tests for FastqFileIP class."""

    def test_fastqfileip_predicts_ip_and_control(self, tmp_path):
        """Test FastqFileIP IP and control detection."""
        f = _write_fastq(tmp_path, "sampleX_IGG_L001_R1_001.fastq.gz")
        fq = FastqFileIP(path=f)
        assert fq.ip == "IGG"
        assert fq.is_control is True
        assert fq.sample_base_without_ip == "sampleX"

    def test_fastqfileip_non_control(self, tmp_path):
        """Test FastqFileIP with non-control IP."""
        f = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        fq = FastqFileIP(path=f)
        assert fq.ip == "H3K4me3"
        assert fq.is_control is False

    def test_fastqfileip_predict_ip_failure(self, tmp_path):
        """Test FastqFileIP IP prediction with malformed filename."""
        f = _write_fastq(tmp_path, "single.fastq.gz")
        fq = FastqFileIP(path=f)
        # Should still work, just use whole name
        assert fq.ip is not None

    def test_fastqfileip_explicit_ip_and_control(self, tmp_path):
        """Test FastqFileIP with explicit IP and control values."""
        f = _write_fastq(tmp_path, "sample_H3K27me3_R1.fastq.gz")
        fq = FastqFileIP(path=f, ip="custom_IP", is_control=True)
        assert fq.ip == "custom_IP"
        assert fq.is_control is True

    def test_fastqfileip_allow_na_validator(self, tmp_path):
        """Test FastqFileIP validator allows NA/NaN values."""
        f = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        
        # IP is auto-predicted from filename, test it works
        fq1 = FastqFileIP(path=f)
        assert fq1.ip == "H3K4me3"
        
        # Test with string IP value
        fq2 = FastqFileIP(path=f, ip="custom_IP")
        assert fq2.ip == "custom_IP"

    def test_fastqfileip_sample_base_without_ip_complex(self, tmp_path):
        """Test sample_base_without_ip with complex suffixes."""
        f = _write_fastq(tmp_path, "sample123_IGG_S5_L001_R1_001.fastq.gz")
        fq = FastqFileIP(path=f)
        # Should remove IP and Illumina suffixes
        assert "IGG" not in fq.sample_base_without_ip
        assert "S5" not in fq.sample_base_without_ip
        assert fq.sample_base_without_ip == "sample123"

    def test_fastqfileip_equality_and_hash(self, tmp_path):
        """Test FastqFileIP equality and hashing."""
        f1 = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        f2 = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        fq1 = FastqFileIP(path=f1)
        fq2 = FastqFileIP(path=f2)
        assert fq1 == fq2
        assert hash(fq1) == hash(fq2)


class TestFastqSet:
    """Tests for FastqSet class."""

    def test_fastqset_single_end(self, tmp_path):
        """Test FastqSet with single-end data."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        fs = FastqSet(r1=r1)
        assert fs.is_paired is False
        assert fs.file_paths == [r1.path]

    def test_fastqset_paired_end(self, tmp_path):
        """Test FastqSet with paired-end data."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fs = FastqSet(r1=r1, r2=r2)
        assert fs.is_paired is True
        assert set(fs.file_paths) == {r1.path, r2.path}

    def test_fastqset_with_explicit_sample_id(self, tmp_path):
        """Test FastqSet with explicit sample_id."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "sample1_R1.fastq.gz"))
        fs = FastqSet(sample_id="sample1", r1=r1)
        assert fs.sample_id == "sample1"

    def test_fastqset_from_fastq_files_single(self, tmp_path):
        """Test FastqSet.from_fastq_files with single-end."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        fs = FastqSet.from_fastq_files([r1])
        assert fs.sample_id == "s1"
        assert fs.is_paired is False

    def test_fastqset_from_fastq_files_paired(self, tmp_path):
        """Test FastqSet.from_fastq_files with paired-end."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fs = FastqSet.from_fastq_files([r1, r2])
        assert fs.sample_id == "s1"
        assert fs.is_paired is True

    def test_fastqset_from_fastq_files_empty(self):
        """Test FastqSet.from_fastq_files with no files."""
        with pytest.raises(ValueError, match="No FASTQ files provided"):
            FastqSet.from_fastq_files([])

    def test_fastqset_from_fastq_files_too_many(self, tmp_path):
        """Test FastqSet.from_fastq_files with too many files."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        r3 = FastqFile(path=_write_fastq(tmp_path, "s1_R3.fastq.gz"))
        with pytest.raises(ValueError, match="Invalid number of FASTQ files"):
            FastqSet.from_fastq_files([r1, r2, r3])


class TestFastqSetIP:
    """Tests for FastqSetIP class."""

    def test_fastqsetip_predicts_antibody(self, tmp_path):
        """Test FastqSetIP antibody prediction."""
        r1 = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        fs = FastqSetIP(sample_id="sample", r1=r1)
        assert fs.ip == "H3K4me3"

    def test_fastqsetip_explicit_antibody(self, tmp_path):
        """Test FastqSetIP with explicit antibody."""
        r1 = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        fs = FastqSetIP(sample_id="sample", r1=r1, ip="custom_AB")
        assert fs.ip == "custom_AB"

    def test_fastqsetip_full_sample_name(self, tmp_path):
        """Test FastqSetIP full_sample_name property."""
        r1 = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K27me3_R1.fastq.gz"))
        fs = FastqSetIP(sample_id="sample", r1=r1)
        assert fs.full_sample_name == "sample_H3K27me3"

    def test_fastqsetip_base_sample_name(self, tmp_path):
        """Test FastqSetIP base_sample_name property."""
        r1 = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        fs = FastqSetIP(sample_id="sample", r1=r1)
        assert fs.base_sample_name == r1.sample_base_without_ip

    def test_fastqsetip_is_control(self, tmp_path):
        """Test FastqSetIP is_control property."""
        r1_control = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        
        fs_control = FastqSetIP(sample_id="sample", r1=r1_control)
        fs_ip = FastqSetIP(sample_id="sample", r1=r1_ip)
        
        assert fs_control.is_control is True
        assert fs_ip.is_control is False


class TestExperimentIP:
    """Tests for ExperimentIP class."""

    def test_experimentip_basic(self, tmp_path):
        """Test ExperimentIP creation."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        assert exp.has_control is True
        assert exp.ip_performed == "H3K4me3"
        assert exp.control_performed == "IGG"

    def test_experimentip_no_control(self, tmp_path):
        """Test ExperimentIP without control."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        
        exp = ExperimentIP(ip=ip_set)
        assert exp.has_control is False
        assert exp.control_fullname is None
        assert exp.control_performed is None

    def test_experimentip_fullnames(self, tmp_path):
        """Test ExperimentIP fullname properties."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K27me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_INPUT_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        assert exp.ip_set_fullname == "sample_H3K27me3"
        assert exp.control_fullname == "sample_INPUT"

    def test_experimentip_fastqs_are_paired(self, tmp_path):
        """Test ExperimentIP fastqs_are_paired property."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r2_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        r2_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R2.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip, r2=r2_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl, r2=r2_ctrl)
        
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        assert exp.fastqs_are_paired is True

    def test_experimentip_fastqs_single_end(self, tmp_path):
        """Test ExperimentIP fastqs_are_paired with single-end."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        assert exp.fastqs_are_paired is False


class TestFastqCollection:
    """Tests for FastqCollection class."""

    def test_fastqcollection_basic(self, tmp_path):
        """Test basic FastqCollection creation."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
        )
        assert fc.assay == Assay.RNA
        assert len(fc.fastq_sets) == 1

    def test_fastqcollection_symlink(self, tmp_path):
        """Test FastqCollection symlink creation."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
        )

        out = tmp_path / "links"
        fc.symlink_fastq_files(out)
        assert (out / "s1_1.fastq.gz").is_symlink()
        assert (out / "s1_2.fastq.gz").is_symlink()

    def test_fastqcollection_fastq_pairs(self, tmp_path):
        """Test FastqCollection fastq_pairs property."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
        )
        pairs = fc.fastq_pairs
        assert "s1" in pairs
        assert len(pairs["s1"]) == 2

    def test_fastqcollection_query(self, tmp_path):
        """Test FastqCollection query method."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1)],
        )
        fs = fc.query("s1")
        assert fs.sample_id == "s1"

    def test_fastqcollection_query_not_found(self, tmp_path):
        """Test FastqCollection query with missing sample."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1)],
        )
        with pytest.raises(ValueError, match="not found"):
            fc.query("missing")

    def test_fastqcollection_is_paired_end(self, tmp_path):
        """Test FastqCollection is_paired_end method."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[Metadata()],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
        )
        df = fc.to_dataframe()
        uid = df.index[0]
        assert fc.is_paired_end(uid) is True

    def test_fastqcollection_from_fastq_files(self, tmp_path):
        """Test FastqCollection.from_fastq_files."""
        r1 = _write_fastq(tmp_path, "s1_R1.fastq.gz")
        r2 = _write_fastq(tmp_path, "s1_R2.fastq.gz")
        fc = FastqCollection.from_fastq_files(
            assay=Assay.RNA,
            files=[r1, r2],
        )
        assert len(fc.fastq_sets) == 1
        assert fc.fastq_sets[0].is_paired is True

    def test_fastqcollection_round_trips_custom_metadata(self, tmp_path):
        r1 = _write_fastq(tmp_path, "sample_R1.fastq.gz")
        r2 = _write_fastq(tmp_path, "sample_R2.fastq.gz")
        df = pd.DataFrame(
            [{"sample_id": "sample", "r1": r1, "r2": r2, "donor": "D1"}]
        )

        collection = FastqCollection.from_dataframe(assay=Assay.RNA, df=df)
        exported = collection.to_dataframe(validate=False)

        assert exported.iloc[0]["donor"] == "D1"

    def test_fastqcollection_from_directory(self, tmp_path):
        """Test FastqCollection.from_directory."""
        _write_fastq(tmp_path, "s1_R1.fastq.gz")
        _write_fastq(tmp_path, "s1_R2.fastq.gz")
        fc = FastqCollection.from_directory(
            assay=Assay.RNA,
            directory=tmp_path,
        )
        assert len(fc.fastq_sets) == 1

    def test_fastqcollection_from_directory_no_files(self, tmp_path):
        """Test FastqCollection.from_directory with no FASTQ files."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        with pytest.raises(FileNotFoundError, match="No FASTQ files found"):
            FastqCollection.from_directory(
                assay=Assay.RNA,
                directory=empty_dir,
            )

    def test_fastqcollection_validate_non_ip_assay(self, tmp_path):
        """Test FastqCollection rejects IP assays."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        with pytest.raises(ValueError, match="requires IP"):
            FastqCollection(
                assay=Assay.CHIP,  # IP assay
                metadata=[],
                fastq_sets=[FastqSet(sample_id="s1", r1=r1)],
            )

    def test_fastqcollection_sample_ids_and_names(self, tmp_path):
        """Test FastqCollection sample_ids and sample_names properties."""
        r1_s1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r1_s2 = FastqFile(path=_write_fastq(tmp_path, "s2_R1.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[
                FastqSet(sample_id="s1", r1=r1_s1),
                FastqSet(sample_id="s2", r1=r1_s2),
            ],
        )
        assert fc.sample_ids == ["s1", "s2"]
        assert fc.sample_names == ["s1", "s2"]

    def test_fastqcollection_primary_file_type(self, tmp_path):
        """Test FastqCollection primary_file_type property."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1)],
        )
        assert fc.primary_file_type == "fastq"

    def test_fastqcollection_fastq_paths(self, tmp_path):
        """Test FastqCollection fastq_paths property."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
        )
        paths = fc.fastq_paths
        assert len(paths) == 2
        assert r1.path in paths
        assert r2.path in paths

    def test_fastqcollection_get_file_paths_variants(self, tmp_path):
        """Test FastqCollection get_file_paths with different kinds."""
        r1 = FastqFile(path=_write_fastq(tmp_path, "s1_R1.fastq.gz"))
        r2 = FastqFile(path=_write_fastq(tmp_path, "s1_R2.fastq.gz"))
        fc = FastqCollection(
            assay=Assay.RNA,
            metadata=[],
            fastq_sets=[FastqSet(sample_id="s1", r1=r1, r2=r2)],
        )
        
        # Test fastq_r1
        r1_paths = fc.get_file_paths("fastq_r1")
        assert len(r1_paths) == 1
        assert r1.path in r1_paths
        
        # Test fastq_r2
        r2_paths = fc.get_file_paths("fastq_r2")
        assert len(r2_paths) == 1
        assert r2.path in r2_paths
        
        # Test invalid kind
        with pytest.raises(ValueError, match="Unsupported file kind"):
            fc.get_file_paths("invalid")


class TestFastqCollectionForIP:
    """Tests for FastqCollectionForIP class."""

    def test_fastqcollectionforip_basic(self, tmp_path):
        """Test basic FastqCollectionForIP creation."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        assert fc.assay == Assay.CHIP
        assert len(fc.experiments) == 1

    def test_fastqcollectionforip_fastq_pairs(self, tmp_path):
        """Test FastqCollectionForIP fastq_pairs property."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r2_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        r2_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R2.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip, r2=r2_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl, r2=r2_ctrl)
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        pairs = fc.fastq_pairs
        assert "sample_H3K4me3" in pairs
        assert "sample_IGG" in pairs
        assert len(pairs["sample_H3K4me3"]) == 2
        assert len(pairs["sample_IGG"]) == 2

    def test_fastqcollectionforip_get_file_paths_fastq(self, tmp_path):
        """Test FastqCollectionForIP get_file_paths for all FASTQs."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        paths = fc.get_file_paths("fastq")
        assert len(paths) == 1
        assert r1_ip.path in paths

    def test_fastqcollectionforip_get_file_paths_r1(self, tmp_path):
        """Test FastqCollectionForIP get_file_paths for R1."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r2_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip, r2=r2_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        paths = fc.get_file_paths("fastq_r1")
        assert len(paths) == 1
        assert r1_ip.path in paths

    def test_fastqcollectionforip_get_file_paths_r2(self, tmp_path):
        """Test FastqCollectionForIP get_file_paths for R2."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r2_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip, r2=r2_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        paths = fc.get_file_paths("fastq_r2")
        assert len(paths) == 1
        assert r2_ip.path in paths

    def test_fastqcollectionforip_get_file_paths_invalid(self, tmp_path):
        """Test FastqCollectionForIP get_file_paths with invalid kind."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        with pytest.raises(ValueError, match="Unsupported file kind"):
            fc.get_file_paths("invalid_kind")

    def test_fastqcollectionforip_get_control_performed(self, tmp_path):
        """Test FastqCollectionForIP get_control_performed method."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        # For IP sample
        assert fc.get_control_performed("sample_H3K4me3") == "sample_IGG"
        # For control sample
        assert fc.get_control_performed("sample_IGG") == "sample_IGG"

    def test_fastqcollectionforip_get_control_performed_no_control(self, tmp_path):
        """Test get_control_performed with no control."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        assert fc.get_control_performed("sample_H3K4me3") is None

    def test_fastqcollectionforip_get_control_performed_not_found(self, tmp_path):
        """Test get_control_performed with missing sample."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        with pytest.raises(KeyError, match="not found"):
            fc.get_control_performed("missing_sample")

    def test_fastqcollectionforip_is_paired_end(self, tmp_path):
        """Test FastqCollectionForIP is_paired_end method."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r2_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip, r2=r2_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        assert fc.is_paired_end("sample_H3K4me3") is True

    def test_fastqcollectionforip_is_paired_end_not_found(self, tmp_path):
        """Test is_paired_end with missing sample."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        with pytest.raises(KeyError, match="not found"):
            fc.is_paired_end("missing_sample")

    def test_fastqcollectionforip_query_ip(self, tmp_path):
        """Test FastqCollectionForIP query for IP sample."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        result = fc.query("sample_H3K4me3")
        assert isinstance(result, FastqSetIP)
        assert result.ip == "H3K4me3"

    def test_fastqcollectionforip_query_full(self, tmp_path):
        """Test FastqCollectionForIP query with full=True."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        result = fc.query("sample_H3K4me3", full=True)
        assert isinstance(result, dict)
        assert "ip" in result
        assert "control" in result

    def test_fastqcollectionforip_query_not_found(self, tmp_path):
        """Test FastqCollectionForIP query with missing sample."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        exp = ExperimentIP(ip=ip_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        with pytest.raises(KeyError, match="not found"):
            fc.query("missing_sample")

    def test_fastqcollectionforip_from_fastq_files(self, tmp_path):
        """Test FastqCollectionForIP.from_fastq_files."""
        r1_ip = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r2_ip = _write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz")
        r1_ctrl = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        r2_ctrl = _write_fastq(tmp_path, "sample_IGG_R2.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip, r2_ip, r1_ctrl, r2_ctrl],
        )
        
        assert len(fc.experiments) == 1
        assert fc.experiments[0].has_control is True
        assert fc.experiments[0].ip.is_paired is True
        assert fc.experiments[0].control.is_paired is True

    def test_fastqcollectionforip_from_fastq_files_no_control(self, tmp_path):
        """Test FastqCollectionForIP.from_fastq_files without control."""
        r1_ip = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r2_ip = _write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip, r2_ip],
        )
        
        assert len(fc.experiments) == 1
        assert fc.experiments[0].has_control is False

    def test_fastqcollectionforip_from_fastq_files_edge_case_only_control(self, tmp_path):
        """Test from_fastq_files with only control files (no IP)."""
        r1_ctrl = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        r2_ctrl = _write_fastq(tmp_path, "sample_IGG_R2.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ctrl, r2_ctrl],
        )
        
        # Should treat control as IP when no IP files found
        assert len(fc.experiments) == 1
        assert fc.experiments[0].has_control is False

    def test_fastqcollectionforip_from_fastq_files_invalid_count(self, tmp_path):
        """Test from_fastq_files with mismatched IP/control counts."""
        r1_ip = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r2_ip = _write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz")
        r1_ctrl1 = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        r2_ctrl1 = _write_fastq(tmp_path, "sample_IGG_R2.fastq.gz")
        r1_ctrl2 = _write_fastq(tmp_path, "sample_INPUT_R1.fastq.gz")
        
        # This should raise an error because control has too many files
        with pytest.raises(ValueError):
            FastqCollectionForIP.from_fastq_files(
                assay=Assay.CHIP,
                files=[r1_ip, r2_ip, r1_ctrl1, r2_ctrl1, r1_ctrl2],
            )


    def test_ip_to_control_pairing(self, tmp_path):
        """Test from_fastq_files with IP to control pairing."""
        r1_ip = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r2_ip = _write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz")
        r1_ctrl1 = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        r2_ctrl1 = _write_fastq(tmp_path, "sample_IGG_R2.fastq.gz")
        r1_ctrl2 = _write_fastq(tmp_path, "sample_INPUT_R1.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip, r2_ip, r1_ctrl1, r2_ctrl1, r1_ctrl2],
            ip_to_control_map={
                "H3K4me3": "IGG",
            }
        )
        assert len(fc.experiments) == 1
        assert fc.experiments[0].has_control is True
    
    def test_ip_to_control_pairing_more_complex(self, tmp_path):
        """Test from_fastq_files with more complex IP to control pairing."""
        r1_ip1 = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r2_ip1 = _write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz")
        r1_ip2 = _write_fastq(tmp_path, "sample_H3K27me3_R1.fastq.gz")
        r2_ip2 = _write_fastq(tmp_path, "sample_H3K27me3_R2.fastq.gz")
        r1_ctrl1 = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        r2_ctrl1 = _write_fastq(tmp_path, "sample_IGG_R2.fastq.gz")
        r1_ctrl2 = _write_fastq(tmp_path, "sample_INPUT_R1.fastq.gz")
        r2_ctrl2 = _write_fastq(tmp_path, "sample_INPUT_R2.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip1, r2_ip1, r1_ip2, r2_ip2, r1_ctrl1, r2_ctrl1, r1_ctrl2, r2_ctrl2],
            ip_to_control_map={
                "H3K4me3": "IGG",
                "H3K27me3": "INPUT",
            }
        )
        assert len(fc.experiments) == 2
        for exp in fc.experiments:
            assert exp.has_control is True
        
    def test_ip_to_control_pairing_more_samples(self, tmp_path):
        """Test from_fastq_files with multiple samples and IP to control pairing."""
        r1_ip1 = _write_fastq(tmp_path, "sample1_H3K4me3_R1.fastq.gz")
        r2_ip1 = _write_fastq(tmp_path, "sample1_H3K4me3_R2.fastq.gz")
        r1_ip2 = _write_fastq(tmp_path, "sample2_H3K27me3_R1.fastq.gz")
        r2_ip2 = _write_fastq(tmp_path, "sample2_H3K27me3_R2.fastq.gz")
        r1_ctrl1 = _write_fastq(tmp_path, "sample1_IGG_R1.fastq.gz")
        r2_ctrl1 = _write_fastq(tmp_path, "sample1_IGG_R2.fastq.gz")
        r1_ctrl2 = _write_fastq(tmp_path, "sample2_INPUT_R1.fastq.gz")
        r2_ctrl2 = _write_fastq(tmp_path, "sample2_INPUT_R2.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip1, r2_ip1, r1_ip2, r2_ip2, r1_ctrl1, r2_ctrl1, r1_ctrl2, r2_ctrl2],
            ip_to_control_map={
                "H3K4me3": "IGG",
                "H3K27me3": "INPUT",
            }
        )
        assert len(fc.experiments) == 2
        for exp in fc.experiments:
            assert exp.has_control is True
        
    
    def test_ip_to_control_pairing_more_samples_one_missing_control(self, tmp_path):
        """Test from_fastq_files with multiple samples where one control is missing."""
        r1_ip1 = _write_fastq(tmp_path, "sample1_H3K4me3_R1.fastq.gz")
        r2_ip1 = _write_fastq(tmp_path, "sample1_H3K4me3_R2.fastq.gz")
        r1_ip2 = _write_fastq(tmp_path, "sample2_H3K27me3_R1.fastq.gz")
        r2_ip2 = _write_fastq(tmp_path, "sample2_H3K27me3_R2.fastq.gz")
        r1_ctrl1 = _write_fastq(tmp_path, "sample1_IGG_R1.fastq.gz")
        r2_ctrl1 = _write_fastq(tmp_path, "sample1_IGG_R2.fastq.gz")
        # Missing control files for sample2
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip1, r2_ip1, r1_ip2, r2_ip2, r1_ctrl1, r2_ctrl1],
            ip_to_control_map={
                "H3K4me3": "IGG",
                "H3K27me3": "INPUT",
            }
        )
        assert len(fc.experiments) == 2
        for exp in fc.experiments:
            if exp.ip.ip == "H3K4me3":
                assert exp.has_control is True
            else:
                assert exp.has_control is False
    
    def test_control_is_broadcasted(self, tmp_path):
        """Test from_fastq_files with single control for multiple IPs."""
        r1_ip1 = _write_fastq(tmp_path, "sample1_H3K4me3_R1.fastq.gz")
        r2_ip1 = _write_fastq(tmp_path, "sample1_H3K4me3_R2.fastq.gz")
        r1_ip2 = _write_fastq(tmp_path, "sample1_H3K27me3_R1.fastq.gz")
        r2_ip2 = _write_fastq(tmp_path, "sample1_H3K27me3_R2.fastq.gz")
        r1_ctrl = _write_fastq(tmp_path, "sample1_IGG_R1.fastq.gz")
        r2_ctrl = _write_fastq(tmp_path, "sample1_IGG_R2.fastq.gz")
        
        fc = FastqCollectionForIP.from_fastq_files(
            assay=Assay.CHIP,
            files=[r1_ip1, r2_ip1, r1_ip2, r2_ip2, r1_ctrl, r2_ctrl],
        )
        assert len(fc.experiments) == 2
        for exp in fc.experiments:
            assert exp.has_control is True
   

    

    def test_fastqcollectionforip_from_directory(self, tmp_path):
        """Test FastqCollectionForIP.from_directory."""
        _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        
        fc = FastqCollectionForIP.from_directory(
            assay=Assay.CHIP,
            directory=tmp_path,
        )
        
        assert len(fc.experiments) == 1
        assert fc.experiments[0].has_control is True

    def test_fastqcollectionforip_to_dataframe(self, tmp_path):
        """Test FastqCollectionForIP.to_dataframe."""
        r1_ip = FastqFileIP(path=_write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz"))
        r1_ctrl = FastqFileIP(path=_write_fastq(tmp_path, "sample_IGG_R1.fastq.gz"))
        
        ip_set = FastqSetIP(sample_id="sample", r1=r1_ip)
        ctrl_set = FastqSetIP(sample_id="sample", r1=r1_ctrl)
        exp = ExperimentIP(ip=ip_set, control=ctrl_set)
        
        fc = FastqCollectionForIP(
            assay=Assay.CHIP,
            experiments=[exp],
            metadata=[Metadata()],
        )
        
        df = fc.to_dataframe()
        assert "sample_id" in df.columns
        assert "ip" in df.columns
        assert "control" in df.columns
        assert "r1" in df.columns
        assert "r1_control" in df.columns

    def test_fastqcollectionforip_from_dataframe(self, tmp_path):
        """Test FastqCollectionForIP.from_dataframe."""
        r1_ip = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r2_ip = _write_fastq(tmp_path, "sample_H3K4me3_R2.fastq.gz")
        r1_ctrl = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        r2_ctrl = _write_fastq(tmp_path, "sample_IGG_R2.fastq.gz")
        
        df = pd.DataFrame([{
            "sample_id": "sample",
            "ip": "H3K4me3",
            "control": "IGG",
            "r1": str(r1_ip),
            "r2": str(r2_ip),
            "r1_control": str(r1_ctrl),
            "r2_control": str(r2_ctrl),
        }])
        
        fc = FastqCollectionForIP.from_dataframe(assay=Assay.CHIP, df=df)
        
        assert len(fc.experiments) == 1
        assert fc.experiments[0].has_control is True
        assert fc.experiments[0].ip.is_paired is True
        assert fc.experiments[0].control.is_paired is True

    def test_fastqcollectionforip_round_trips_custom_metadata(self, tmp_path):
        r1_ip = _write_fastq(tmp_path, "sample_H3K4me3_R1.fastq.gz")
        r1_ctrl = _write_fastq(tmp_path, "sample_IGG_R1.fastq.gz")
        df = pd.DataFrame(
            [{
                "sample_id": "sample",
                "ip": "H3K4me3",
                "control": "IGG",
                "r1": r1_ip,
                "r1_control": r1_ctrl,
                "donor": "D1",
            }]
        )

        collection = FastqCollectionForIP.from_dataframe(assay=Assay.CHIP, df=df)
        exported = collection.to_dataframe()

        assert exported.iloc[0]["donor"] == "D1"


# =============================================================================
# BAM MODULE TESTS
# =============================================================================

class TestBamFile:
    """Tests for BamFile class."""

    def test_create_bam_file(self, tmp_path):
        """Test creating a BamFile."""
        bam_path = tmp_path / "sample1.bam"
        bam_path.touch()
        
        bam = BamFile(path=bam_path)
        assert bam.path.is_absolute()
        assert bam.sample_id == "sample1"

    def test_bam_file_not_found(self, tmp_path):
        """Test BamFile with missing file."""
        with pytest.raises(FileNotFoundError, match="BAM file not found"):
            BamFile(path=tmp_path / "missing.bam")


class TestBamCollection:
    """Tests for BamCollection class."""

    def test_bam_collection_from_files(self, tmp_path):
        """Test BamCollection.from_files."""
        bam1 = tmp_path / "sample1.bam"
        bam2 = tmp_path / "sample2.bam"
        bam1.touch()
        bam2.touch()
        
        collection = BamCollection.from_files(
            assay=Assay.ATAC,
            files=[bam1, bam2]
        )
        
        assert len(collection.bam_files) == 2
        assert collection.assay == Assay.ATAC

    def test_bam_collection_from_directory(self, tmp_path):
        """Test BamCollection.from_directory."""
        bam1 = tmp_path / "sample1.bam"
        bam2 = tmp_path / "sample2.bam"
        bam1.touch()
        bam2.touch()
        
        collection = BamCollection.from_directory(
            assay=Assay.ATAC,
            directory=tmp_path
        )
        
        assert len(collection.bam_files) == 2


# =============================================================================
# BIGWIG MODULE TESTS
# =============================================================================

class TestBigWigFile:
    """Tests for BigWigFile class."""

    def test_create_bigwig_file(self, tmp_path):
        """Test creating a BigWigFile."""
        bw_path = tmp_path / "sample1.bigWig"
        bw_path.touch()
        
        bw = BigWigFile(path=bw_path)
        assert bw.path.is_absolute()
        assert bw.stem == "sample1"

    def test_bigwig_file_not_found(self, tmp_path):
        """Test BigWigFile with missing file."""
        with pytest.raises(FileNotFoundError, match="bigWig file not found"):
            BigWigFile(path=tmp_path / "missing.bigWig")

    def test_infer_sample_id_rna_stranded(self, tmp_path):
        """Test inferring sample ID from RNA stranded files."""
        bw_plus = tmp_path / "sample1_plus.bigWig"
        bw_minus = tmp_path / "sample1_minus.bigWig"
        bw_plus.touch()
        bw_minus.touch()
        
        bw_p = BigWigFile(path=bw_plus)
        bw_m = BigWigFile(path=bw_minus)
        
        assert bw_p.infer_sample_id(Assay.RNA) == "sample1"
        assert bw_m.infer_sample_id(Assay.RNA) == "sample1"

    def test_strand_detection(self, tmp_path):
        """Test strand detection for RNA files."""
        bw_plus = tmp_path / "sample1_plus.bigWig"
        bw_plus.touch()
        bw = BigWigFile(path=bw_plus)
        
        assert bw.is_strand_specific(Assay.RNA) is True
        assert bw.strand(Assay.RNA) == "plus"
        assert bw.is_strand_specific(Assay.CHIP) is False


class TestBigWigCollection:
    """Tests for BigWigCollection class."""

    def test_bigwig_collection_chip(self, tmp_path):
        """Test BigWigCollection for ChIP assay."""
        bw1 = tmp_path / "sample1.bigWig"
        bw2 = tmp_path / "sample2.bigWig"
        bw1.touch()
        bw2.touch()
        
        collection = BigWigCollection.from_files(
            assay=Assay.CHIP,
            files=[bw1, bw2]
        )
        
        assert collection.assay == Assay.CHIP
        assert len(collection.bigwig_files) == 2
        assert collection.primary_file_type == "bigwig"

    def test_bigwig_collection_rna_stranded(self, tmp_path):
        """Test BigWigCollection for RNA with stranded files."""
        bw_plus = tmp_path / "sample1_plus.bigWig"
        bw_minus = tmp_path / "sample1_minus.bigWig"
        bw_plus.touch()
        bw_minus.touch()
        
        collection = BigWigCollection.from_directory(
            assay=Assay.RNA,
            directory=tmp_path
        )
        
        # Both files are detected, sample ID appears twice (plus and minus)
        assert len(collection.sample_ids) == 2
        assert all(sid == "sample1" for sid in collection.sample_ids)


# =============================================================================
# GROUPING MODULE TESTS
# =============================================================================

class TestSampleGrouping:
    """Tests for sample grouping classes."""

    def test_samplegroup_basic(self):
        """Test basic SampleGroup creation."""
        sg = SampleGroup(name="group1", samples=["s1", "s2", "s3"])
        assert sg.name == "group1"
        assert len(sg.samples) == 3

    def test_samplegroups_from_dataframe(self):
        """Test SampleGroups.from_dataframe."""
        df = pd.DataFrame(
            {
                "scaling_group": ["g1", "g1", "g2"],
            },
            index=["s1", "s2", "s3"],
        )

        groups = SampleGroups.from_dataframe(df, subset_column="scaling_group")
        assert len(groups) == 2
        names = {g.name for g in groups.groups}
        assert names == {"g1", "g2"}
        assert set(groups.get_samples("g1")) == {"s1", "s2"}

    def test_samplegroupings_container(self):
        """Test SampleGroupings container functionality."""
        g1 = SampleGroups(groups=[SampleGroup(name="g", samples=["a", "b"])])
        sg = SampleGroupings()
        sg.add_grouping("consensus", g1)
        assert "consensus" in sg
        assert sg.get_grouping("consensus").groups[0].name == "g"


# =============================================================================
# HELPERS MODULE TESTS
# =============================================================================

class TestHelpers:
    """Tests for helper functions."""

    def test_get_sample_collection_fastq(self, tmp_path):
        """Test get_sample_collection with FASTQ metadata."""
        csv_path = tmp_path / "metadata.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "r1", "r2"])
            w.writerow([
                "s1",
                str(_write_fastq(tmp_path, "s1_R1.fastq.gz")),
                str(_write_fastq(tmp_path, "s1_R2.fastq.gz")),
            ])

        sc = get_sample_collection(Assay.RNA, csv_path)
        assert isinstance(sc, CollectionLike)
        paths = sc.get_file_paths("fastq")
        assert len(paths) == 2

    def test_get_sample_collection_with_ip(self, tmp_path):
        """Test get_sample_collection with IP metadata (FastqCollectionForIP)."""
        csv_path = tmp_path / "metadata.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "ip", "r1"])
            w.writerow(["s1", "H3K4me3", str(_write_fastq(tmp_path, "s1_R1.fastq.gz"))])

        sc = get_sample_collection(Assay.CHIP, csv_path)
        assert isinstance(sc, CollectionLike)
        assert sc.assay == Assay.CHIP

    def test_get_sample_collection_missing_file(self, tmp_path):
        """Test get_sample_collection with missing file."""
        with pytest.raises(FileNotFoundError):
            get_sample_collection(Assay.RNA, tmp_path / "missing.csv")

    def test_get_sample_collection_unknown_format(self, tmp_path):
        """Test get_sample_collection with unknown format."""
        csv_path = tmp_path / "metadata.csv"
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "unknown_column"])
            w.writerow(["s1", "value"])

        with pytest.raises(RuntimeError, match="Could not determine the type of collection"):
            get_sample_collection(Assay.RNA, csv_path)

    def test_get_sample_collection_bam(self, tmp_path):
        """Test get_sample_collection with BAM metadata."""
        csv_path = tmp_path / "metadata.csv"
        bam_file = tmp_path / "sample.bam"
        bam_file.touch()
        
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "bam"])
            w.writerow(["s1", str(bam_file)])

        sc = get_sample_collection(Assay.ATAC, csv_path)
        assert isinstance(sc, BamCollection)
        assert len(sc.bam_files) == 1

    def test_get_sample_collection_bigwig(self, tmp_path):
        """Test get_sample_collection with BigWig metadata."""
        csv_path = tmp_path / "metadata.csv"
        bw_file = tmp_path / "sample.bigWig"
        bw_file.touch()
        
        with csv_path.open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "bigwig"])
            w.writerow(["s1", str(bw_file)])

        sc = get_sample_collection(Assay.ATAC, csv_path)
        assert isinstance(sc, BigWigCollection)
        assert len(sc.bigwig_files) == 1


# =============================================================================
# INTERFACES MODULE TESTS
# =============================================================================

class TestInterfaces:
    """Tests for interface protocols."""

    def test_ensure_collection_valid(self, tmp_path):
        """Test ensure_collection with valid collection."""
        bam = tmp_path / "sample.bam"
        bam.touch()
        collection = BamCollection.from_files(assay=Assay.ATAC, files=[bam])
        
        result = ensure_collection(collection)
        assert result == collection

    def test_ensure_collection_invalid(self):
        """Test ensure_collection with invalid object."""
        with pytest.raises(TypeError, match="does not implement required collection interface"):
            ensure_collection("not a collection")

    def test_detect_file_type_fastq(self, tmp_path):
        """Test detect_file_type with FASTQ files."""
        files = [
            tmp_path / "s1.fastq.gz",
            tmp_path / "s2.fastq.gz",
            tmp_path / "s3.fq.gz",
        ]
        for f in files:
            f.touch()
        
        file_type = detect_file_type(files)
        assert file_type == "fastq"

    def test_detect_file_type_bam(self, tmp_path):
        """Test detect_file_type with BAM files."""
        files = [tmp_path / "s1.bam", tmp_path / "s2.bam"]
        for f in files:
            f.touch()
        
        file_type = detect_file_type(files)
        assert file_type == "bam"

    def test_detect_file_type_bigwig(self, tmp_path):
        """Test detect_file_type with BigWig files."""
        files = [tmp_path / "s1.bigWig", tmp_path / "s2.bw"]
        for f in files:
            f.touch()
        
        file_type = detect_file_type(files)
        assert file_type == "bigwig"

    def test_detect_file_type_mixed(self, tmp_path):
        """Test detect_file_type with mixed files (no clear winner)."""
        files = [tmp_path / "s1.fastq.gz", tmp_path / "s2.bam"]
        for f in files:
            f.touch()
        
        file_type = detect_file_type(files)
        # With a tie, should return None
        assert file_type is None

    def test_detect_file_type_unknown(self, tmp_path):
        """Test detect_file_type with unknown files."""
        files = [tmp_path / "s1.txt", tmp_path / "s2.pdf"]
        for f in files:
            f.touch()
        
        file_type = detect_file_type(files)
        assert file_type is None


# =============================================================================
# VALIDATION MODULE TESTS
# =============================================================================

class TestViewpointsFile:
    """Tests for ViewpointsFile validation schema."""

    def test_valid_viewpoints_file(self):
        """Test valid viewpoints file validation."""
        data = {
            "Chromosome": ["chr1", "chr2", "chr3"],
            "Start": [100, 200, 300],
            "End": [200, 300, 400],
            "Name": ["viewpoint1", "viewpoint2", "viewpoint3"],
            "Strand": ["+", "-", "+"],
            "Score": [10.0, 20.0, 30.0],
        }
        df = pd.DataFrame(data)

        validated = ViewpointsFile.validate(df)
        assert validated is not None
        assert len(validated) == 3

    def test_viewpoints_without_optional_columns(self):
        """Test viewpoints file without optional columns."""
        data = {
            "Chromosome": ["chr1", "chr2"],
            "Start": [100, 200],
            "End": [200, 300],
            "Name": ["viewpoint1", "viewpoint2"],
        }
        df = pd.DataFrame(data)

        validated = ViewpointsFile.validate(df)
        assert validated is not None

    def test_viewpoint_names_invalid(self):
        """Test viewpoint names with spaces are invalid."""
        data = {
            "Chromosome": ["chr1"],
            "Start": [100],
            "End": [200],
            "Name": ["my viewpoint"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_viewpoint_names"):
            ViewpointsFile.validate(df)


class TestDesignDataFrame:
    """Tests for DesignDataFrame validation schema."""

    def test_valid_design_dataframe(self):
        """Test valid design dataframe validation."""
        data = {
            "sample_id": ["sample1", "sample2", "sample3"],
            "r1": ["read1.fq", "read2.fq", "read3.fq"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None
        assert len(validated) == 3

    def test_design_with_optional_columns(self):
        """Test design dataframe with all optional columns."""
        data = {
            "sample_id": ["sample1", "sample2"],
            "scaling_group": ["batch1", "batch1"],
            "consensus_group": ["group1", "group2"],
            "deseq2": [0, 1],
            "r1": ["read1.fq", "read2.fq"],
        }
        df = pd.DataFrame(data)

        validated = DesignDataFrame.validate(df)
        assert validated is not None

    def test_unique_sample_ids(self):
        """Test sample IDs must be unique."""
        data = {
            "sample_id": ["sample1", "sample1"],
            "r1": ["read1.fq", "read2.fq"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_sample_name_is_unique"):
            DesignDataFrame.validate(df)

    def test_sample_id_with_spaces_invalid(self):
        """Test sample IDs with spaces are invalid."""
        data = {
            "sample_id": ["sample 1"],
            "r1": ["read1.fq"],
        }
        df = pd.DataFrame(data)

        with pytest.raises(SchemaError, match="check_sample_id"):
            DesignDataFrame.validate(df)
