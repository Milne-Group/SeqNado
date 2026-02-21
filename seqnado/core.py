from enum import Enum

from pydantic import BaseModel, field_validator

# =============================================================================
# CONSTANTS
# =============================================================================
NONE_VALUES = [None, "None", "none", "null", "Null", "NULL", ".", "", "NA"]

# ============================================================
# ENUMS
# ============================================================


class Assay(Enum):
    """Supported sequencing assay types."""

    RNA = "RNA"
    ATAC = "ATAC"
    SNP = "SNP"
    CHIP = "ChIP"
    CAT = "CUT&TAG"
    METH = "METH"
    MCC = "MCC"
    CRISPR = "CRISPR"
    MULTIOMICS = "MULTIOMICS"

    @classmethod
    def all_assays(cls):
        """Return all supported assays."""
        return list(cls)

    @property
    def clean_name(self):
        """Return a short name for the assay."""
        short_names = {
            self.RNA: "rna",
            self.ATAC: "atac",
            self.SNP: "snp",
            self.CHIP: "chip",
            self.CAT: "cat",
            self.METH: "meth",
            self.MCC: "mcc",
            self.CRISPR: "crispr",
            self.MULTIOMICS: "multiomics",
        }

        if self in short_names:
            return short_names[self]
        else:
            raise ValueError(f"Unknown assay type: {self}")

    @classmethod
    def from_clean_name(cls, clean_name):
        """Return the assay type from a short name."""
        for assay in cls:
            if assay.clean_name == clean_name:
                return assay
        raise ValueError(f"Unknown clean name: {clean_name}")

    @classmethod
    def all_assay_clean_names(cls):
        """Return a list of all clean names for assays."""
        return [assay.clean_name for assay in cls]

    @classmethod
    def non_multiomics_assays(cls):
        """Return assays that are not multiomics."""
        return [assay for assay in cls if assay != cls.MULTIOMICS]

    @classmethod
    def non_ip_assays(cls):
        """Return assays that don't require IP (immunoprecipitation)."""
        ip_assays = {cls.CHIP, cls.CAT}
        return [assay for assay in cls if assay not in ip_assays]

    @classmethod
    def ip_assays(cls):
        """Return assays that require IP (immunoprecipitation)."""
        return [cls.CHIP, cls.CAT]


AssaysWithPeakCalling = (Assay.ATAC, Assay.CHIP, Assay.CAT, Assay.MCC)
AssaysWithHeatmaps = (Assay.ATAC, Assay.CHIP, Assay.RNA, Assay.CAT)
AssaysWithSpikein = (Assay.ATAC, Assay.CHIP, Assay.CAT, Assay.RNA)


class FileType(Enum):
    """Supported file types."""

    BAM = "bam"
    BED = "bed"
    BIGWIG = "bigwig"
    FASTQ = "fastq"
    FASTA = "fasta"
    TXT = "txt"
    TSV = "tsv"
    XLSX = "xlsx"
    TAG_DIRECTORY = "homer_tag_directory"


class PileupMethod(Enum):
    """Methods for creating pileup files."""

    DEEPTOOLS = "deeptools"
    HOMER = "homer"
    BAMNADO = "bamnado"


class DataScalingTechnique(Enum):
    """Methods for scaling genomic data."""

    UNSCALED = "unscaled"
    CSAW = "csaw"
    CPM = "cpm"
    RPKM = "rpkm"
    SPIKEIN = "spikein"
    MERGED = "merged"


class PeakCallingMethod(Enum):
    """Methods for calling peaks."""

    MACS2 = "macs2"
    MACS3 = "macs3"
    HOMER = "homer"
    LANCEOTRON = "lanceotron"
    SEACR = "seacr"
    LANCEOTRON_MCC = "lanceotron-mcc"


class MotifMethod(Enum):
    """Methods for motif analysis."""

    HOMER = "homer"
    MEME = "meme"


class PCRDuplicateHandling(Enum):
    """Methods for handling PCR duplicates."""

    REMOVE = "remove"
    MARK = "mark"
    NONE = "keep"


class PCRDuplicateTool(Enum):
    PICARD = "picard"
    SAMTOOLS = "samtools"
    NONE = "None"


class SpikeInMethod(Enum):
    """Methods for spike-in normalization."""

    ORLANDO = "orlando"
    WITH_INPUT = "with_input"
    DESEQ2 = "deseq2"
    EDGER = "edger"


class SNPCallingMethod(Enum):
    """Methods for SNP calling."""

    BCFTOOLS = "bcftools"
    DEEPVARIANT = "deepvariant"


class QuantificationMethod(Enum):
    """Methods for quantification."""

    FEATURE_COUNTS = "feature_counts"
    SALMON = "salmon"


class MethylationMethod(Enum):
    """Methods for methylation calling."""

    TAPS = "taps"
    BISULFITE = "bisulfite"


class Molecule(Enum):
    rna_total = "total RNA"
    rna_polya = "polyA RNA"
    rna_cytoplasmic = "cytoplasmic RNA"
    rna_nuclear = "nuclear RNA"
    dna_genomic = "genomic DNA"
    protein = "protein"
    other = "other"


class Organism(Enum):
    """Supported organisms."""

    HUMAN = "Homo sapiens"
    MOUSE = "Mus musculus"
    RAT = "Rattus norvegicus"
    ZEBRAFISH = "Danio rerio"
    DROSOPHILA = "Drosophila melanogaster"
    C_ELEGANS = "Caenorhabditis elegans"
    UNKNOWN = "Unknown"


class LibraryType(Enum):
    """Supported library types."""

    SINGLE = "single-end"
    PAIRED = "paired-end"


class GenomicCoordinate(BaseModel):
    """Configuration for genomic coordinates."""

    chromosome: str
    start: int
    end: int

    @field_validator("start", "end")
    def validate_coordinates(cls, v: int) -> int:
        if v < 0:
            raise ValueError("Genomic coordinates must be non-negative.")
        return v

    # Check that end is greater than start
    @field_validator("end")
    def validate_end_greater_than_start(cls, v: int, info) -> int:
        # Only validate if start is available (i.e., it passed validation)
        if "start" in info.data and v < info.data["start"]:
            raise ValueError("End coordinate must be greater than start coordinate.")
        return v

    @classmethod
    def from_string(cls, coord_str: str) -> "GenomicCoordinate":
        """
        Create a GenomicCoordinate instance from a string representation.
        """
        chromosome, positions = coord_str.split(":")
        start, end = map(int, positions.split("-"))
        return cls(chromosome=chromosome, start=start, end=end)
