import shlex
from typing import Optional, Literal, Annotated
from enum import Enum
from pathlib import Path

from pydantic import BaseModel, Field, field_validator, model_validator, BeforeValidator
from pydantic.functional_serializers import PlainSerializer

from seqnado import Assay
from typing import get_origin, get_args, Union


def none_str_to_none(v):
    """Convert string 'None' to actual None value."""
    if isinstance(v, str) and v.strip().lower() == "none":
        return None
    return v


# =============================================================================
# Core Configuration Classes
# =============================================================================

class CommandLineArguments(BaseModel):
    """Base class for CLI options with validation and filtering capabilities."""
    
    value: str = Field(default="", description="CLI options string")
    exclude: set[str] | None = Field(
        None, description="Options to exclude from the final command"
    )
    include: set[str] | None = Field(
        None, description="Options to include in the final command"
    )

    def model_post_init(self, context):
        # Initialize sets if None
        if self.exclude is None:
            self.exclude = set()
        if self.include is None:
            self.include = set()


    @field_validator("value", mode="before")
    @classmethod
    def validate_options_syntax(cls, v: str) -> str:
        if not isinstance(v, str):
            raise ValueError("Value must be a string")
        try:
            shlex.split(v)
        except ValueError as e:
            raise ValueError(f"Invalid CLI option string: {e}")
        if any(c in v for c in ["\n", "\r", "\x00"]):
            raise ValueError("Value string contains unsafe characters")
        return v

    @model_validator(mode="after")
    def ensure_leading_space(self) -> "CommandLineArguments":
        if self.value and not self.value.startswith(" "):
            self.value = " " + self.value
        return self

    def _filter_excluded_options(self, options_str: str) -> str:
        """
        Filter an options string by applying include (allow-list) and exclude (deny-list)

        Behavior:
        - If `include` is non-empty, only tokens that match an include pattern are kept.
          A match is either exact equality or startswith pattern + '='. If a token is
          an option (e.g. '--opt') that consumes a following value (non-dashed token),
          the value is kept alongside the option.
                - After include filtering (if any), exclude patterns are applied and matching
                    tokens are removed. If the excluded option is provided with a separate value
                    (e.g. "--opt value"), the very next token is also removed. The "--opt=value"
                    form is matched via startswith("--opt=") and removed as a single token.
        - Exclude patterns override include: if the same pattern is present in both,
          the token will be removed.
        """
        if not (self.include or self.exclude):
            return options_str
        try:
            tokens = shlex.split(options_str.strip())

            def matches(tok: str, pattern: str) -> bool:
                return tok == pattern or tok.startswith(pattern + "=")

            # Apply excludes first to the token list
            if self.exclude:
                filtered = []
                i = 0
                while i < len(tokens):
                    tok = tokens[i]
                    excluded = False
                    for pat in self.exclude:
                        if matches(tok, pat):
                            excluded = True
                            # Drop the following token if present to handle separate value form
                            if tok == pat and i + 1 < len(tokens):
                                i += 1
                            break
                    if not excluded:
                        filtered.append(tok)
                    i += 1
                tokens = filtered

            # Ensure include patterns are present; append missing ones.
            if self.include:
                def present(pat: str) -> bool:
                    for t in tokens:
                        if matches(t, pat):
                            return True
                    return False

                for pat in self.include:
                    if not present(pat):
                        # append the pattern as-is; if user intended a value they should
                        # provide it via the pattern (e.g. '--opt=value')
                        tokens.append(pat)

            return " " + shlex.join(tokens) if tokens else ""
        except ValueError:
            return options_str

    # --- Dynamic update helpers for include / exclude sets -----------------
    def add_include(self, *patterns: str) -> "CommandLineArguments":
        """Add one or more include patterns dynamically and return self for chaining."""
        for p in patterns:
            if isinstance(p, str) and p:
                self.include.add(p)
        return self

    def remove_include(self, *patterns: str) -> "CommandLineArguments":
        """Remove one or more include patterns dynamically and return self for chaining."""
        for p in patterns:
            self.include.discard(p)
        return self

    def clear_include(self) -> "CommandLineArguments":
        """Clear all include patterns."""
        self.include.clear()
        return self

    def set_include(self, patterns: set[str]) -> "CommandLineArguments":
        """Replace the include set with `patterns` (must be an iterable of strings)."""
        self.include = set(patterns) if patterns is not None else set()
        return self

    def add_exclude(self, *patterns: str) -> "CommandLineArguments":
        """Add one or more exclude patterns dynamically and return self for chaining."""
        for p in patterns:
            if isinstance(p, str) and p:
                self.exclude.add(p)
        return self

    def remove_exclude(self, *patterns: str) -> "CommandLineArguments":
        """Remove one or more exclude patterns dynamically and return self for chaining."""
        for p in patterns:
            self.exclude.discard(p)
        return self

    def clear_exclude(self) -> "CommandLineArguments":
        """Clear all exclude patterns."""
        self.exclude.clear()
        return self

    def set_exclude(self, patterns: set[str]) -> "CommandLineArguments":
        """Replace the exclude set with `patterns` (must be an iterable of strings)."""
        self.exclude = set(patterns) if patterns is not None else set()
        return self

    @property
    def option_string_filtered(self) -> str:
        return self._filter_excluded_options(self.value)

    @property
    def option_string_raw(self) -> str:
        return self.value
    
    def __str__(self) -> str:
        return self._filter_excluded_options(self.value).strip()


def serialize_options_to_string(options: CommandLineArguments) -> str:
    """Custom serializer that returns only the filtered option string."""
    return options.option_string_filtered.strip()


# Create the annotated Options type with custom serialization
CommandLineArgumentsType = Annotated[
    CommandLineArguments,
    PlainSerializer(serialize_options_to_string, return_type=str)
]


class ToolConfig(BaseModel):
    """Base configuration for a tool with threads and CLI options."""
    
    threads: int = Field(default=1, ge=1, description="Number of threads to use")
    command_line_arguments: CommandLineArgumentsType = Field(default_factory=CommandLineArguments, description="CLI options")

    @field_validator("command_line_arguments", mode="before")
    @classmethod
    def validate_command_line_arguments(cls, v):
        if isinstance(v, CommandLineArguments):
            return v
        if isinstance(v, str):
            return CommandLineArguments(value=v)
        if isinstance(v, dict):
            return CommandLineArguments(**v)
        raise TypeError("command_line_arguments must be a string, dict, or CommandLineArguments instance")


# =============================================================================
# Enums
# =============================================================================

class CutadaptMode(Enum):
    DEFAULT = "default"
    CRISPR = "crispr"


class MacsMode(Enum):
    ATAC = "atac"
    CHIP_BROAD = "chip_broad"
    CHIP_TF = "chip_tf"
    GENERIC = "generic"


# =============================================================================
# Simple Tool Configurations (single ToolConfig)
# =============================================================================

class Cutadapt(ToolConfig):
    """Cutadapt adapter trimming tool configuration.

    When mode='crispr', adapters are auto-detected from FASTQ files.
    Otherwise, specify adapters in command_line_arguments.
    """

    mode: CutadaptMode = Field(default=CutadaptMode.CRISPR, description="Cutadapt mode: 'crispr' for auto-detection, 'default' for manual")

    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 4)

        # Set default command_line_arguments based on mode
        mode = data.get("mode", CutadaptMode.CRISPR)
        if isinstance(mode, str):
            mode = CutadaptMode(mode.lower())

        if mode == CutadaptMode.CRISPR:
            # For CRISPR mode, don't include adapter specs - they'll be auto-detected
            # Only include parameters that affect trimming behavior
            data.setdefault(
                "command_line_arguments",
                CommandLineArguments(value="--cut 0 -l 20 -m 20 -M 20 --discard-untrimmed"),
            )
        else:
            # For default mode, use empty defaults (user should specify adapters)
            data.setdefault(
                "command_line_arguments",
                CommandLineArguments(value=""),
            )

        return data


class LanceotronMCC(ToolConfig):
    """Lanceotron MCC peak calling tool configuration."""
    
    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 8)
        data.setdefault(
            "command_line_arguments",
            CommandLineArguments(
                value="--max-peak-size 2000 --no-mcc-filter-enriched-regions --no-mcc-filter-background"
            ),
        )
        return data


class Methyldackel(ToolConfig):
    """Methyldackel methylation analysis tool configuration."""
    
    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 8)
        data.setdefault("command_line_arguments", CommandLineArguments(value=""))
        return data


# =============================================================================
# Complex Tool Configurations (multiple sub-tools)
# =============================================================================

class Bowtie2(BaseModel):
    """Bowtie2 alignment tool configuration."""
    
    align: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--very-sensitive")
        ),
        description="Alignment configuration"
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--threads {threads} --quiet")
        ),
        description="Index building configuration"
    )


class Samtools(BaseModel):
    """Samtools suite configuration."""
    
    sort: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="-@ {threads} -m 2G")
        ),
        description="Sort configuration"
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="")),
        description="Index configuration"
    )
    view: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="-bS")),
        description="View configuration"
    )
    merge: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="")),
        description="Merge configuration"
    )


class Deeptools(BaseModel):
    """Deeptools suite configuration."""
    
    bam_coverage: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(
                value="--binSize 10 --normalizeUsing RPKM", 
                exclude={"-e", "--extendReads"}
            )
        ),
        description="BAM coverage configuration"
    )
    compute_matrix: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="")),
        description="Compute matrix configuration"
    )
    plot_heatmap: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            command_line_arguments=CommandLineArguments(
                value="--colorMap RdYlBu --whatToShow 'heatmap and colorbar' --boxAroundHeatmaps no"
            ),
        ),
        description="Plot heatmap configuration"
    )



class Homer(BaseModel):
    """HOMER motif analysis suite configuration."""

    make_tag_directory: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Tag directory creation configuration"
    )
    make_bigwig: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="BigWig creation configuration"
    )
    annotate_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Peak annotation configuration"
    )
    find_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Peak finding configuration"
    )
    find_motifs_genome: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="Motif finding configuration"
    )
    known_motif_db: Annotated[Optional[str], BeforeValidator(none_str_to_none)] = Field(
        default=None,
        description="Path to known motif database file"
    )


class Meme(BaseModel):
    """MEME Suite motif analysis configuration."""

    meme_chip: ToolConfig = Field(
        default_factory=lambda: ToolConfig(command_line_arguments=CommandLineArguments(value="")),
        description="MEME-ChIP motif analysis configuration"
    )
    known_motif_db: Annotated[Optional[str], BeforeValidator(none_str_to_none)] = Field(
        default=None,
        description="Path to known motif database file"
    )


class Star(BaseModel):
    """STAR RNA-seq aligner configuration."""
    
    align: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            command_line_arguments=CommandLineArguments(
                value="--quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --outSAMattributes Standard"
            ),
        ),
        description="Alignment configuration"
    )
    index: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8,
            command_line_arguments=CommandLineArguments(value="--runThreadN {threads} --genomeSAindexNbases 14"),
        ),
        description="Index building configuration"
    )


class BcfTools(BaseModel):
    """BCFtools variant calling suite configuration."""
    
    call: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="-c -v -f GQ,DP")
        ),
        description="Variant calling configuration"
    )
    consensus: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="-s -c")),
        description="Consensus calling configuration"
    )

class Salmon(BaseModel):
    """Salmon quantification tool configuration."""

    quant: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="")
        ),
        description="Quantification configuration"
    )


# =============================================================================
# Single-purpose Tool Configurations
# =============================================================================

class Bamnado(BaseModel):
    """Bamnado coverage analysis tool configuration."""
    
    bam_coverage: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--bin-size 10 --norm-method raw")
        ),
        description="BAM coverage analysis configuration"
    )
    bigwig_compare: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--bin-size 10")
        ),
        description="BigWig comparison configuration"
    )
    bigwig_aggregate: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="--bin-size 10")
        ),
        description="BigWig aggregation configuration"
    )


class Subread(BaseModel):
    """Subread feature counting tool configuration."""
    
    feature_counts: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=16, 
            command_line_arguments=CommandLineArguments(value="")
        ),
        description="Feature counting configuration"
    )


class MageckCount(ToolConfig):
    """MAGeCK count configuration."""

    norm_method: Literal["none", "median", "total", "control"] = Field(
        default="median",
        description="Normalization method for read counts"
    )
    pdf_report: bool = Field(
        default=False,
        description="Generate PDF quality control report"
    )

    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 8)
        data.setdefault("command_line_arguments", CommandLineArguments(value=""))
        return data


class MageckTest(ToolConfig):
    """MAGeCK test (RRA) configuration."""

    control_samples: Optional[str] = Field(
        default=None,
        description="Comma-separated list of control sample labels"
    )
    treatment_samples: Optional[str] = Field(
        default=None,
        description="Comma-separated list of treatment sample labels"
    )

    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 4)
        data.setdefault("command_line_arguments", CommandLineArguments(value=""))
        return data


class MageckMLE(ToolConfig):
    """MAGeCK MLE configuration."""

    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 8)
        data.setdefault("command_line_arguments", CommandLineArguments(value=""))
        return data


class Mageck(BaseModel):
    """MAGeCK CRISPR screen analysis tool configuration."""

    count: MageckCount = Field(
        default_factory=MageckCount,
        description="MAGeCK count configuration for guide RNA quantification"
    )
    test: Optional[MageckTest] = Field(
        default=None,
        description="MAGeCK test (RRA) configuration for comparing treatment vs control"
    )
    mle: MageckMLE = Field(
        default_factory=MageckMLE,
        description="MAGeCK MLE configuration for complex experimental designs"
    )
    sgrna_library: Optional[str] = Field(
        default=None,
        description="Path to sgRNA library file (.txt or .csv format)"
    )
    design_matrix: Optional[str] = Field(
        default=None,
        description="Path to design matrix file for MAGeCK MLE analysis"
    )
    control_sgrna: Optional[str] = Field(
        default=None,
        description="Comma-separated list of control sgRNA IDs or file path (prefix with --control-sgrna)"
    )
    control_gene: Optional[str] = Field(
        default=None,
        description="Comma-separated list of control gene names or file path (prefix with --control-gene)"
    )


class Lanceotron(BaseModel):
    """Lanceotron peak calling tool configuration."""
    
    call_peaks: ToolConfig = Field(
        default_factory=lambda: ToolConfig(threads=8, command_line_arguments=CommandLineArguments(value="-c 0.5")),
        description="Peak calling configuration"
    )


class Picard(BaseModel):
    """Picard tools configuration."""
    
    mark_duplicates: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=8, command_line_arguments=CommandLineArguments(value="REMOVE_DUPLICATES=true")
        ),
        description="Mark duplicates configuration"
    )


class Trimgalore(BaseModel):
    """Trim Galore adapter trimming tool configuration."""
    
    trim: ToolConfig = Field(
        default_factory=lambda: ToolConfig(
            threads=4, command_line_arguments=CommandLineArguments(value="--2colour 20")
        ),
        description="Trimming configuration"
    )


# =============================================================================
# Special Configuration Classes
# =============================================================================

class Macs(BaseModel):
    """MACS peak calling tool configuration with mode-specific defaults."""
    
    version: Literal["2", "3"] = Field(default="3", description="MACS version to use")
    mode: MacsMode = Field(default=MacsMode.GENERIC, description="Analysis mode")
    call_peaks: Optional[ToolConfig] = Field(default=None, description="Peak calling configuration")

    def _get_options_for_mode(self) -> str:
        """Get default options based on the analysis mode."""
        match self.mode:
            case MacsMode.ATAC:
                return "--nomodel --shift -100 --extsize 200"
            case MacsMode.CHIP_BROAD:
                return "--broad --nomodel --shift -100 --extsize 200"
            case MacsMode.CHIP_TF:
                return "--nomodel --shift -100 --extsize 200"
            case _:
                return "-f BAMPE"

    @model_validator(mode="after")
    def set_default_call_peak(self) -> "Macs":
        if self.call_peaks is None:
            self.call_peaks = ToolConfig(
                command_line_arguments=CommandLineArguments(value=self._get_options_for_mode())
            )
        return self


class Seacr(BaseModel):
    """SEACR peak calling tool configuration (non-ToolConfig based)."""
    
    threshold: float = Field(default=0.01, ge=0, le=1, description="Peak calling threshold")
    normalization: Optional[str] = Field(default="non", description="Normalization method")
    stringency: Literal["stringent", "relaxed"] = Field(
        default="stringent", description="Peak calling stringency"
    )


class FastqScreen(ToolConfig):
    """FastqScreen tool configuration."""

    config: str = Field(default="", description="Path to FastqScreen configuration file")

    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 4)
        data.setdefault("config", str(Path.home() / ".config/seqnado/fastq_screen.conf"))
        data.setdefault("command_line_arguments", CommandLineArguments(value="--subset 10000"))
        return data
    
    @field_validator("config", mode="before")
    @classmethod
    def validate_config_path(cls, v):
        if v is None or v == "None":
            return ""
        if isinstance(v, Path):
            return str(v)
        return str(v) if v else ""

class FastQC(ToolConfig):
    """FastQC quality control tool configuration."""
    
    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 4)
        data.setdefault("command_line_arguments", CommandLineArguments(value="--noextract"))
        return data

class Qualimap(ToolConfig):
    """Qualimap quality control tool configuration."""
    
    @model_validator(mode="before")
    @classmethod
    def set_defaults(cls, data: dict) -> dict:
        data = dict(data)
        data.setdefault("threads", 4)
        data.setdefault("command_line_arguments", CommandLineArguments(value="--java-mem-size=4G"))
        return data


# =============================================================================
# Main Configuration Class
# =============================================================================

def get_assay_specific_tools(assay: Assay) -> list[type[BaseModel]]:
    """Get the list of tool classes appropriate for a given assay type."""

    qc_tools = [FastqScreen, FastQC, Qualimap]

    generic_dna_tools = [
        Bowtie2,
        Bamnado,
        Deeptools,
        Homer,
        Meme,
        Lanceotron,
        Macs,
        Picard,
        Samtools,
        Trimgalore,
        Subread,
    ]

    tools = {
        Assay.ATAC: qc_tools + generic_dna_tools,
        Assay.CHIP: qc_tools + generic_dna_tools,
        Assay.CAT: qc_tools + generic_dna_tools + [Seacr],
        Assay.RNA: qc_tools + [Star, Deeptools, Samtools, Trimgalore, Salmon, Subread, Homer, Bamnado],
        Assay.SNP: qc_tools + [Bowtie2, Trimgalore, Samtools, BcfTools],
        Assay.METH: qc_tools + [Bowtie2, Trimgalore, Methyldackel, Samtools, Picard],
        Assay.CRISPR: qc_tools + [Qualimap, Cutadapt, Bowtie2, Subread, Samtools, Mageck],
        Assay.MCC: qc_tools + [Trimgalore, Bowtie2, Samtools, Deeptools, LanceotronMCC, Bamnado],
    }
    return tools.get(assay, [])


class ThirdPartyToolsConfig(BaseModel):
    """Configuration for all third-party bioinformatics tools."""
    
    # Quality control tools
    fastq_screen: Annotated[Optional[FastqScreen], BeforeValidator(none_str_to_none)] = Field(default=None, description="FastqScreen configuration")
    fastqc: Annotated[Optional[FastQC], BeforeValidator(none_str_to_none)] = Field(
        default=None, description="FastQC quality control configuration"
    )
    qualimap: Annotated[Optional[Qualimap], BeforeValidator(none_str_to_none)] = Field(
        default=None, description="Qualimap quality control configuration"
    )

    # Alignment tools
    bowtie2: Annotated[Optional[Bowtie2], BeforeValidator(none_str_to_none)] = Field(default=None, description="Bowtie2 aligner configuration")
    star: Annotated[Optional[Star], BeforeValidator(none_str_to_none)] = Field(default=None, description="STAR RNA-seq aligner configuration")

    # Processing tools
    samtools: Annotated[Optional[Samtools], BeforeValidator(none_str_to_none)] = Field(default=None, description="Samtools suite configuration")
    picard: Annotated[Optional[Picard], BeforeValidator(none_str_to_none)] = Field(default=None, description="Picard tools configuration")
    
    # Trimming tools
    cutadapt: Annotated[Optional[Cutadapt], BeforeValidator(none_str_to_none)] = Field(default=None, description="Cutadapt trimming configuration")
    trim_galore: Annotated[Optional[Trimgalore], BeforeValidator(none_str_to_none)] = Field(default=None, description="Trim Galore configuration")
    
    # Peak calling tools
    macs: Annotated[Optional[Macs], BeforeValidator(none_str_to_none)] = Field(default=None, description="MACS peak caller configuration")
    lanceotron: Annotated[Optional[Lanceotron], BeforeValidator(none_str_to_none)] = Field(default=None, description="Lanceotron peak caller configuration")
    lanceotron_mcc: Annotated[Optional[LanceotronMCC], BeforeValidator(none_str_to_none)] = Field(default=None, description="Lanceotron MCC configuration")
    seacr: Annotated[Optional[Seacr], BeforeValidator(none_str_to_none)] = Field(default=None, description="SEACR peak caller configuration")
    
    # Analysis tools
    deeptools: Annotated[Optional[Deeptools], BeforeValidator(none_str_to_none)] = Field(default=None, description="Deeptools suite configuration")
    homer: Annotated[Optional[Homer], BeforeValidator(none_str_to_none)] = Field(default=None, description="HOMER suite configuration")
    meme: Annotated[Optional[Meme], BeforeValidator(none_str_to_none)] = Field(default=None, description="MEME Suite motif analysis configuration")
    bamnado: Annotated[Optional[Bamnado], BeforeValidator(none_str_to_none)] = Field(default=None, description="Bamnado coverage analysis configuration")
    subread: Annotated[Optional[Subread], BeforeValidator(none_str_to_none)] = Field(default=None, description="Subread feature counting configuration")
    salmon: Annotated[Optional[Salmon], BeforeValidator(none_str_to_none)] = Field(default=None, description="Salmon quantification configuration")
    methyldackel: Annotated[Optional[Methyldackel], BeforeValidator(none_str_to_none)] = Field(default=None, description="Methyldackel methylation analysis configuration")
    bcftools: Annotated[Optional[BcfTools], BeforeValidator(none_str_to_none)] = Field(default=None, description="BCFtools variant calling suite configuration")
    mageck: Annotated[Optional[Mageck], BeforeValidator(none_str_to_none)] = Field(default=None, description="MAGeCK CRISPR screen analysis configuration")
    
    @classmethod
    def _class_to_field_map(cls) -> dict[type[BaseModel], str]:
        """Map tool classes to the corresponding field name on this model."""
        mapping: dict[type[BaseModel], str] = {}
        for fname, f in cls.model_fields.items():
            ann = f.annotation
            origin = get_origin(ann)
            if origin is Union:
                # e.g. Optional[FastqScreen] -> Union[FastqScreen, NoneType]
                for a in get_args(ann):
                    if isinstance(a, type) and issubclass(a, BaseModel):
                        mapping[a] = fname
            elif isinstance(ann, type) and issubclass(ann, BaseModel):
                mapping[ann] = fname
        return mapping

    @classmethod
    def for_assay(cls, assay: Assay, **overrides) -> "ThirdPartyToolsConfig":
        assay_tools = get_assay_specific_tools(assay)
        if not assay_tools:
            raise ValueError(f"No tools configured for assay {assay.value}")

        class_to_field = cls._class_to_field_map()

        defaults: dict[str, BaseModel] = {}
        for tool_class in assay_tools:
            # exact match on the declared field’s expected type
            field_name = class_to_field.get(tool_class)
            if field_name is None:
                # Optional: try subclass match if your tool classes have inheritance
                for k, v in class_to_field.items():
                    if issubclass(tool_class, k):
                        field_name = v
                        break
            if field_name is None:
                raise ValueError(
                    f"No field on {cls.__name__} is annotated for tool class {tool_class.__name__}"
                )
            defaults[field_name] = tool_class()
            
        
        defaults.update(overrides)
        return cls(**defaults)

    def get_configured_tools(self) -> dict[str, BaseModel]:
        """Get a dictionary of all configured (non-None) tools."""
        return {
            name: tool for name, tool in self.model_dump(exclude_none=True).items()
            if tool is not None
        }

    def has_tool(self, tool_name: str) -> bool:
        """Check if a specific tool is configured."""
        return getattr(self, tool_name, None) is not None