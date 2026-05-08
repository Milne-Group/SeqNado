from enum import Enum
from pathlib import Path
from typing import Annotated, Union 

from pydantic import (
    BaseModel,
    BeforeValidator,
    Field,
    computed_field,
    field_validator,
    model_validator,
)

from seqnado import Assay, PeakCallingMethod, SpikeInMethod

from .configs import (
    BigwigConfig,
    GenomeConfig,
    MCCConfig,
    MethylationConfig,
    PCRDuplicatesConfig,
    PeakCallingConfig,
    PlottingConfig,
    ProjectConfig,
    QCConfig,
    RNAQuantificationConfig,
    SNPCallingConfig,
    SpikeInConfig,
    UCSCHubConfig,
    none_str_to_none,
)
from .mixins import (
    CommonComputedFieldsMixin,
    MethylationMixin,
    PeakCallingMixin,
    SNPCallingMixin,
)
from .third_party_tools import Seacr, ThirdPartyToolsConfig


class BaseAssayConfig(CommonComputedFieldsMixin):
    """Base configuration for all assays."""

    bigwigs: BigwigConfig | None = None
    plotting: PlottingConfig | None = None
    ucsc_hub: UCSCHubConfig | None = None

    # Boolean flags for optional features
    create_geo_submission_files: bool = False

    # Lifted from ChIPAssayConfig / CATAssayConfig / RNAAssayConfig
    spikein: Annotated[SpikeInConfig | None, BeforeValidator(none_str_to_none)] = None


class ATACAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to ATAC-seq assays."""

    tn5_shift: bool = False
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class ChIPAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to ChIP-seq assays."""

    spikein: Annotated[SpikeInConfig | None, BeforeValidator(none_str_to_none)] = None
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class CATAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to CAT-seq assays."""

    spikein: Annotated[SpikeInConfig | None, BeforeValidator(none_str_to_none)] = None
    tn5_shift: bool = False
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False


class RNAAssayConfig(BaseAssayConfig):
    """Configuration specific to RNA-seq assays."""

    spikein: Annotated[SpikeInConfig | None, BeforeValidator(none_str_to_none)] = None
    rna_quantification: RNAQuantificationConfig | None = None
    create_heatmaps: bool = False

    @field_validator("spikein")
    @classmethod
    def validate_spikein_methods(cls, v):
        """Filter out incompatible spike-in methods for RNA-seq.
        
        WITH_INPUT requires paired input/control samples, which is a ChIP-seq concept.
        RNA-seq supports only ORLANDO, DESEQ2, and EDGER.
        """
        if v is not None and SpikeInMethod.WITH_INPUT in v.method:
            from seqnado.utils import warn_once
            warn_once(
                "The 'with_input' spike-in method is not compatible with RNA-seq and will be skipped. "
                "RNA-seq supports: 'orlando', 'deseq2', 'edger'. "
                "WITH_INPUT requires paired input/control samples (ChIP-seq concept)."
            )
            # Filter out the incompatible method
            v.method = [m for m in v.method if m != SpikeInMethod.WITH_INPUT]
        return v


class SNPAssayConfig(BaseAssayConfig, SNPCallingMixin):
    """Configuration specific to SNP calling assays."""

    snp_calling: SNPCallingConfig | None = None
    snp_database: str | None = None
    create_heatmaps: bool = False


class MCCAssayConfig(BaseAssayConfig, PeakCallingMixin):
    """Configuration specific to MCC (Capture-C) assays."""

    mcc: MCCConfig | None = None
    peak_calling: PeakCallingConfig | None = None
    create_heatmaps: bool = False

class MethylationAssayConfig(BaseAssayConfig, MethylationMixin):
    """Configuration specific to methylation assays."""

    # METH uses methylation-specific bigwig generation (TAPS/WGBS via methyldackel)
    # not standard pileup methods, so bigwigs config is not applicable
    bigwigs: None = None
    methylation: MethylationConfig | None = None
    ucsc_hub: None
    create_heatmaps: bool = False


class CRISPRAssayConfig(BaseAssayConfig):
    """Configuration specific to CRISPR assays."""

    # CRISPR-specific options
    ucsc_hub: None = None
    create_heatmaps: bool = False
    use_mageck: bool = False

# Union type for all assay-specific configurations
AssaySpecificConfig = Union[
    ATACAssayConfig,
    ChIPAssayConfig,
    CATAssayConfig,
    RNAAssayConfig,
    SNPAssayConfig,
    MCCAssayConfig,
    MethylationAssayConfig,
    CRISPRAssayConfig,
]

# Class constant for assay config mapping to reduce redundancy
ASSAY_CONFIG_MAP = {
    Assay.ATAC: ATACAssayConfig,
    Assay.CHIP: ChIPAssayConfig,
    Assay.CAT: CATAssayConfig,
    Assay.RNA: RNAAssayConfig,
    Assay.SNP: SNPAssayConfig,
    Assay.MCC: MCCAssayConfig,
    Assay.METH: MethylationAssayConfig,
    Assay.CRISPR: CRISPRAssayConfig,
}


class AssayConfig(Enum):
    ATAC = ATACAssayConfig
    CHIP = ChIPAssayConfig
    CAT = CATAssayConfig
    RNA = RNAAssayConfig
    SNP = SNPAssayConfig
    MCC = MCCAssayConfig
    METH = MethylationAssayConfig
    CRISPR = CRISPRAssayConfig


class SeqnadoConfig(BaseModel):
    """Configuration for the SeqNado workflow."""

    assay: Assay
    project: ProjectConfig
    genome: GenomeConfig
    metadata: Path
    qc: QCConfig = Field(default_factory=QCConfig)
    pcr_duplicates: PCRDuplicatesConfig = Field(default_factory=PCRDuplicatesConfig)
    assay_config: AssaySpecificConfig | None = None
    third_party_tools: ThirdPartyToolsConfig | None = Field(None, description="Configuration for third-party tools.")

    # If no third_party_tools config is provided, use defaults
    @model_validator(mode="before")
    def set_default_third_party_tools(cls, values):
        if "third_party_tools" not in values or values["third_party_tools"] is None:
            assay = values.get("assay")
            # Normalize assay to Assay enum if it's a string
            if isinstance(assay, str):
                assay = Assay(assay)
            values["third_party_tools"] = ThirdPartyToolsConfig.for_assay(assay)

        return values

    @model_validator(mode="before")
    def set_default_pcr_duplicates(cls, values):
        """Set default PCR duplicate handling based on assay type."""
        from seqnado import PCRDuplicateHandling

        if "pcr_duplicates" not in values or values["pcr_duplicates"] is None:
            assay = values.get("assay")
            # Normalize assay to Assay enum if it's a string
            if isinstance(assay, str):
                assay = Assay(assay)
            # Default to REMOVE for ATAC, ChIP, CAT, SNP, and METH; KEEP for RNA
            if assay in [Assay.ATAC, Assay.CHIP, Assay.CAT, Assay.SNP, Assay.METH]:
                values["pcr_duplicates"] = PCRDuplicatesConfig(strategy=PCRDuplicateHandling.REMOVE)
            else:
                values["pcr_duplicates"] = PCRDuplicatesConfig(strategy=PCRDuplicateHandling.NONE)

        return values

    @model_validator(mode="before")
    def migrate_remove_blacklist_to_qc(cls, values):
        """Support legacy top-level remove_blacklist by moving it under qc."""
        if not isinstance(values, dict):
            return values

        if "remove_blacklist" in values:
            qc_values = dict(values.get("qc") or {})
            qc_values.setdefault("remove_blacklist", values.pop("remove_blacklist"))
            values["qc"] = qc_values

        return values

    @model_validator(mode="after")
    def sync_peak_caller_tool_defaults(self) -> "SeqnadoConfig":
        """Ensure peak-caller-specific tool configs exist when those methods are selected."""
        peak_calling = getattr(self.assay_config, "peak_calling", None)
        methods = getattr(peak_calling, "method", None) or []

        if self.third_party_tools is None:
            self.third_party_tools = ThirdPartyToolsConfig.for_assay(self.assay)

        if PeakCallingMethod.SEACR in methods and self.third_party_tools.seacr is None:
            self.third_party_tools.seacr = Seacr()

        return self

    @model_validator(mode="after")
    def validate_qc_remove_blacklist(self) -> "SeqnadoConfig":
        """Require a genome blacklist when blacklist removal is enabled."""
        if self.qc.remove_blacklist and not self.genome.blacklist:
            raise ValueError(
                "qc.remove_blacklist can only be True if genome blacklist is provided."
            )
        return self

    @classmethod
    def from_yaml(cls, path: Path) -> "SeqnadoConfig":
        """Load configuration from a YAML file."""
        import yaml

        with open(path, "r") as f:
            data = yaml.safe_load(f)

        return cls(**data)

    @computed_field
    @property
    def organism(self) -> str:
        """Return the organism (string) from the genome configuration."""
        return self.genome.organism

    @computed_field
    @property
    def shift_for_tn5_insertion(self) -> bool:
        """Return the Tn5 shift configuration for the specified assay."""
        return hasattr(self.assay_config, "tn5_shift") and self.assay_config.tn5_shift
    
    @computed_field
    @property
    def mcc_viewpoints(self) -> str:
        """Return the MCC viewpoints file path."""
        if self.assay_config and hasattr(self.assay_config, "mcc"):
            return str(self.assay_config.mcc.viewpoints)
        return ""

    @field_validator("assay_config", mode="before")
    def validate_assay_config_matches_assay(cls, v, info):
        """Ensure the assay_config type matches the specified assay."""
        if v is None:
            return v

        assay = info.data.get("assay")
        if assay is None:
            return v

        expected_config_class = ASSAY_CONFIG_MAP.get(assay)
        if expected_config_class and not isinstance(v, expected_config_class):
            if isinstance(v, dict):
                # Try to create the appropriate config from dict
                return expected_config_class(**v)
            else:
                raise ValueError(
                    f"assay_config must be of type {expected_config_class.__name__} for assay {assay.value}"
                )

        return v

    @classmethod
    def create_assay_config(cls, assay: Assay, **kwargs) -> AssaySpecificConfig:
        """Create the appropriate assay config for the given assay type."""
        config_class = cls._ASSAY_CONFIG_MAP.get(assay)
        if config_class is None:
            raise ValueError(
                f"No configuration class available for assay {assay.value}"
            )

        return config_class(**kwargs)
