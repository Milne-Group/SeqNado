"""Core enums, constants, and utility functions for SeqNado design module."""

from pathlib import Path 
import re
from enum import Enum
from typing import Optional, Union, Callable, Self, TYPE_CHECKING, Any
from abc import ABC, abstractmethod
from pydantic import BaseModel, Field, computed_field, field_validator, model_validator, ValidationInfo
from seqnado import Assay, Organism
import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================

ILLUMINA_FILENAME_PATTERNS = {
    r"_S\d+_": "_",
    r"_L00\d_": "_",
    r"_R?[12](_001)?$": "_",
    r"__": "_",
    r"_$": "",
}

INPUT_CONTROL_SUBSTRINGS = ["input", "mock", "igg", "control"]


# =============================================================================
# Models
# =============================================================================

class Metadata(BaseModel):
    """
    Metadata for samples. Optional fields can be set to None.

    Field relationships:
    - `consensus_group`: Groups samples for BAM merging (e.g., replicates)
    - `condition`: Groups samples for bigwig comparisons when perform_comparisons=true
    - `group`/`deseq2`: Used for differential expression analysis (RNA-seq)
    - `scaling_group`: Groups samples for normalization factor calculation

    These fields can be used independently or in combination.
    """
    assay: Assay | None = Field(
        default=None,
        description="Assay type, should be one of the Assay enum values"
    )
    consensus_group: str | None = Field(
        default=None,
        description="Grouping variable for merging samples into consensus, defaults to None if not specified"
    )
    scaling_group: str = Field(
        default="default",
        description="Grouping variable for scaling samples, defaults to 'default' if not specified"
    )
    group: str | None = Field(
        default=None,
        description="Experimental group name (e.g., control, treated, WT, KO), used for DESeq2 analysis"
    )
    deseq2: int | None = Field(
        default=None,
        description="DESeq2 binary encoding: 0 for control/reference group, 1 for treatment/comparison group"
    )
    condition: str | None = Field(
        default=None,
        description="Biological condition or experimental group (e.g., control, treated, WT, KO). When perform_comparisons is enabled, samples with the same condition are aggregated into condition-level bigwigs, and pairwise subtractions are generated across conditions. Requires at least 2 unique condition values."
    )

    @field_validator("group", "deseq2", mode="before")
    @classmethod
    def coerce_none_values(cls, v, info: ValidationInfo):
        import numpy as np
        import pandas as pd
        """Convert none-like values to None for optional fields."""
        # Define a list of values that should be treated as 'None'
        # This includes None, pd.NA, np.nan, and common string representations of None

        none_vals = [
            None,
            "None",
            "none",
            "null",
            "Null",
            "NULL",
            ".",
            "",
            "NA",
        ]
        
        # Check for pandas NA
        if v is pd.NA:
            return None
            
        # Check for numpy nan
        if isinstance(v, float) and np.isnan(v):
            return None
            
        # Check for string/literal None values
        if v in none_vals:
            return None
        
        # For deseq2, check if validation is required
        # Only validate deseq2 for RNA assays when spike-in is configured
        if info.field_name == "deseq2":
            validate_deseq2 = info.context.get('validate_deseq2', False) if info.context else False
            validate_assay = info.context.get('assay', None) if info.context else None
            
            # Only enforce deseq2 validation for RNA assays with spike-in
            if validate_deseq2 and validate_assay == Assay.RNA and v is None:
                raise ValueError("deseq2 field is required for RNA assay when spikein is configured")
            
        return v

    @model_validator(mode='after')
    def set_mcc_defaults(self) -> Self:
        """Set default consensus_group for MCC assay."""
        if self.assay == Assay.MCC and self.consensus_group is None:
            self.consensus_group = 'default'
        return self



#==============================================================================
# Base class for collections
#==============================================================================
class BaseCollection(BaseModel):
    """
    Base class for all design types providing common functionality
    that is *not* tied to FASTQ files.
    """

    assay: Assay
    metadata: list[Metadata]

    # --------- Common metadata utilities (format-agnostic) ---------

    @classmethod
    def _build_metadata(
        cls,
        sample_name: str,
        metadata: Callable[[str], Metadata] | Metadata | None,
        assay: Assay,
    ) -> Metadata:
        """Build metadata for a sample ensuring assay is always set.

        Rules:
        - If callable: call with sample_name to get Metadata.
        - If Metadata instance: use directly.
        - If None: create default.
        In all cases force metadata.assay to the provided assay if it's None or different.
        """
        if callable(metadata):
            md = metadata(sample_name)
        elif isinstance(metadata, Metadata):
            md = metadata
        else:
            md = Metadata()
        # Always stamp assay (requirement: always include assay in metadata)
        if md.assay != assay:
            md.assay = assay
        return md

    @classmethod
    def _prepare_metadata_for_directory(
        cls,
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> Callable[[str], Metadata] | Metadata:
        """Prepare metadata parameter for directory-based construction."""
        if not callable(metadata) and not isinstance(metadata, Metadata):
            metadata = Metadata(**kwargs)
        return metadata

    # ----------------------- I/O scaffolding -----------------------

    @classmethod
    def from_csv(cls, path: str | Path, *args, **kwargs) -> Self:
        """Build a collection from a CSV file by delegating to from_dataframe."""
        df = pd.read_csv(path)
        # Extract validate_deseq2 and assay_for_validation if present in kwargs for passing to from_dataframe
        validate_deseq2 = kwargs.pop('validate_deseq2', False)
        assay_for_validation = kwargs.pop('assay_for_validation', None)
        return cls.from_dataframe(df=df, validate_deseq2=validate_deseq2, assay_for_validation=assay_for_validation, *args, **kwargs)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame) -> Self:  # type: ignore[override]
        """Build a collection from a pandas DataFrame."""
        raise NotImplementedError("Subclasses must implement from_dataframe")

    def to_dataframe(self) -> pd.DataFrame:
        """Export the design to a pandas DataFrame."""
        raise NotImplementedError("Subclasses must implement to_dataframe")

    # --------------------- Common abstract API ---------------------

    @property
    def sample_names(self) -> list[str]:
        raise NotImplementedError("Subclasses must implement sample_names")



# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def predict_organism(genome: str) -> Organism:
    """Predict organism from genome identifier."""
    if "hg" in genome:
        return Organism.HUMAN
    elif "mm" in genome:
        return Organism.MOUSE
    else:
        return Organism.UNKNOWN


def is_valid_path(path: str | Path | None) -> bool:
    """Check if a path is valid and exists."""
    if path is None:
        return False
    
    try:
        p = Path(path)
        return p.exists() and str(p) not in ["-", ".", "", "None"]
    except (TypeError, ValueError):
        return False


# =============================================================================
# STRING/FILENAME PROCESSING FUNCTIONS
# =============================================================================

def clean_sample_name(name: str, patterns: dict = ILLUMINA_FILENAME_PATTERNS) -> str:
    """Clean sample name using regex patterns."""
    for pattern, replacement in patterns.items():
        name = re.sub(pattern, replacement, name)
    return name


def extract_read_number(filename: str) -> Optional[int]:
    """Extract read number from FASTQ filename."""
    regex_std_illumina_paired = re.compile(r".*_R?([12])(_001)?")
    
    match = regex_std_illumina_paired.match(filename)
    if match:
        return int(match.group(1))
    return None


def is_control_sample(ip_name: str) -> bool:
    """Check if sample is a control based on IP name."""
    return any(substring in ip_name.lower() for substring in INPUT_CONTROL_SUBSTRINGS)


