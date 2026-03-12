"""FASTQ file handling classes for genomic sequencing workflows."""

from __future__ import annotations

import re
from collections import defaultdict
from itertools import chain
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd
from loguru import logger
from pandera.typing import DataFrame
from pydantic import BaseModel, Field, computed_field, field_validator

from .core import (
    Assay,
    BaseCollection,
    Metadata,
    clean_sample_name,
    extract_read_number,
    is_control_sample,
    is_valid_path,
)
from .validation import DesignDataFrame


# =============================================================================
# BASE FASTQ CLASSES
# =============================================================================


class FastqFile(BaseModel):
    """Represents a single FASTQ file with metadata extraction."""

    path: Path
    use_resolved_name: bool = False

    def model_post_init(self, __context: dict[str, any] | None) -> None:
        """Validate file path after initialization."""
        if self.use_resolved_name:
            self.path = self.path.resolve()
        else:
            self.path = self.path.absolute()

        if not is_valid_path(self.path):
            raise FileNotFoundError(f"FASTQ file not found: {self.path}")

    @property
    def stem(self) -> str:
        """Get file stem without .gz extension."""
        return Path(str(self.path).removesuffix(".gz")).stem

    @computed_field
    @property
    def filename_base(self) -> str:
        """Extract base filename identifier by removing common suffixes."""
        name = self.stem
        if name.endswith("_001"):
            name = name.removesuffix("_001")
        return name

    @computed_field
    @property
    def sample_base(self) -> str:
        """Extract base sample name by cleaning Illumina patterns."""
        return clean_sample_name(self.filename_base)

    @computed_field
    @property
    def read_number(self) -> int | None:
        """Extract read number (1 or 2) from filename."""
        read_num = extract_read_number(self.filename_base)
        if read_num is None:
            logger.debug(
                f"Could not extract read number from '{self.filename_base}', "
                f"assuming single-end sequencing"
            )
        return read_num

    @computed_field
    @property
    def is_paired(self) -> bool:
        """Check if file is part of a paired-end experiment."""
        return self.read_number is not None

    @computed_field
    @property
    def is_lane(self) -> bool:
        """Check if file is lane-split."""
        return "_L00" in self.filename_base

    def __lt__(self, other: FastqFile) -> bool:
        """Compare FastqFile objects by path."""
        return self.path < other.path

    def __gt__(self, other: FastqFile) -> bool:
        """Compare FastqFile objects by path."""
        return self.path > other.path

    def __eq__(self, other: object) -> bool:
        """Check equality of FastqFile objects by path."""
        if not isinstance(other, FastqFile):
            return NotImplemented
        return self.path == other.path

    def __hash__(self) -> int:
        """Make FastqFile hashable based on path."""
        return hash(self.path)


# =============================================================================
# IP EXPERIMENT CLASSES
# =============================================================================


class FastqFileIP(FastqFile):
    """FASTQ file for IP experiments with antibody information."""

    ip: str = Field(default=None, description="IP antibody performed on the sample")
    is_control: bool = Field(default=None, description="Whether sample is a control")

    def model_post_init(self, __context: dict[str, any] | None) -> None:
        """Initialize and predict IP information."""
        super().model_post_init(__context)

        if self.ip is None:
            self.ip = self._predict_ip()

        if self.is_control is None:
            self.is_control = self._predict_is_control()

    def _predict_ip(self) -> str | None:
        """Predict IP antibody from sample name."""
        try:
            return self.sample_base.split("_")[-1]
        except IndexError:
            logger.warning(f"Could not predict IP for {self.sample_base}")
            return None

    def _predict_is_control(self) -> bool:
        """Predict if sample is a control."""
        if not self.ip:
            return False
        return is_control_sample(self.ip)

    @computed_field
    @property
    def sample_base_without_ip(self) -> str:
        """Get sample base name without antibody suffix."""
        if not self.ip:
            return self.sample_base

        # Remove IP and common Illumina suffixes
        pattern = rf"(_{re.escape(self.ip)})?(_S\d+)?(_L00\d)?(_R?[12])?(_001)?"
        return re.sub(pattern, "", self.filename_base)

    @field_validator("ip")
    @classmethod
    def allow_na_or_nan(cls, v: str | None | float) -> str | None | float:
        """Allow None, pd.NA, or np.nan values for IP."""
        if v is None or v is pd.NA or (isinstance(v, float) and np.isnan(v)):
            return v
        if not isinstance(v, str):
            raise ValueError("ip must be a string, None, pd.NA, or np.nan")
        return v


class ExperimentIP(BaseModel):
    ip: FastqSetIP
    control: FastqSetIP | None = None

    @property
    def has_control(self) -> bool:
        """Check if the experiment has a control sample."""
        return self.control is not None

    @property
    def ip_set_fullname(self) -> str:
        """Get the full sample name for the IP sample including the antibody used."""
        return self.ip.full_sample_name

    @property
    def control_fullname(self) -> str:
        """Get the full sample name for the control sample including the control performed."""
        return self.control.full_sample_name if self.control else None

    @property
    def ip_performed(self) -> str:
        """Get the antibody used for the IP sample."""
        return self.ip.ip

    @property
    def control_performed(self) -> str:
        """Get the antibody used for the control sample."""
        return self.control.ip if self.control else None

    @property
    def fastqs_are_paired(self) -> bool:
        """Check if both IP and control samples are paired-end."""
        ip = self.ip.is_paired
        control = self.control.is_paired if self.control else True
        return ip and control


# =============================================================================
# FASTQ SET CLASSES
# =============================================================================


class FastqSet(BaseModel):
    """Represents a set of FASTQ files for a single sample."""

    sample_id: str = Field(default=None, description="Base sample identifier")
    r1: FastqFile
    r2: FastqFile | None = None

    @property
    def file_paths(self) -> list[Path]:
        """Get list of FASTQ file paths."""
        paths = [self.r1.path]
        if self.is_paired and self.r2:
            paths.append(self.r2.path)
        return paths

    @property
    def is_paired(self) -> bool:
        """Check if sample has paired-end reads."""
        return self.r2 is not None and self.r2.path.exists()

    @classmethod
    def from_fastq_files(cls, fq: list[FastqFile], **kwargs) -> FastqSet:
        """Create FastqSet from list of FastqFile objects."""
        if not fq:
            raise ValueError("No FASTQ files provided")

        sample_id = fq[0].sample_base

        match len(fq):
            case 1:
                return cls(sample_id=sample_id, r1=fq[0], **kwargs)
            case 2:
                return cls(sample_id=sample_id, r1=fq[0], r2=fq[1], **kwargs)
            case _:
                raise ValueError(
                    f"Invalid number of FASTQ files for {sample_id}: {len(fq)}. "
                    f"Expected 1 or 2 files."
                )


class FastqSetIP(FastqSet):
    """FASTQ set for IP experiments."""

    r1: FastqFileIP
    r2: FastqFileIP | None = None
    ip: str = Field(default=None, description="IP performed")

    def model_post_init(self, __context: dict[str, any] | None) -> None:
        """Initialize and predict information."""
        super().model_post_init(__context)

        if self.ip is None:
            self.ip = self._predict_ip()

    def _predict_ip(self) -> str:
        """Predict the target IP from the r1 file."""
        return self.r1.ip

    @property
    def full_sample_name(self) -> str:
        """Get complete sample name including IP/control information."""
        return f"{self.sample_id}_{self.ip}"

    @property
    def base_sample_name(self) -> str:
        """Get base sample name without IP information."""
        return self.r1.sample_base_without_ip

    @property
    def is_control(self) -> bool:
        """Check if this is a control sample."""
        return self.r1.is_control


class BaseFastqCollection(BaseCollection):
    """
    FASTQ-specific mixin/base providing helpers and expectations
    for collections that manage FASTQ files.
    """

    # -------- FASTQ discovery (FASTQ-specific) --------

    @classmethod
    def _discover_files(
        cls,
        directory: str | Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
    ) -> list[Path]:
        """Discover FASTQ files in a directory using glob patterns."""
        dir_path = Path(directory)
        files = list(chain.from_iterable(dir_path.glob(p) for p in glob_patterns))
        if not files:
            raise FileNotFoundError(f"No FASTQ files found in {dir_path}")
        return files

    # -------- FASTQ-facing API surface --------

    @property
    def fastq_pairs(self) -> dict[str, list[Path]]:
        """Map sample name -> list of FASTQ paths (e.g., [R1, R2])."""
        raise NotImplementedError("Subclasses must implement fastq_pairs")

    @property
    def fastq_paths(self) -> list[Path]:
        """Flattened list of FASTQ paths across all samples."""
        return list(chain.from_iterable(self.fastq_pairs.values()))

    def symlink_fastq_files(self, output_dir: str | Path) -> None:
        """Symlink FASTQ files to a specified output directory.

        Args:
            output_dir: Directory to create symlinks in.
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        df = self.to_dataframe()
        for uid, row in df.iterrows():
            # For single-end data, don't add read number suffix
            if pd.isna(row["r2"]):
                symlink_path = output_dir / f"{uid}.fastq.gz"
                if not symlink_path.exists():
                    symlink_path.symlink_to(Path(row["r1"]).resolve().absolute())
            else:
                # For paired-end data, add _1 and _2 suffixes
                symlink_path_r1 = output_dir / f"{uid}_1.fastq.gz"
                if not symlink_path_r1.exists():
                    symlink_path_r1.symlink_to(Path(row["r1"]).resolve().absolute())
                symlink_path_r2 = output_dir / f"{uid}_2.fastq.gz"
                if not symlink_path_r2.exists():
                    symlink_path_r2.symlink_to(Path(row["r2"]).resolve().absolute())


class FastqCollection(BaseFastqCollection):
    """
    Represents a collection of sequencing samples (FASTQ files) grouped into named sets,
    with optional per-sample metadata.

    Attributes:
        fastq_sets: List of FastqSet objects (paired or single-end samples).
        metadata:  List of Metadata objects corresponding one-to-one with fastq_sets.
    """

    fastq_sets: list[FastqSet]

    @property
    def primary_file_type(self) -> str:
        return "fastq"

    @field_validator("assay")
    @classmethod
    def validate_non_ip_assay(cls, v: Assay) -> Assay:
        """Ensure the assay doesn't require IP (immunoprecipitation)."""
        if v in Assay.ip_assays():
            raise ValueError(
                f"Assay '{v.value}' requires IP and should use IPSampleCollection instead"
            )
        return v

    @property
    def sample_ids(self) -> list[str]:
        """
        Returns all sample IDs in the design.
        """
        return [fs.sample_id for fs in self.fastq_sets]

    @property
    def sample_names(self) -> list[str]:
        """
        Returns all sample names in the design.
        """
        return self.sample_ids

    @property
    def fastq_paths(self) -> list[Path]:
        """
        Flattens all R1/R2 file paths into a single list.
        """
        return [
            path
            for fs in self.fastq_sets
            for path in ([fs.r1.path] + ([fs.r2.path] if fs.r2 else []))
        ]

    def get_file_paths(self, kind: str | None = None) -> list[Path]:
        kind = kind or "fastq"
        match kind:
            case "fastq":
                return self.fastq_paths
            case "fastq_r1":
                return [fs.r1.path for fs in self.fastq_sets]
            case "fastq_r2":
                return [fs.r2.path for fs in self.fastq_sets if fs.r2]
            case _:
                raise ValueError(f"Unsupported file kind '{kind}' for SampleCollection")

    @property
    def fastq_pairs(self) -> dict[str, list[Path]]:
        """
        Returns a dictionary mapping sample names to their FASTQ file paths.
        """
        return {
            fs.sample_id: [fs.r1.path, fs.r2.path] if fs.r2 else [fs.r1.path]
            for fs in self.fastq_sets
        }

    def query(self, sample_name: str) -> FastqSet:
        """
        Retrieve the FastqSet by its sample name.

        Raises:
            ValueError: If sample_name not found.
        """
        try:
            return next(fs for fs in self.fastq_sets if fs.sample_id == sample_name)
        except StopIteration:
            raise ValueError(f"Sample '{sample_name}' not found in SampleCollection")

    def is_paired_end(self, uid: str) -> bool:
        """
        Check if the given sample ID is paired-end.
        """
        return self.to_dataframe().loc[uid, "r2"] is not None

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **fastqset_kwargs: Any,
    ) -> FastqCollection:
        """
        Build a SampleCollection by scanning a list of FASTQ paths:

        1. Convert raw paths to FastqFile.
        2. Group by `sample_base` and sort by read_number.
        3. Create FastqSet (single- or paired-end) for each sample.
        4. Generate Metadata via `metadata(sample_name)`, or default.

        Args:
            files: Iterable of file paths (strings or Path).
            metadata:
                - Callable(sample_name) → Metadata to customize per-sample metadata.
                - Single Metadata instance applied to all.
                - None → defaults to Metadata().
            fastqset_kwargs: Extra fields forwarded to FastqSet constructor.
        """
        # Convert and sort
        fq_files = [FastqFile(path=Path(f)) for f in files]
        fq_files.sort(key=lambda x: (x.sample_base, x.read_number))

        # Group by sample_stem
        groups: dict[str, list[FastqFile]] = defaultdict(list)
        for fq in fq_files:
            groups[fq.sample_base].append(fq)

        _fastq_sets: list[FastqSet] = []
        _metadata: list[Metadata] = []
        for sample, fqs in groups.items():
            # Build FastqSet
            if len(fqs) == 1:
                fs = FastqSet(sample_id=sample, r1=fqs[0], **fastqset_kwargs)
            elif len(fqs) == 2:
                fs = FastqSet(sample_id=sample, r1=fqs[0], r2=fqs[1], **fastqset_kwargs)
            else:
                raise ValueError(
                    f"Unexpected number of FASTQ files for '{sample}': {len(fqs)}"
                )
            _fastq_sets.append(fs)

            # Build Metadata using base class method
            _metadata.append(cls._build_metadata(sample, metadata, assay))

        return cls(assay=assay, fastq_sets=_fastq_sets, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> FastqCollection:
        """
        Recursively scan a directory for FASTQ files and build a SampleCollection.

        Args:
            directory: Root path to search.
            glob_patterns: Filename patterns to include.
            metadata: Callable(sample_name) → Metadata or single Metadata instance.
            **kwargs: Extra fields converted directly to a shared Metadata.
        """
        files = cls._discover_files(directory, glob_patterns)
        metadata = cls._prepare_metadata_for_directory(metadata, **kwargs)
        return cls.from_fastq_files(assay=assay, files=files, metadata=metadata)

    def to_dataframe(self, validate: bool = True) -> pd.DataFrame:
        """
        Export the design to a pandas DataFrame, validated by DataFrameDesign.

        Columns: sample_name, r1, r2, plus all metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []

        if self.metadata:
            for fs, md in zip(self.fastq_sets, self.metadata):
                row: dict[str, Any] = {
                    "sample_id": fs.sample_id,
                    "r1": fs.r1.path,
                    "r2": fs.r2.path if fs.r2 else None,
                    "uid": f"{fs.sample_id}",
                }
                metadata_dict = md.model_dump(exclude_none=True)
                # Convert Assay enum to string value for schema validation
                if "assay" in metadata_dict and hasattr(metadata_dict["assay"], "value"):
                    metadata_dict["assay"] = metadata_dict["assay"].value
                row.update(metadata_dict)
                rows.append(row)
        else:
            for fs in self.fastq_sets:
                row = {
                    "sample_id": fs.sample_id,
                    "r1": fs.r1.path,
                    "r2": fs.r2.path if fs.r2 else None,
                    "uid": f"{fs.sample_id}",
                }
                rows.append(row)

        if not rows:
            # Return empty DataFrame with expected columns
            df = pd.DataFrame(columns=["sample_id", "r1", "r2", "uid"]).set_index("uid")
        else:
            df = pd.DataFrame(rows).sort_values("sample_id").set_index("uid")

        # Define column order: critical columns first (assay, sample info, files), then other metadata
        core_cols = ["assay", "sample_id", "r1", "r2"]
        metadata_cols = [col for col in df.columns if col not in core_cols]
        ordered_cols = core_cols + sorted(metadata_cols)
        df = df[[col for col in ordered_cols if col in df.columns]]

        if validate:
            return DataFrame[DesignDataFrame](df)
        else:
            return df

    @classmethod
    def from_dataframe(
        cls, assay: Assay, df: pd.DataFrame, validate_deseq2: bool = False, assay_for_validation: Assay | None = None, **fastqset_kwargs: Any
    ) -> FastqCollection:
        """
        Build a SampleCollection from a DataFrame, validated by DataFrameDesign.

        Expects columns: sample_name, r1, r2, plus any metadata fields.
        
        Args:
            assay: The assay type
            df: DataFrame with sample metadata
            validate_deseq2: If True, require deseq2 field to be non-null (for RNA assays)
            assay_for_validation: Assay type to check in validation context
            **fastqset_kwargs: Additional kwargs for FastqSet
        """
        df = DesignDataFrame.validate(df)
        fastq_sets: list[FastqSet] = []
        metadata: list[Metadata] = []
        metadata_fields = set(Metadata.model_fields.keys())
        
        # Use provided assay_for_validation or fall back to assay
        validation_assay = assay_for_validation or assay

        for rec in df.to_dict(orient="records"):
            # Build FastqSet
            r2_path = rec.get("r2")
            fs = FastqSet(
                sample_id=rec["sample_id"],
                r1=FastqFile(path=rec["r1"]),
                r2=FastqFile(path=r2_path) if pd.notna(r2_path) else None,
                **fastqset_kwargs,
            )
            fastq_sets.append(fs)

            # Collect metadata with validation context
            meta_fields = {k: rec.get(k) for k in metadata_fields if k in rec}
            metadata.append(Metadata.model_validate(meta_fields, context={'validate_deseq2': validate_deseq2, 'assay': validation_assay}))

        return cls(assay=assay, fastq_sets=fastq_sets, metadata=metadata)


def create_experiments_from_fastqs(
    ip_list: list[FastqFileIP],
    ctrl_list: list[FastqFileIP],
    name: str,
    **exp_kwargs: Any,
) -> ExperimentIP:
    """
    Create an ExperimentIP from lists of IP and control FastqFileIP objects.
    Args:
        ip_list: List of FastqFileIP for the IP sample (1 or 2 files).
        ctrl_list: List of FastqFileIP for the control sample (0, 1, or 2 files).
        name: Base sample name for the experiment.
        exp_kwargs: Additional fields for ExperimentIP.
    Returns:
        ExperimentIP object with populated IP and control FastqSetIP.
    """

    n_ip_fastqs = len(ip_list)
    n_control_fastqs = len(ctrl_list)

    match (n_ip_fastqs, n_control_fastqs):
        case 0, _:
            raise ValueError(f"No IP FASTQ files found for sample '{name}'")
        case ((1 | 2), 0):
            ip_set = FastqSetIP(
                sample_id=name,
                r1=ip_list[0],
                r2=ip_list[1] if len(ip_list) > 1 else None,
            )
            ctrl_set = None
        case ((1 | 2), (1 | 2)):
            ip_set = FastqSetIP(
                sample_id=name,
                r1=ip_list[0],
                r2=ip_list[1] if len(ip_list) > 1 else None,
            )
            ctrl_set = FastqSetIP(
                sample_id=name,
                r1=ctrl_list[0],
                r2=ctrl_list[1] if len(ctrl_list) > 1 else None,
            )

        case _:
            raise ValueError(
                f"Unexpected number of FASTQ files for sample '{name}': "
                f"IP files={n_ip_fastqs}, control files={n_control_fastqs}"
            )
    return ExperimentIP(ip=ip_set, control=ctrl_set, **exp_kwargs)


class FastqCollectionForIP(BaseFastqCollection):
    """
    Represents an IP (e.g., ChIP/CAT) experiment design, consisting of paired IP/control FastqSetIP objects and per-experiment metadata.
        assay (Assay): Assay type (e.g., ChIP, CAT).
        experiments (list[ExperimentIP]): List of IPExperiment objects, each containing .ip and optional .control FastqSetIP.
        metadata (list[Metadata]): List of Metadata objects matching the experiments list.
        sample_names (list[str]): All unique IP and control set names, sorted.
        ips_performed (list[str]): Unique IP antibodies/proteins used across experiments.
        controls_performed (list[str]): Unique control antibodies used across experiments.
        query(sample_name, full=False):
            Scan a directory for IP/control FASTQ files and build an IPSampleCollection.
        to_dataframe():
    Example:
        >>> collection = IPSampleCollection.from_directory(
        ...     assay=Assay.CHIP,
        ...     directory="/data/fastq",
        ...     glob_patterns=["*.fastq.gz"]
        ... )
        >>> df = collection.to_dataframe()
    """

    experiments: list[ExperimentIP]

    @property
    def primary_file_type(self) -> str:
        return "fastq"

    @field_validator("assay")
    @classmethod
    def validate_assay(cls, assay: Assay) -> Assay:
        if assay not in Assay.ip_assays():
            raise ValueError(
                f"Assay '{assay.value}' should use `SampleCollection` instead"
            )
        return assay

    @property
    def sample_ids(self) -> list[str]:
        """
        All unique IP and control set IDs, sorted.
        """
        ip_names = [exp.ip_set_fullname for exp in self.experiments]
        ctrl_names = [
            exp.control_fullname for exp in self.experiments if exp.has_control
        ]
        return sorted(set(ip_names + ctrl_names))

    @property
    def sample_names(self) -> list[str]:
        """
        All unique IP and control set names, sorted.
        """
        return self.sample_ids

    @property
    def ip_sample_names(self) -> list[str]:
        """
        Only IP sample names (excludes control samples).
        Use this for peak calling to avoid calling peaks on input controls.
        """
        return sorted({exp.ip_set_fullname for exp in self.experiments})

    @property
    def ips_performed(self) -> list[str]:
        """
        Unique IP antibodies/proteins used across experiments.
        """
        return list({exp.ip_performed for exp in self.experiments})

    @property
    def controls_performed(self) -> list[str]:
        """
        Unique control antibodies used across experiments.
        """
        return list(
            {exp.control_performed for exp in self.experiments if exp.has_control}
        )

    @property
    def fastq_pairs(self) -> dict[str, list[Path]]:
        """
        Returns a dictionary mapping sample names to their FASTQ file paths.
        """
        files = dict()
        for exp in self.experiments:
            files[exp.ip_set_fullname] = (
                [exp.ip.r1.path, exp.ip.r2.path] if exp.ip.r2 else [exp.ip.r1.path]
            )
            if exp.has_control:
                files[exp.control_fullname] = (
                    [exp.control.r1.path, exp.control.r2.path]
                    if exp.control.r2
                    else [exp.control.r1.path]
                )
        return files

    def get_file_paths(self, kind: str | None = None) -> list[Path]:
        kind = kind or "fastq"
        all_fastqs = list({p for paths in self.fastq_pairs.values() for p in paths})
        match kind:
            case "fastq":
                return sorted(all_fastqs)
            case "fastq_r1":
                return sorted(
                    {exp.ip.r1.path for exp in self.experiments}
                    | {exp.control.r1.path for exp in self.experiments if exp.control}
                )
            case "fastq_r2":
                return sorted(
                    {
                        p
                        for exp in self.experiments
                        for p in [
                            exp.ip.r2.path if exp.ip.r2 else None,
                            exp.control.r2.path
                            if exp.control and exp.control.r2
                            else None,
                        ]
                        if p is not None
                    }
                )
            case _:
                raise ValueError(
                    f"Unsupported file kind '{kind}' for SampleCollectionForIP"
                )

    def get_control_performed(self, uid: str) -> str:
        """
        Get the control name for the given sample (IP or control).
        Returns the control fullname associated with this sample.
        """
        # Try to find in experiments
        for exp in self.experiments:
            # If it's an IP sample, return its control
            if exp.ip_set_fullname == uid:
                return exp.control_fullname if exp.has_control else None
            # If it's a control sample, return itself
            if exp.has_control and exp.control_fullname == uid:
                return exp.control_fullname

        # Fallback to DataFrame lookup
        df = self.to_dataframe()
        if uid in df.index:
            return df.loc[uid, "control"]

        raise KeyError(f"Sample '{uid}' not found in FastqCollectionForIP")

    def is_paired_end(self, uid: str) -> bool:
        """
        Check if the given sample ID is paired-end.
        For IP samples, checks the IP fastq r2.
        For control samples, checks the control fastq r2.
        """
        # Try to find the sample in experiments directly
        for exp in self.experiments:
            # Check if it's the IP sample
            if exp.ip_set_fullname == uid:
                return exp.ip.is_paired
            # Check if it's the control sample
            if exp.has_control and exp.control_fullname == uid:
                return exp.control.is_paired

        # Fallback: try DataFrame lookup (for backwards compatibility)
        df = self.to_dataframe()
        if uid in df.index:
            return df.loc[uid, "r2"] is not None

        raise KeyError(f"Sample '{uid}' not found in FastqCollectionForIP")

    def query(
        self,
        sample_name: str,
        full: bool = False,
    ) -> FastqSetIP | dict[str, FastqSetIP]:
        """
        Retrieve IP (and optionally control) FastqSetIP by sample name.

        Args:
            sample_name: IP or control set fullname.
            full: If True, return dict {"ip": ..., "control": ...}; else return the matched set only.

        Raises:
            KeyError: If no match.
        """
        for exp in self.experiments:
            if exp.ip_set_fullname == sample_name:
                out = {"ip": exp.ip, "control": exp.control}
                return out if full else exp.ip
            if exp.has_control and exp.control_fullname == sample_name:
                out = {"ip": exp.ip, "control": exp.control}
                return out if full else exp.control
        raise KeyError(f"Sample '{sample_name}' not found in IPSampleCollection")

    @classmethod
    def from_fastq_files(
        cls,
        assay: Assay,
        files: Iterable[str | Path],
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        ip_to_control_map: dict[str, str] | None = None,
        **exp_kwargs: Any,
    ) -> FastqCollectionForIP:
        """
        Build IPSampleCollection from a mixture of IP/control FASTQ paths.

        Groups by FastqFileIP.sample_base_without_ip and read number.

        Args:
            assay: The assay type (must be an IP assay).
            files: Iterable of file paths (strings or Path).
            metadata:
                - Callable(sample_name) → Metadata to customize per-sample metadata.
                - Single Metadata instance applied to all.
                - None → defaults to Metadata().
            ip_to_control_map: Optional mapping from IP antibody names to control antibody names.
            exp_kwargs: Extra fields forwarded to IPExperiment constructor.
        """
        # Convert and sort
        fastqs: list[FastqFileIP] = [FastqFileIP(path=Path(f)) for f in files]

        # Bucket by sample_base_without_ip
        buckets: dict[str, dict[str, list[FastqFileIP]]] = defaultdict(
            lambda: {"ip": [], "control": []}
        )
        for f in fastqs:
            key = f.sample_base_without_ip
            side = "control" if f.is_control else "ip"
            buckets[key][side].append(f)

        experiments: list[ExperimentIP] = []
        _metadata: list[Metadata] = []
        for name, sides in buckets.items():
            # Sort by read_number
            ip_list = sorted(sides["ip"], key=lambda x: x.read_number or 0)
            ctrl_list = sorted(sides["control"], key=lambda x: x.read_number or 0)

            if not ip_list and not ctrl_list:
                raise ValueError(f"No valid FASTQ files found for sample '{name}'")
            elif not ip_list and ctrl_list:
                ip_list = ctrl_list
                ctrl_list = []

            # Build ExperimentIP
            n_ip_fastqs = len(ip_list)
            n_ctrl_fastqs = len(ctrl_list)

            if n_ip_fastqs <= 2 and n_ctrl_fastqs <= 2:
                # Simple case, only two files exist from the same condition, just pair them up if the controls exist
                # if not create empty control
                exp = create_experiments_from_fastqs(ip_list, ctrl_list, name, **exp_kwargs)
                experiments.append(exp)
                 # Build Metadata using base class method
                _metadata.append(cls._build_metadata(name, metadata, assay))
            else:
                # Complex case, more than 2 files exist for either IP or control.
                # Likely same condition by muliple antibodies
                ip_groups: dict[str, list[FastqFileIP]] = defaultdict(list)
                ctrl_groups: dict[str, list[FastqFileIP]] = defaultdict(list)
                for f in ip_list:
                    ip_groups[f.ip].append(f)
                
                # If we have a mapping from IP to control, use it to group controls
                if ip_to_control_map:
                    for f in ctrl_list:
                        for ip_ab, ctrl_ab in ip_to_control_map.items():
                            if f.ip == ctrl_ab:
                                ctrl_groups[ip_ab].append(f)
                
                # We don't have a mapping, this is fine if we have a single control.
                # In this case we will assign the same control to all IPs
                # If multiple controls exist without a mapping, raise an error, we can't disambiguate
                # Will warn the user to provide a mapping
                elif len(ctrl_list) > 2:
                    raise ValueError(
                        f"Multiple control FASTQ files found for sample '{name}' without a mapping from IP to control."
                        f" Please provide an IP to control mapping to allow for correct assignment of controls to IPs."
                        f" Hint: IP antibodies {set([f.ip for f in ip_list])} control antibodies: {set([f.ip for f in ctrl_list])}"
                        f" If using seqnado CLI, provide the mapping via the --ip-to-control flag 'IP1:CTRL1,IP2:CTRL2"
                        f" e.g., --ip-to-control {ip_list[0].ip}:{ctrl_list[0].ip}"
                    )
                else:
                    for f in ctrl_list:
                        for ip_ab in ip_groups.keys():
                            ctrl_groups[ip_ab].append(f)
                
                for ip_ab, ip_files in ip_groups.items():
                    ctrl_files = ctrl_groups.get(ip_ab, [])
                    exp = create_experiments_from_fastqs(ip_files, ctrl_files, name, **exp_kwargs)
                    experiments.append(exp)
                     # Build Metadata using base class method
                    _metadata.append(cls._build_metadata(name, metadata, assay))


        return cls(assay=assay, experiments=experiments, metadata=_metadata)

    @classmethod
    def from_directory(
        cls,
        assay: Assay,
        directory: str | Path,
        glob_patterns: Iterable[str] = ("*.fq", "*.fq.gz", "*.fastq", "*.fastq.gz"),
        metadata: Callable[[str], Metadata] | Metadata | None = None,
        **kwargs: Any,
    ) -> FastqCollectionForIP:
        """
        Scan a directory for IP/control FASTQ files and build a IPSampleCollection.

        Args:
            directory: Root path to search.
            glob_patterns: Filename patterns to include.
            metadata: Callable(sample_name) → Metadata or single Metadata instance.
            **kwargs: Extra fields converted directly to a shared Metadata.
        """
        files = cls._discover_files(directory, glob_patterns)
        metadata = cls._prepare_metadata_for_directory(metadata, **kwargs)
        return cls.from_fastq_files(assay=assay, files=files, metadata=metadata)

    def to_dataframe(self) -> pd.DataFrame:
        """
        Export IP design to pandas DataFrame, validated by DataFrameDesignIP.

        Columns: sample_name, ip, control, ip_r1, ip_r2, control_r1, control_r2, plus metadata fields.
        """
        import pandas as pd

        rows: list[dict[str, Any]] = []
        for exp, md in zip(self.experiments, self.metadata):
            base = exp.ip.sample_id
            row: dict[str, Any] = {
                "sample_id": base,
                # IP/control labels
                "ip": exp.ip_performed,
                "control": exp.control_performed if exp.has_control else None,
                "uid": f"{base}_{exp.ip_performed}",
                # IP reads
                "r1": exp.ip.r1.path,
                "r2": exp.ip.r2.path if exp.ip.r2 else None,
                # Control reads
                "r1_control": exp.control.r1.path if exp.control else None,
                "r2_control": exp.control.r2.path
                if exp.control and exp.control.r2
                else None,
            }
            row.update(md.model_dump(exclude_none=True, mode="json"))
            rows.append(row)

        df = pd.DataFrame(rows).sort_values("sample_id").set_index("uid")

        # Define column order: critical columns first (assay, sample info, IP info, files), then other metadata
        core_cols = [
            "assay",
            "sample_id",
            "ip",
            "control",
            "r1",
            "r2",
            "r1_control",
            "r2_control",
        ]
        metadata_cols = [col for col in df.columns if col not in core_cols]
        ordered_cols = core_cols + sorted(metadata_cols)
        df = df[[col for col in ordered_cols if col in df.columns]]

        return DataFrame[DesignDataFrame](df)

    @classmethod
    def from_dataframe(
        cls, assay: Assay, df: pd.DataFrame, validate_deseq2: bool = False, assay_for_validation: Assay | None = None, **exp_kwargs: Any
    ) -> FastqCollectionForIP:
        """
        Build a FastqCollectionForIP from a DataFrame.

        Expects columns: sample_id, r1, r2 (optional), r1_control (optional),
        r2_control (optional), ip, control (optional), plus any metadata fields.
        
        Args:
            assay: The assay type
            df: DataFrame with sample metadata
            validate_deseq2: If True, require deseq2 field to be non-null (for RNA assays)
            assay_for_validation: Assay type to check in validation context
            **exp_kwargs: Additional kwargs for ExperimentIP
        """
        df = DesignDataFrame.validate(df)
        experiments: list[ExperimentIP] = []
        metadata: list[Metadata] = []
        metadata_fields = set(Metadata.model_fields.keys())
        
        # Use provided assay_for_validation or fall back to assay
        validation_assay = assay_for_validation or assay

        for rec in df.to_dict(orient="records"):
            # Build IP FastqSetIP
            r2_path = rec.get("r2")
            ip_set = FastqSetIP(
                sample_id=rec["sample_id"],
                r1=FastqFileIP(path=rec["r1"]),
                r2=FastqFileIP(path=r2_path) if pd.notna(r2_path) else None,
                ip=rec.get("ip"),
            )

            # Build control FastqSetIP if present
            r1_ctrl = rec.get("r1_control")
            r2_ctrl = rec.get("r2_control")
            control_set = None
            if pd.notna(r1_ctrl):
                control_set = FastqSetIP(
                    sample_id=rec["sample_id"],
                    r1=FastqFileIP(path=r1_ctrl),
                    r2=FastqFileIP(path=r2_ctrl) if pd.notna(r2_ctrl) else None,
                    ip=rec.get("control"),
                )

            experiments.append(
                ExperimentIP(ip=ip_set, control=control_set, **exp_kwargs)
            )

            # Collect metadata with validation context
            meta_fields = {k: rec.get(k) for k in metadata_fields if k in rec}
            metadata.append(Metadata.model_validate(meta_fields, context={'validate_deseq2': validate_deseq2, 'assay': validation_assay}))
        return cls(assay=assay, experiments=experiments, metadata=metadata)

    def symlink_fastq_files(self, output_dir):
        # Link the R1 and R2 files
        super().symlink_fastq_files(output_dir)

        # Additionally, link control files with appropriate names
        output_dir = Path(output_dir)
        df = self.to_dataframe()
        for _, row in df.iterrows():
            # paired-end control data needs both R1 and R2
            if pd.notna(row["r1_control"]) and not pd.isna(row["r2_control"]):
                symlink_path_ctrl_r1 = output_dir / f"{row['sample_id']}_{row['control']}_1.fastq.gz"
                symlink_path_ctrl_r2 = output_dir / f"{row['sample_id']}_{row['control']}_2.fastq.gz"

                if not symlink_path_ctrl_r1.exists():
                    symlink_path_ctrl_r1.symlink_to(
                        Path(row["r1_control"]).resolve().absolute()
                    )
                if not symlink_path_ctrl_r2.exists():
                    symlink_path_ctrl_r2.symlink_to(
                        Path(row["r2_control"]).resolve().absolute()
                    )
            # single-end control data needs only R1 - no _[12] suffix
            elif pd.notna(row["r1_control"]):
                symlink_path_ctrl_r1 = output_dir / f"{row['sample_id']}_{row['control']}.fastq.gz"