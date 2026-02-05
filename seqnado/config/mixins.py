from pathlib import Path
from typing import Optional, Any

from pydantic import BaseModel, field_validator, model_validator, ValidationInfo
from seqnado import Assay


# -----------------------------------------------------------------------------
# Path validation helpers
# -----------------------------------------------------------------------------

class PathValidatorMixin:
    @staticmethod
    def validate_path_exists(
        v: Path | str | None,
        field_name: str = "path",
        info: ValidationInfo | None = None,
    ) -> Path | str | None:
        if v is None:
            return v

        if info and info.context and info.context.get("skip_path_validation", False):
            return v

        path = Path(v) if isinstance(v, str) else v
        if not path.exists():
            raise ValueError(f"{field_name} {v} does not exist.")

        return v

    @staticmethod
    def validate_required_when(
        v: str | None,
        condition: bool,
        field_name: str,
    ) -> str | None:
        if condition and not v:
            raise ValueError(f"{field_name} must be provided when condition is met.")
        return v


# -----------------------------------------------------------------------------
# Common computed toggles
# -----------------------------------------------------------------------------

class CommonComputedFieldsMixin(BaseModel):
    create_bigwigs: Optional[bool] = None
    plot_with_plotnado: Optional[bool] = None
    create_dataset: Optional[bool] = None
    create_ucsc_hub: Optional[bool] = None
    has_spikein: Optional[bool] = None
    create_quantification_files: bool = False

    create_geo_submission_files: bool = False

    @field_validator("create_bigwigs", mode="before")
    def compute_create_bigwigs(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        if issubclass(cls, (MethylationMixin, SNPCallingMixin)):
            return False

        return info.data.get("bigwigs") is not None

    @field_validator("plot_with_plotnado", mode="before")
    def compute_plot_with_plotnado(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("plotting") is not None

    @field_validator("create_dataset", mode="before")
    def compute_create_dataset(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("dataset_for_ml") is not None

    @field_validator("create_ucsc_hub", mode="before")
    def compute_create_ucsc_hub(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("ucsc_hub") is not None

    @field_validator("has_spikein", mode="before")
    def compute_has_spikein(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("spikein") is not None

    @field_validator("create_geo_submission_files", mode="before")
    def validate_geo_submission_files(
        cls,
        v: Any,
    ) -> bool:
        return bool(v)
    
    @field_validator("create_quantification_files", mode="before")
    def validate_create_quantification_files(
        cls,
        v: Any,
        info: ValidationInfo,
    ) -> bool:
        
        if v is not None:
            return bool(v)

        # Check to see if we need to create quantification files
        # For RNA-seq this is always True
        # For other assays, this will be left to False unless explicitly set to True
        if info.data.get('assay') == Assay.RNA:
            return True

       


    @model_validator(mode="after")
    def compute_missing_computed_fields(cls, model: "CommonComputedFieldsMixin") -> "CommonComputedFieldsMixin":
        # Populate any computed booleans left as None based on other config fields
        if model.create_bigwigs is None:
            if isinstance(model, (MethylationMixin, SNPCallingMixin)):
                model.create_bigwigs = False
            else:
                model.create_bigwigs = getattr(model, "bigwigs", None) is not None

        if model.plot_with_plotnado is None:
            model.plot_with_plotnado = getattr(model, "plotting", None) is not None

        if model.create_dataset is None:
            model.create_dataset = getattr(model, "dataset_for_ml", None) is not None

        if model.create_ucsc_hub is None:
            model.create_ucsc_hub = getattr(model, "ucsc_hub", None) is not None

        if model.has_spikein is None:
            model.has_spikein = getattr(model, "spikein", None) is not None

        return model


# -----------------------------------------------------------------------------
# Assay-specific mixins
# -----------------------------------------------------------------------------

class PeakCallingMixin(BaseModel):
    call_peaks: Optional[bool] = None

    @field_validator("call_peaks", mode="before")
    def compute_call_peaks(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("peak_calling") is not None

    @model_validator(mode="after")
    def compute_missing_call_peaks(cls, model: "PeakCallingMixin") -> "PeakCallingMixin":
        if model.call_peaks is None:
            model.call_peaks = getattr(model, "peak_calling", None) is not None
        return model


class SNPCallingMixin(BaseModel):
    call_snps: Optional[bool] = None

    @field_validator("call_snps", mode="before")
    def compute_call_snps(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("snp_calling") is not None

    @model_validator(mode="after")
    def compute_missing_call_snps(cls, model: "SNPCallingMixin") -> "SNPCallingMixin":
        if model.call_snps is None:
            model.call_snps = getattr(model, "snp_calling", None) is not None
        return model


class MethylationMixin(BaseModel):
    call_methylation: Optional[bool] = None

    @field_validator("call_methylation", mode="before")
    def compute_call_methylation(
        cls,
        v: Optional[bool],
        info: ValidationInfo,
    ) -> Optional[bool]:
        if v is not None:
            return v

        return info.data.get("methylation") is not None

    @model_validator(mode="after")
    def compute_missing_call_methylation(cls, model: "MethylationMixin") -> "MethylationMixin":
        if model.call_methylation is None:
            model.call_methylation = getattr(model, "methylation", None) is not None
        return model
