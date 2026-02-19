"""Helper functions for Snakemake workflow rules."""

from .bam import get_bam_files_for_consensus, get_bam_split
from .common import (
    define_memory_requested,
    define_time_requested,
    format_deeptools_options,
    get_alignment_input,
    get_fastq_paths,
    get_group_for_sample,
    get_scale_method,
)
from .crispr import get_cutadapt_adapter_args
from .geo import get_files_for_symlink, get_symlinked_files
from .mcc import (
    extract_viewpoints,
    get_mcc_bam_files_for_merge,
    get_n_cis_scaling_factor,
    identify_extracted_bam_files,
    redefine_viewpoints,
    viewpoint_to_grouped_viewpoint,
)
from .multiomics import get_assay_all_inputs
from .normalization import get_norm_factor_spikein, get_scaling_factor
from .peaks import (
    correct_macs_options,
    get_control_file,
    get_lanceotron_call_peaks_threshold,
)
from .qc import format_frip_enrichment_options, format_qualimap_options
from .quant import get_bams_to_count

__all__ = [
    # Common
    "define_memory_requested",
    "define_time_requested",
    "get_group_for_sample",
    "get_scale_method",
    "get_fastq_paths",
    "get_alignment_input",
    "format_deeptools_options",
    # Peaks
    "get_control_file",
    "correct_macs_options",
    "get_lanceotron_call_peaks_threshold",
    # GEO
    "get_files_for_symlink",
    "get_symlinked_files",
    # MCC
    "get_n_cis_scaling_factor",
    "get_mcc_bam_files_for_merge",
    "identify_extracted_bam_files",
    "redefine_viewpoints",
    "extract_viewpoints",
    "viewpoint_to_grouped_viewpoint",
    # CRISPR
    "get_cutadapt_adapter_args",
    # QC
    "format_qualimap_options",
    "format_frip_enrichment_options",
    # Normalization
    "get_norm_factor_spikein",
    "get_scaling_factor",
    # BAM
    "get_bam_files_for_consensus",
    "get_bam_split",
    # Multiomics
    "get_assay_all_inputs",
    # Quant
    "get_bams_to_count",
]
