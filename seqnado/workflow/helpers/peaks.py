"""Helper functions for peak calling workflows."""

import re

from seqnado import DataScalingTechnique, FileType
from seqnado.config.third_party_tools import CommandLineArguments
from seqnado.inputs import FastqCollection, FastqCollectionForIP


def get_control_file(wildcards, file_type: FileType, INPUT_FILES, OUTPUT_DIR):
    """
    Get the control file for a given treatment sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample' and 'treatment'.
        file_type (FileType): The type of file to retrieve (e.g., FileType.BAM, FileType.BIGWIG, FileType.TAG_DIRECTORY).
        INPUT_FILES: The input files collection (FastqCollection or FastqCollectionForIP).
        OUTPUT_DIR: The output directory path.

    Returns:
        str or None: The path to the control file if it exists, otherwise None.
    """

    if isinstance(INPUT_FILES, FastqCollection):
        return "NO_CONTROL_PRESENT_FOR_FILE"  # Dummy value to prevent the rule from ever being able to run if not an IP
    elif isinstance(INPUT_FILES, FastqCollectionForIP):
        search_term = wildcards.sample_id
        control_name = INPUT_FILES.get_control_performed(search_term)
        if not control_name:
            return "NO_CONTROL_PRESENT_FOR_FILE"

    # Match control with correct filetype
    match file_type:
        case FileType.BAM:
            return OUTPUT_DIR + f"/aligned/{control_name}.bam"
        case FileType.BIGWIG:
            return (
                OUTPUT_DIR
                + f"/bigwigs/deeptools/{DataScalingTechnique.PER_SAMPLE.value}/{control_name}.bigWig"
            )
        case FileType.TAG_DIRECTORY:
            return OUTPUT_DIR + f"/tag_dirs/{control_name}"
        case _:
            raise ValueError(
                f"Unsupported file type '{file_type}' for control file retrieval"
            )


def correct_macs_options(
    wildcards, options: CommandLineArguments, INPUT_FILES, SAMPLE_GROUPINGS
):
    """
    Correct MACS options based on whether the input is paired-end or single-end.

    Args:
        wildcards: Snakemake wildcards object.
        options (CommandLineArguments): The MACS command line arguments to correct.
        INPUT_FILES: The input files collection.
        SAMPLE_GROUPINGS: The sample groupings object.

    Returns:
        CommandLineArguments: Corrected command line arguments.
    """
    # Correct the options based on whether the input is paired or not
    if hasattr(wildcards, "group"):
        sample = (
            SAMPLE_GROUPINGS.get_grouping("consensus")
            .get_group(wildcards.group)
            .samples
        )
        if not all(INPUT_FILES.is_paired_end(sample_name) for sample_name in sample):
            options = CommandLineArguments(value=options.value, exclude={"-f"})
    else:
        sample_id = wildcards.sample_id
        is_paired = INPUT_FILES.is_paired_end(sample_id)
        if not is_paired:
            options = CommandLineArguments(value=options.value, exclude={"-f"})
    return options


def get_lanceotron_call_peaks_threshold(wildcards, CONFIG):
    """
    Extract the Lanceotron peak calling threshold from config options.

    Args:
        wildcards: Snakemake wildcards object.
        CONFIG: The configuration object.

    Returns:
        str: The threshold value (default: "0.5").
    """
    options = str(CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments)
    threshold_pattern = re.compile(r"\-c\s*(\d+\.?\d*)")
    match = threshold_pattern.search(options)
    if match:
        return match.group(1)
    else:
        return "0.5"  # Default threshold
