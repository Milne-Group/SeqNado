"""Helper functions for QC workflows."""

from seqnado import Assay
from seqnado.config.third_party_tools import CommandLineArguments


def format_qualimap_options(wildcards, options: CommandLineArguments, INPUT_FILES, ASSAY) -> str:
    """
    Format the command line options for qualimap based on the input files and parameters.

    Mainly this removes the paired-end options if single-ended and also adds the correct
    options depending on assay type.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
        options (CommandLineArguments): The qualimap command line arguments to format.
        INPUT_FILES: The input files collection.
        ASSAY: The assay type.

    Returns:
        str: Formatted command line options string.
    """
    is_paired = INPUT_FILES.is_paired_end(wildcards.sample)
    match (is_paired, ASSAY):
        case (True, Assay.RNA):
            options = CommandLineArguments(value=options.value, include={"--paired", "--sorted"})
        case (True, _):
            options = CommandLineArguments(
                value=options.value, include={"--collect-overlap-pairs"}
            )
        case (False, Assay.RNA):
            options = CommandLineArguments(value=options.value, exclude={"--paired", "--sorted"})
        case (False, _):
            options = CommandLineArguments(
                value=options.value, exclude={"--collect-overlap-pairs"}
            )
    
    # Also remove any memory specification as this is handled elsewhere
    options = CommandLineArguments(value=options.value, exclude={"--java-mem-size"})
    return str(options)


def format_frip_enrichment_options(wildcards, options: CommandLineArguments, INPUT_FILES):
    """
    Format FRIP enrichment options based on whether input is paired-end.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
        options (CommandLineArguments): The command line arguments to format.
        INPUT_FILES: The input files collection.

    Returns:
        str: Formatted command line options string.
    """
    is_paired = INPUT_FILES.is_paired_end(wildcards.sample)
    if is_paired:
        options = CommandLineArguments(value=options.value, include={"--paired"})
    else:
        options = CommandLineArguments(value=options.value, exclude={"--paired"})
    return str(options)
