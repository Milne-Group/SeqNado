"""Helper functions for quantification workflows."""

from seqnado.config import SeqnadoConfig
from seqnado.inputs.fastq import FastqCollection, FastqCollectionForIP
from seqnado.inputs.grouping import SampleGroup, SampleGroupings
from pathlib import Path
from seqnado.config.third_party_tools import CommandLineArguments


def get_bams_to_count(
    config: SeqnadoConfig, sample_names: list[str], output_dir: str
) -> list[str]:
    """
    Get list of BAM files to use for counting.

    Returns spike-in filtered BAMs if spike-in is configured, otherwise returns
    regular aligned BAMs.

    Args:
        config (SeqNadoConfig): The SeqNado configuration object.
        sample_names (list): List of sample names to get BAM files for.
        output_dir (str): The base output directory.

    Returns:
        list: List of BAM file paths to count.
    """
    if config.assay_config.has_spikein:
        return [
            output_dir + f"/aligned/spikein/filtered/{sample}.bam"
            for sample in sample_names
        ]
    return [output_dir + f"/aligned/{sample}.bam" for sample in sample_names]


def get_bams_to_count_grouped(
    config: SeqnadoConfig,
    sample_groupings: SampleGroupings,
    output_dir: str,
    group_name: str,
    grouping: str = "consensus",
) -> list[str]:
    """
    Get list of BAM files to use for counting.

    Returns spike-in filtered BAMs if spike-in is configured, otherwise returns
    regular aligned BAMs.

    Args:
        config (SeqnadoConfig): The SeqNado configuration object.
        sample_group (SampleGroup): The sample grouping object.
        output_dir (str): The base output directory.

    Returns:
        list: List of BAM file paths to count.
    """

    try:
        sample_group = sample_groupings.get_grouping(grouping).get_group(group_name)
        if config.assay_config.has_spikein:
            return [
                output_dir + f"/aligned/spikein/filtered/{sample}.bam"
                for sample in sample_group.samples
            ]
        return [
            output_dir + f"/aligned/{sample}.bam" for sample in sample_group.samples
        ]
    except KeyError:
        return []  # Return empty list if group or grouping not found


def correct_featurecounts_options_for_peaks(
    options: CommandLineArguments,
    input_files: FastqCollection | FastqCollectionForIP,
    sample_groupings: SampleGroupings,
    group_name: str,
) -> str:
    """
    Correct featureCounts options for peak counting.

    Args:
        options (CommandLineArguments): The original featureCounts options.
    Returns:
        str: The corrected featureCounts options for peak counting.
    """

    # We don't have a type so remove any options that specify a type:
    options = options.add_exclude("-t", "--type")

    # Check that all samples are paired end, if not remove the -p option:
    samples = (
            sample_groupings.get_grouping("consensus")
            .get_group(group_name)
            .samples
        )
    
    if not all(
        input_files.is_paired_end(sample_name) for sample_name in samples
    ):
        options = options.add_exclude("-p", "--pairedEnd")
    else:
        options = options.add_include("-p")

    return str(options)
