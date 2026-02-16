"""Helper functions for quantification workflows."""
from seqnado.config import SeqnadoConfig
from seqnado.inputs.grouping import SampleGroup, SampleGroupings
from pathlib import Path


def get_bams_to_count(config: SeqnadoConfig, sample_names: list[str], output_dir: str) -> list[str]:
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



def get_bams_to_count_grouped(config: SeqnadoConfig, sample_groupings: SampleGroupings, output_dir: str, group_name: str, grouping: str = "consensus") -> list[str]:
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
        return [output_dir + f"/aligned/{sample}.bam" for sample in sample_group.samples]
    except KeyError:
        return []  # Return empty list if group or grouping not found
