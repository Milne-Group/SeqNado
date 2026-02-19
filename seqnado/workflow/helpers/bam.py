"""Helper functions for BAM file operations."""


def get_bam_files_for_consensus(wildcards, SAMPLE_GROUPINGS, OUTPUT_DIR):
    """
    Get BAM files for merging based on sample names in consensus group.

    Args:
        wildcards: Snakemake wildcards object containing 'group'.
        SAMPLE_GROUPINGS: The sample groupings object.
        OUTPUT_DIR: The output directory path.

    Returns:
        list: List of BAM file paths to merge.
    """
    groups = SAMPLE_GROUPINGS.get_grouping("consensus").get_group(wildcards.group)
    sample_names = groups.samples
    bam_files = [OUTPUT_DIR + f"/aligned/{sample}.bam" for sample in sample_names]
    return bam_files


def get_bam_split(wildcards, checkpoints):
    """
    Get split BAM file from methylation checkpoint.

    Args:
        wildcards: Snakemake wildcards object.
        checkpoints: Snakemake checkpoints object.

    Returns:
        Path to the split BAM file.
    """
    checkpoint_output = checkpoints.methylation_bam_splits.get(
        sample=wildcards.sample, genome=wildcards.genome
    ).output
    return checkpoint_output.bam
