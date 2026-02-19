"""Common helper functions used across workflow rules."""

from typing import Dict, List, Union

from seqnado import DataScalingTechnique


def define_memory_requested(
    attempts: int = 1, initial_value: int = 1, scale: float = 1
) -> str:
    """
    Define the memory requested for the job,
    returns a string like "4G" avoiding decimals for qualimap.
    """
    mem_value = int(initial_value) * 2 ** (int(attempts) - 1)
    mem_value = int(mem_value * float(scale))
    return f"{mem_value}G"


def define_time_requested(
    attempts: int = 1, initial_value: int = 1, scale: float = 1
) -> str:
    """
    Define the time requested for the job.

    Base time is 1 hour.
    """
    time = int(initial_value) * 2 ** (int(attempts) - 1)
    time = time * float(scale)
    return f"{time}h"


def get_group_for_sample(
    wildcards, design, strip: str = ""
):
    """
    Get the group for a sample from the design.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
        design: FastqCollection or FastqCollectionForIP object.
        strip: Optional string to strip from sample name.

    Returns:
        str: The group name for the sample.

    Raises:
        KeyError: If sample not found in normalisation groups.
    """
    from seqnado.inputs import SampleGroups

    df = design.to_dataframe()
    scaling_groups = SampleGroups.from_dataframe(df, subset_column="scaling_group")

    sample_name = wildcards.sample.strip(strip)
    try:
        group = scaling_groups.sample_to_group[sample_name]
        return group
    except KeyError:
        raise KeyError(f"Sample {wildcards.sample} not found in normalisation groups.")


def get_scale_method(config: Dict) -> List[str]:
    """
    Returns the scale method based on the config.
    """

    method = [DataScalingTechnique.UNSCALED]

    if config.get("spikein"):
        method.append(DataScalingTechnique.SPIKEIN)
    elif config.get("scale"):
        method.append(DataScalingTechnique.CSAW)
    return [m.value for m in method]


def get_fastq_paths(
    wildcards, output_dir: str, directory: str, paired: bool = True
) -> Union[dict, str]:
    """
    Get fastq file paths for a given directory and pairing.

    Args:
        wildcards: Snakemake wildcards object
        output_dir: Base output directory
        directory: Subdirectory ('fastqs' for raw, 'trimmed' for trimmed)
        paired: Whether to return paired-end (dict) or single-end (str) paths

    Returns:
        dict with 'fq1' and 'fq2' keys for paired-end, or str for single-end
    """
    if paired:
        return {
            "fq1": f"{output_dir}/{directory}/{wildcards.sample}_1.fastq.gz",
            "fq2": f"{output_dir}/{directory}/{wildcards.sample}_2.fastq.gz",
        }
    else:
        return f"{output_dir}/{directory}/{wildcards.sample}.fastq.gz"


def get_alignment_input(
    wildcards, output_dir: str, config, paired: bool = True
) -> Union[dict, str]:
    """Return input paths based on trimming config.

    Args:
        wildcards: Snakemake wildcards object
        output_dir: Base output directory
        config: Configuration object
        paired: Whether to return paired-end (dict) or single-end (str) paths

    Returns:
        dict with 'fq1' and 'fq2' keys for paired-end, or str for single-end
    """

    trimming_enabled = getattr(config.qc, "trim_fastq", True)
    directory = "trimmed" if trimming_enabled else "fastqs"
    return get_fastq_paths(wildcards, output_dir, directory, paired=paired)


def format_deeptools_options(wildcards, options, input_files, sample_groupings=None, raw_counts=False):
    """
    Format the command line options for deeptools based on the input files and parameters.

    Mainly this removes the extend reads option if single ended.

    Args:
        wildcards: Snakemake wildcards object
        options: CommandLineArguments or string with options
        input_files: FastqCollection object to check if samples are paired-end
        sample_groupings: Optional SampleGroupings object for grouped samples

    Returns:
        Formatted options string
    """
    from seqnado.config.third_party_tools import CommandLineArguments

    # Convert string to CommandLineArguments if needed
    if isinstance(options, str):
        options = CommandLineArguments(value=options)

    try:
        if hasattr(wildcards, "group"):
            if sample_groupings is None:
                raise ValueError(
                    "sample_groupings required when wildcards has 'group' attribute"
                )
            sample = (
                sample_groupings.get_grouping("consensus")
                .get_group(wildcards.group)
                .samples
            )
            if not all(
                input_files.is_paired_end(sample_name) for sample_name in sample
            ):
                options = CommandLineArguments(
                    value=str(options),
                    exclude={"--extendReads", "-e", "--samFlagInclude 3"},
                )
        else:
            search_term = f"{wildcards.sample}"
            is_paired = input_files.is_paired_end(search_term)
            if not is_paired:
                options = CommandLineArguments(
                    value=str(options),
                    exclude={"--extendReads", "-e", "--samFlagInclude 3"},
                )
    except KeyError:
        pass

    if raw_counts:
        options = CommandLineArguments(
            value=str(options),
            exclude={"--scaleFactor", "--normalizeUsing", "--effectiveGenomeSize"},
        )

    return str(options)
