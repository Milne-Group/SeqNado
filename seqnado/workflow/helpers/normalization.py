"""Helper functions for normalization workflows."""

import json
import pandas as pd


def get_norm_factor_spikein(wildcards, OUTPUT_DIR, CONFIG, negative=False):
    """
    Get normalization factor from spike-in data.

    Args:
        wildcards: Snakemake wildcards object containing 'sample' and 'spikein_method'.
        OUTPUT_DIR: The output directory path.
        CONFIG: The configuration object.
        negative (bool): If True, return negative scaling factor. Default False.

    Returns:
        float: The normalization scaling factor.
    """
    method = wildcards.spikein_method
    norm_file = OUTPUT_DIR + f"/resources/{method}/normalisation_factors.json"

    with open(norm_file, "r") as f:
        norm_factors = json.load(f)

    scale_factor = float(norm_factors.get(wildcards.sample, 1.0))

    if negative:
        scale_factor = -scale_factor

    return scale_factor


def get_scaling_factor(wildcards, scaling_file):
    """
    Get scaling factor from a TSV file for a specific sample.

    Args:
        wildcards: Snakemake wildcards object containing 'sample'.
        scaling_file: Path to the TSV file containing scaling factors.

    Returns:
        float: The scaling factor for the sample.
    """
    df = pd.read_csv(scaling_file, sep="\t")
    scale = df.loc[df["sample"] == wildcards.sample, "scale_factor"].values[0]

    return float(scale)


def get_scaling_factor_for_merged_group(wildcards, sample_groupings, output_dir):
    """
    Average CSAW scaling factors for all samples in a merged consensus group.

    Args:
        wildcards: Snakemake wildcards object containing 'group'.
        sample_groupings: SampleGroupings object with consensus grouping.
        output_dir: The output directory path.

    Returns:
        float: The averaged scaling factor for the group.
    """
    samples = sample_groupings.get_grouping("consensus").get_group(wildcards.group).samples
    tsv_path = f"{output_dir}/resources/{wildcards.group}_scaling_factors.tsv"
    df = pd.read_csv(tsv_path, sep="\t")

    factors = [
        float(df.loc[df["sample"] == s, "scale_factor"].values[0])
        for s in samples
        if not df.loc[df["sample"] == s, "scale_factor"].empty
    ]

    return sum(factors) / len(factors) if factors else 1.0


def get_norm_factor_spikein_for_merged_group(wildcards, sample_groupings, output_dir, negative=False):
    """
    Average spike-in normalisation factors for all samples in a merged consensus group.

    Args:
        wildcards: Snakemake wildcards object containing 'group' and 'spikein_method'.
        sample_groupings: SampleGroupings object with consensus grouping.
        output_dir: The output directory path.
        negative (bool): If True, return negative averaged factor. Default False.

    Returns:
        float: The averaged (and optionally negated) normalisation factor.
    """
    samples = sample_groupings.get_grouping("consensus").get_group(wildcards.group).samples
    method = wildcards.spikein_method
    norm_file = f"{output_dir}/resources/{method}/normalisation_factors.json"

    with open(norm_file, "r") as f:
        norm_factors = json.load(f)

    factors = [float(norm_factors.get(s, 1.0)) for s in samples]
    avg = sum(factors) / len(factors) if factors else 1.0

    return -avg if negative else avg
