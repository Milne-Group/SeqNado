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
