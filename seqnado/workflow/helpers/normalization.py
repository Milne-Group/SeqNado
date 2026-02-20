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
    Compute the CSAW scaling factor for a merged consensus-group bigwig.

    The merged BAM is the physical concatenation of all per-sample BAMs, so its
    total library size is sum(L_i).  The per-sample scale factors are
    s_i = mean_L / L_i, which means L_i = mean_L / s_i.  The correct factor for
    the merged BAM is:

        scale_merged = mean_L / sum(L_i) = 1 / sum(1/s_i)

    This ensures the merged bigwig is on the same scale as individual normalized
    bigwigs (i.e. represents the average replicate, not n times the signal).
    Note: arithmetic mean of s_i is incorrect and inflates the merged signal.

    Args:
        wildcards: Snakemake wildcards object containing 'group'.
        sample_groupings: SampleGroupings object with consensus grouping.
        output_dir: The output directory path.

    Returns:
        float: The scaling factor for the merged group BAM.
    """
    samples = sample_groupings.get_grouping("consensus").get_group(wildcards.group).samples
    tsv_path = f"{output_dir}/resources/{wildcards.group}_scaling_factors.tsv"
    df = pd.read_csv(tsv_path, sep="\t")

    factors = [
        float(df.loc[df["sample"] == s, "scale_factor"].values[0])
        for s in samples
        if not df.loc[df["sample"] == s, "scale_factor"].empty
    ]

    return 1.0 / sum(1.0 / f for f in factors) if factors else 1.0


def get_norm_factor_spikein_for_merged_group(wildcards, sample_groupings, output_dir, negative=False):
    """
    Compute a spike-in normalisation factor for a merged consensus-group bigwig.

    The merged BAM is the physical concatenation of all per-sample reference-genome
    BAMs (spike-in reads have already been split out).  To match a "merge-then-split"
    calculation the factor must be derived from the *pooled* spike-in counts rather
    than from the arithmetic mean of per-sample factors (which gives incorrect results
    when samples differ in sequencing depth).

    Method-specific behaviour
    -------------------------
    orlando
        Exact pooled factor: 1e6 / sum(spikein_reads_i).
        Reads the raw bam_split stats for each sample.

    with_input
        True pooled formula, equivalent to applying the with_input calculation
        directly to the merged BAM.  The normalisation table written by
        calculate_spikein_norm_factors.py already contains per-sample IP and
        paired-input counts, so we sum them and apply:
            factor = sum(S_ctrl_i) * 1e7 / (sum(S_ip_i) * sum(R_ctrl_i))
        This is mathematically identical to "merge-then-calculate".

    deseq2 / edgeR
        Factors are derived from statistical models fitted to count matrices rather
        than directly from spike-in read totals.  Arithmetic mean is used as a
        fallback; re-running the model on pooled data would be required for a fully
        correct answer.

    Args:
        wildcards: Snakemake wildcards object containing 'group' and 'spikein_method'.
        sample_groupings: SampleGroupings object with consensus grouping.
        output_dir: The output directory path.
        negative (bool): If True, negate the returned factor. Default False.

    Returns:
        float: The pooled (or weighted) normalisation factor for the group.
    """
    samples = sample_groupings.get_grouping("consensus").get_group(wildcards.group).samples
    method = wildcards.spikein_method
    norm_file = f"{output_dir}/resources/{method}/normalisation_factors.json"

    with open(norm_file, "r") as f:
        norm_factors = json.load(f)

    if method == "orlando":
        # Sum spike-in reads across all samples and apply the orlando formula.
        total_spikein = 0
        for sample in samples:
            stats_path = f"{output_dir}/aligned/spikein/{sample}_stats.tsv"
            df = pd.read_csv(stats_path, sep="\t")
            total_spikein += int(df["spikein_reads"].values[0])
        factor = (1e6 / total_spikein) if total_spikein > 0 else 1.0

    elif method == "with_input":
        # True pooled formula: sum S_ctrl, S_ip, R_ctrl across group samples and
        # apply the same expression used per-sample:
        #   factor = sum(S_ctrl_i) * 1e7 / (sum(S_ip_i) * sum(R_ctrl_i))
        norm_table_path = f"{output_dir}/resources/with_input/normalisation_factors.tsv"
        norm_table = pd.read_csv(norm_table_path, sep="\t")
        group_rows = norm_table[norm_table["sample"].isin(samples)]
        total_s_ctrl = group_rows["spikein_reads_control"].sum()
        total_s_ip = group_rows["spikein_reads_ip"].sum()
        total_r_ctrl = group_rows["reference_reads_control"].sum()
        factor = (
            (total_s_ctrl * 1e7) / (total_s_ip * total_r_ctrl)
            if total_s_ip > 0 and total_r_ctrl > 0
            else 1.0
        )

    else:
        # deseq2, edgeR: arithmetic mean fallback.
        factors = [float(norm_factors.get(s, 1.0)) for s in samples]
        factor = sum(factors) / len(factors) if factors else 1.0

    return -factor if negative else factor
