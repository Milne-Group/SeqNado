"""Calculate CSAW-style scaling factors from genomic bin read counts.

Reads featureCounts output from genomic bins and calculates
library-size normalization factors for each sample within a scaling group.
Output is a TSV file with 'sample' and 'scale_factor' columns.
"""

import sys
import pandas as pd
from pathlib import Path

log_file = snakemake.log[0]  # noqa: F821

with open(log_file, "w") as log:
    def _log(msg):
        log.write(msg + "\n")
        log.flush()

    try:
        counts_file = snakemake.input.counts  # noqa: F821
        design_file = snakemake.input.design  # noqa: F821
        group = snakemake.wildcards.group  # noqa: F821
        output_file = snakemake.output.scaling_factors  # noqa: F821

        _log(f"Calculating CSAW scaling factors for group: {group}")
        _log(f"Counts file: {counts_file}")
        _log(f"Design file: {design_file}")

        # Parse featureCounts output (skip comment lines starting with #)
        counts = pd.read_csv(counts_file, sep="\t", comment="#")

        # featureCounts columns: Geneid, Chr, Start, End, Strand, Length, {bam_paths}...
        meta_cols = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
        sample_cols = [c for c in counts.columns if c not in meta_cols]

        # Extract sample names from BAM file paths (strip directory and .bam extension)
        sample_name_map = {col: Path(col).stem for col in sample_cols}
        counts_data = counts[sample_cols].rename(columns=sample_name_map)

        # Calculate library sizes (total counts per sample)
        lib_sizes = counts_data.sum(axis=0)
        _log(f"Library sizes:\n{lib_sizes.to_string()}")

        # Read design metadata
        design = pd.read_csv(design_file)

        # Reconstruct sample names (uid = sample_id + "_" + ip for IP-based assays)
        if "uid" not in design.columns:
            if "sample_id" in design.columns and "ip" in design.columns:
                design["uid"] = design["sample_id"] + "_" + design["ip"]
            elif "sample_id" in design.columns:
                design["uid"] = design["sample_id"]
            else:
                raise ValueError("Cannot determine sample names from metadata.csv: no 'uid', 'sample_id', or 'ip' columns found.")

        # Get samples for this group.
        # The same {group}_scaling_factors.tsv pattern is used for:
        #   - Individual sample scaling (group = scaling_group value)
        #   - Merged bigwig scaling (group = consensus_group value)
        # Try scaling_group first, then consensus_group.
        group_samples = []
        if "scaling_group" in design.columns:
            mask = design["scaling_group"] == group
            if mask.any():
                group_samples = design.loc[mask, "uid"].tolist()
                _log(f"Found group '{group}' in scaling_group column")

        if not group_samples and "consensus_group" in design.columns:
            mask = design["consensus_group"] == group
            if mask.any():
                group_samples = design.loc[mask, "uid"].tolist()
                _log(f"Found group '{group}' in consensus_group column")

        if not group_samples:
            _log(f"Group '{group}' not found in scaling_group or consensus_group; using all samples")
            group_samples = design["uid"].tolist()

        _log(f"Samples in group '{group}': {group_samples}")

        # Filter to samples present in the counts data
        available_samples = [s for s in group_samples if s in lib_sizes.index]
        if not available_samples:
            raise ValueError(
                f"No samples from group '{group}' found in featureCounts output. "
                f"Expected: {group_samples}, Available: {list(lib_sizes.index)}"
            )

        lib_sizes_group = lib_sizes[available_samples]
        _log(f"Library sizes for group '{group}':\n{lib_sizes_group.to_string()}")

        # Calculate scale factors: normalize to equal effective depth within group.
        # scale_factor = mean_lib_size / lib_size
        # Samples with larger libraries get smaller scale factors (scaled down),
        # samples with smaller libraries get larger scale factors (scaled up).
        mean_lib_size = lib_sizes_group.mean()
        scale_factors = mean_lib_size / lib_sizes_group

        _log(f"Scale factors for group '{group}':\n{scale_factors.to_string()}")

        # Output TSV with 'sample' and 'scale_factor' columns
        result = pd.DataFrame({
            "sample": scale_factors.index,
            "scale_factor": scale_factors.values,
        })
        result.to_csv(output_file, sep="\t", index=False)
        _log(f"Scaling factors written to {output_file}")

    except Exception as e:
        _log(f"ERROR: {e}")
        raise
