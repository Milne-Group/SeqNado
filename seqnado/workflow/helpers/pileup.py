from snakemake.io import expand
from seqnado.inputs import SampleGroupings
from seqnado import Assay


def get_condition_input_bigwigs(wildcards, pileup_method, spikein_method=None, output_dir: str | None = None, sample_groupings: SampleGroupings | None = None, assay: Assay | None = None, strand: str | None = None):
    """
    Return per-sample bigwig paths for a condition group.

    For RNA-seq, handles stranded bigwigs (_plus, _minus) if strand is specified.

    Args:
        wildcards: Snakemake wildcards (must have 'condition')
        pileup_method: The pileup method (e.g., "bamnado", "deeptools")
        spikein_method: Optional spike-in method; if provided, returns spikein-normalized paths
        output_dir: Output directory path
        sample_groupings: Sample groupings object
        assay: Assay type (needed for strand detection)
        strand: Strand specifier ("plus" or "minus") for RNA-seq; None for non-RNA or unstranded

    Returns:
        List of bigwig file paths for all samples in the condition group
    """
    samples = sample_groupings.get_grouping("condition").get_samples(wildcards.condition)

    # Build suffix for stranded bigwigs (RNA-seq only)
    suffix = f"_{strand}" if strand else ""

    if spikein_method:
        return expand(
            output_dir + f"/bigwigs/{pileup_method}/spikein/{spikein_method}/{{sample}}{suffix}.bigWig",
            sample=samples,
        )
    return expand(
        output_dir + f"/bigwigs/{pileup_method}/unscaled/{{sample}}{suffix}.bigWig",
        sample=samples,
    )