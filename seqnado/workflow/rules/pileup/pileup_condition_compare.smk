"""
Condition-based bigwig averaging and comparison rules for non-MCC assays.

These rules aggregate per-sample bigwigs to condition-level means and generate
pairwise subtractions. They are gated on:
  - CONFIG.assay_config.bigwigs.perform_comparisons == True
  - len(SAMPLE_GROUPINGS.get_grouping("condition").group_names) >= 2
  - CONFIG.third_party_tools.bamnado is not None

For RNA-seq, stranded bigwigs (_plus, _minus) are handled with strand wildcards.
"""

import itertools
from seqnado import Assay
from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
from seqnado.workflow.helpers.pileup import get_condition_input_bigwigs


# Gate all rules on: perform_comparisons == True AND len(condition_groups) >= 2 AND bamnado is configured
_PERFORM_COMPARISONS = (
    CONFIG.assay_config.bigwigs
    and getattr(CONFIG.assay_config.bigwigs, "perform_comparisons", False)
    and CONFIG.third_party_tools.bamnado is not None
)

_CONDITION_GROUPS = (
    SAMPLE_GROUPINGS.get_grouping("condition").group_names
    if _PERFORM_COMPARISONS
    else []
)

_PILEUP_METHODS = (
    CONFIG.assay_config.bigwigs.pileup_method
    if (CONFIG.assay_config.bigwigs and CONFIG.assay_config.bigwigs.pileup_method)
    else []
)

# Check if this is an RNA assay (which has stranded bigwigs)
_IS_RNA = ASSAY == Assay.RNA
_STRANDS = ["plus", "minus"] if _IS_RNA else [""]


if _PERFORM_COMPARISONS and len(_CONDITION_GROUPS) >= 2:

    rule make_bigwigs_aggregated:
        """
        Aggregate per-sample unscaled bigwigs to condition-level means.
        For RNA-seq, handles stranded bigwigs (plus/minus).
        """
        input:
            bigwigs=lambda wildcards: get_condition_input_bigwigs(
                wildcards,
                pileup_method=wildcards.pileup_method,
                spikein_method=None,
                output_dir=OUTPUT_DIR,
                sample_groupings=SAMPLE_GROUPINGS,
                strand=wildcards.strand if _IS_RNA else None,
            ),
        output:
            bigwig=OUTPUT_DIR + "/bigwigs/{pileup_method}/aggregated/{condition}{strand}.bigWig",
        params:
            options=str(
                CONFIG.third_party_tools.bamnado.bigwig_aggregate.command_line_arguments
            ),
        log:
            OUTPUT_DIR + "/logs/bigwig/{pileup_method}_aggregated_{condition}{strand}.log",
        benchmark:
            OUTPUT_DIR + "/.benchmark/bigwig/{pileup_method}_aggregated_{condition}{strand}.tsv",
        container:
            "docker://ghcr.io/alsmith151/bamnado:latest",
        message:
            "Aggregating bigwigs for condition {wildcards.condition} using {wildcards.pileup_method}",
        wildcard_constraints:
            pileup_method="|".join([m.value for m in _PILEUP_METHODS]) if _PILEUP_METHODS else ".*",
            condition="|".join(_CONDITION_GROUPS) if _CONDITION_GROUPS else ".*",
            strand="_plus|_minus|" if _IS_RNA else "",
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        shell:
            """
            bamnado \
            bigwig-aggregate \
            --bigwigs {input.bigwigs} \
            -o {output.bigwig} \
            -m mean \
            {params.options} \
            > {log} 2>&1
            """

    rule make_bigwigs_subtraction:
        """
        Generate pairwise subtraction bigwigs from condition-level aggregates (unscaled source).
        For RNA-seq, handles stranded bigwigs (plus/minus).
        """
        input:
            bw1=OUTPUT_DIR + "/bigwigs/{pileup_method}/aggregated/{condition1}{strand}.bigWig",
            bw2=OUTPUT_DIR + "/bigwigs/{pileup_method}/aggregated/{condition2}{strand}.bigWig",
        output:
            bigwig=OUTPUT_DIR + "/bigwigs/{pileup_method}/subtraction/{condition1}_vs_{condition2}{strand}.bigWig",
        params:
            options=str(
                CONFIG.third_party_tools.bamnado.bigwig_compare.command_line_arguments
            ),
        log:
            OUTPUT_DIR + "/logs/bigwig/{pileup_method}_subtraction_{condition1}_vs_{condition2}{strand}.log",
        benchmark:
            OUTPUT_DIR + "/.benchmark/bigwig/{pileup_method}_subtraction_{condition1}_vs_{condition2}{strand}.tsv",
        container:
            "docker://ghcr.io/alsmith151/bamnado:latest",
        message:
            "Subtracting {wildcards.condition2} from {wildcards.condition1} bigwigs using {wildcards.pileup_method}",
        wildcard_constraints:
            pileup_method="|".join([m.value for m in _PILEUP_METHODS]) if _PILEUP_METHODS else ".*",
            condition1="|".join(_CONDITION_GROUPS) if _CONDITION_GROUPS else ".*",
            condition2="|".join(_CONDITION_GROUPS) if _CONDITION_GROUPS else ".*",
            strand="_plus|_minus|" if _IS_RNA else "",
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        shell:
            """
            bamnado \
            bigwig-compare \
            --bw1 {input.bw1} \
            --bw2 {input.bw2} \
            -o {output.bigwig} \
            -c subtraction \
            {params.options} \
            > {log} 2>&1
            """

    # Spike-in normalized rules (only if has_spikein is configured)
    _HAS_SPIKEIN = getattr(CONFIG.assay_config, "has_spikein", False)
    _SPIKEIN_METHODS = (
        [s.value for s in CONFIG.assay_config.spike_in.method]
        if (_HAS_SPIKEIN and CONFIG.assay_config.spike_in and CONFIG.assay_config.spike_in.method)
        else []
    )

    if _HAS_SPIKEIN and _SPIKEIN_METHODS:

        rule make_bigwigs_spikein_aggregated:
            """
            Aggregate per-sample spike-in normalized bigwigs to condition-level means.
            For RNA-seq, handles stranded bigwigs (plus/minus).
            """
            input:
                bigwigs=lambda wildcards: get_condition_input_bigwigs(
                    wildcards,
                    pileup_method=wildcards.pileup_method,
                    spikein_method=wildcards.spikein_method,
                    output_dir=OUTPUT_DIR,
                    sample_groupings=SAMPLE_GROUPINGS,
                    strand=wildcards.strand if _IS_RNA else None,
                ),
            output:
                bigwig=OUTPUT_DIR + "/bigwigs/{pileup_method}/spikein/{spikein_method}/aggregated/{condition}{strand}.bigWig",
            params:
                options=str(
                    CONFIG.third_party_tools.bamnado.bigwig_aggregate.command_line_arguments
                ),
            log:
                OUTPUT_DIR + "/logs/bigwig/{pileup_method}_spikein_{spikein_method}_aggregated_{condition}{strand}.log",
            benchmark:
                OUTPUT_DIR + "/.benchmark/bigwig/{pileup_method}_spikein_{spikein_method}_aggregated_{condition}{strand}.tsv",
            container:
                "docker://ghcr.io/alsmith151/bamnado:latest",
            message:
                "Aggregating spike-in normalized bigwigs for condition {wildcards.condition} ({wildcards.spikein_method}) using {wildcards.pileup_method}",
            wildcard_constraints:
                pileup_method="|".join([m.value for m in _PILEUP_METHODS]) if _PILEUP_METHODS else ".*",
                spikein_method="|".join(_SPIKEIN_METHODS) if _SPIKEIN_METHODS else ".*",
                condition="|".join(_CONDITION_GROUPS) if _CONDITION_GROUPS else ".*",
                strand="_plus|_minus|" if _IS_RNA else "",
            resources:
                mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
                runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            shell:
                """
                bamnado \
                bigwig-aggregate \
                --bigwigs {input.bigwigs} \
                -o {output.bigwig} \
                -m mean \
                {params.options} \
                > {log} 2>&1
                """

        rule make_bigwigs_spikein_subtraction:
            """
            Generate pairwise subtraction bigwigs from spike-in normalized condition aggregates.
            For RNA-seq, handles stranded bigwigs (plus/minus).
            """
            input:
                bw1=OUTPUT_DIR + "/bigwigs/{pileup_method}/spikein/{spikein_method}/aggregated/{condition1}{strand}.bigWig",
                bw2=OUTPUT_DIR + "/bigwigs/{pileup_method}/spikein/{spikein_method}/aggregated/{condition2}{strand}.bigWig",
            output:
                bigwig=OUTPUT_DIR + "/bigwigs/{pileup_method}/spikein/{spikein_method}/subtraction/{condition1}_vs_{condition2}{strand}.bigWig",
            params:
                options=str(
                    CONFIG.third_party_tools.bamnado.bigwig_compare.command_line_arguments
                ),
            log:
                OUTPUT_DIR + "/logs/bigwig/{pileup_method}_spikein_{spikein_method}_subtraction_{condition1}_vs_{condition2}{strand}.log",
            benchmark:
                OUTPUT_DIR + "/.benchmark/bigwig/{pileup_method}_spikein_{spikein_method}_subtraction_{condition1}_vs_{condition2}{strand}.tsv",
            container:
                "docker://ghcr.io/alsmith151/bamnado:latest",
            message:
                "Subtracting spike-in normalized {wildcards.condition2} from {wildcards.condition1} bigwigs ({wildcards.spikein_method}) using {wildcards.pileup_method}",
            wildcard_constraints:
                pileup_method="|".join([m.value for m in _PILEUP_METHODS]) if _PILEUP_METHODS else ".*",
                spikein_method="|".join(_SPIKEIN_METHODS) if _SPIKEIN_METHODS else ".*",
                condition1="|".join(_CONDITION_GROUPS) if _CONDITION_GROUPS else ".*",
                condition2="|".join(_CONDITION_GROUPS) if _CONDITION_GROUPS else ".*",
                strand="_plus|_minus|" if _IS_RNA else "",
            resources:
                mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
                runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            shell:
                """
                bamnado \
                bigwig-compare \
                --bw1 {input.bw1} \
                --bw2 {input.bw2} \
                -o {output.bigwig} \
                -c subtraction \
                {params.options} \
                > {log} 2>&1
                """
