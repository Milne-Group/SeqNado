from seqnado.workflow.helpers.common import (
    define_memory_requested,
    define_time_requested,
    format_deeptools_options,
)
from seqnado.workflow.helpers.normalization import (
    get_scaling_factor_for_merged_group,
    get_norm_factor_spikein_for_merged_group,
)


rule make_bigwigs_deeptools_merged_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/csaw/{group}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        options=lambda wc: format_deeptools_options(
            wc,
            str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
            INPUT_FILES,
            SAMPLE_GROUPINGS,
            raw_counts=True,
        ),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/csaw/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/csaw/{group}.tsv",
    message: "Making CSAW-scaled merged bigWig with deeptools for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1
        """


rule make_bigwigs_deeptools_merged_spikein:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=False),
        options=lambda wc: format_deeptools_options(
            wc,
            str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
            INPUT_FILES,
            SAMPLE_GROUPINGS,
            raw_counts=True,
        ),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}.tsv",
    message: "Making spike-in normalized merged bigWig with deeptools for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1
        """


rule make_bigwigs_bamnado_merged_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/csaw/{group}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/csaw/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/bamnado/merged/csaw/{group}.tsv",
    message: "Making CSAW-scaled merged bigWig with bamnado for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} --scale-factor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule make_bigwigs_bamnado_merged_spikein:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=False),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}.tsv",
    message: "Making spike-in normalized merged bigWig with bamnado for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} --scale-factor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """
