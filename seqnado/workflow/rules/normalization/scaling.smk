from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


rule feature_counts_genome:
    input:
        bam=expand(OUTPUT_DIR + "/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand(OUTPUT_DIR + "/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=OUTPUT_DIR + "/resources/genomic_bins.saf",
    output:
        counts=OUTPUT_DIR + "/resources/binned_counts/read_counts.tsv",
    params:
        options=lambda wc: str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
        paired="-p --countReadPairs" if INPUT_FILES.is_paired_end(SAMPLE_NAMES[0]) else "",
    threads:
        CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/resources/binned_counts/featurecounts.log"
    benchmark:
        OUTPUT_DIR + "/.benchmark/resources/binned_counts/featurecounts.tsv"
    message:
        "Calculating feature counts for all samples on genomic bins"
    shell: """
    featureCounts \
    -F SAF \
    -a {input.annotation} \
    -T {threads} \
    --donotsort \
    {params.paired} \
    {params.options} \
    -o {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """


rule calculate_csaw_scaling_factors:
    input:
        counts=OUTPUT_DIR + "/resources/binned_counts/read_counts.tsv",
        design=OUTPUT_DIR + "/metadata.csv",
    output:
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/normalization/csaw/{group}_scaling_factors.log"
    benchmark:
        OUTPUT_DIR + "/.benchmark/normalization/csaw/{group}_scaling_factors.tsv"
    message:
        "Calculating CSAW scaling factors for group {wildcards.group}"
    script:
        "../../scripts/calculate_csaw_scaling_factors.py"
