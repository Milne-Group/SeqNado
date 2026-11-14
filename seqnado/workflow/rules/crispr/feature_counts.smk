rule feature_counts:
    input:
        bam=expand(OUTPUT_DIR + "/aligned/{sample}.bam", sample=SAMPLE_NAMES),
        bai=expand(OUTPUT_DIR + "/aligned/{sample}.bam.bai", sample=SAMPLE_NAMES),
        annotation=str(CONFIG.genome.gtf),
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
        paired= "-p --countReadPairs" if INPUT_FILES.is_paired_end(SAMPLE_NAMES[0]) else "",
        annotation_format= "SAF" if CONFIG.genome.gtf.suffix == ".saf" else "GTF",
    threads: CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/readcounts/featurecounts/featurecounts.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/featurecounts/featurecounts.tsv",
    message: "Counting features using featureCounts",
    shell: """
    featureCounts \
    -a {input.annotation} \
    -F {params.annotation_format} \
    -T {threads} \
    {params.options} \
    {params.paired} \
    --donotsort \
    -o {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """