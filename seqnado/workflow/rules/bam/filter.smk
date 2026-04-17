from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested, get_read_count_flags


rule bam_filter:
    input:
        bam=OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam.bai",
    output:
        bam=OUTPUT_DIR + "/aligned/filtered/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/filtered/{sample}.bam.bai",
    params:
        options=str(CONFIG.third_party_tools.samtools.view.command_line_arguments),
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        count_flags=lambda wildcards: get_read_count_flags(wildcards, INPUT_FILES),
    threads: CONFIG.third_party_tools.samtools.view.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_filter.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_filter.tsv",
    message: "Filtering aligned BAM for sample {wildcards.sample} using samtools",
    shell: f"""
    before=$(samtools view -c {{params.count_flags}} {{input.bam}}) &&
    samtools view -@ {{threads}} -h -b {{input.bam}} {{params.options}} > {{output.bam}} 2>> {{log}} &&
    samtools index {{output.bam}} >> {{log}} 2>&1 &&
    after=$(samtools view -c {{params.count_flags}} {{output.bam}}) &&
    {emit_read_logs("Filter", "{wildcards.sample}", "{params.read_log}")}
    """
