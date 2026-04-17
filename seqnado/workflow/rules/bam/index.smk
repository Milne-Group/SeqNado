from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


rule bam_sort:
    input:
        bam=OUTPUT_DIR + "/aligned/raw/{sample}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/sorted/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=16, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    threads: CONFIG.third_party_tools.samtools.sort.threads
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_sort.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_sort.tsv",
    message: "Sorting aligned BAM for sample {wildcards.sample} using samtools",
    shell: f"""
    samtools sort {{input.bam}} -@ {{threads}} -o {{output.bam}} -m 900M >> {{log}} 2>&1 
    """

rule bam_sort_by_qname:
    input:
        bam=OUTPUT_DIR + "/aligned/sorted/{sample}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/sorted_by_qname/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    threads: CONFIG.third_party_tools.samtools.sort.threads
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_filter_qname.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_filter_qname.tsv",
    message: "Sorting aligned BAM by QNAME for sample {wildcards.sample} using samtools",
    shell: f"""
    samtools sort -n {{input.bam}} -@ {{threads}} -o {{output.bam}} -m 900M >> {{log}} 2>&1 
    """


rule bam_index:
    input:
        bam=OUTPUT_DIR + "/aligned/sorted/{sample}.bam",
    output:
        bai=temp(OUTPUT_DIR + "/aligned/sorted/{sample}.bam.bai"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_index.tsv",
    message: "Indexing aligned BAM for sample {wildcards.sample} using samtools",
    shell: f"""
    samtools index -@ {{threads}} -b {{input.bam}} >> {{log}} 2>&1
    """


rule bam_move_to_final_location:
    input:
        bam=OUTPUT_DIR + "/aligned/filtered/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/filtered/{sample}.bam.bai",
    output:
        bam=OUTPUT_DIR + "/aligned/{sample,[A-Za-z\\d\\-_]+}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample,[A-Za-z\\d\\-_]+}.bam.bai",
    params:
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        log_entity="{sample}",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_final.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_final.tsv",
    message: "Moving final BAM for sample {wildcards.sample} to final location",
    shell: f"""
    before=$(samtools view -c {{input.bam}}) &&
    cp {{input.bam}} {{output.bam}} >> {{log}} 2>&1 &&
    cp {{input.bai}} {{output.bai}} >> {{log}} 2>&1 &&
    after=$(samtools view -c {{output.bam}}) &&
    {emit_read_logs("Finalise", "{params.log_entity}", "{params.read_log}")}
    """
