from seqnado.workflow.helpers.common import (
    define_memory_requested, 
    define_time_requested,
    get_alignment_input,
)


rule align_paired:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq1"],
        fq2=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq2"],
    output:
        bam_dir=temp(directory(OUTPUT_DIR + "/aligned/star/{sample}.star/")),
        log_out=temp(OUTPUT_DIR + "/aligned/star/{sample}.star/Log.final.out"),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/star/{sample}.star/",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_paired",
    threads: CONFIG.third_party_tools.star.align.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=35, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using STAR",
    shell: f"""
    STAR \
    --genomeDir {{params.index}} \
    --readFilesIn {{input.fq1}} {{input.fq2}} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN {{threads}} \
    --outSAMattrRGline ID:{{wildcards.sample}} SM:{{wildcards.sample}} \
    --outFileNamePrefix {{params.prefix}} \
    {{params.options}} \
    > {{log}} 2>&1 &&
    mapped=$(samtools view -c {{output.bam_dir}}/Aligned.sortedByCoord.out.bam) &&
    before=0 &&
    after="$mapped" &&
    {emit_read_logs("Aligned ({params.rule_label})", "{wildcards.sample}", "{params.read_log}")}
    """


rule align_single:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=False),
    output:
        bam_dir=temp(directory(OUTPUT_DIR + "/aligned/star/{sample}.star/")),
        log_out=temp(OUTPUT_DIR + "/aligned/star/{sample}.star/Log.final.out"),
    params:
        index=str(CONFIG.genome.index.prefix),
        options=str(CONFIG.third_party_tools.star.align.command_line_arguments),
        prefix=OUTPUT_DIR + "/aligned/star/{sample}.star/",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_single",
    threads: CONFIG.third_party_tools.star.align.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=35, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using STAR",
    shell: f"""
    STAR \
    --genomeDir {{params.index}} \
    --readFilesIn {{input.fq1}} \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN {{threads}} \
    --outSAMattrRGline ID:{{wildcards.sample}} SM:{{wildcards.sample}} \
    --outFileNamePrefix {{params.prefix}} \
    {{params.options}} \
    > {{log}} 2>&1 &&
    mapped=$(samtools view -c {{output.bam_dir}}/Aligned.sortedByCoord.out.bam) &&
    before=0 &&
    after="$mapped" &&
    {emit_read_logs("Aligned ({params.rule_label})", "{wildcards.sample}", "{params.read_log}")}
    """


rule rename_aligned:
    input:
        bam_dir=OUTPUT_DIR + "/aligned/star/{sample}.star/",
        log_out=OUTPUT_DIR + "/aligned/star/{sample}.star/Log.final.out",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    params:
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/rename_aligned/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/rename_aligned/{sample}.tsv",
    message: "Renaming aligned BAM for sample {wildcards.sample} to standard format",
    shell: f"""
    before=$(samtools view -c {{input.bam_dir}}/Aligned.sortedByCoord.out.bam) &&
    mv {{input.bam_dir}}/Aligned.sortedByCoord.out.bam {{output.bam}} >> {{log}} 2>&1 &&
    after=$(samtools view -c {{output.bam}}) &&
    {emit_read_logs("Rename Aligned BAM", "{wildcards.sample}", "{params.read_log}")}
    """
