from seqnado.workflow.helpers.common import (
    define_time_requested, 
    define_memory_requested,
    get_alignment_input,
)


rule align_paired:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq1"],
        fq2=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq2"],
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_paired",
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using Bowtie2",
    shell: f"""
    bowtie2 \
        -p {{threads}} \
        -x {{params.index}} \
        -1 {{input.fq1}} \
        -2 {{input.fq2}} \
        {{params.rg}} \
        {{params.options}} \
        2> {{log}} \
    | samtools view -bS - > {{output.bam}} &&
    mapped=$(samtools view -c {{output.bam}}) &&
    before=0 &&
    after="$mapped" &&
    {emit_read_logs("Aligned ({params.rule_label})", "{wildcards.sample}", "{params.read_log}")}
    """


rule align_single:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=False),
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
        rg="--rg-id {sample} --rg SM:{sample}",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_single",
    threads: CONFIG.third_party_tools.bowtie2.align.threads,
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning reads for sample {wildcards.sample} using Bowtie2",
    shell: f"""
    bowtie2 \
        -p {{threads}} \
        -x {{params.index}} \
        -U {{input.fq1}} \
        {{params.rg}} \
        {{params.options}} \
        2> {{log}} \
    | samtools view -bS - > {{output.bam}} &&
    mapped=$(samtools view -c {{output.bam}}) &&
    before=0 &&
    after="$mapped" &&
    {emit_read_logs("Aligned ({params.rule_label})", "{wildcards.sample}", "{params.read_log}")}
    """


ruleorder: align_paired > align_single
