from seqnado.workflow.helpers.common import (
    define_time_requested,
    define_memory_requested,
    get_alignment_input,
)


rule align_crispr_paired:
    input:
        fq1=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq1"],
        fq2=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=True)["fq2"],
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
    output:
        bam_unsorted=temp(OUTPUT_DIR + "/aligned/{sample}.unsorted.bam"),
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: CONFIG.third_party_tools.bowtie2.align.threads
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning CRISPR sample {wildcards.sample} using Bowtie2 (paired-end)",
    shell: """
    bowtie2 -p {threads} -x {params.index} -1 {input.fq1} -2 {input.fq2} {params.options} 2> {log} |
    samtools view -bS - > {output.bam_unsorted} &&
    samtools sort -@ {threads} {output.bam_unsorted} -o {output.bam}
    """


rule align_crispr_single:
    input:
        fq=lambda wildcards: get_alignment_input(wildcards, OUTPUT_DIR, CONFIG, paired=False),
    params:
        index=CONFIG.genome.index.prefix,
        options=str(CONFIG.third_party_tools.bowtie2.align.command_line_arguments),
    output:
        bam_unsorted=temp(OUTPUT_DIR + "/aligned/{sample}.unsorted.bam"),
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: CONFIG.third_party_tools.bowtie2.align.threads
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/align/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}.tsv",
    message: "Aligning CRISPR sample {wildcards.sample} using Bowtie2 (single-end)",
    shell: """
    bowtie2 -p {threads} -x {params.index} -U {input.fq} {params.options} 2> {log} |
    samtools view -bS - > {output.bam_unsorted} &&
    samtools sort -@ {threads} {output.bam_unsorted} -o {output.bam}
    """


rule index_crispr_bam:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    output:
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    threads: CONFIG.third_party_tools.samtools.index.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/align/{sample}_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/align/{sample}_index.tsv",
    message: "Indexing BAM file for CRISPR sample {wildcards.sample}",
    shell: """
    samtools index {input.bam} {output.bai} > {log} 2>&1
    """
    


ruleorder: align_crispr_paired > align_crispr_single
