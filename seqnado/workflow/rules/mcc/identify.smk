from seqnado.workflow.helpers.mcc import identify_extracted_bam_files, redefine_viewpoints

use rule sort_bam_by_qname as sort_genomic_aligned_reads with:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bam=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_qname_sorted.bam"),
        read_log=temp(OUTPUT_DIR + "/qc/mcc/{sample}_qname_sort.tsv"),
    threads: CONFIG.third_party_tools.samtools.sort.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/mcc/sort_genomic_aligned_reads/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/mcc/sort_genomic_aligned_reads/{sample}.tsv",
        

rule identify_viewpoint_reads:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_qname_sorted.bam",
    output:
        bam=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_annotated.bam"),
    params:
        output_dir=OUTPUT_DIR + "/mcc/{sample}/reporters/raw/",
    threads: 1
    resources:
        mem="1GB",
    container: "docker://ghcr.io/alsmith151/mccnado:latest"
    log: OUTPUT_DIR + "/logs/identify_viewpoint_reads/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/identify_viewpoint_reads/{sample}.tsv",
    message: "Identifying viewpoint reads for sample {wildcards.sample}",
    shell:
        """
        mccnado annotate-bam {input.bam} {output.bam} > {log} 2>&1
        """

rule deduplicate_bam_file:
    input: 
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_annotated.bam",
    output:
        bam=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_deduplicated.bam"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/mccnado:latest"
    log: OUTPUT_DIR + "/logs/deduplicate_bam/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deduplicate_bam/{sample}.tsv",
    message: "Deduplicating BAM file for sample {wildcards.sample}",
    shell: """
    mccnado deduplicate-bam {input.bam} {output.bam} > {log} 2>&1
    """

        
use rule sort_bam as sort_bam_viewpoints with:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_deduplicated.bam",
    output:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
        read_log=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}_read_log.txt",
    log: OUTPUT_DIR + "/logs/sort_bam_viewpoints/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/sort_bam_viewpoints/{sample}.tsv",
    message: "Sorting BAM file for viewpoints for sample {wildcards.sample}",

use rule index_bam as index_bam_viewpoints with:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
    output:
        bai=temp(OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam.bai"),
    log: OUTPUT_DIR + "/logs/index_bam_viewpoints/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/index_bam_viewpoints/{sample}.tsv",
    message: "Indexing BAM file for viewpoints for sample {wildcards.sample}",





