# Combine bam files for merge
from seqnado.workflow.helpers.bam import get_bam_files_for_consensus


rule bam_merge:
    input:
        bams=lambda wc: get_bam_files_for_consensus(wc, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS, OUTPUT_DIR=OUTPUT_DIR),
    output:
        temp(OUTPUT_DIR + "/aligned/merged/{group}.bam"),
    threads: CONFIG.third_party_tools.samtools.merge.threads
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/merge_bam/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/merge_bam/{group}.tsv",
    message: "Merging BAM files for group {wildcards.group} using samtools",
    shell: """
    echo "Input BAM files: {input.bams}" > {log} 2>&1
    samtools merge {output} {input} -@ {threads} >> {log} 2>&1
    """


use rule bam_index as bam_index_consensus with:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        bai=temp(OUTPUT_DIR + "/aligned/merged/{group}.bam.bai"),
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    threads: 8
    log: OUTPUT_DIR + "/logs/merge_bam/{group}_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/merge_bam/{group}_index.tsv",
    message: "Indexing merged BAM for group {wildcards.group} using samtools"