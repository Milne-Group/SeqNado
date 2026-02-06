from seqnado.workflow.helpers.mcc import (
    get_n_cis_scaling_factor,
    get_mcc_bam_files_for_merge,
)
from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested 

rule merge_mcc_bams:
    input:
        bams=lambda wc: get_mcc_bam_files_for_merge(
            wc, SAMPLE_GROUPINGS=SAMPLE_GROUPINGS, OUTPUT_DIR=OUTPUT_DIR
        ),
    output:
        bam=temp(OUTPUT_DIR + "/mcc/{group}/{group}.bam"),
    threads: CONFIG.third_party_tools.samtools.merge.threads
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/merge_bam/{group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/merge_bam/{group}.tsv",
    message:
        "Merging BAM files for group {wildcards.group}",
    shell:
        """
    samtools merge -o {output.bam} {input.bams} -@ {threads}
    """


use rule index_bam as index_bam_merged with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
    output:
        bai=OUTPUT_DIR + "/mcc/{group}/{group}.bam.bai",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
    log:
        OUTPUT_DIR + "/logs/index_bam_merged/{group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/index_bam_merged/{group}.tsv",
    message:
        "Indexing merged BAM for group {wildcards.group}"