
###################################################
# Spike-in alignment and processing for DNA assays #
###################################################

use rule align_paired as align_paired_spikein with:
    params:
        options="--no-mixed --no-discordant",
        index=CONFIG.genome.index.prefix,
        rg="--rg-id {sample} --rg SM:{sample}",
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_paired_spikein",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/raw/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/aligned_spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned_spikein/{sample}.tsv",
    message: "Aligning spike-in reads for sample {wildcards.sample} using Bowtie2"


use rule align_single as align_single_spikein with:
    params:
        read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
        rule_label="align_single_spikein",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/raw/{sample}.bam"),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/aligned_spikein/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned_spikein/{sample}.tsv",
    message: "Aligning spike-in reads for sample {wildcards.sample} using Bowtie2"
