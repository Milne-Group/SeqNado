from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested, get_read_count_flags
from seqnado import PCRDuplicateHandling, PCRDuplicateTool

picard_config = getattr(getattr(CONFIG, "third_party_tools", None), "picard", None)
picard_mark_duplicates = getattr(picard_config, "mark_duplicates", None)
samtools_config = getattr(getattr(CONFIG, "third_party_tools", None), "samtools", None)
samtools_sort = getattr(samtools_config, "sort", None)

if (
    CONFIG.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE
    and CONFIG.pcr_duplicates.tool == PCRDuplicateTool.SAMTOOLS
):
    rule bam_remove_duplicates_using_samtools:
        input:
            bam=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai"),
        params:
            read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
            count_flags=lambda wildcards: get_read_count_flags(wildcards, INPUT_FILES),
        threads: samtools_sort.threads if samtools_sort is not None else 8
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_remove_duplicates.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_remove_duplicates.tsv",
        message: "Removing duplicates from aligned BAM for sample {wildcards.sample} using samtools",
        shell: f"""
        echo "Removing duplicates with samtools" > {{log}} 2>&1 &&
        before=$(samtools view -c {{params.count_flags}} {{input.bam}}) &&
        samtools rmdup -@ {{threads}} {{input.bam}} {{output.bam}} >> {{log}} 2>&1 &&
        samtools index {{output.bam}} >> {{log}} 2>&1 &&
        after=$(samtools view -c {{params.count_flags}} {{output.bam}}) &&
        {emit_read_logs("Remove Duplicates", "{wildcards.sample}", "{params.read_log}")}
        """

elif CONFIG.pcr_duplicates.strategy == PCRDuplicateHandling.REMOVE and picard_mark_duplicates is not None:
    rule bam_remove_duplicates_using_picard:
        input:
            bam=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai"),
            metrics=OUTPUT_DIR + "/qc/library_complexity/{sample}.metrics",
        threads: picard_mark_duplicates.threads
        params:
            options=str(picard_mark_duplicates.command_line_arguments),
            read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
            count_flags=lambda wildcards: get_read_count_flags(wildcards, INPUT_FILES),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_remove_duplicates.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_remove_duplicates.tsv",
        message: "Removing duplicates from aligned BAM for sample {wildcards.sample} using Picard",
        shell: f"""
        echo "Removing duplicates with Picard MarkDuplicates" > {{log}} 2>&1 &&
        before=$(samtools view -c {{params.count_flags}} {{input.bam}}) &&
        picard MarkDuplicates I={{input.bam}} O={{output.bam}} M={{output.metrics}} CREATE_INDEX=true {{params.options}} >> {{log}} 2>&1 &&
        mv {OUTPUT_DIR}/aligned/duplicates_removed/{{wildcards.sample}}.bai {{output.bai}} &&
        after=$(samtools view -c {{params.count_flags}} {{output.bam}}) &&
        {emit_read_logs("Remove Duplicates", "{wildcards.sample}", "{params.read_log}")}
        """

else:
    rule bam_ignore_duplicates:
        input:
            bam=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai"),
        params:
            read_log=read_log_shared_path(OUTPUT_DIR, "{sample}"),
            count_flags=lambda wildcards: get_read_count_flags(wildcards, INPUT_FILES),
        threads: 8
        resources:
            mem="500MB",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_remove_duplicates.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_remove_duplicates.tsv",
        message: "Skipping duplicate removal for sample {wildcards.sample}",
        shell: f"""
        echo "Skipping duplicate removal" > {{log}} 2>&1 &&
        before=$(samtools view -c {{params.count_flags}} {{input.bam}}) &&
        cp {{input.bam}} {{output.bam}} >> {{log}} 2>&1 &&
        cp {{input.bai}} {{output.bai}} >> {{log}} 2>&1 &&
        after=$(samtools view -c {{params.count_flags}} {{output.bam}}) &&
        {emit_read_logs("Remove Duplicates", "{wildcards.sample}", "{params.read_log}")}
        """
