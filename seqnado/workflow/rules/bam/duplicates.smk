from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
from seqnado import PCRDuplicateTool

if CONFIG.pcr_duplicates.tool == PCRDuplicateTool.PICARD:
    rule bam_remove_duplicates_using_picard:
        input:
            bam=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai"),
            metrics=OUTPUT_DIR + "/qc/library_complexity/{sample}.metrics",
            read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_remove_duplicates.tsv"),
        threads: 8
        params:
            options=str(CONFIG.third_party_tools.picard.duplicate_removal.command_line_arguments),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_remove_duplicates.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_remove_duplicates.tsv",
        message: "Removing duplicates from aligned BAM for sample {wildcards.sample} using Picard",
        shell: """
        picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.metrics} --REMOVE_DUPLICATES true --CREATE_INDEX true {params.options} 2> {log} &&
        mv seqnado_output/aligned/duplicates_removed/{wildcards.sample}.bai {output.bai} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """

elif CONFIG.pcr_duplicates.tool == PCRDuplicateTool.SAMTOOLS:
    rule bam_remove_duplicates_using_samtools:
        input:
            bam=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai"),
            read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_remove_duplicates.tsv"),
        threads: CONFIG.third_party_tools.samtools.duplicate_removal.threads
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=5, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_remove_duplicates.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_remove_duplicates.tsv",
        message: "Removing duplicates from aligned BAM for sample {wildcards.sample} using samtools",
        shell: """
        samtools rmdup -@ {threads} {input.bam} {output.bam} &&
        samtools index {output.bam} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """
else:
    rule bam_ignore_duplicates:
        input:
            bam=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/blacklist_regions_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai"),
            read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_remove_duplicates.tsv"),
        threads: 8
        resources:
            mem="500MB",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_remove_duplicates.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_remove_duplicates.tsv",
        message: "Skipping duplicate removal for sample {wildcards.sample}",
        shell: """
        mv {input.bam} {output.bam} &&
        mv {input.bai} {output.bai} &&
        echo -e "duplicate removal\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
        """
