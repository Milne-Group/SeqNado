from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

if CONFIG.shift_for_tn5_insertion:
    rule shift_atac_alignments:
        input:
            bam=OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai",
        output:
            tmp=temp(OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam.tmp"),
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        threads: 1
        container: "docker://ghcr.io/alsmith151/bamnado:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_atac_shift.log",
        message: "Shifting ATAC-seq alignments for sample {wildcards.sample} using bamnado",
        shell: """
        bamnado modify --input {input.bam} --output {output.tmp} --tn5-shift 2>&1 | tee {log}
        """

    rule sort_and_index_shifted_bam:
        input:
            tmp=OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam.tmp",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam"),
            bai=temp(OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"),
            read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_atac_shift.tsv"),
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_atac_shift.tsv",
        message: "Sorting and indexing shifted ATAC-seq alignments for sample {wildcards.sample}",
        shell: """
        samtools sort {input.tmp} -@ {threads} -o {output.bam} &&
        samtools index {output.bam} &&
        echo -e "ATAC shift\t$(samtools view -c {output.bam})" > {output.read_log}
        """

else:
    rule move_bam_to_temp_location:
        input:
            bam=OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam",
            bai=OUTPUT_DIR + "/aligned/duplicates_removed/{sample}.bam.bai",
        output:
            bam=temp(OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam"),
            bai=temp(
                OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam.bai"
            ),
            read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_atac_shift.tsv"),
        threads: 1
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_atac_shift.log",
        benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_atac_shift.tsv",
        message: "Skipping ATAC-seq shift for sample {wildcards.sample}",
        shell: """
        mv {input.bam} {output.bam} &&
        mv {input.bam}.bai {output.bai} &&
        echo -e "ATAC shift\t$(samtools view -c {output.bam})" >> {output.read_log}
        """
