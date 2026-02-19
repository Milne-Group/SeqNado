from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


rule bam_filter:
    input:
        bam=OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/shifted_for_tn5_insertion/{sample}.bam.bai",
    output:
        bam=OUTPUT_DIR + "/aligned/filtered/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/filtered/{sample}.bam.bai",
        read_log=temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_filter.tsv"),
    params:
        options=str(CONFIG.third_party_tools.samtools.view.command_line_arguments),
    threads: CONFIG.third_party_tools.samtools.view.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_filter.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_filter.tsv",
    message: "Filtering aligned BAM for sample {wildcards.sample} using samtools",
    shell: """
    samtools view -@ {threads} -h -b {input.bam} {params.options} > {output.bam} &&
    samtools index {output.bam} &&
    echo -e "filtering\t$(samtools view -c {output.bam})" >> {output.read_log} 2>&1 | tee -a {log}
    """
