from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


rule trimgalore_paired:
    # Trim reads using trimgalore
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
    output:
        trimmed1=temp(OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp(OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz"),
    threads: CONFIG.third_party_tools.trim_galore.trim.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    params:
        options=str(CONFIG.third_party_tools.trim_galore.trim.command_line_arguments),
        trim_dir=OUTPUT_DIR + "/trimmed",
    log: OUTPUT_DIR + "/logs/trimming/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/trimming/{sample}.tsv",
    message: "Trimming reads for sample {wildcards.sample} using Trim Galore",
    shell:
        """
        trim_galore --cores {threads} {params.options} --basename {wildcards.sample} --paired --output_dir {params.trim_dir} {input.fq1} {input.fq2} >> {log} 2>&1 &&
        mv {params.trim_dir}/{wildcards.sample}_val_1.fq.gz {output.trimmed1} &&
        mv {params.trim_dir}/{wildcards.sample}_val_2.fq.gz {output.trimmed2}
        """


rule trimgalore_single:
    # Trim reads using trimgalore
    input:
        fq=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
    output:
        trimmed=temp(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz"),
    threads: CONFIG.third_party_tools.trim_galore.trim.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    params:
        options=str(CONFIG.third_party_tools.trim_galore.trim.command_line_arguments),
        trim_dir=OUTPUT_DIR + "/trimmed",
    log: OUTPUT_DIR + "/logs/trimming/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/trimming/{sample}.tsv",
    message: "Trimming reads for sample {wildcards.sample} using Trim Galore",
    shell:
        """
        trim_galore --cores {threads} {params.options} --basename {wildcards.sample} --output_dir {params.trim_dir} {input.fq} >> {log} 2>&1 &&
        mv {params.trim_dir}/{wildcards.sample}_trimmed.fq.gz {output.trimmed}
        """


ruleorder: trimgalore_paired > trimgalore_single
