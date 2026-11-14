from seqnado.workflow.helpers.crispr import get_cutadapt_adapter_args


rule crispr_trimming_paired:
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
        adapters=OUTPUT_DIR + "/resources/{sample}_adapters.json",
    output:
        trimmed1=temp(OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz"),
        trimmed2=temp(OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz"),
    params:
        options=lambda wc: get_cutadapt_adapter_args(wc, CONFIG=CONFIG, OUTPUT_DIR=OUTPUT_DIR),
        trim_dir=OUTPUT_DIR + "/trimmed",
    threads: CONFIG.third_party_tools.cutadapt.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/trimming/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/trimming/{sample}.tsv",
    message: "Trimming adapters from CRISPR sample {wildcards.sample} using cutadapt (paired-end)",
    shell: """
    cutadapt {params.options} -o {output.trimmed1} -p {output.trimmed2} {input.fq1} {input.fq2} > {log} 2>&1
    """


rule crispr_trimming_single:
    input:
        fq=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
        adapters=OUTPUT_DIR + "/resources/{sample}_adapters.json",
    output:
        trimmed=temp(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz"),
    params:
        options=lambda wc: get_cutadapt_adapter_args(wc, CONFIG=CONFIG, OUTPUT_DIR=OUTPUT_DIR),
        trim_dir=OUTPUT_DIR + "/trimmed",
    threads: CONFIG.third_party_tools.cutadapt.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/trimming/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/trimming/{sample}.tsv",
    message: "Trimming adapters from CRISPR sample {wildcards.sample} using cutadapt (single-end)",
    shell: """
    cutadapt {params.options} -o {output.trimmed} {input.fq} > {log} 2>&1
    """


ruleorder: crispr_trimming_paired > crispr_trimming_single