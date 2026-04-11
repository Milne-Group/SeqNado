CUTADAPT_MODE = CONFIG.third_party_tools.cutadapt.mode if CONFIG.third_party_tools.cutadapt else None
CUTADAPT_MODE_VALUE = CUTADAPT_MODE.value if CUTADAPT_MODE else None
AUTO_DETECT_ADAPTERS = CUTADAPT_MODE_VALUE == 'crispr'


rule detect_crispr_adapters_paired:
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
    output:
        adapters=OUTPUT_DIR + "/resources/{sample}_adapters.json",
    params:
        n_reads=10000,
        min_length=6,
        max_length=100,
        min_frequency=0.7,
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/adapter_detection/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/adapter_detection/{sample}.tsv",
    message: "Detecting adapters for CRISPR sample {wildcards.sample} (paired-end)",
    script:
        "../../scripts/detect_crispr_adapters.py"


rule detect_crispr_adapters_single:
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
    output:
        adapters=OUTPUT_DIR + "/resources/{sample}_adapters.json",
    params:
        n_reads=10000,
        min_length=6,
        max_length=100,
        min_frequency=0.7,
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/adapter_detection/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/adapter_detection/{sample}.tsv",
    message: "Detecting adapters for CRISPR sample {wildcards.sample} (single-end)",
    script:
        "../../scripts/detect_crispr_adapters.py"


ruleorder: detect_crispr_adapters_paired > detect_crispr_adapters_single