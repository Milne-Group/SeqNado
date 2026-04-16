from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
import shutil

localrules: fastq_screen_config


rule fastq_screen_config:
    input:
        conf=CONFIG.third_party_tools.fastq_screen.config,
    output:
        conf=temp(OUTPUT_DIR + "/qc/fastq_screen/fastq_screen.conf"),
    log: OUTPUT_DIR + "/logs/fastq_screen/fastq_screen_config.log",
    benchmark: OUTPUT_DIR + "/.benchmark/fastq_screen/fastq_screen_config.tsv",
    message: "Copying fastq_screen configuration file to output directory",
    run:
        shutil.copy(input.conf, output.conf)


rule fastq_screen_paired:
    input:
        fq=OUTPUT_DIR + "/fastqs/{sample}_{read}.fastq.gz",
        conf=OUTPUT_DIR + "/qc/fastq_screen/fastq_screen.conf",
    output:
        fq_screen=temp(OUTPUT_DIR + "/qc/fastq_screen/{sample}_{read}_screen.html"),
        fq_screen_png=temp(OUTPUT_DIR + "/qc/fastq_screen/{sample}_{read}_screen.png"),
        fq_screen_txt=OUTPUT_DIR + "/qc/fastq_screen/{sample}_{read}_screen.txt",
    params:
        outdir=temp(OUTPUT_DIR + "/qc/fastq_screen"),
        options=str(CONFIG.third_party_tools.fastq_screen.command_line_arguments),
    threads: CONFIG.third_party_tools.fastq_screen.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/fastq_screen/{sample}_{read}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/fastq_screen/{sample}_{read}.tsv",
    message: "Running fastq_screen for sample {wildcards.sample} read {wildcards.read}",
    shell: """ 
    fastq_screen --conf {input.conf} {params.options} --aligner bowtie2 --threads {threads} {input.fq}  --outdir {params.outdir} > {log} 2>&1 
    """


use rule fastq_screen_paired as fastq_screen_single with:
    input:
        fq=OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
        conf=OUTPUT_DIR + "/qc/fastq_screen/fastq_screen.conf",
    output:
        fq_screen=temp(OUTPUT_DIR + "/qc/fastq_screen/{sample}_screen.html"),
        fq_screen_png=temp(OUTPUT_DIR + "/qc/fastq_screen/{sample}_screen.png"),
        fq_screen_txt=temp(OUTPUT_DIR + "/qc/fastq_screen/{sample}_screen.txt"),
    log: OUTPUT_DIR + "/logs/fastq_screen/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/fastq_screen/{sample}.tsv",
    message: "Running fastq_screen for sample {wildcards.sample}",


ruleorder: fastq_screen_paired > fastq_screen_single 
