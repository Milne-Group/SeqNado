from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested
from seqnado.workflow.helpers.bam import get_bam_split

REF_GENOME=CONFIG.assay_config.methylation.reference_genome if CONFIG.assay_config.methylation and CONFIG.assay_config.methylation.reference_genome else CONFIG.genome.name
SPIKEIN_GENOMES=CONFIG.assay_config.methylation.spikein_genomes if CONFIG.assay_config.methylation and CONFIG.assay_config.methylation.spikein_genomes else []
GENOMES=[REF_GENOME]+SPIKEIN_GENOMES

checkpoint methylation_bam_splits:
    input:
        OUTPUT_DIR + "/aligned/{sample}.bam"
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_{genome}.bam"),
        bai=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_{genome}.bam.bai")
    params:
        ref_genome=REF_GENOME
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/methylation/bam_splits/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation_bam_splits/{sample}_{genome}.tsv"
    message: "Splitting BAM for sample {wildcards.sample} into genome {wildcards.genome}"
    shell: """
    if [[ "{wildcards.genome}" == "{params.ref_genome}" ]]; then
        samtools view -h {input} | awk '{{if($0 ~ /^@/ || $3 ~ /^chr/) print}}' | samtools view -b -o {output.bam} 2> {log}
    else
        samtools view -h {input} | awk '{{if($0 ~ /^@/ || $3 ~ /^{wildcards.genome}/) print}}' | samtools view -b -o {output.bam} 2> {log}
    fi
    samtools index {output.bam}
    """


rule methyldackel_bias:
    input:
        bam=lambda wc: get_bam_split(wc, checkpoints)
    output:
        bias=OUTPUT_DIR + "/methylation/methyldackel/bias/{sample}_{genome}.txt",
    params:
        fasta=CONFIG.genome.fasta,
        prefix=OUTPUT_DIR + "/methylation/methyldackel/bias/{sample}_{genome}"
    threads: CONFIG.third_party_tools.methyldackel.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container:"oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:OUTPUT_DIR + "/logs/methylation/methyldackel/bias/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/methyldackel/bias/{sample}_{genome}.tsv"
    message: "Running MethylDackel bias for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    MethylDackel mbias -@ {threads} --txt {params.fasta} {input.bam} {params.prefix} > {output.bias} 2> {log}
    """


rule calculate_conversion:
    input:
        expand(OUTPUT_DIR + "/methylation/methyldackel/bias/{sample}_{genome}.txt", sample=SAMPLE_NAMES, genome=GENOMES)
    output:
        conversion=OUTPUT_DIR + "/methylation/methylation_conversion.tsv",
        plot=OUTPUT_DIR + "/methylation/methylation_conversion.png"
    params:
        assay=str(CONFIG.assay_config.methylation.method) if CONFIG.assay_config.methylation and hasattr(CONFIG.assay_config.methylation, 'method') else None,
        method=str(CONFIG.assay_config.methylation.method) if CONFIG.assay_config.methylation and hasattr(CONFIG.assay_config.methylation, 'method') else None
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/methylation/conversion.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/calculate_conversion.tsv"
    message: "Calculating methylation conversion rates across all samples"
    script: "../../scripts/methylation_conversion.py"


rule methyldackel_extract:
    input:
        bam=rules.methylation_bam_splits.output.bam
    output:
        bdg=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}_CpG.bedGraph"
    params:
        fasta=CONFIG.genome.fasta,
        options=CONFIG.third_party_tools.methyldackel.command_line_arguments,
        prefix=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}"
    threads: CONFIG.third_party_tools.methyldackel.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/methylation/methyldackel/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/methyldackel/{sample}_{genome}.tsv"
    message: "Running MethylDackel extract for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    MethylDackel extract -@ {threads} {params.options} -o {params.prefix} {params.fasta} {input.bam} > {log} 2>&1
    """


rule taps_inverted:
    input:
        bdg=rules.methyldackel_extract.output.bdg
    output:
        taps=OUTPUT_DIR + "/methylation/methyldackel/{sample}_{genome}_CpG_inverted.bedGraph"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/methylation/taps_inverted/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/taps_inverted/{sample}_{genome}.tsv"
    message: "Converting to TAPS methylation for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    awk -v OFS="\t" '{{print $1, $2, $3, (100-$4), $5, $6}}' {input.bdg} > {output.taps} 2> {log}
    rm {input.bdg}
    """

localrules:
    calculate_conversion