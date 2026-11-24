from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested
from seqnado.workflow.helpers.bam import get_bam_split
from seqnado.workflow.helpers.spikein import reference_contig_regex, spikein_contig_regex

REF_GENOME=CONFIG.assay_config.methylation.reference_genome if CONFIG.assay_config.methylation and CONFIG.assay_config.methylation.reference_genome else CONFIG.genome.name
SPIKEIN_GENOMES=CONFIG.assay_config.methylation.spikein_genomes if CONFIG.assay_config.methylation and CONFIG.assay_config.methylation.spikein_genomes else []
GENOMES=[REF_GENOME]+SPIKEIN_GENOMES
METH_REF_CONTIG_REGEX=reference_contig_regex(
    CONFIG.assay_config.methylation.endogenous_contig_prefix
    if CONFIG.assay_config.methylation
    else None
)
# All spike-in contigs; excluded from the reference split.
METH_ALL_SPIKEIN_REGEX=spikein_contig_regex(SPIKEIN_GENOMES) if SPIKEIN_GENOMES else ""


checkpoint methylation_bam_splits:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_{genome}.bam"),
        bai=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_{genome}.bam.bai")
    params:
        script=workflow.source_path("../../scripts/split_bam_by_contig.sh"),
        ref_regex=METH_REF_CONTIG_REGEX,
        # Reference: exclude every spike-in. Spike-in: select just this genome.
        exo_regex=lambda wc: METH_ALL_SPIKEIN_REGEX if wc.genome == REF_GENOME else spikein_contig_regex([wc.genome]),
        out_flag=lambda wc: "--ref-out" if wc.genome == REF_GENOME else "--exo-out",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/methylation/bam_splits/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation_bam_splits/{sample}_{genome}.tsv"
    message: "Splitting BAM for sample {wildcards.sample} into genome {wildcards.genome}"
    shell: """
        bash {params.script} --bam {input.bam} --ref-regex {params.ref_regex:q} --exo-regex {params.exo_regex:q} {params.out_flag} {output.bam} --on-empty-ref error --on-empty-exo warn --index --label {wildcards.sample}_{wildcards.genome} --log {log}
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
    container:CONTAINERS.SEQNADO_PIPELINE
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
        assay=CONFIG.assay_config.methylation.method.value if CONFIG.assay_config.methylation else None,
        method=CONFIG.assay_config.methylation.method.value if CONFIG.assay_config.methylation else None
    container: CONTAINERS.SEQNADO_PIPELINE
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
    container: CONTAINERS.SEQNADO_PIPELINE
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
    log: OUTPUT_DIR + "/logs/methylation/taps_inverted/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/taps_inverted/{sample}_{genome}.tsv"
    message: "Converting to TAPS methylation for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    mkdir -p $(dirname {output.taps}) &&
    mkdir -p $(dirname {log}) &&
    awk -v OFS="\t" '{{print $1, $2, $3, (100-$4), $5, $6}}' {input.bdg} > {output.taps} 2> {log}
    """

rule bigwigs_meth_taps:
    input:
        bdg=rules.taps_inverted.output.taps
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/taps/{sample}_{genome}.bigWig",
        bdg=temp(OUTPUT_DIR + "/bigwigs/taps/{sample}_{genome}.bedGraph")
    params:
        chrom_sizes=CONFIG.genome.chromosome_sizes,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/methylation/taps/bigwig/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/taps/bigwig/{sample}_{genome}.tsv"
    message: "Converting to bigWig for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    awk -v OFS="\t" '{{print $1, $2, $3, $4}}' {input.bdg} > {output.bdg} 2> {log}
    bedGraphToBigWig {output.bdg} {params.chrom_sizes} {output.bigwig} > {log} 2>&1
    """

rule bigwigs_meth_wgbs:
    input:
        bdg=rules.methyldackel_extract.output.bdg
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/wgbs/{sample}_{genome}.bigWig",
        bdg=temp(OUTPUT_DIR + "/bigwigs/wgbs/{sample}_{genome}.bedGraph")
    params:
        chrom_sizes=CONFIG.genome.chromosome_sizes,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/methylation/wgbs/bigwig/{sample}_{genome}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/methylation/wgbs/bigwig/{sample}_{genome}.tsv"
    message: "Converting to bigWig for sample {wildcards.sample} and genome {wildcards.genome}"
    shell: """
    awk -v OFS="\t" '{{print $1, $2, $3, $4}}' {input.bdg} > {output.bdg} 2> {log}
    bedGraphToBigWig {output.bdg} {params.chrom_sizes} {output.bigwig} > {log} 2>&1
    """
