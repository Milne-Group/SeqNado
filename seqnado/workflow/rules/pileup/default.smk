import re
from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested, format_deeptools_options
from seqnado.config.third_party_tools import CommandLineArguments

include: "common.smk"

rule homer_make_bigwigs:
    input:
        homer_tag_directory=OUTPUT_DIR + "/tag_dirs/{sample}",
    output:
        homer_bigwig=OUTPUT_DIR + "/bigwigs/homer/unscaled/{sample}.bigWig",
    params:
        genome_name=CONFIG.genome.name,
        genome_chrom_sizes=CONFIG.genome.chromosome_sizes,
        options=str(CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments),
        outdir=OUTPUT_DIR + "/bigwigs/homer/",
        temp_bw=lambda wc, output: output.homer_bigwig.replace(
            ".bigWig", ".ucsc.bigWig"
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/homer/makebigwigs_{sample}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/makebigwigs_{sample}.tsv"
    message:
        "Making bigWig with HOMER for sample {wildcards.sample}"
    shell:
        """
    makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} {params.options} > {log} 2>&1 &&
    mv {params.outdir}/{wildcards.sample}.ucsc.bigWig {output.homer_bigwig}
    """

rule deeptools_make_bigwigs:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options(
            wildcards,
            str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
            INPUT_FILES,
            SAMPLE_GROUPINGS,
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/pileups/deeptools/unscaled/{sample}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/deeptools/makebigwigs_{sample}.tsv"
    message:
        "Making bigWig with deeptools for sample {wildcards.sample}"
    shell:
        """
    bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule bamnado_bam_coverage:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/unscaled/{sample}.bigWig",
    params:
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads
    container:
        "docker://ghcr.io/alsmith151/bamnado:latest"
    log:
        OUTPUT_DIR + "/logs/pileups/bamnado/{sample}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/bamnado/makebigwigs_{sample}.tsv"
    message:
        "Making bigWig with bamnado for sample {wildcards.sample}"
    shell:
        """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule fragment_bedgraph:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        filtered=temp(OUTPUT_DIR + "/bedgraphs/{sample}.filtered.bam"),
        sort=temp(OUTPUT_DIR + "/bedgraphs/{sample}.sorted.bam"),
        bed=temp(OUTPUT_DIR + "/bedgraphs/{sample}.bed"),
        bed_log=temp(OUTPUT_DIR + "/logs/bedgraphs/{sample}_bamtobed.log"),
        fragments=temp(OUTPUT_DIR + "/bedgraphs/{sample}.fragments.bed"),
        bdg=temp(OUTPUT_DIR + "/bedgraphs/{sample}.bedGraph"),
    params:
        chromosome_sizes=CONFIG.genome.chromosome_sizes,
    threads: 16
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=12, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/bedgraphs/{sample}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/bedgraphs/fragment_bedgraph_{sample}.tsv"
    message:
        "Generating fragment bedGraph for sample {wildcards.sample}"
    shell:
        """
    echo "Generating fragment bedGraph for sample {wildcards.sample}" 2>&1 | tee {log}

    samtools view -@ {threads} -q 30 -f 2 -h {input.bam} | grep -v chrM > {output.filtered} 2>&1 | tee -a {log}
    samtools sort -@ {threads} -m 900M -o {output.sort} -T {output.sort}.tmp {output.filtered} 2>&1 | tee -a {log}
    bedtools bamtobed -bedpe -i {output.sort} > {output.bed} 2>&1 | tee -a {output.bed_log}
    awk '$1==$4 && $6-$2 < 1000' {output.bed} > {output.fragments}.temp 2>&1 | tee -a {log}
    awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $6}}' {output.fragments}.temp | sort -k1,1 -k2,2n -k3,3n > {output.fragments} 2>&1 | tee -a {log}
    bedtools genomecov -bg -i {output.fragments} -g {params.chromosome_sizes} > {output.bdg} 2>&1 | tee -a {log}
    rm -f {output.fragments}.temp
    """
