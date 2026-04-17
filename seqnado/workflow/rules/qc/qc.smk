import os

from seqnado import Assay, AssaysWithPeakCalling, QuantificationMethod
from seqnado.config.third_party_tools import CommandLineArguments
from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested
from seqnado.outputs.core import SeqNadoReportFiles

from seqnado.workflow.helpers.qc import format_qualimap_options, format_frip_enrichment_options
##############################################
#                   FastQC                   #
############################################## 


rule fastqc_raw_paired:
    input:
        fq1=OUTPUT_DIR + "/fastqs/{sample}_1.fastq.gz",
        fq2=OUTPUT_DIR + "/fastqs/{sample}_2.fastq.gz",
    output:
        html1=OUTPUT_DIR + "/qc/fastqc_raw/{sample}_1_fastqc.html",
        html2=OUTPUT_DIR + "/qc/fastqc_raw/{sample}_2_fastqc.html",
        zip1=OUTPUT_DIR + "/qc/fastqc_raw/{sample}_1_fastqc.zip",
        zip2=OUTPUT_DIR + "/qc/fastqc_raw/{sample}_2_fastqc.zip",
    params:
        extra="--quiet",
        output_dir=OUTPUT_DIR + "/qc/fastqc_raw/",
    threads: CONFIG.third_party_tools.fastqc.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/fastqc_raw/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/fastqc_raw/{sample}.tsv",
    message: "Running FastQC on raw FASTQ files for sample {wildcards.sample}",
    shell:
        """
        fastqc -o {params.output_dir} {input.fq1} {input.fq2} > {log} 2>&1
        """


rule fastqc_raw_single:
    input:
        OUTPUT_DIR + "/fastqs/{sample}.fastq.gz",
    output:
        html=OUTPUT_DIR + "/qc/fastqc_raw/{sample}_fastqc.html",
        zip=OUTPUT_DIR + "/qc/fastqc_raw/{sample}_fastqc.zip",
    params:
        extra="--quiet",
        output_dir=OUTPUT_DIR + "/qc/fastqc_raw/",
    threads: CONFIG.third_party_tools.fastqc.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/fastqc_raw/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/fastqc_raw/{sample}.tsv",
    message: "Running FastQC on raw FASTQ files for sample {wildcards.sample}",
    shell:
        """
        fastqc -o {params.output_dir} {input} > {log} 2>&1
        """

##############################################
#                  Qualimap                  #
##############################################


rule qualimap_bamqc:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    output:
        html=OUTPUT_DIR + "/qc/qualimap_bamqc/{sample}/qualimapReport.html",
    params:
        output_dir=OUTPUT_DIR + "/qc/qualimap_bamqc/{sample}/",
        options=lambda wc: format_qualimap_options(wc, CONFIG.third_party_tools.qualimap.command_line_arguments, INPUT_FILES=INPUT_FILES, ASSAY=ASSAY)
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/qualimap_bamqc/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/qualimap_bamqc/{sample}.tsv",
    message: "Running Qualimap BAMQC for sample {wildcards.sample}",
    shell:
        """
        qualimap --java-mem-size={resources.mem} bamqc \
        {params.options} \
        -nt {threads} \
        -bam {input.bam} \
        -outdir {params.output_dir} \
        > {log} 2>&1
        """


rule qualimap_rnaseq:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    output:
        html=OUTPUT_DIR + "/qc/qualimap_rnaseq/{sample}/qualimapReport.html",
    params:
        output_dir=OUTPUT_DIR + "/qc/qualimap_rnaseq/{sample}/",
        annotation=CONFIG.genome.gtf,
        options=lambda wc: format_qualimap_options(wc, CONFIG.third_party_tools.qualimap.command_line_arguments, INPUT_FILES=INPUT_FILES, ASSAY=ASSAY)
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 16
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/qualimap_rnaseq/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/qualimap_rnaseq/{sample}.tsv",
    message: "Running Qualimap RNA-seq for sample {wildcards.sample}",
    shell:
        """
        qualimap --java-mem-size={resources.mem} rnaseq \
        {params.options} \
        -bam {input.bam} \
        -gtf {params.annotation} \
        -outdir {params.output_dir} \
        > {log} 2>&1    
        """

##############################################
#            BAM filtering stats             #
##############################################


rule bam_stats:
    input: 
        sort=OUTPUT_DIR + "/qc/alignment_post_process/{sample}_sort.tsv",
        blacklist=OUTPUT_DIR + "/qc/alignment_post_process/{sample}_blacklist.tsv",
        remove_duplicates=OUTPUT_DIR + "/qc/alignment_post_process/{sample}_remove_duplicates.tsv",
        atac_shift=OUTPUT_DIR + "/qc/alignment_post_process/{sample}_atac_shift.tsv",
        filtered=OUTPUT_DIR + "/qc/alignment_post_process/{sample}_filter.tsv",
        final=OUTPUT_DIR + "/qc/alignment_post_process/{sample}_final.tsv",
    output: temp(OUTPUT_DIR + "/qc/alignment_post_process/{sample}_alignment_stats.tsv")
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_post_process/{sample}_alignment_stats.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_post_process/{sample}_alignment_stats.tsv",
    message: "Compiling alignment post-processing stats for sample {wildcards.sample}",
    shell:
        """
        cat {input.sort} {input.blacklist} {input.remove_duplicates} {input.atac_shift} {input.filtered} {input.final} > {output}
        """


rule prepare_stats_report:
    input:
        expand(
            OUTPUT_DIR + "/qc/alignment_post_process/{sample}_alignment_stats.tsv",
            sample=SAMPLE_NAMES,
        ),
    output:
        OUTPUT_DIR + "/qc/alignment_stats.tsv",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/alignment_stats.log",
    benchmark: OUTPUT_DIR + "/.benchmark/alignment_stats.tsv",
    message: "Generating alignment stats report for all samples",
    script:
        "../../scripts/alignment_stats.py"


##############################################
#               Frip Enrichment              #
##############################################


rule frip_enrichment:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        peak=OUTPUT_DIR + "/peaks/{directory}/{sample}.bed",
    output:
        pdf=OUTPUT_DIR + "/qc/frip_enrichment/{directory}/{sample}_frip.pdf",
        frip_count=OUTPUT_DIR + "/qc/frip_enrichment/{directory}/{sample}_frip.txt",
    params:
        options=lambda wc: format_frip_enrichment_options(wc, CommandLineArguments(), INPUT_FILES=INPUT_FILES),
    threads: 16
    log: OUTPUT_DIR + "/logs/frip_enrichment/{directory}/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/frip_enrichment/{directory}/{sample}.tsv",
    message: "Running FRiP enrichment for sample {wildcards.sample}",
    shell:
        """
        plotEnrichment -p {threads} \
        --bamfiles {input.bam} \
        --BED {input.peak} \
        --outRawCounts {output.frip_count} \
        --plotFile {output.pdf} \
        > {log} 2>&1
        """

##############################################
#                  MultiQC                   #
##############################################

multiqc_input_files = SeqNadoReportFiles(
    assay=ASSAY,
    samples=INPUT_FILES,
    config=CONFIG,
    sample_groupings=SAMPLE_GROUPINGS,
    output_dir=OUTPUT_DIR,
).gather_input_files


rule seqnado_report:
    input:
        multiqc_input_files,
    output:
        report = OUTPUT_DIR + f"/seqnado_report_{ASSAY.clean_name}.html",
    params:
        multiqc_config = "/opt/seqnado/multiqc_config.yaml",
        output_dir = OUTPUT_DIR,
        assay_name = ASSAY.clean_name,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/multiqc/seqnado_report.log",
    benchmark: OUTPUT_DIR + "/.benchmark/multiqc/seqnado_report.tsv",
    message: "Generating SeqNado report",
    shell:
        """
        multiqc -o {params.output_dir} {params.output_dir} \
        --config {params.multiqc_config} \
        --filename "seqnado_report_{params.assay_name}.html" \
        --no-data-dir \
        --force > {log} 2>&1
        """


ruleorder: fastqc_raw_paired > fastqc_raw_single > seqnado_report
