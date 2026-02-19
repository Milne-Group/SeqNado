import os
import re

from seqnado.workflow.helpers.geo import load_geo_samples
from seqnado.workflow.helpers.common import (
    define_memory_requested,
    define_time_requested,
)

SCALE_RESOURCES = float(os.environ.get("SCALE_RESOURCES", "1"))
OUT_DIR = "seqnado_output"

rule geo_prefetch:
    """
    Download SRA data using prefetch
    """
    output:
        sra=temp(config["geo_outdir"] + "/sra_cache/{srr}/{srr}.sra")
    params:
        outdir=config["geo_outdir"] + "/sra_cache",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    message:
        "Prefetching SRA data for {wildcards.srr}"
    log:
        OUT_DIR + "/logs/geo_download/geo_prefetch/{srr}.log"
    benchmark:
        OUT_DIR + "/.benchmark/geo_prefetch/{srr}.benchmark.tsv"
    container:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    shell:
        """
        exec &> {log}
        echo "Downloading {wildcards.srr}"
        prefetch {wildcards.srr} --max-size 50G -O {params.outdir}
        echo "Done prefetching {wildcards.srr}"
        """


rule geo_fastq_dump_paired:
    """
    Extract paired-end FASTQ files from prefetched SRA data
    """
    wildcard_constraints:
        sample_name="|".join(re.escape(s) for s in config.get("geo_samples_paired", {}).keys()),
    input:
        sra=lambda wc: "{outdir}/sra_cache/{srr}/{srr}.sra".format(
            outdir=config["geo_outdir"],
            srr=config.get("geo_samples_paired", {}).get(wc.sample_name, {}).get("srr", "INVALID"),
        ),
    output:
        r1=config["geo_outdir"] + "/{sample_name}_R1.fastq",
        r2=config["geo_outdir"] + "/{sample_name}_R2.fastq",
    params:
        srr=lambda wc: config.get("geo_samples_paired", {}).get(wc.sample_name, {}).get("srr", ""),
        outdir=config["geo_outdir"],
    threads: 8
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    message:
        "Extracting paired-end FASTQ files for {wildcards.sample_name} ({params.srr})"
    log:
        OUT_DIR + "/logs/geo_download/geo_fastq_dump/{sample_name}.log"
    benchmark:
        OUT_DIR + "/.benchmark/geo_fastq_dump/{sample_name}.benchmark.tsv"
    container:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    shell:
        """
        exec &> {log}

        echo "Extracting FASTQ files from {params.srr}"
        fasterq-dump {input.sra} --split-3 --skip-technical --outdir {params.outdir} -e {threads}

        echo "Found split paired-end files"
        mv {params.outdir}/{params.srr}_1.fastq {output.r1}
        mv {params.outdir}/{params.srr}_2.fastq {output.r2}
        rm -f {params.outdir}/{params.srr}.fastq

        echo "Done extracting {wildcards.sample_name}"
        """


rule geo_fastq_dump_single:
    """
    Extract single-end FASTQ files from prefetched SRA data
    """
    wildcard_constraints:
        sample_name="|".join(re.escape(s) for s in config.get("geo_samples_single", {}).keys()),
    input:
        sra=lambda wc: "{outdir}/sra_cache/{srr}/{srr}.sra".format(
            outdir=config["geo_outdir"],
            srr=config.get("geo_samples_single", {}).get(wc.sample_name, {}).get("srr", "INVALID"),
        ),
    output:
        r1=config["geo_outdir"] + "/{sample_name}.fastq",
    params:
        srr=lambda wc: config.get("geo_samples_single", {}).get(wc.sample_name, {}).get("srr", ""),
        outdir=config["geo_outdir"],
    threads: 8
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    message:
        "Extracting single-end FASTQ files for {wildcards.sample_name} ({params.srr})"
    log:
        OUT_DIR + "/logs/geo_download/geo_fastq_dump/{sample_name}.log"
    benchmark:
        OUT_DIR + "/.benchmark/geo_fastq_dump/{sample_name}.benchmark.tsv"
    container:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    shell:
        """
        exec &> {log}

        echo "Extracting FASTQ files from {params.srr}"
        fasterq-dump {input.sra} --skip-technical --outdir {params.outdir} -e {threads}

        echo "Renaming file"
        mv {params.outdir}/{params.srr}.fastq {output.r1}

        echo "Done extracting {wildcards.sample_name}"
        """

rule geo_compress_fastq_files:
    input:
        config["geo_outdir"] + "/{filename}.fastq"
    output:
        config["geo_outdir"] + "/{filename}.fastq.gz"
    threads: 32
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    message:
        "Compressing {wildcards.filename}.fastq"
    log:
        OUT_DIR + "/logs/geo_download/geo_compress/{filename}.log"
    benchmark:
        OUT_DIR + "/.benchmark/geo_compress/{filename}.benchmark.tsv"
    shell:
        """
        exec &> {log}
        if command -v pigz >/dev/null 2>&1; then
            pigz -p {threads} {input}
        else
            gzip -f {input}
        fi
        """

rule geo_download_all:
    """
    Download all GEO samples specified in metadata file
    """
    message:
        "Downloading all GEO samples"
    input:
        paired=lambda wildcards: expand(
            "{outdir}/{sample}_R1.fastq.gz",
            outdir=config.get("geo_outdir", "geo_data"),
            sample=config.get("geo_samples_paired", {}).keys()
        ) + expand(
            "{outdir}/{sample}_R2.fastq.gz",
            outdir=config.get("geo_outdir", "geo_data"),
            sample=config.get("geo_samples_paired", {}).keys()
        ),
        single=lambda wildcards: expand(
            "{outdir}/{sample}.fastq.gz",
            outdir=config.get("geo_outdir", "geo_data"),
            sample=config.get("geo_samples_single", {}).keys()
        )
    output:
        touch(OUT_DIR + "/logs/geo_download/geo_download_complete.txt")
    log:
        OUT_DIR + "/logs/geo_download/geo_download_all.log"
    benchmark:
        OUT_DIR + "/.benchmark/geo_download_all.benchmark.tsv"
    shell:
        """
        touch {output}
        echo "GEO download complete." > {log}
        """

ruleorder: geo_fastq_dump_paired > geo_fastq_dump_single
