"""
Rules for downloading GEO/SRA data
"""

import pandas as pd


# Load sample metadata from TSV file
def load_geo_samples(metadata_file):
    """Load GEO sample metadata from TSV file"""
    samples = pd.read_csv(metadata_file, sep="\t")
    samples_dict = {}
    for _, row in samples.iterrows():
        srr = row["run_accession"]
        sample_name = f"{row['library_name']}-{row['sample_title']}"
        samples_dict[sample_name] = {
            "srr": srr,
            "gsm": row["library_name"],
            "sample": row["sample_title"]
        }
    return samples_dict


# Prioritize paired rule over single rule to avoid ambiguity
ruleorder: geo_download_paired > geo_download_single


rule geo_download_paired:
    """
    Download paired-end FASTQ files from GEO/SRA with retry logic
    """
    output:
        r1=config["geo_outdir"] + "/{sample_name}_R1.fastq.gz",
        r2=config["geo_outdir"] + "/{sample_name}_R2.fastq.gz"
    params:
        srr=lambda wildcards: config["geo_samples_paired"][wildcards.sample_name]["srr"],
        outdir=config["geo_outdir"],
        max_retries=5
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240  # 4 hours in minutes
    log:
        "logs/geo_download/{sample_name}.log"
    container:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    shell:
        """
        exec &> {log}
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        echo "Downloading paired-end data: {params.srr}"
        
        # Retry prefetch up to {params.max_retries} times with increasing wait
        MAX_RETRIES={params.max_retries}
        for attempt in $(seq 1 $MAX_RETRIES); do
            echo "prefetch attempt $attempt of $MAX_RETRIES"
            prefetch {params.srr} --max-size 50G && break
            if [ $attempt -lt $MAX_RETRIES ]; then
                echo "prefetch failed, retrying in $((attempt * 60))s..."
                sleep $((attempt * 60))
            else
                echo "prefetch failed after $MAX_RETRIES attempts"
                exit 1
            fi
        done
        
        echo "Extracting FASTQ files"
        fasterq-dump {params.srr} --split-files --threads {threads}
        
        # Rename files
        mv {params.srr}_1.fastq {wildcards.sample_name}_R1.fastq
        mv {params.srr}_2.fastq {wildcards.sample_name}_R2.fastq
        
        # Compress
        echo "Compressing FASTQ files"
        pigz -p {threads} {wildcards.sample_name}_R1.fastq
        pigz -p {threads} {wildcards.sample_name}_R2.fastq
        
        # Remove SRA files to save space
        rm -rf {params.srr}.sra
        
        echo "Done downloading {wildcards.sample_name}"
        """


rule geo_download_single:
    """
    Download single-end FASTQ files from GEO/SRA with retry logic
    """
    output:
        r1=config["geo_outdir"] + "/{sample_name}.fastq.gz"
    params:
        srr=lambda wildcards: config["geo_samples_single"][wildcards.sample_name]["srr"],
        outdir=config["geo_outdir"],
        max_retries=5
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240  # 4 hours in minutes
    log:
        "logs/geo_download/{sample_name}.log"
    container:
        "docker://quay.io/biocontainers/sra-tools:3.0.10--h9f5acd7_0"
    shell:
        """
        exec &> {log}
        
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        echo "Downloading single-end data: {params.srr}"
        
        # Retry prefetch up to {params.max_retries} times with increasing wait
        MAX_RETRIES={params.max_retries}
        for attempt in $(seq 1 $MAX_RETRIES); do
            echo "prefetch attempt $attempt of $MAX_RETRIES"
            prefetch {params.srr} --max-size 50G && break
            if [ $attempt -lt $MAX_RETRIES ]; then
                echo "prefetch failed, retrying in $((attempt * 60))s..."
                sleep $((attempt * 60))
            else
                echo "prefetch failed after $MAX_RETRIES attempts"
                exit 1
            fi
        done
        
        echo "Extracting FASTQ file"
        fasterq-dump {params.srr} --threads {threads}
        
        # Rename file (handle both _1.fastq and .fastq outputs)
        if [ -f {params.srr}_1.fastq ]; then
            mv {params.srr}_1.fastq {wildcards.sample_name}.fastq
        elif [ -f {params.srr}.fastq ]; then
            mv {params.srr}.fastq {wildcards.sample_name}.fastq
        fi
        
        # Compress
        echo "Compressing FASTQ file"
        pigz -p {threads} {wildcards.sample_name}.fastq
        
        # Remove SRA files to save space
        rm -rf {params.srr}.sra
        
        echo "Done downloading {wildcards.sample_name}"
        """


rule geo_download_all:
    """
    Download all GEO samples specified in metadata file
    """
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
        touch("logs/geo_download_complete.txt")

