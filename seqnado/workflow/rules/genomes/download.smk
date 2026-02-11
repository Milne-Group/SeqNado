################################
# Genome Download Rules
# Downloads FASTA, GTF, chrom.sizes, and blacklist files from UCSC/ENCODE
################################

BLACKLIST_BASE_URL = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists"
PATH = "http://hgdownload.soe.ucsc.edu/goldenPath/"

rule download_fasta:
    output:
        fasta=OUTPUT_DIR + "/{genome}/sequence/{genome}.fa",
    params:
        url=lambda wc: f"{PATH}{wc.genome}/bigZips/{wc.genome}.fa.gz",
        gz=lambda wc: OUTPUT_DIR + f"/{wc.genome}/sequence/{wc.genome}.fa.gz",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_fasta.log",
    message:
        "Downloading FASTA for {wildcards.genome} from UCSC"
    shell:
        """
        wget --no-hsts -q {params.url} -O {params.gz} 2>&1 | tee {log} &&
        gunzip {params.gz} 2>&1 | tee -a {log}
        """


rule download_chrom_sizes:
    output:
        chrom_sizes=OUTPUT_DIR + "/{genome}/sequence/{genome}.chrom.sizes",
    params:
        url=lambda wc: f"{PATH}{wc.genome}/bigZips/{wc.genome}.chrom.sizes",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_chrom_sizes.log",
    message:
        "Downloading chromosome sizes for {wildcards.genome} from UCSC"
    shell:
        """
        wget --no-hsts -q {params.url} -O {output.chrom_sizes} 2>&1 | tee {log}
        """


rule download_gtf:
    output:
        gtf=OUTPUT_DIR + "/{genome}/genes/{genome}.ncbiRefSeq.gtf",
    params:
        url=lambda wc: f"{PATH}{wc.genome}/bigZips/genes/{wc.genome}.ncbiRefSeq.gtf.gz",
        gz=lambda wc: OUTPUT_DIR + f"/{wc.genome}/genes/{wc.genome}.ncbiRefSeq.gtf.gz",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_gtf.log",
    message:
        "Downloading GTF for {wildcards.genome} from UCSC"
    shell:
        """
        wget --no-hsts -q {params.url} -O {params.gz} 2>&1 | tee {log} &&
        gunzip {params.gz} 2>&1 | tee -a {log}
        """


rule download_blacklist:
    output:
        blacklist=OUTPUT_DIR + "/{genome}/{genome}-blacklist.bed.gz",
    params:
        url=lambda wc: f"{BLACKLIST_BASE_URL}/{wc.genome}-blacklist.v2.bed.gz",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_blacklist.log",
    message:
        "Downloading blacklist for {wildcards.genome}"
    shell:
        """
        wget --no-hsts -q {params.url} -O {output.blacklist} 2>{log} || {{
            echo "No blacklist available for {wildcards.genome}, creating empty placeholder" >> {log}
            rm -f {output.blacklist}
            touch {output.blacklist}
        }}
        """
