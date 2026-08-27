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
        custom=lambda wc: CUSTOM_FASTA.get(wc.genome, ""),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_fasta.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/download/{genome}_fasta.tsv",
    message:
        "Preparing FASTA for {wildcards.genome}"
    shell:
        """
        if [ -n "{params.custom}" ]; then
            case "{params.custom}" in
                *.gz) zcat "{params.custom}" > {output.fasta} 2> {log} ;;
                *) cp "{params.custom}" {output.fasta} 2> {log} ;;
            esac
        else
            if ! wget --no-hsts -q {params.url} -O {params.gz} 2> {log}; then
                {{
                    echo "ERROR: FASTA for genome '{wildcards.genome}' not found on UCSC.";
                    echo "  Tried: {params.url}";
                    echo "  If '{wildcards.genome}' is not a UCSC assembly, rerun with --fasta /path/to/your.fa (and --gtf if you need a STAR index).";
                }} | tee -a {log} >&2
                rm -f {params.gz}
                exit 1
            fi
            gunzip {params.gz} 2>&1 | tee -a {log}
        fi
        """


rule download_chrom_sizes:
    input:
        fai=lambda wc: (
            OUTPUT_DIR + f"/{wc.genome}/sequence/{wc.genome}.fa.fai"
            if wc.genome in CUSTOM_FASTA
            else []
        ),
    output:
        chrom_sizes=OUTPUT_DIR + "/{genome}/sequence/{genome}.chrom.sizes",
    params:
        url=lambda wc: f"{PATH}{wc.genome}/bigZips/{wc.genome}.chrom.sizes",
        is_custom=lambda wc: wc.genome in CUSTOM_FASTA,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_chrom_sizes.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/download/{genome}_chrom_sizes.tsv",
    message:
        "Preparing chromosome sizes for {wildcards.genome}"
    shell:
        """
        if [ "{params.is_custom}" = "True" ]; then
            cut -f1,2 {input.fai} > {output.chrom_sizes} 2> {log}
        else
            if ! wget --no-hsts -q {params.url} -O {output.chrom_sizes} 2> {log}; then
                {{
                    echo "ERROR: Chromosome sizes for genome '{wildcards.genome}' not found on UCSC.";
                    echo "  Tried: {params.url}";
                    echo "  If '{wildcards.genome}' is not a UCSC assembly, rerun with --fasta /path/to/your.fa (chrom.sizes is then derived automatically).";
                }} | tee -a {log} >&2
                rm -f {output.chrom_sizes}
                exit 1
            fi
        fi
        """


rule download_gtf:
    output:
        gtf=OUTPUT_DIR + "/{genome}/genes/{genome}.ncbiRefSeq.gtf",
    params:
        url=lambda wc: f"{PATH}{wc.genome}/bigZips/genes/{wc.genome}.ncbiRefSeq.gtf.gz",
        gz=lambda wc: OUTPUT_DIR + f"/{wc.genome}/genes/{wc.genome}.ncbiRefSeq.gtf.gz",
        custom=lambda wc: CUSTOM_GTF.get(wc.genome, ""),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    retries: 2
    log:
        OUTPUT_DIR + "/logs/download/{genome}_gtf.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/download/{genome}_gtf.tsv",
    message:
        "Preparing GTF for {wildcards.genome}"
    shell:
        """
        if [ -n "{params.custom}" ]; then
            case "{params.custom}" in
                *.gz) zcat "{params.custom}" > {output.gtf} 2> {log} ;;
                *) cp "{params.custom}" {output.gtf} 2> {log} ;;
            esac
        else
            if ! wget --no-hsts -q {params.url} -O {params.gz} 2> {log}; then
                {{
                    echo "ERROR: GTF for genome '{wildcards.genome}' not found on UCSC.";
                    echo "  Tried: {params.url}";
                    echo "  If '{wildcards.genome}' is not a UCSC assembly, rerun with --gtf /path/to/your.gtf (only needed if you want a STAR index).";
                }} | tee -a {log} >&2
                rm -f {params.gz}
                exit 1
            fi
            gunzip {params.gz} 2>&1 | tee -a {log}
        fi
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
    benchmark:
        OUTPUT_DIR + "/.benchmark/download/{genome}_blacklist.tsv",
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
