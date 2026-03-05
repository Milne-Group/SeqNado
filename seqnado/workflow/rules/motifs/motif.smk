


rule get_fasta:
    input:
        peaks=OUTPUT_DIR + "/peaks/{method}/{sample}.bed",
    output:
        fasta=temp(OUTPUT_DIR + "/peaks/{method}/{sample}.fa"),
    params:
        genome=CONFIG.genome.fasta,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/motifs/meme/{method}/fasta/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/motifs/meme/{method}/fasta/{sample}.tsv",
    message: "Extracting FASTA sequences for peaks of sample {wildcards.sample}",
    shell: """
    bedtools getfasta -fullHeader -fi {params.genome} -bed {input.peaks} -fo {output.fasta}
    """


rule motif_meme_chip:
    input:
        fasta=OUTPUT_DIR + "/peaks/{method}/{sample}.fa",
    output:
        touch=temp(OUTPUT_DIR + "/motifs/meme/{method}/{sample}/.completed"),
    params:
        meme_dir=OUTPUT_DIR + "/motifs/meme/{method}/{sample}/",
        meme_chip_params=CONFIG.third_party_tools.meme.meme_chip.command_line_arguments.value if CONFIG.third_party_tools.meme else "",
        meme_chip_db=CONFIG.third_party_tools.meme.known_motif_db if CONFIG.third_party_tools.meme and CONFIG.third_party_tools.meme.known_motif_db else "",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads: 8
    container: "docker://quay.io/biocontainers/meme:5.5.9--pl5321h1ca524f_0"
    log: OUTPUT_DIR + "/logs/motifs/meme/{method}/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/motifs/meme/{method}/{sample}.tsv",
    message: "Running MEME-ChIP motif analysis for sample {wildcards.sample}"
    shell: """
    if [ -s {input.fasta} ] && grep -q "^>" {input.fasta}; then
        CMD="meme-chip -meme-p {threads} -oc {params.meme_dir}"
        [ -n "{params.meme_chip_db}" ] && CMD="$CMD -db {params.meme_chip_db}"
        $CMD{params.meme_chip_params} {input.fasta} > {log} 2>&1
    else
        echo "No sequences found in {input.fasta}. Skipping MEME-ChIP analysis." > {log}
        mkdir -p {params.meme_dir}
    fi
    touch {output.touch}
    """


rule motif_homer:
    input:
        peaks=OUTPUT_DIR + "/peaks/{method}/{sample}.bed",
    output:
        touch=temp(OUTPUT_DIR + "/motifs/homer/{method}/{sample}/.completed"),
    params:
        genome=CONFIG.genome.fasta,
        homer_dir=OUTPUT_DIR + "/motifs/homer/{method}/{sample}/",
        homer_params=CONFIG.third_party_tools.homer.find_motifs_genome.command_line_arguments.value,
        preparsed_dir="/tmp/homer_preparsed_{sample}_{method}",
        known_db=CONFIG.third_party_tools.homer.known_motif_db if CONFIG.third_party_tools.homer.known_motif_db else "",
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/motifs/homer/{method}/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/motifs/homer/{method}/{sample}.tsv",
    message: "Running HOMER motif analysis for sample {wildcards.sample}"
    shell: """
    if [ -s {input.peaks} ] && [ $(wc -l < {input.peaks}) -gt 0 ]; then
        mkdir -p {params.preparsed_dir}
        TMPDIR=$(mktemp -d)
        trap "rm -rf $TMPDIR" EXIT
        CMD="findMotifsGenome.pl {input.peaks} {params.genome} {params.homer_dir} -preparsedDir {params.preparsed_dir}"
        [ -n "{params.known_db}" ] && CMD="$CMD -mknown {params.known_db}"
        TMPDIR=$TMPDIR $CMD{params.homer_params} > {log} 2>&1
    else
        echo "No peaks found in {input.peaks}. Skipping HOMER analysis." > {log}
        mkdir -p {params.homer_dir}
    fi
    touch {output.touch}
    """


ruleorder: motif_meme_chip > motif_homer
