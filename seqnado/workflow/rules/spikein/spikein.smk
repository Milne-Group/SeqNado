
use rule bam_sort as bam_sort_spikein with:
    input:
        bam=OUTPUT_DIR + "/aligned/spikein/raw/{sample}.bam",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/sorted/{sample}.bam"),
        read_log=temp(OUTPUT_DIR + "/aligned/spikein/sorted/{sample}_read.log"),
    log: OUTPUT_DIR + "/logs/aligned_spikein/{sample}_sort.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned_spikein/{sample}_sort.tsv",
    message: "Sorting spike-in aligned BAM for sample {wildcards.sample} using samtools"


use rule bam_index as bam_index_spikein with:
    input:
        bam=rules.bam_sort_spikein.output.bam,
    output:
        bai=temp(OUTPUT_DIR + "/aligned/spikein/sorted/{sample}.bam.bai"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/aligned_spikein/{sample}_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned_spikein/{sample}_index.tsv",
    message: "Indexing spike-in aligned BAM for sample {wildcards.sample} using samtools"


rule filter_bam_spikein:
    input:
        bam=rules.bam_sort_spikein.output.bam,
    output:
        bam=temp(OUTPUT_DIR + "/aligned/spikein/filtered/{sample}.bam"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/aligned_spikein/{sample}_filter.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned_spikein/{sample}_filter.tsv",
    message: "Filtering spike-in aligned BAM for sample {wildcards.sample} using samtools",
    shell: """
    samtools view -b -F 260 -@ 8 {input.bam} > {output.bam} &&
    echo 'Filtered bam number of mapped reads:' > {log} 2>&1 &&
    samtools view -c {output.bam} >> {log} 2>&1
    """


use rule bam_index as bam_index_spikein_filtered with:
    input:
        bam=rules.filter_bam_spikein.output.bam,
    output:
        bai=temp(OUTPUT_DIR + "/aligned/spikein/filtered/{sample}.bam.bai"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/aligned_spikein/{sample}_filter_index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned_spikein/{sample}_filter_index.tsv",
    message: "Indexing filtered spike-in aligned BAM for sample {wildcards.sample} using samtools"


rule bam_split:
    input:
        bam=rules.filter_bam_spikein.output.bam,
        bai=rules.bam_index_spikein_filtered.output.bai,
    output:
        ref_bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_ref.bam"),
        exo_bam=temp(OUTPUT_DIR + "/aligned/spikein/{sample}_exo.bam"),
        stats=OUTPUT_DIR + "/aligned/spikein/{sample}_stats.tsv",
    params:
        genome_prefix=CONFIG.assay_config.spikein.endogenous_genome,
        exo_prefix=CONFIG.assay_config.spikein.exogenous_genome,
        prefix=OUTPUT_DIR + "/aligned/spikein/{sample}",
        map_qual=30,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bam_split/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bam_split/{sample}.tsv",
    message: "Splitting BAM into reference and spike-in for sample {wildcards.sample} using samtools",
    shell: """
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^chr/) print}}' | samtools view -b -o {output.ref_bam} &&
        samtools view -h {input.bam} | awk '{{if($0 ~ /^@/ || $3 ~ /^{params.exo_prefix}/) print}}' | samtools view -b -o {output.exo_bam} &&
        echo -e "sample\treference_reads\tspikein_reads" > {output.stats} &&
        echo -e "{wildcards.sample}\t$(samtools view -c {output.ref_bam})\t$(samtools view -c {output.exo_bam})" >> {output.stats}
        """


rule bam_move_ref:
    input:
        bam=rules.bam_split.output.ref_bam,
    output:
        bam=temp(OUTPUT_DIR + "/aligned/raw/{sample}.bam"),
    container: None
    log: OUTPUT_DIR + "/logs/bam_move_ref/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bam_move_ref/{sample}.tsv",
    message: "Moving reference BAM for sample {wildcards.sample} to raw aligned directory",
    shell: """
    mv {input.bam} {output.bam}
    """

