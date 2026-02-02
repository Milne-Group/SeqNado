


rule viewpoints_to_fasta:
    input:
        bed=str(CONFIG.assay_config.mcc.viewpoints),
        genome=str(CONFIG.genome.fasta),
    output:
        fasta=OUTPUT_DIR + "/resources/viewpoints.fa",
    log: OUTPUT_DIR + "/logs/bed_to_fasta/viewpoints.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bed_to_fasta/viewpoints.tsv",
    message: "Converting viewpoints BED to FASTA",
    shell: """
    bedtools getfasta -fi {input.genome} -bed {input.bed} -fo {output.fasta} -name 2> {log} &&
    cat {output.fasta} | sed -E 's/:+/-/g' > {output.fasta}.tmp &&
    mv {output.fasta}.tmp {output.fasta}
    """


rule fasta_index:
    input:
        fasta=OUTPUT_DIR + "/resources/viewpoints.fa",
    output:
        index=OUTPUT_DIR + "/resources/viewpoints.fa.fai",
    log: OUTPUT_DIR + "/logs/bed_to_fasta/index.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bed_to_fasta/index.tsv",
    message: "Indexing viewpoints FASTA",
    shell: """
    samtools faidx {input.fasta} -o {output.index}
    """


rule exclusion_regions:
    input:
        bed=str(CONFIG.assay_config.mcc.viewpoints),
    output:
        bed=OUTPUT_DIR + "/resources/exclusion_regions.bed"
    params:
        genome=CONFIG.genome.chromosome_sizes,
        exclusion_zone=CONFIG.assay_config.mcc.exclusion_zone
    log: OUTPUT_DIR + "/logs/exclusion_regions.log"
    benchmark: OUTPUT_DIR + "/.benchmark/exclusion_regions.tsv",
    message: "Creating exclusion regions BED file",
    shell: """
    bedtools slop -i {input.bed} -g {params.genome}  -b {params.exclusion_zone} > {output.bed}
    """


rule minimap2_to_viewpoints:
    input:
        fq=OUTPUT_DIR + "/flashed/{sample}/{sample}.extendedFrags.fastq.gz",
        viewpoints=OUTPUT_DIR + "/resources/viewpoints.fa",
    output:
        bam=temp(OUTPUT_DIR + "/aligned/aligned_to_viewpoints/{sample}.bam"),
    threads: 4
    resources:
        mem="4GB",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/aligned/aligned_to_viewpoints/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/aligned/aligned_to_viewpoints/{sample}.tsv",
    message: "Aligning reads to viewpoints for sample {wildcards.sample}",
    shell: """
    minimap2 -x sr -a -k 8 -w 1 --cs=long {input.viewpoints} {input.fq} 2> {log} |
    samtools view -h -b -o {output.bam} 2>> {log} &&
    samtools sort -@ 4 -o {output.bam}.sorted {output.bam} 2>> {log} &&
    mv {output.bam}.sorted {output.bam} &&
    samtools index {output.bam}
    """ 
