include: 'common.smk'

rule deeptools_make_bigwigs_rna_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample}_plus.bigWig",
    params:
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/deeptools/unscaled/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/makebigwigs_{sample}_plus.tsv",
    message: "Making plus strand bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage {params.options} -p {threads} --filterRNAstrand forward -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule deeptools_make_bigwigs_rna_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/unscaled/{sample}_minus.bigWig",
    params:
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/deeptools/unscaled/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/makebigwigs_{sample}_minus.tsv",
    message: "Making minus strand bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} --filterRNAstrand reverse --scaleFactor -1 > {log} 2>&1
    """

rule homer_make_bigwigs_rna:
    input:
        homer_tag_directory=OUTPUT_DIR + "/tag_dirs/{sample}",
    output:
        bw_plus=OUTPUT_DIR + "/bigwigs/homer/unscaled/{sample}_plus.bigWig",
        bw_minus=OUTPUT_DIR + "/bigwigs/homer/unscaled/{sample}_minus.bigWig",
    params:
        genome_name=CONFIG.genome.name,
        genome_chrom_sizes=CONFIG.genome.chromosome_sizes,
        options=str(CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments),
        outdir=OUTPUT_DIR + "/bigwigs/homer/unscaled/",
        temp_bw_plus=lambda wc, output: output.bw_plus.replace(
            "_plus.bigWig", "pos.ucsc.bigWig"
        ),
        temp_bw_minus=lambda wc, output: output.bw_minus.replace(
            "_minus.bigWig", "neg.ucsc.bigWig"
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
    makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} -strand separate {params.options} > {log} 2>&1 &&
    mv {params.temp_bw_plus} {output.bw_plus} &&
    mv {params.temp_bw_minus} {output.bw_minus}
    """


rule bamnado_bam_coverage_rna_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/{sample}_plus.bigWig",
    params:
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/makebigwigs_{sample}_plus.tsv",
    message: "Making plus strand bigWig with bamnado for sample {wildcards.sample}"
    shell: """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward > {log} 2>&1
    """

rule bamnado_bam_coverage_rna_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/{sample}_minus.bigWig",
    params:
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources: 
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/makebigwigs_{sample}_minus.tsv",
    message: "Making minus strand bigWig with bamnado for sample {wildcards.sample}"
    shell: """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse --scale-factor -1 > {log} 2>&1
    """


ruleorder: deeptools_make_bigwigs_rna_plus > deeptools_make_bigwigs_rna_minus