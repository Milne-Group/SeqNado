from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

# Merged, stranded bigwigs for RNA-seq

rule make_bigwigs_deeptools_rna_merged_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/unscaled/{group}_plus.bigWig",
    params:
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/{group}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/merged/{group}_plus.tsv",
    message: "Making plus strand merged bigWig with deeptools for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage {params.options} -p {threads} --filterRNAstrand forward -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule make_bigwigs_deeptools_rna_merged_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/unscaled/{group}_minus.bigWig",
    params:
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/{group}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/merged/{group}_minus.tsv",
    message: "Making minus strand merged bigWig with deeptools for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage {params.options} -p {threads} --filterRNAstrand reverse --scaleFactor -1 -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule make_bigwigs_bamnado_rna_merged_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/unscaled/{group}_plus.bigWig",
    params:
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/{group}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/merged/{group}_plus.tsv",
    message: "Making plus strand merged bigWig with bamnado for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward > {log} 2>&1
        """


rule make_bigwigs_bamnado_rna_merged_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/unscaled/{group}_minus.bigWig",
    params:
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/{group}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/merged/{group}_minus.tsv",
    message: "Making minus strand merged bigWig with bamnado for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse > {log} 2>&1
        """


rule make_bigwigs_homer_rna_merged:
    input:
        homer_tag_directory=OUTPUT_DIR + "/tag_dirs/merged/unscaled/{group}",
    output:
        bw_plus=OUTPUT_DIR + "/bigwigs/homer/merged/unscaled/{group}_plus.bigWig",
        bw_minus=OUTPUT_DIR + "/bigwigs/homer/merged/unscaled/{group}_minus.bigWig",
    params:
        genome_name=CONFIG.genome.name,
        genome_chrom_sizes=CONFIG.genome.chromosome_sizes,
        options=str(CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments),
        outdir=OUTPUT_DIR + "/bigwigs/homer/merged/unscaled/",
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
        OUTPUT_DIR + "/logs/homer/makebigwigs/merged/{group}_plus.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/makebigwigs/merged/{group}_plus.tsv"
    message:
        "Making plus strand merged bigWig with HOMER for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} -strand separate {params.options} > {log} 2>&1 &&
        mv {params.temp_bw_plus} {output.bw_plus} &&
        mv {params.temp_bw_minus} {output.bw_minus}
        """
