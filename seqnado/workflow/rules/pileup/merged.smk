from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested


# Pileup for grouped sample


rule bamnado_make_bigwigs_merged:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/{group}.bigWig",
    params:
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/bamnado/merged/{group}.tsv",
    message: "Making bigWig with bamnado for merged sample {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule deeptools_make_bigwigs_merged:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/{group}.bigWig",
    params:
       options=lambda wildcards: format_deeptools_options(
           wildcards,
           str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
           INPUT_FILES,
           SAMPLE_GROUPINGS,
       ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.deeptools.bam_coverage.threads,
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/{group}.tsv",
    message: "Making bigWig with deeptools for merged sample {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage {params.options} -p {threads} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """
        

rule homer_make_tag_directory_merged:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        homer_tag_directory=directory(OUTPUT_DIR + "/tag_dirs/merged/{group}"),
    params:
        options=str(
            CONFIG.third_party_tools.homer.make_tag_directory.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/homer/maketagdirectory_merged_{group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/maketagdirectory_merged_{group}.tsv"
    message:
        "Making tag directory with HOMER for merged sample {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1
        """


rule homer_make_bigwigs_merged:
    input:
        homer_tag_directory=OUTPUT_DIR + "/tag_dirs/merged/{group}",
    output:
        homer_bigwig=OUTPUT_DIR + "/bigwigs/homer/merged/{group}.bigWig",
    params:
        genome_name=CONFIG.genome.name,
        genome_chrom_sizes=CONFIG.genome.chromosome_sizes,
        options=str(CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments),
        outdir=OUTPUT_DIR + "/bigwigs/homer/merged/",
        temp_bw=lambda wc, output: output.homer_bigwig.replace(
            ".bigWig", ".ucsc.bigWig"
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
        OUTPUT_DIR + "/logs/homer/makebigwigs_merged_{group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/makebigwigs_merged_{group}.tsv"
    message:
        "Making bigWig with HOMER for merged sample {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} {params.options} > {log} 2>&1 &&
        mv {params.temp_bw} {output.homer_bigwig}
        """
