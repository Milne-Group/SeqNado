from seqnado.workflow.helpers.common import (
    define_memory_requested,
    define_time_requested,
    format_deeptools_options,
)
from seqnado.workflow.helpers.normalization import (
    get_scaling_factor_for_merged_group,
    get_norm_factor_spikein_for_merged_group,
)


rule deeptools_make_bigwigs_merged_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/csaw/merged/{group}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        options=lambda wc: format_deeptools_options(
            wc,
            str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
            INPUT_FILES,
            SAMPLE_GROUPINGS,
            raw_counts=True,
        ),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/csaw/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/csaw/merged/{group}.tsv",
    message: "Making CSAW-scaled merged bigWig with deeptools for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1
        """


rule deeptools_make_bigwigs_merged_spikein:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{spikein_method}/merged/{group}.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=False),
        options=lambda wc: format_deeptools_options(
            wc,
            str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
            INPUT_FILES,
            SAMPLE_GROUPINGS,
            raw_counts=True,
        ),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{spikein_method}/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{spikein_method}/merged/{group}.tsv",
    message: "Making spike-in normalized merged bigWig with deeptools for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1
        """


rule bamnado_make_bigwigs_merged_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/csaw/merged/{group}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/csaw/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/bamnado/csaw/merged/{group}.tsv",
    message: "Making CSAW-scaled merged bigWig with bamnado for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} --scale-factor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule bamnado_make_bigwigs_merged_spikein:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/spikein/{spikein_method}/merged/{group}.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=False),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/spikein/{spikein_method}/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/bamnado/spikein/{spikein_method}/merged/{group}.tsv",
    message: "Making spike-in normalized merged bigWig with bamnado for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} --scale-factor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
        """


rule homer_make_bigwigs_merged_scale:
    input:
        homer_tag_directory=OUTPUT_DIR + "/tag_dirs/merged/{group}",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        homer_bigwig=OUTPUT_DIR + "/bigwigs/homer/csaw/merged/{group}.bigWig",
    params:
        genome_name=CONFIG.genome.name,
        genome_chrom_sizes=CONFIG.genome.chromosome_sizes,
        options=str(CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments),
        outdir=OUTPUT_DIR + "/bigwigs/homer/csaw/merged/",
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        temp_bw=lambda wc, output: output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/homer/csaw/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/homer/csaw/merged/{group}.tsv",
    message: "Making CSAW-scaled merged bigWig with HOMER for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} {params.options} -scale {params.scale} > {log} 2>&1 &&
        mv {params.temp_bw} {output.homer_bigwig}
        """


rule homer_make_bigwigs_merged_spikein:
    input:
        homer_tag_directory=OUTPUT_DIR + "/tag_dirs/merged/{group}",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        homer_bigwig=OUTPUT_DIR + "/bigwigs/homer/spikein/{spikein_method}/merged/{group}.bigWig",
    params:
        genome_name=CONFIG.genome.name,
        genome_chrom_sizes=CONFIG.genome.chromosome_sizes,
        options=str(CONFIG.third_party_tools.homer.make_bigwig.command_line_arguments),
        outdir=OUTPUT_DIR + "/bigwigs/homer/spikein/{spikein_method}/merged/",
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=False),
        temp_bw=lambda wc, output: output.homer_bigwig.replace(".bigWig", ".ucsc.bigWig"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/homer/spikein/{spikein_method}/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/homer/spikein/{spikein_method}/merged/{group}.tsv",
    message: "Making spike-in normalized merged bigWig with HOMER for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        makeBigWig.pl {input.homer_tag_directory} {params.genome_name} -chromSizes {params.genome_chrom_sizes} -url INSERT_URL -webdir {params.outdir} {params.options} -scale {params.scale} > {log} 2>&1 &&
        mv {params.temp_bw} {output.homer_bigwig}
        """
