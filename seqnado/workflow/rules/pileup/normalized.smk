from seqnado.workflow.helpers.common import (
    define_memory_requested,
    define_time_requested,
    format_deeptools_options,
    get_group_for_sample,
)
from seqnado.workflow.helpers.normalization import (
    get_norm_factor_spikein,
    get_scaling_factor,
)


rule deeptools_make_bigwigs_scale:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/csaw/{sample}.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
        ),
        options=lambda wc: format_deeptools_options(wc, str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments), INPUT_FILES, SAMPLE_GROUPINGS, raw_counts=True),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/scaled/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/scaled/{sample}.tsv",
    message: "Making scaled bigWig with deeptools for sample {wildcards.sample}"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} > {log} 2>&1
        """


use rule deeptools_make_bigwigs_scale as deeptools_make_bigwigs_spikein with:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{spikein_method}/{sample}.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options(wildcards, str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments), INPUT_FILES, SAMPLE_GROUPINGS, raw_counts=True),
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=False),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{spikein_method}/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{spikein_method}/{sample}.tsv",
    message: "Making spike-in normalized bigWig with deeptools for sample {wildcards.sample} using {wildcards.spikein_method}"


rule deeptools_make_bigwigs_rna_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{spikein_method}/{sample}_plus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options(wildcards, str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments), INPUT_FILES, SAMPLE_GROUPINGS, raw_counts=True),
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=False),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{spikein_method}/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{spikein_method}/{sample}_plus.tsv",
    message: "Making plus strand spike-in normalized bigWig with deeptools for sample {wildcards.sample} using {wildcards.spikein_method}"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand forward > {log} 2>&1
        """


rule deeptools_make_bigwigs_rna_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{spikein_method}/{sample}_minus.bigWig",
    params:
        options=lambda wildcards: format_deeptools_options(wildcards, str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments), INPUT_FILES, SAMPLE_GROUPINGS, raw_counts=True),
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=True),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{spikein_method}/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/pileups/deeptools/spikein/{spikein_method}/{sample}_minus.tsv",
    message: "Making minus strand spike-in normalized bigWig with deeptools for sample {wildcards.sample} using {wildcards.spikein_method}"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} -p {threads} --scaleFactor {params.scale} {params.options} --filterRNAstrand reverse > {log} 2>&1
        """