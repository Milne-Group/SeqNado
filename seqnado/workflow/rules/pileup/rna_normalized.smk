from seqnado.workflow.helpers.common import (
    define_time_requested,
    define_memory_requested,
    get_group_for_sample,
)
from seqnado.workflow.helpers.normalization import (
    get_norm_factor_spikein,
    get_scaling_factor,
)

# Spike-in and CSAW normalized, stranded bigwigs for individual RNA samples

rule make_bigwigs_deeptools_rna_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{spikein_method}/{sample}_plus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=False),
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{spikein_method}/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/spikein/{spikein_method}/{sample}_plus.tsv",
    message: "Making plus strand spike-in normalized bigWig with deeptools for sample {wildcards.sample} using {wildcards.spikein_method}"
    shell: """
    bamCoverage {params.options} -p {threads} --filterRNAstrand forward --scaleFactor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule make_bigwigs_deeptools_rna_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/spikein/{spikein_method}/{sample}_minus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=True),
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/deeptools/spikein/{spikein_method}/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/spikein/{spikein_method}/{sample}_minus.tsv",
    message: "Making minus strand spike-in normalized bigWig with deeptools for sample {wildcards.sample} using {wildcards.spikein_method}"
    shell: """
    bamCoverage {params.options} -p {threads} --filterRNAstrand reverse --scaleFactor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule make_bigwigs_bamnado_rna_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/spikein/{spikein_method}/{sample}_plus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=False),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/spikein/{spikein_method}/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/spikein/{spikein_method}/{sample}_plus.tsv",
    message: "Making plus strand spike-in normalized bigWig with bamnado for sample {wildcards.sample} using {wildcards.spikein_method}"
    shell: """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward --scale-factor {params.scale} > {log} 2>&1
    """


rule make_bigwigs_bamnado_rna_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/spikein/{spikein_method}/{sample}_minus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein(wc, OUTPUT_DIR, CONFIG, negative=True),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources: 
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/spikein/{spikein_method}/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/spikein/{spikein_method}/{sample}_minus.tsv",
    message: "Making minus strand spike-in normalized bigWig with bamnado for sample {wildcards.sample} using {wildcards.spikein_method}"
    shell: """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse > {log} 2>&1
    """


# CSAW (library-size) normalized, stranded bigwigs for individual RNA samples

rule make_bigwigs_deeptools_rna_csaw_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/csaw/{sample}_plus.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
        ),
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/deeptools/csaw/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/csaw/{sample}_plus.tsv",
    message: "Making plus strand CSAW-scaled bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage {params.options} -p {threads} --filterRNAstrand forward --scaleFactor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule make_bigwigs_deeptools_rna_csaw_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/csaw/{sample}_minus.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
        ),
        options=str(CONFIG.third_party_tools.deeptools.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.deeptools.bam_coverage.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/pileups/deeptools/csaw/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/deeptools/csaw/{sample}_minus.tsv",
    message: "Making minus strand CSAW-scaled bigWig with deeptools for sample {wildcards.sample}"
    shell: """
    bamCoverage {params.options} -p {threads} --filterRNAstrand reverse --scaleFactor {params.scale} -b {input.bam} -o {output.bigwig} > {log} 2>&1
    """


rule make_bigwigs_bamnado_rna_csaw_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/csaw/{sample}_plus.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
        ),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/csaw/{sample}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/csaw/{sample}_plus.tsv",
    message: "Making plus strand CSAW-scaled bigWig with bamnado for sample {wildcards.sample}"
    shell: """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward --scale-factor {params.scale} > {log} 2>&1
    """


rule make_bigwigs_bamnado_rna_csaw_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
        scaling_factors=lambda wc: OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/csaw/{sample}_minus.bigWig",
    params:
        scale=lambda wc: get_scaling_factor(
            wc,
            OUTPUT_DIR + f"/resources/{get_group_for_sample(wc, INPUT_FILES)}_scaling_factors.tsv",
        ),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/pileups/bamnado/csaw/{sample}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/csaw/{sample}_minus.tsv",
    message: "Making minus strand CSAW-scaled bigWig with bamnado for sample {wildcards.sample}"
    shell: """
    export RAYON_NUM_THREADS={threads}
    bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse --scale-factor {params.scale} > {log} 2>&1
    """

ruleorder: make_bigwigs_deeptools_rna_spikein_plus > make_bigwigs_deeptools_rna_spikein_minus
ruleorder: make_bigwigs_bamnado_rna_spikein_plus > make_bigwigs_bamnado_rna_spikein_minus
ruleorder: make_bigwigs_deeptools_rna_csaw_plus > make_bigwigs_deeptools_rna_csaw_minus
ruleorder: make_bigwigs_bamnado_rna_csaw_plus > make_bigwigs_bamnado_rna_csaw_minus
