from seqnado.workflow.helpers.common import (
    define_time_requested,
    define_memory_requested,
    format_deeptools_options,
    get_group_for_sample,
)
from seqnado.workflow.helpers.normalization import (
    get_norm_factor_spikein,
    get_norm_factor_spikein_for_merged_group,
    get_scaling_factor_for_merged_group,
)

rule make_bigwigs_deeptools_rna_merged_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}_plus.bigWig",
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
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}_plus.tsv",
    message: "Making plus strand spike-in normalized merged bigWig with deeptools for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} --filterRNAstrand forward > {log} 2>&1
        """


rule make_bigwigs_deeptools_rna_merged_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}_minus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=True),
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
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/spikein/{spikein_method}/{group}_minus.tsv",
    message: "Making minus strand spike-in normalized merged bigWig with deeptools for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} --filterRNAstrand reverse > {log} 2>&1
        """


rule make_bigwigs_bamnado_rna_merged_spikein_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}_plus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=False),
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/merged/spikein/{spikein_method}/{group}_plus.tsv",
    message: "Making plus strand spike-in normalized merged bigWig with bamnado for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward --scale-factor {params.scale} > {log} 2>&1
        """


rule make_bigwigs_bamnado_rna_merged_spikein_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{spikein_method}/normalisation_factors.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}_minus.bigWig",
    params:
        scale=lambda wc: get_norm_factor_spikein_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR, negative=True),
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=6, attempts=attempt, scale=SCALE_RESOURCES),
    threads:
        CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/spikein/{spikein_method}/{group}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/merged/spikein/{spikein_method}/{group}_minus.tsv",
    message: "Making minus strand spike-in normalized merged bigWig with bamnado for group {wildcards.group} using {wildcards.spikein_method}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse > {log} 2>&1
        """


# CSAW (library-size) normalized, stranded bigwigs for merged RNA samples

rule make_bigwigs_deeptools_rna_merged_csaw_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/csaw/{group}_plus.bigWig",
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
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/csaw/{group}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/csaw/{group}_plus.tsv",
    message: "Making plus strand CSAW-scaled merged bigWig with deeptools for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} --filterRNAstrand forward > {log} 2>&1
        """


rule make_bigwigs_deeptools_rna_merged_csaw_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/csaw/{group}_minus.bigWig",
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
    log: OUTPUT_DIR + "/logs/bigwigs/deeptools/merged/csaw/{group}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bigwigs/deeptools/merged/csaw/{group}_minus.tsv",
    message: "Making minus strand CSAW-scaled merged bigWig with deeptools for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bigwig} --scaleFactor {params.scale} -p {threads} {params.options} --filterRNAstrand reverse > {log} 2>&1
        """


rule make_bigwigs_bamnado_rna_merged_csaw_plus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/csaw/{group}_plus.bigWig",
    params:
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/csaw/{group}_plus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/merged/csaw/{group}_plus.tsv",
    message: "Making plus strand CSAW-scaled merged bigWig with bamnado for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand forward --scale-factor {params.scale} > {log} 2>&1
        """


rule make_bigwigs_bamnado_rna_merged_csaw_minus:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
        bai=OUTPUT_DIR + "/aligned/merged/{group}.bam.bai",
        scaling_factors=OUTPUT_DIR + "/resources/{group}_scaling_factors.tsv",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/bamnado/merged/csaw/{group}_minus.bigWig",
    params:
        scale=lambda wc: get_scaling_factor_for_merged_group(wc, SAMPLE_GROUPINGS, OUTPUT_DIR),
        options=str(CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments),
    threads: CONFIG.third_party_tools.bamnado.bam_coverage.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/alsmith151/bamnado:latest"
    log: OUTPUT_DIR + "/logs/bigwigs/bamnado/merged/csaw/{group}_minus.log",
    benchmark: OUTPUT_DIR + "/.benchmark/bamnado/merged/csaw/{group}_minus.tsv",
    message: "Making minus strand CSAW-scaled merged bigWig with bamnado for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        export RAYON_NUM_THREADS={threads}
        bamnado bam-coverage {params.options} -b {input.bam} -o {output.bigwig} --strand reverse --scale-factor {params.scale} > {log} 2>&1
        """

ruleorder: make_bigwigs_deeptools_rna_merged_plus > make_bigwigs_deeptools_rna_merged_minus
ruleorder: make_bigwigs_deeptools_rna_merged_csaw_plus > make_bigwigs_deeptools_rna_merged_csaw_minus
ruleorder: make_bigwigs_bamnado_rna_merged_csaw_plus > make_bigwigs_bamnado_rna_merged_csaw_minus

