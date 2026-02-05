import itertools
from seqnado.workflow.helpers.mcc import (
    get_n_cis_scaling_factor,
    get_mcc_bam_files_for_merge,
)


rule make_bigwigs_mcc_replicates:
    input:
        bam=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam",
        bai=OUTPUT_DIR + "/mcc/replicates/{sample}/{sample}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
        cis_or_trans_stats=OUTPUT_DIR
        + "/resources/replicates/{sample}_ligation_stats.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/replicates/{sample}_{viewpoint_group}.bigWig",
    params:
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc, OUTPUT_DIR=OUTPUT_DIR),
    log:
        OUTPUT_DIR + "/logs/bigwig/{sample}_{viewpoint_group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/bigwig/{sample}_{viewpoint_group}.tsv",
    container:
        "docker://ghcr.io/alsmith151/bamnado:latest",
    message:
        "Generating bigWig for MCC sample {wildcards.sample} and viewpoint group {wildcards.viewpoint_group}",
    threads: 8,
    shell:
        """
    export RAYON_NUM_THREADS={threads} &&
    bamnado \
    bam-coverage \
    -b {input.bam} \
    -o {output.bigwig} \
    --scale-factor {params.scale_factor} \
    --blacklisted-locations {input.excluded_regions} \
    --min-mapq 0 \
    --filter-tag VP \
    --filter-tag-value {wildcards.viewpoint_group} \
    {params.options} > {log} 2>&1
    """


use rule make_bigwigs_mcc_replicates as make_bigwigs_mcc_grouped_norm with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
        bai=OUTPUT_DIR + "/mcc/{group}/{group}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
        cis_or_trans_stats=OUTPUT_DIR + "/resources/{group}_ligation_stats.json",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/n_cis/{group}_{viewpoint_group}.bigWig",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
        viewpoint_group="|".join(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    params:
        scale_factor=lambda wc: get_n_cis_scaling_factor(wc, OUTPUT_DIR=OUTPUT_DIR),
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    log:
        OUTPUT_DIR + "/logs/bigwig/{group}_{viewpoint_group}_n_cis.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/bigwig/{group}_{viewpoint_group}_n_cis.tsv",
    message:
        "Generating n_cis normalized bigWig for MCC group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}",
    container:
        "docker://ghcr.io/alsmith151/bamnado:latest",


use rule make_bigwigs_mcc_replicates as make_bigwigs_mcc_grouped_raw with:
    input:
        bam=OUTPUT_DIR + "/mcc/{group}/{group}.bam",
        bai=OUTPUT_DIR + "/mcc/{group}/{group}.bam.bai",
        excluded_regions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
    output:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/unscaled/{group}_{viewpoint_group}.bigWig",
    params:
        scale_factor=1,
        options=str(
            CONFIG.third_party_tools.bamnado.bam_coverage.command_line_arguments
        ),
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
        viewpoint_group="|".join(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    log:
        OUTPUT_DIR + "/logs/bigwig/{group}_{viewpoint_group}_unscaled.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/bigwig/{group}_{viewpoint_group}_unscaled.tsv",
    container:
        "docker://ghcr.io/alsmith151/bamnado:latest",
    message:
        "Generating unscaled bigWig for MCC group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}"


rule confirm_bigwigs_generated:
    input:
        replicates=expand(
            OUTPUT_DIR + "/bigwigs/mcc/replicates/{sample}_{viewpoint_group}.bigWig",
            sample=SAMPLE_NAMES,
            viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values(),
        ) if CONFIG.assay_config.mcc.create_replicate_bigwigs else [],
        n_cis_consensus=expand(
            OUTPUT_DIR + "/bigwigs/mcc/n_cis/{group}_{viewpoint_group}.bigWig",
            group=SAMPLE_GROUPINGS.get_grouping("consensus").group_names,
            viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values(),
        ),
        unscaled_consensus=expand(
            OUTPUT_DIR + "/bigwigs/mcc/unscaled/{group}_{viewpoint_group}.bigWig",
            group=SAMPLE_GROUPINGS.get_grouping("consensus").group_names,
            viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values(),
        ),
        subtractions=[
            str(Path(OUTPUT_DIR) / f"bigwigs/mcc/subtractions/{group1}_vs_{group2}_{viewpoint_group}.bigWig")
            for group1, group2 in itertools.product(
                SAMPLE_GROUPINGS.get_grouping("condition").group_names,
                SAMPLE_GROUPINGS.get_grouping("condition").group_names,
            )
            if group1 != group2
            for viewpoint_group in VIEWPOINT_TO_GROUPED_VIEWPOINT.values()
        ] if len(
            SAMPLE_GROUPINGS.get_grouping("condition").group_names
        ) >= 2 else [],
    output:
        touch(OUTPUT_DIR + "/bigwigs/mcc/.mcc_bigwigs_generated.txt"),
    message:
        "Confirming all MCC bigWigs have been generated"
