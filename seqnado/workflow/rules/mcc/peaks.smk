rule call_mcc_peaks:
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/{group}_{viewpoint_group}.bigWig",
        viewpoints=str(CONFIG.assay_config.mcc.viewpoints),
        exclusions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed",
        scores=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed.tsv",
        unresolved=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed.unresolved.bed",
        provenance=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed.json",
    params:
        options=str(CONFIG.third_party_tools.lanceotron_mcc.command_line_arguments),
    threads: CONFIG.third_party_tools.lanceotron_mcc.threads
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
        viewpoint_group="|".join(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    container:
        CONTAINERS.MCC_PEAKS
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=8, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=1, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/call_mcc_peaks/{group}_{viewpoint_group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/call_mcc_peaks/{group}_{viewpoint_group}.tsv"
    message:
        "Calling MCC peaks for group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}"
    shell:
        """
        mcc-peaks call-bigwig-peaks \
            --bigwigs {input.bigwig:q} \
            --viewpoint-bed {input.viewpoints:q} \
            --viewpoint {wildcards.viewpoint_group:q} \
            --exclude-bed {input.exclusions:q} \
            --threads {threads} \
            --output {output.peaks:q} \
            {params.options} > {log:q} 2>&1
        """


rule confirm_peaks_generated:
    input:
        expand(
            OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed",
            group=SAMPLE_GROUPINGS.get_grouping("consensus").group_names,
            viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values(),
        ),
    output:
        touch(temp(OUTPUT_DIR + "/peaks/mcc/.mcc_peaks_called.txt")),
    message:
        "Confirming MCC peaks have been called for all groups and viewpoint groups"
    log:
        OUTPUT_DIR + "/logs/mcc/mcc_peaks_called.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/mcc/mcc_peaks_called.tsv",
    shell:
        """
        echo "MCC peak calling completed successfully." > {output}
        echo "MCC peak calling completed successfully." > {log}
        """
