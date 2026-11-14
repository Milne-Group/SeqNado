from seqnado.workflow.helpers.mcc import (
    get_mcc_peak_methods,
    get_mcc_bigwigs_for_peak_calling,
)

MCC_PEAK_METHODS_SELECTED = get_mcc_peak_methods(CONFIG)
_MCC_PEAKS_TOOL = CONFIG.third_party_tools.mcc_peaks


# Legacy LanceOtron-derived heuristic (`peak_calling.method: lanceotron-mcc`).
# Whole-chromosome GPU calling on the grouped raw bigWig.
rule call_mcc_peaks:
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/{group}_{viewpoint_group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed",
    params:
        options=str(CONFIG.third_party_tools.lanceotron_mcc.command_line_arguments),
        image=CONTAINERS.LANCEOTRON_MCC,
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
        viewpoint_group="|".join(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    container:
        None
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=8, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=1, attempts=attempt, scale=SCALE_RESOURCES
        ),
        gpu=1,
    log:
        OUTPUT_DIR + "/logs/call_mcc_peaks/{group}_{viewpoint_group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/call_mcc_peaks/{group}_{viewpoint_group}.tsv"
    message:
        "Calling MCC peaks for group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}"
    shell:
        """
        # If apptainer is available, run inside the container (with NV GPU support).
        if command -v apptainer >/dev/null 2>&1; then
            echo "INFO: apptainer found — running inside container" >&2
            apptainer exec --nv {params.image} \
                mcc-peaks call-legacy-peaks \
                --bigwig {input.bigwig} \
                --output {output.peaks} \
                {params.options} > {log} 2>&1 || exit $?
        else
            echo "INFO: apptainer not found — attempting to run mcc-peaks from current environment" >&2
            # Fall back to running the CLI directly (e.g., from a conda env). Fail with clear error if not available.
            if ! command -v mcc-peaks >/dev/null 2>&1; then
                echo "ERROR: mcc-peaks not found in PATH. Activate the conda environment that provides it or install it." >&2
                exit 127
            fi
            mcc-peaks call-legacy-peaks \
                --bigwig {input.bigwig} \
                --output {output.peaks} \
                {params.options} > {log} 2>&1 || exit $?
        fi
        """


# Distance-aware caller (`peak_calling.method: mcc-peaks`). CPU only. Uses each
# replicate bigWig when `mcc.create_replicate_bigwigs` is set (majority support),
# otherwise the grouped raw bigWig.
rule call_mcc_peaks_bigwig:
    input:
        bigwigs=lambda wc: get_mcc_bigwigs_for_peak_calling(
            wc,
            SAMPLE_GROUPINGS=SAMPLE_GROUPINGS,
            OUTPUT_DIR=OUTPUT_DIR,
            use_replicates=CONFIG.assay_config.mcc.create_replicate_bigwigs,
        ),
        viewpoints=str(CONFIG.assay_config.mcc.viewpoints),
        exclusions=OUTPUT_DIR + "/resources/exclusion_regions.bed",
    output:
        peaks=OUTPUT_DIR + "/peaks/mcc-peaks/{group}_{viewpoint_group}.bed",
    params:
        options=str(_MCC_PEAKS_TOOL.command_line_arguments) if _MCC_PEAKS_TOOL else "",
    threads: _MCC_PEAKS_TOOL.threads if _MCC_PEAKS_TOOL else 4
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
        viewpoint_group="|".join(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    container:
        CONTAINERS.LANCEOTRON_MCC
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=1, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/call_mcc_peaks_bigwig/{group}_{viewpoint_group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/call_mcc_peaks_bigwig/{group}_{viewpoint_group}.tsv"
    message:
        "Calling distance-aware MCC peaks for group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}"
    shell:
        """
        mcc-peaks call-bigwig-peaks \
            --bigwigs {input.bigwigs:q} \
            --viewpoint-bed {input.viewpoints:q} \
            --viewpoint {wildcards.viewpoint_group:q} \
            --exclude-bed {input.exclusions:q} \
            --threads {threads} \
            --output {output.peaks:q} \
            {params.options} > {log} 2>&1
        """


rule confirm_peaks_generated:
    input:
        expand(
            OUTPUT_DIR + "/peaks/{method}/{group}_{viewpoint_group}.bed",
            method=MCC_PEAK_METHODS_SELECTED,
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
