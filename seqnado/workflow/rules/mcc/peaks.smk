rule call_mcc_peaks: 
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/mcc/{group}_{viewpoint_group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron-mcc/{group}_{viewpoint_group}.bed",
    params:
        options=str(CONFIG.third_party_tools.lanceotron_mcc.command_line_arguments),
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
            apptainer exec --nv library://asmith151/lanceotron/lanceotron-mcc:latest \
                lanceotron-mcc call-mcc-peaks \
                --bigwig {input.bigwig} \
                --outfile {output.peaks} \
                {params.options} > {log} 2>&1 || exit $?
        else
            echo "INFO: apptainer not found — attempting to run lanceotron-mcc from current environment" >&2
            # Fall back to running the CLI directly (e.g., from a conda env). Fail with clear error if not available.
            if ! command -v lanceotron-mcc >/dev/null 2>&1; then
                echo "ERROR: lanceotron-mcc not found in PATH. Activate the conda environment that provides it or install it." >&2
                exit 127
            fi
            lanceotron-mcc call-mcc-peaks \
                --bigwig {input.bigwig} \
                --outfile {output.peaks} \
                {params.options} > {log} 2>&1 || exit $?
        fi
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
