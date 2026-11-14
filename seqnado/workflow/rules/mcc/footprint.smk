from seqnado.workflow.helpers.mcc import get_mcc_footprint_peaks

_MCC_FOOTPRINT_TOOL = CONFIG.third_party_tools.mcc_footprint


# Base-pair ligation-junction footprints (`mcc.footprint: true`) from the grouped
# MCCNado pairs. Reporter regions are the peaks of the first MCC peak calling
# method; without peak calling the capture oligos are profiled. SeqNado's
# exclusion BED is deliberately not passed: it covers the viewpoints and would
# drop the capture-side junctions being profiled.
rule mcc_footprint:
    input:
        pairs=OUTPUT_DIR + "/mcc/{group}/ligation_junctions/{viewpoint_group}.pairs.gz",
        viewpoints=str(CONFIG.assay_config.mcc.viewpoints),
        chrom_sizes=str(CONFIG.genome.chromosome_sizes),
        peaks=lambda wc: get_mcc_footprint_peaks(wc, config=CONFIG, OUTPUT_DIR=OUTPUT_DIR),
    output:
        manifest=OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}.manifest.tsv",
        profiles=OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}.profiles.tsv.gz",
        contacts=OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}.contacts.tsv",
        regions=OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}.regions.tsv",
        metrics=OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}.metrics.json",
    params:
        prefix=OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}",
        peaks=lambda wc, input: " ".join(f"--peaks '{p}'" for p in input.peaks),
        options=str(_MCC_FOOTPRINT_TOOL.command_line_arguments) if _MCC_FOOTPRINT_TOOL else "",
    threads: _MCC_FOOTPRINT_TOOL.threads if _MCC_FOOTPRINT_TOOL else 1
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names),
        viewpoint_group="|".join(VIEWPOINT_TO_GROUPED_VIEWPOINT.values()),
    container:
        CONTAINERS.LANCEOTRON_MCC
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        # Single-threaded pass over the pairs file: ~16 min for a 92 MB pairs.gz.
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    log:
        OUTPUT_DIR + "/logs/mcc_footprint/{group}_{viewpoint_group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/mcc_footprint/{group}_{viewpoint_group}.tsv"
    message:
        "Profiling ligation-junction footprints for group {wildcards.group} and viewpoint group {wildcards.viewpoint_group}"
    shell:
        """
        printf 'sample\\tviewpoint\\tpairs\\n%s\\t%s\\t%s\\n' \\
            {wildcards.group:q} {wildcards.viewpoint_group:q} "$(realpath {input.pairs:q})" \\
            > {output.manifest:q}
        mcc-peaks footprint \\
            --manifest {output.manifest:q} \\
            --viewpoint {wildcards.viewpoint_group:q} \\
            --viewpoint-bed {input.viewpoints:q} \\
            --chrom-sizes {input.chrom_sizes:q} \\
            --output {params.prefix:q} \\
            {params.peaks} {params.options} > {log} 2>&1
        """


rule confirm_footprints_generated:
    input:
        expand(
            OUTPUT_DIR + "/mcc/footprints/{group}/{viewpoint_group}.metrics.json",
            group=SAMPLE_GROUPINGS.get_grouping("consensus").group_names,
            viewpoint_group=VIEWPOINT_TO_GROUPED_VIEWPOINT.values(),
        ),
    output:
        touch(temp(OUTPUT_DIR + "/mcc/footprints/.mcc_footprints_generated.txt")),
    message:
        "Confirming MCC junction footprints have been generated"
    log:
        OUTPUT_DIR + "/logs/mcc/mcc_footprints_generated.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/mcc/mcc_footprints_generated.tsv",
    shell:
        """
        echo "MCC junction footprinting completed successfully." > {log}
        """
