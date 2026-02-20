


rule call_peaks_lanceotron_no_input_consensus:
    input:
        bigwig=OUTPUT_DIR + "/bigwigs/deeptools/merged/unscaled/{group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/lanceotron/merged/{group}.bed",
        ltron_peaks=temp(OUTPUT_DIR + "/peaks/lanceotron/merged/{group}_L-tron.bed"),
    threads:
        CONFIG.third_party_tools.lanceotron.call_peaks.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=16, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=6, attempts=attempt, scale=SCALE_RESOURCES
        ),
    params:
        outdir=OUTPUT_DIR + "/peaks/lanceotron/merged",
        options=str(CONFIG.third_party_tools.lanceotron.call_peaks.command_line_arguments)
    container:
        "oras://ghcr.io/alsmith151/seqnado_ml_cpu:latest"
    log: OUTPUT_DIR + "/logs/lanceotron/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/lanceotron/merged/{group}.tsv",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    message: "Calling peaks with LanceOtron for merged group {wildcards.group}"
    shell: """
    lanceotron callPeaks {input.bigwig} -f {params.outdir} --skipheader  {params.options} > {log} 2>&1 &&
    cat {output.ltron_peaks} | cut -f 1-3 > {output.peaks}
    """

use rule call_peaks_macs2_no_input as call_peaks_macs2_no_input_consensus with:
    input:
        treatment = OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        peaks = OUTPUT_DIR + "/peaks/macs2/merged/{group}.bed",
    log: OUTPUT_DIR + "/logs/macs2/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/macs2/merged/{group}.tsv",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    message: "Calling peaks with MACS2 for merged group {wildcards.group}"


use rule call_peaks_macs3_no_input as call_peaks_macs3_no_input_consensus with:
    input:
        treatment = OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        peaks = OUTPUT_DIR + "/peaks/macs3/merged/{group}.bed",
    log: OUTPUT_DIR + "/logs/macs3/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/macs3/merged/{group}.tsv",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    message: "Calling peaks with MACS3 for merged group {wildcards.group}"


use rule call_peaks_homer_no_input as call_peaks_homer_no_input_consensus with:
    input:
        treatment = OUTPUT_DIR + "/tag_dirs/merged/unscaled/{group}",
    output:
        peaks = OUTPUT_DIR + "/peaks/homer/merged/{group}.bed",
    log: OUTPUT_DIR + "/logs/homer/merged/{group}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/homer/merged/{group}.tsv",
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    message: "Calling peaks with HOMER for merged group {wildcards.group}"


rule call_peaks_seacr_consensus:
    input:
        treatment=OUTPUT_DIR + "/bigwigs/deeptools/merged/unscaled/{group}.bigWig",
    output:
        peaks=OUTPUT_DIR + "/peaks/seacr/merged/{group}.bed",
        temp_peaks=temp(OUTPUT_DIR + "/peaks/seacr/merged/{group}.stringent.bed"),
    params:
        threshold=lambda wildcards: CONFIG.third_party_tools.seacr.threshold,
        normalization=lambda wildcards: CONFIG.third_party_tools.seacr.normalization,
        stringency=lambda wildcards: CONFIG.third_party_tools.seacr.stringency,
        basename=lambda wc, output: output.peaks.replace(".bed", ""),
    threads: 1
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
        OUTPUT_DIR + "/logs/seacr/merged/{group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/seacr/merged/{group}.tsv"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    message:
        "Calling peaks with SEACR for merged group {wildcards.group}"
    shell:
        """
    if ! SEACR.sh {input.treatment} {params.threshold} {params.normalization} {params.stringency} {params.basename} > {log} 2>&1; then
        touch {output.peaks}
    else
        awk 'BEGIN{{OFS="\\t"}} {{print $1, $2, $3}}' {output.temp_peaks} > {output.peaks} 2>> {log}
    fi
    """

