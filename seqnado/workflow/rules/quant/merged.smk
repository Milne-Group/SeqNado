# Define rules for generating merged SAF files and counting reads using featureCounts
# Requires a rule to generate the merged peak files beforehand
from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
from seqnado.workflow.helpers.quant import get_bams_to_count_grouped, correct_featurecounts_options_for_peaks


rule merged_saf:
    input:
        peaks=OUTPUT_DIR + "/peaks/lanceotron/merged/{group}.bed",
    output:
        saf=temp(OUTPUT_DIR + "/readcounts/feature_counts/{group}.saf"),
    threads: 1,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/readcounts/feature_counts/{group}_saf.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/feature_counts/{group}_saf.tsv",
    message: "Generating SAF file from merged peaks for group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names)    
    shell: """
    awk 'BEGIN{{OFS="\t"}}{{print $1":"$2"-"$3,$1,$2,$3,"*"}}' {input.peaks} > {output.saf}
    """


rule merged_counts:
    input:
        bam=lambda wildcards: get_bams_to_count_grouped(CONFIG, SAMPLE_GROUPINGS, OUTPUT_DIR, group_name=wildcards.group),
        bai=lambda wildcards: expand("{bam}.bai", bam=get_bams_to_count_grouped(CONFIG, SAMPLE_GROUPINGS, OUTPUT_DIR, group_name=wildcards.group)),
        saf=rules.merged_saf.output.saf,
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/{group}_counts.tsv",
    params:
        options=lambda wc: correct_featurecounts_options_for_peaks(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments, input_files=INPUT_FILES, sample_groupings=SAMPLE_GROUPINGS, group_name=wc.group),
    threads: CONFIG.third_party_tools.subread.feature_counts.threads,
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/readcounts/feature_counts/{group}_counts.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/feature_counts/{group}_counts.tsv",
    message: "Running featureCounts to quantify reads for merged peaks in group {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping("consensus").group_names)
    shell: """
    featureCounts \
    -a {input.saf} \
    -F SAF \
    -T {threads} \
    --donotsort \
    {params.options} \
    -o {output.counts} \
    {input.bam} > {log} 2>&1
    """

localrules:
    merged_saf