from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
from seqnado.workflow.helpers.quant import get_bams_to_count


rule feature_counts:
    input:
        bam=get_bams_to_count(config=CONFIG, sample_names=SAMPLE_NAMES, output_dir=OUTPUT_DIR),
        bai=[bam.replace(".bam", ".bam.bai") for bam in get_bams_to_count(config=CONFIG, sample_names=SAMPLE_NAMES, output_dir=OUTPUT_DIR)],
        annotation=str(CONFIG.genome.gtf),
    output:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
    params:
        options=str(CONFIG.third_party_tools.subread.feature_counts.command_line_arguments),
        annotation_format= "SAF" if CONFIG.genome.gtf.suffix == ".saf" else "GTF",
        paired= "-p --countReadPairs" if INPUT_FILES.is_paired_end(SAMPLE_NAMES[0]) else "",
        stranded="-s " + str(getattr(CONFIG.assay_config.rna_quantification, "strandedness", 0)),
    threads: 
        CONFIG.third_party_tools.subread.feature_counts.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=3, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log: OUTPUT_DIR + "/logs/readcounts/featurecounts/featurecounts.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/featurecounts/featurecounts.tsv",
    message: "Running featureCounts to quantify reads for all samples"
    shell: """
    # Determine annotation format from file extension
    annotation="{input.annotation}"
    if [[ "$annotation" == *.saf ]]; then
        format_flag="-F SAF"
    else
        format_flag="-F GTF"
    fi
    
    featureCounts \
    $format_flag \
    -a {input.annotation} \
    -F {params.annotation_format} \
    -T {threads} \
    --donotsort \
    {params.paired} \
    {params.stranded} \
    {params.options} \
    -o {output.counts} \
    {input.bam} \
    > {log} 2>&1
    """
