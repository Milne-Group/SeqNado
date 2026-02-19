from seqnado import SpikeInMethod


rule calculate_normalisation_factors_orlando:
    input:
        lambda wc: expand(
            rules.bam_split.output.stats,
            sample=SAMPLE_NAMES,
        ),
    output:
        normalisation_table=OUTPUT_DIR + "/resources/orlando/normalisation_factors.tsv",
        normalisation_factors=OUTPUT_DIR + "/resources/orlando/normalisation_factors.json",
    log: OUTPUT_DIR + "/logs/normalisation_factors.log",
    benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
    message: "Calculating normalisation factors using Orlando method"
    script:
        "../../scripts/calculate_spikein_norm_orlando.py"


rule calculate_normalisation_factors_input:
    input:
        lambda wc: expand(
            rules.bam_split.output.stats,
            sample=SAMPLE_NAMES,
        ),
        design=OUTPUT_DIR + "/metadata.csv",
    output:
        normalisation_table=OUTPUT_DIR + "/resources/with_input/normalisation_factors.tsv",
        normalisation_factors=OUTPUT_DIR + "/resources/with_input/normalisation_factors.json",
    log: OUTPUT_DIR + "/logs/normalisation_factors.log",
    benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
    message: "Calculating normalisation factors with input control"
    script:
        "../../scripts/calculate_spikein_norm_factors.py"


rule calculate_normalisation_factors_deseq2:
    input:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
        design=OUTPUT_DIR + "/metadata.csv",
    output:
        normalisation_table=OUTPUT_DIR + "/resources/deseq2/normalisation_factors.tsv",
        normalisation_factors=OUTPUT_DIR + "/resources/deseq2/normalisation_factors.json",
    params:
        spikein_genes=CONFIG.assay_config.spikein.control_genes or [],
        bam_dir=OUTPUT_DIR + "/aligned",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/normalisation_factors.log",
    benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
    message: "Calculating DESeq2 spike-in normalisation factors"
    script:
        "../../scripts/calculate_spikein_norm_factors_deseq.R"


rule calculate_normalisation_factors_edger:
    input:
        counts=OUTPUT_DIR + "/readcounts/feature_counts/read_counts.tsv",
        design=OUTPUT_DIR + "/metadata.csv",
    output:
        normalisation_table=OUTPUT_DIR + "/resources/edger/normalisation_factors.tsv",
        normalisation_factors=OUTPUT_DIR + "/resources/edger/normalisation_factors.json",
    params:
        spikein_genes=CONFIG.assay_config.spikein.control_genes or [],
        bam_dir=OUTPUT_DIR + "/aligned",
        outdir=OUTPUT_DIR + "/resources/edger",
        outfile="normalisation_factors.json",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    log: OUTPUT_DIR + "/logs/normalisation_factors.log",
    benchmark: OUTPUT_DIR + "/.benchmark/normalisation_factors.tsv",
    message: "Calculating edgeR spike-in normalisation factors"
    script:
        "../../scripts/calculate_spikein_norm_factors_edger.R"
