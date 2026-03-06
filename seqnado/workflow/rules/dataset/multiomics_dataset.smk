from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

SCALE_RESOURCES = 1

MULTI_BAM_FILES = multiomics_builder.dataset_bam_files
CHROMOSOME_SIZES = CONFIGS_PER_ASSAY[ASSAYS[0]].genome.chromosome_sizes

rule multiomics_make_dataset:
    """Create a dataset from bam files using QuantNado."""
    input:
        bam_files=MULTI_BAM_FILES,
    output:
        dataset=directory(OUTPUT_DIR + "/multiomics/dataset.zarr"),
    params:
        chromosome_sizes=CHROMOSOME_SIZES,
        dataset=OUTPUT_DIR + "/multiomics/dataset",
    threads: 8
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest",
    log: OUTPUT_DIR + "/logs/multiomics/dataset/quantnado.log",
    benchmark: OUTPUT_DIR + "/.benchmark/multiomics/dataset/quantnado.tsv",
    message: "Making dataset from bam files using QuantNado."
    shell: """
    quantnado create-dataset  \
    --output {params.dataset} \
    --chromsizes {params.chromosome_sizes} \
    --max-workers {threads} \
    --log-file {log} \
    --verbose \
    --overwrite \
    --resume \
    {input.bam_files}
    """
