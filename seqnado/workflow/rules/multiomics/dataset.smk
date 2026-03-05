from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

SCALE_RESOURCES = 1
BAM_FILES = multiomics_builder.dataset_bam_files

rule make_dataset:
    """Create a dataset from bam files using QuantNado."""
    input:
        OUTPUT_DIR + "/multiomics_summary.txt",
        bam_files=BAM_FILES,
    output:
        dataset=OUTPUT_DIR + "/multiomics/dataset.zarr",
    params:
        chromosome_sizes=lambda wildcards: LOADED_CONFIGS[ASSAYS[0].clean_name]["genome"]["chromosome_sizes"],
        dataset=OUTPUT_DIR + "/multiomics/dataset",
    threads: 8
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest",
    log: OUTPUT_DIR + "/multiomics/logs/make_dataset.log",
    benchmark: OUTPUT_DIR + "/multiomics/.benchmark/make_dataset.tsv",
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