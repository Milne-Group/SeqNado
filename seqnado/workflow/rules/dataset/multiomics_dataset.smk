import os

from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested

SCALE_RESOURCES = float(os.environ.get("SCALE_RESOURCES", "1"))


def _get_multiomics_samples():
    samples = []
    for output_files in OUTPUTS_PER_ASSAY.values():
        assay_samples = output_files.ip_sample_names or output_files.sample_names
        samples.extend(assay_samples)
    return sorted(set(samples))


MULTIOMICS_SAMPLES = _get_multiomics_samples()


def _get_multiomics_store_paths():
    stores = []
    for assay, output_files in OUTPUTS_PER_ASSAY.items():
        assay_samples = output_files.ip_sample_names or output_files.sample_names
        stores.extend(
            [
                f"{OUTPUT_DIR}/{assay.clean_name}/dataset/{sample_id}.zarr"
                for sample_id in assay_samples
            ]
        )
    return stores

rule combine_multiomics_dataset:
    input: 
        stores = lambda wildcards: _get_multiomics_store_paths(),
    output: 
        dataset=directory(OUTPUT_DIR + "/dataset/dataset.zarr"),
    threads: 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=32, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest"
    log: OUTPUT_DIR + "/logs/dataset/dataset_combine.log"
    benchmark: OUTPUT_DIR + "/.benchmark/dataset/dataset_combine.tsv"
    message: "Combining multi-omics datasets using QuantNado."
    shell: """
    quantnado dataset combine \
    --stores {input.stores} \
    --output {output.dataset} \
    --overwrite \
    --log-file {log}
    """


rule finalize_multiomics_dataset:
    input:
        OUTPUT_DIR + "/dataset/dataset.zarr"
    output:
        touch(OUTPUT_DIR + "/dataset/.complete")
    message:
        "Finalizing multiomics dataset at {output}."
