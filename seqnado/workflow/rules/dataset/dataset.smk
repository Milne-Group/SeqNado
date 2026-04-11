from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested

DATASET_SAMPLE_NAMES = OUTPUT.ip_sample_names or OUTPUT.sample_names

if ASSAY.clean_name == "snp":
    rule create_dataset_snp:
        input:
            vcf_file=OUTPUT_DIR + "/variant/{sample_id}.vcf.gz",
        output:
            dataset=temp(directory(OUTPUT_DIR + "/dataset/{sample_id}.zarr")),
        params:
            assay=ASSAY.value,
            output_dir=OUTPUT_DIR + "/dataset",
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(
                initial_value=32, attempts=attempt, scale=SCALE_RESOURCES
            ),
            runtime=lambda wildcards, attempt: define_time_requested(
                initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
            ),
        container: "docker://ghcr.io/milne-group/quantnado-ci:latest"
        log: OUTPUT_DIR + "/logs/dataset/{sample_id}.log"
        benchmark: OUTPUT_DIR + "/.benchmark/dataset/{sample_id}.tsv"
        message: "Creating dataset for sample {wildcards.sample_id} using QuantNado."
        shell: """
        quantnado dataset create \
        --sample {wildcards.sample_id} \
        --output-dir {params.output_dir} \
        --overwrite \
        --log-file {log} \
        --assay "{params.assay}" \
        --vcf-file {input.vcf_file}
        """


elif ASSAY.clean_name == "meth":
    METH_BDG_FILES = OUTPUT.meth_file_map
    rule create_dataset_meth:
        input:
            bam=OUTPUT_DIR + "/aligned/{sample_id}.bam",
            bai=OUTPUT_DIR + "/aligned/{sample_id}.bam.bai",
            bdg=lambda wc: METH_BDG_FILES[wc.sample_id]
        output:
            dataset=temp(directory(OUTPUT_DIR + "/dataset/{sample_id}.zarr")),
        params:
            assay=ASSAY.value,
            output_dir=OUTPUT_DIR + "/dataset",
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(
                initial_value=32, attempts=attempt, scale=SCALE_RESOURCES
            ),
            runtime=lambda wildcards, attempt: define_time_requested(
                initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
            ),
        container: "docker://ghcr.io/milne-group/quantnado-ci:latest"
        log: OUTPUT_DIR + "/logs/dataset/{sample_id}.log"
        benchmark: OUTPUT_DIR + "/.benchmark/dataset/{sample_id}.tsv"
        message: "Creating dataset for sample {wildcards.sample_id} using QuantNado."
        shell: """
        quantnado dataset create \
        --sample {wildcards.sample_id} \
        --output-dir {params.output_dir} \
        --overwrite \
        --log-file {log} \
        --assay "{params.assay}" \
        --bam-file {input.bam} \
        --methylation_file {input.bdg}
        """

elif ASSAY.clean_name == "chip" or ASSAY.clean_name == "cat":
    def _get_ip_for_sample(wildcards):
        sample = INPUT_FILES.query(wildcards.sample_id)
        return sample.ip

    rule create_dataset_ip:
        input:
            bam=OUTPUT_DIR + "/aligned/{sample_id}.bam",
            bai=OUTPUT_DIR + "/aligned/{sample_id}.bam.bai",
        output:
            dataset=temp(directory(OUTPUT_DIR + "/dataset/{sample_id}.zarr")),
        params:
            assay=ASSAY.value,
            output_dir=OUTPUT_DIR + "/dataset",
            ip=_get_ip_for_sample,
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(
                initial_value=32, attempts=attempt, scale=SCALE_RESOURCES
            ),
            runtime=lambda wildcards, attempt: define_time_requested(
                initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
            ),
        container: "docker://ghcr.io/milne-group/quantnado-ci:latest"
        log: OUTPUT_DIR + "/logs/dataset/{sample_id}.log"
        benchmark: OUTPUT_DIR + "/.benchmark/dataset/{sample_id}.tsv"
        message: "Creating dataset for sample {wildcards.sample_id} using QuantNado."
        shell: """
        echo "Using IP {params.ip} for sample {wildcards.sample_id}" >> {log}
        quantnado dataset create \
        --sample {wildcards.sample_id} \
        --output-dir {params.output_dir} \
        --overwrite \
        --log-file {log} \
        --assay "{params.assay}" \
        --bam-file {input.bam} \
        --ip {params.ip}
        """

else:
    rule create_dataset:
        input:
            bam=OUTPUT_DIR + "/aligned/{sample_id}.bam",
            bai=OUTPUT_DIR + "/aligned/{sample_id}.bam.bai",
        output:
            dataset=temp(directory(OUTPUT_DIR + "/dataset/{sample_id}.zarr")),
        params:
            assay=ASSAY.value,
            output_dir=OUTPUT_DIR + "/dataset",
        threads: 1
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(
                initial_value=32, attempts=attempt, scale=SCALE_RESOURCES
            ),
            runtime=lambda wildcards, attempt: define_time_requested(
                initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
            ),
        container: "docker://ghcr.io/milne-group/quantnado-ci:latest"
        log: OUTPUT_DIR + "/logs/dataset/{sample_id}.log"
        benchmark: OUTPUT_DIR + "/.benchmark/dataset/{sample_id}.tsv"
        message: "Creating dataset for sample {wildcards.sample_id} using QuantNado."
        shell: """
        quantnado dataset create \
        --sample {wildcards.sample_id} \
        --output-dir {params.output_dir} \
        --overwrite \
        --log-file {log} \
        --assay "{params.assay}" \
        --bam-file {input.bam}
        """


rule combine_dataset:
    input:
        stores=lambda wildcards: expand(
            OUTPUT_DIR + "/dataset/{sample_id}.zarr",
            sample_id=DATASET_SAMPLE_NAMES,
        )
    output:
        dataset=directory(OUTPUT_DIR + "/dataset.zarr")
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
    message: "Combining single-assay datasets using QuantNado."
    shell: """
    quantnado dataset combine \
    --stores {input.stores} \
    --output {output.dataset} \
    --overwrite \
    --log-file {log}
    """


rule finalize_dataset:
    input:
        OUTPUT_DIR + "/dataset.zarr"
    output:
        touch(OUTPUT_DIR + "/dataset.zarr.complete")
    log:
        OUTPUT_DIR + "/logs/dataset/dataset_finalize.log"
    benchmark:
        OUTPUT_DIR + "/.benchmark/dataset/dataset_finalize.tsv"
    message:
        "Finalizing assay dataset at {output}."
