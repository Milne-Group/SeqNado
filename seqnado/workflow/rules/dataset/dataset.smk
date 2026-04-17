import re

from seqnado.workflow.helpers.common import define_memory_requested, define_time_requested

DATASET_SAMPLE_NAMES = OUTPUT.ip_sample_names or OUTPUT.sample_names
DATASET_PATH = f"{OUTPUT_DIR}/dataset/{CONFIG.project.date}_{CONFIG.project.name}.zarr"
DATASET_SAMPLE_PATTERN = "|".join(re.escape(sample) for sample in DATASET_SAMPLE_NAMES)


def _get_dataset_args(wildcards):
    assay = ASSAY.clean_name
    sample_id = wildcards.sample_id
    bam = OUTPUT_DIR + f"/aligned/{sample_id}.bam"
    bai = OUTPUT_DIR + f"/aligned/{sample_id}.bam.bai"

    if assay == "snp":
        return {
            "bam": [],
            "bai": [],
            "vcf": OUTPUT_DIR + f"/variant/{sample_id}.vcf.gz",
            "bdg": [],
            "primary_input_flag": f"--vcf-file {OUTPUT_DIR}/variant/{sample_id}.vcf.gz",
            "extra_args": "",
        }

    if assay in ("chip", "cat"):
        return {
            "bam": bam,
            "bai": bai,
            "vcf": [],
            "bdg": [],
            "primary_input_flag": f"--bam-file {bam}",
            "extra_args": f"--ip {INPUT_FILES.query(sample_id).ip}",
        }

    if assay == "rna":
        strandedness = (
            CONFIG.assay_config.rna_quantification.strandedness
            if CONFIG.assay_config.rna_quantification
            else 0
        )
        return {
            "bam": bam,
            "bai": bai,
            "vcf": [],
            "bdg": [],
            "primary_input_flag": f"--bam-file {bam}",
            "extra_args": f"--stranded {strandedness}",
        }

    if assay == "meth":
        return {
            "bam": bam,
            "bai": bai,
            "vcf": [],
            "bdg": OUTPUT.meth_file_map[sample_id],
            "primary_input_flag": f"--bam-file {bam}",
            "extra_args": f"--methylation_file {OUTPUT.meth_file_map[sample_id]}",
        }

    if assay == "atac":
        return {
            "bam": bam,
            "bai": bai,
            "vcf": [],
            "bdg": [],
            "primary_input_flag": f"--bam-file {bam}",
            "extra_args": "",
        }

    raise ValueError(f"Unsupported assay type: {assay}")


rule dataset_create:
    input:
        bam=lambda wc: _get_dataset_args(wc)["bam"],
        bai=lambda wc: _get_dataset_args(wc)["bai"],
        vcf=lambda wc: _get_dataset_args(wc)["vcf"],
        bdg=lambda wc: _get_dataset_args(wc)["bdg"],
    output:
        dataset=temp(directory(OUTPUT_DIR + "/dataset/{sample_id}.zarr")),
    wildcard_constraints:
        sample_id=DATASET_SAMPLE_PATTERN,
    params:
        assay=ASSAY.value,
        output_dir=OUTPUT_DIR + "/dataset",
        primary_input_flag=lambda wc: _get_dataset_args(wc)["primary_input_flag"],
        extra_args=lambda wc: _get_dataset_args(wc)["extra_args"],
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
    {params.primary_input_flag} \
    {params.extra_args}
    """


rule dataset_combine:
    input:
        stores=lambda wildcards: expand(
            OUTPUT_DIR + "/dataset/{sample_id}.zarr",
            sample_id=DATASET_SAMPLE_NAMES,
        )
    output:
        dataset=directory(DATASET_PATH)
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
