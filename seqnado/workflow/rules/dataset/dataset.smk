from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

SCALE_RESOURCES = 1

BAM_FILES = OUTPUT.bam_files
VCF_FILES = OUTPUT.vcf_files
BEDGRAPH_FILES = OUTPUT.bedgraph_files

CHROMOSOME_SIZES = CONFIG.genome.chromosome_sizes

rule make_dataset:
    """Create a dataset from bam files using QuantNado."""
    input:
        bam_files=BAM_FILES,
        bedgraph_files=BEDGRAPH_FILES,
        vcf_files=VCF_FILES,
    output:
        dataset=directory(OUTPUT_DIR + "/dataset.zarr"),
    params:
        chromosome_sizes=CHROMOSOME_SIZES,
        dataset=OUTPUT_DIR + "/dataset",
    threads: 16
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest",
    log: OUTPUT_DIR + "/logs/dataset/quantnado.log",
    benchmark: OUTPUT_DIR + "/.benchmark/dataset/quantnado.tsv",
    message: "Making dataset from bam files using QuantNado."
    shell: """
    quantnado create-dataset  \
    --output {params.dataset} \
    --bam {input.bam_files} \
    --bedgraph {input.bedgraph_files} \
    --vcf {input.vcf_files} \
    --chromsizes {params.chromosome_sizes} \
    --max-workers {threads} \
    --log-file {log} \
    --verbose \
    --overwrite \
    --resume \
    {input.bam_files}
    """
    
                                                                                                                                                                                                                                                                                                        