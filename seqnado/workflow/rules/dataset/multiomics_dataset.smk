from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

SCALE_RESOURCES = 1

BAM_FILES = multiomics_builder.dataset_bam_files
BAM_SAMPLE_NAMES = multiomics_builder.dataset_bam_sample_names
VCF_FILES = multiomics_builder.dataset_vcf_files
VCF_SAMPLE_NAMES = multiomics_builder.dataset_vcf_sample_names
METHYLATION_FILES = multiomics_builder.dataset_bedgraph_files
METHYLATION_SAMPLE_NAMES = multiomics_builder.dataset_methylation_sample_names
STRANDED_CONFIG = multiomics_builder.dataset_stranded_config
CHROMOSOME_SIZES = CONFIGS_PER_ASSAY[ASSAYS[0]].genome.chromosome_sizes

rule multiomics_make_dataset:
    """Create a dataset from bam files using QuantNado."""
    input:
        bam_files=BAM_FILES,
        vcf_files=VCF_FILES,
        methylation_files=METHYLATION_FILES,
    output:
        dataset=directory(OUTPUT_DIR + "/multiomics/dataset.zarr"),
    params:
        dataset=OUTPUT_DIR + "/multiomics/dataset.zarr",
        bam_sample_names=BAM_SAMPLE_NAMES,
        vcf_sample_names=VCF_SAMPLE_NAMES,
        methylation_sample_names=METHYLATION_SAMPLE_NAMES,
        stranded_config=STRANDED_CONFIG,
        chromosome_sizes=CHROMOSOME_SIZES,
    threads: 8
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=64, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest",
    log: OUTPUT_DIR + "/multiomics/logs/dataset/quantnado.log",
    benchmark: OUTPUT_DIR + "/.benchmark/multiomics/dataset/quantnado.tsv",
    message: "Making dataset from bam files using QuantNado."
    run:
        cmd_basic = [
            "quantnado create-dataset",
            f"--output {params.dataset}",
            f"--chromsizes {params.chromosome_sizes}",
            f"--max-workers {threads}",
            f"--log-file {log}",
            "--chr-workers 1",
            "--local-staging",
            "--staging-dir $TMPDIR",
            "--construction-compression fast",
            "--verbose",
            "--overwrite",
        ]

        # Add BAM files (quantnado --bam takes a single comma-separated list)
        if input.bam_files:
            cmd_basic.append(f"--bam {','.join(input.bam_files)}")
        if params.bam_sample_names:
            cmd_basic.append(f"--bam-sample-names {','.join(params.bam_sample_names)}")

        # Add Methylation files
        if input.methylation_files:
            cmd_basic.append(f"--methylation {','.join(input.methylation_files)}")
            if params.methylation_sample_names:
                cmd_basic.append(f"--methylation-sample-names {','.join(params.methylation_sample_names)}")

        # Add VCF files
        if input.vcf_files:
            cmd_basic.append(f"--vcf {','.join(input.vcf_files)}")
            if params.vcf_sample_names:
                cmd_basic.append(f"--vcf-sample-names {','.join(params.vcf_sample_names)}")

        # Add --stranded for RNA samples with non-zero strandedness
        if params.stranded_config:
            import json
            json_str = json.dumps(params.stranded_config).replace("{", "{{").replace("}", "}}")
            cmd_basic.append(f"--stranded '{json_str}'")

        cmd = " ".join(cmd_basic)
        shell(cmd)
