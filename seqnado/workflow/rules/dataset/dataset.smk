from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested
import json

SCALE_RESOURCES = 1

BAM_FILES = OUTPUT.bam_files
VCF_FILES = OUTPUT.vcf_files
BEDGRAPH_FILES = OUTPUT.bedgraph_files

CHROMOSOME_SIZES = CONFIG.genome.chromosome_sizes

# VCF files are named {sample}.vcf.gz — Path.stem gives {sample}.vcf, not {sample}
DATASET_VCF_SAMPLE_NAMES = OUTPUT.sample_names if VCF_FILES else []

# Methylation files are named {sample}_{genome}_CpG.bedGraph — stem is not the sample name
DATASET_METHYLATION_SAMPLE_NAMES = OUTPUT.sample_names if BEDGRAPH_FILES else []

# Strandedness: only applies to RNA assay with non-zero strandedness
_STRAND_MAP = {1: "F", 2: "R"}
_strandedness = getattr(
    getattr(CONFIG.assay_config, "rna_quantification", None), "strandedness", 0
)
_library_type = _STRAND_MAP.get(_strandedness)
_bam_sample_names = OUTPUT.ip_sample_names or OUTPUT.sample_names
DATASET_STRANDED_CONFIG = (
    {s: _library_type for s in _bam_sample_names} if _library_type else {}
)


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
        dataset=OUTPUT_DIR + "/dataset.zarr",
        vcf_sample_names=DATASET_VCF_SAMPLE_NAMES,
        methylation_sample_names=DATASET_METHYLATION_SAMPLE_NAMES,
        stranded_config=DATASET_STRANDED_CONFIG,
    threads: 16
    resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest",
    log: OUTPUT_DIR + "/logs/dataset/quantnado.log",
    benchmark: OUTPUT_DIR + "/.benchmark/dataset/quantnado.tsv",
    message: "Making dataset from bam files using QuantNado."
    run:
        # Build command with repeated flags for each file
        cmd_parts = ["quantnado create-dataset", f"--output {params.dataset}"]

        # Add --bam for each BAM file
        for bam_file in input.bam_files:
            cmd_parts.append(f"--bam {bam_file}")

        # Add --methylation for each bedgraph file, with sample name overrides
        if input.bedgraph_files:
            for bedgraph_file in input.bedgraph_files:
                cmd_parts.append(f"--methylation {bedgraph_file}")
            if params.methylation_sample_names:
                cmd_parts.append(f"--methylation-sample-names {','.join(params.methylation_sample_names)}")

        # Add --vcf for each VCF file, with sample name overrides
        if input.vcf_files:
            for vcf_file in input.vcf_files:
                cmd_parts.append(f"--vcf {vcf_file}")
            if params.vcf_sample_names:
                cmd_parts.append(f"--vcf-sample-names {','.join(params.vcf_sample_names)}")

        # Add --stranded for RNA samples with non-zero strandedness
        if params.stranded_config:
            cmd_parts.append(f"--stranded '{json.dumps(params.stranded_config)}'")

        # Add remaining options
        cmd_parts.extend([
            f"--chromsizes {params.chromosome_sizes}",
            f"--max-workers {threads}",
            "--chr-workers 1",
            "--local-staging",
            "--staging-dir $TMPDIR",
            "--construction-compression fast",
            f"--log-file {log}",
            "--verbose",
            "--overwrite"
        ])

        cmd = " ".join(cmd_parts)
        shell(cmd)
