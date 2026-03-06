from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

SCALE_RESOURCES = 1

MULTI_BAM_FILES = multiomics_builder.dataset_bam_files
MULTI_VCF_FILES = multiomics_builder.dataset_vcf_files
MULTI_BEDGRAPH_FILES = multiomics_builder.dataset_bedgraph_files
CHROMOSOME_SIZES = CONFIGS_PER_ASSAY[ASSAYS[0]].genome.chromosome_sizes

rule multiomics_make_dataset:
    """Create a dataset from bam files using QuantNado."""
    input:
        summary=OUTPUT_DIR + "/multiomics_summary.txt",
    output:
        dataset=directory(OUTPUT_DIR + "/multiomics/dataset.zarr"),
    params:
        chromosome_sizes=CHROMOSOME_SIZES,
        dataset=OUTPUT_DIR + "/multiomics/dataset.zarr",
        bam_files=MULTI_BAM_FILES,
        vcf_files=MULTI_VCF_FILES,
        bedgraph_files=MULTI_BEDGRAPH_FILES,
    threads: 8
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "docker://ghcr.io/milne-group/quantnado-ci:latest",
    log: OUTPUT_DIR + "/logs/multiomics/dataset/quantnado.log",
    benchmark: OUTPUT_DIR + "/.benchmark/multiomics/dataset/quantnado.tsv",
    message: "Making dataset from bam files using QuantNado."
    run:
        # Build command with repeated flags for each file
        cmd_parts = ["quantnado create-dataset", f"--output {params.dataset}"]
        
        # Add --bam for each BAM file
        for bam_file in params.bam_files:
            cmd_parts.append(f"--bam {bam_file}")
        
        # Add --bedgraph for each BEDGRAPH file
        if params.bedgraph_files:
            for bedgraph_file in params.bedgraph_files:
                cmd_parts.append(f"--bedgraph {bedgraph_file}")
        
        # Add --vcf for each VCF file
        if params.vcf_files:
            for vcf_file in params.vcf_files:
                cmd_parts.append(f"--vcf {vcf_file}")
        
        # Add remaining options
        cmd_parts.extend([
            f"--chromsizes {params.chromosome_sizes}",
            f"--max-workers {threads}",
            f"--log-file {log[0]}",
            "--verbose",
            "--overwrite",
            "--resume"
        ])
        
        cmd = " ".join(cmd_parts)
        shell(cmd)
