from seqnado.workflow.helpers.common import define_time_requested, define_memory_requested

SCALE_RESOURCES = 1

MULTI_BAM_FILES = multiomics_builder.dataset_bam_files
MULTI_VCF_FILES = multiomics_builder.dataset_vcf_files
MULTI_VCF_SAMPLE_NAMES = multiomics_builder.dataset_vcf_sample_names
MULTI_BEDGRAPH_FILES = multiomics_builder.dataset_bedgraph_files
MULTI_METHYLATION_SAMPLE_NAMES = multiomics_builder.dataset_methylation_sample_names
MULTI_STRANDED_CONFIG = multiomics_builder.dataset_stranded_config
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
        vcf_sample_names=MULTI_VCF_SAMPLE_NAMES,
        bedgraph_files=MULTI_BEDGRAPH_FILES,
        methylation_sample_names=MULTI_METHYLATION_SAMPLE_NAMES,
        stranded_config=MULTI_STRANDED_CONFIG,
    threads: 8
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=64, attempts=attempt, scale=SCALE_RESOURCES),
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
        
        # Add --methylation for each bedgraph file, with sample name overrides
        if params.bedgraph_files:
            for bedgraph_file in params.bedgraph_files:
                cmd_parts.append(f"--methylation {bedgraph_file}")
            if params.methylation_sample_names:
                cmd_parts.append(f"--methylation-sample-names {','.join(params.methylation_sample_names)}")

        # Add --vcf for each VCF file, with sample name overrides
        if params.vcf_files:
            for vcf_file in params.vcf_files:
                cmd_parts.append(f"--vcf {vcf_file}")
            if params.vcf_sample_names:
                cmd_parts.append(f"--vcf-sample-names {','.join(params.vcf_sample_names)}")
        
        # Add --stranded for RNA samples with non-zero strandedness
        if params.stranded_config:
            import json
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
