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
        dataset=OUTPUT_DIR + "/dataset.zarr",
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
        
        # Add --bedgraph for each BEDGRAPH file
        if input.bedgraph_files:
            for bedgraph_file in input.bedgraph_files:
                cmd_parts.append(f"--bedgraph {bedgraph_file}")
        
        # Add --vcf for each VCF file
        if input.vcf_files:
            for vcf_file in input.vcf_files:
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
    
                                                                                                                                                                                                                                                                                                        