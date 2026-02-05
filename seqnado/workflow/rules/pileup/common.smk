rule homer_make_tag_directory:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
    output:
        homer_tag_directory=directory(OUTPUT_DIR + "/tag_dirs/{sample}"),
    wildcard_constraints:
        sample=r"(?!merged/).*",
    params:
        options=str(
            CONFIG.third_party_tools.homer.make_tag_directory.command_line_arguments
        ),
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(
            initial_value=4, attempts=attempt, scale=SCALE_RESOURCES
        ),
        runtime=lambda wildcards, attempt: define_time_requested(
            initial_value=2, attempts=attempt, scale=SCALE_RESOURCES
        ),
    container:
        "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:
        OUTPUT_DIR + "/logs/homer/maketagdirectory_{sample}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/maketagdirectory_{sample}.tsv"
    message:
        "Making tag directory with HOMER for sample {wildcards.sample}"
    shell:
        """
    makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1
    """