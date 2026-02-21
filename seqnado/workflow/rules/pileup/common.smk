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
        ) if CONFIG.third_party_tools.homer else "",
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

rule homer_make_tag_directory_merged:
    input:
        bam=OUTPUT_DIR + "/aligned/merged/{group}.bam",
    output:
        homer_tag_directory=directory(OUTPUT_DIR + "/tag_dirs/merged/unscaled/{group}"),
    params:
        options=str(
            CONFIG.third_party_tools.homer.make_tag_directory.command_line_arguments
        ) if CONFIG.third_party_tools.homer else "",
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
        OUTPUT_DIR + "/logs/homer/maketagdirectory/merged/{group}.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/homer/maketagdirectory/merged/{group}.tsv"
    message:
        "Making tag directory with HOMER for merged sample {wildcards.group}"
    wildcard_constraints:
        group="|".join(SAMPLE_GROUPINGS.get_grouping('consensus').group_names),
    shell:
        """
        makeTagDirectory {output.homer_tag_directory} {input.bam} {params.options} > {log} 2>&1
        """



ruleorder: homer_make_tag_directory > homer_make_tag_directory_merged
