rule create_mageck_design_matrix:
    input:
        metadata=CONFIG.metadata,
    output:
        design_matrix=OUTPUT_DIR + "/readcounts/mageck/design_matrix.txt",
    params:
        design_column="deseq2",
    log:
        OUTPUT_DIR + "/logs/readcounts/mageck/create_design_matrix.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/readcounts/mageck/create_design_matrix.tsv",
    message:
        "Creating MAGeCK design matrix from metadata",
    script:
        "../scripts/create_mageck_design_matrix.py"


rule mageck_count:
    input:
        fastqs=lambda wildcards: (
            expand(OUTPUT_DIR + "/trimmed/{sample}_1.fastq.gz", sample=SAMPLE_NAMES) +
            expand(OUTPUT_DIR + "/trimmed/{sample}_2.fastq.gz", sample=SAMPLE_NAMES)
            if any(INPUT_FILES.is_paired_end(s) for s in SAMPLE_NAMES)
            else expand(OUTPUT_DIR + "/trimmed/{sample}.fastq.gz", sample=SAMPLE_NAMES)
        ),
    output:
        counts=OUTPUT_DIR + "/readcounts/mageck/mageck_count.count.txt",
        normalized=OUTPUT_DIR + "/readcounts/mageck/mageck_count.count_normalized.txt",
        summary=OUTPUT_DIR + "/readcounts/mageck/mageck_count.countsummary.txt",
    params:
        options=str(CONFIG.third_party_tools.mageck.count.command_line_arguments),
        sgrna_library=CONFIG.third_party_tools.mageck.sgrna_library if CONFIG.third_party_tools.mageck.sgrna_library else "",
        output_prefix=OUTPUT_DIR + "/readcounts/mageck/mageck_count",
        sample_labels=",".join(SAMPLE_NAMES),
        fastq_1=lambda wildcards: " ".join([OUTPUT_DIR + f"/trimmed/{s}_1.fastq.gz" for s in SAMPLE_NAMES]),
        fastq_2=lambda wildcards: (
            " ".join([OUTPUT_DIR + f"/trimmed/{s}_2.fastq.gz" for s in SAMPLE_NAMES])
            if any(INPUT_FILES.is_paired_end(s) for s in SAMPLE_NAMES)
            else ""
        ),
        fastq_single=lambda wildcards: (
            " ".join([OUTPUT_DIR + f"/trimmed/{s}.fastq.gz" for s in SAMPLE_NAMES])
            if not any(INPUT_FILES.is_paired_end(s) for s in SAMPLE_NAMES)
            else ""
        ),
        is_paired=lambda wildcards: any(INPUT_FILES.is_paired_end(s) for s in SAMPLE_NAMES),
        control_sgrna=CONFIG.third_party_tools.mageck.control_sgrna if CONFIG.third_party_tools.mageck.control_sgrna else "",
        control_gene=CONFIG.third_party_tools.mageck.control_gene if CONFIG.third_party_tools.mageck.control_gene else "",
        norm_method=CONFIG.third_party_tools.mageck.count.norm_method,
        pdf_report="--pdf-report" if CONFIG.third_party_tools.mageck.count.pdf_report else "",
    threads: CONFIG.third_party_tools.mageck.count.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.MAGECK
    log: OUTPUT_DIR + "/logs/readcounts/mageck/mageck_count.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/mageck/mageck_count.tsv",
    message: "Counting guide RNAs using MAGeCK",
    shell: """
    if [ "{params.is_paired}" = "True" ]; then
        mageck count \
        --fastq {params.fastq_1} \
        --fastq-2 {params.fastq_2} \
        -l {params.sgrna_library} \
        -n {params.output_prefix} \
        --sample-label {params.sample_labels} \
        --norm-method {params.norm_method} \
        {params.pdf_report} \
        {params.control_sgrna} \
        {params.control_gene} \
        {params.options}
    else
        mageck count \
        --fastq {params.fastq_single} \
        -l {params.sgrna_library} \
        -n {params.output_prefix} \
        --sample-label {params.sample_labels} \
        --norm-method {params.norm_method} \
        {params.pdf_report} \
        {params.control_sgrna} \
        {params.control_gene} \
        {params.options}
    fi > {log} 2>&1
    """


rule mageck_test:
    input:
        counts=OUTPUT_DIR + "/readcounts/mageck/mageck_count.count.txt",
    output:
        gene_summary=OUTPUT_DIR + "/readcounts/mageck/mageck_test.gene_summary.txt",
        sgrna_summary=OUTPUT_DIR + "/readcounts/mageck/mageck_test.sgrna_summary.txt",
    params:
        options=str(CONFIG.third_party_tools.mageck.test.command_line_arguments) if CONFIG.third_party_tools.mageck.test else "",
        output_prefix=OUTPUT_DIR + "/readcounts/mageck/mageck_test",
        control_samples=CONFIG.third_party_tools.mageck.test.control_samples if CONFIG.third_party_tools.mageck.test and CONFIG.third_party_tools.mageck.test.control_samples else "",
        treatment_samples=CONFIG.third_party_tools.mageck.test.treatment_samples if CONFIG.third_party_tools.mageck.test and CONFIG.third_party_tools.mageck.test.treatment_samples else "",
        control_sgrna=CONFIG.third_party_tools.mageck.control_sgrna if CONFIG.third_party_tools.mageck.control_sgrna else "",
        control_gene=CONFIG.third_party_tools.mageck.control_gene if CONFIG.third_party_tools.mageck.control_gene else "",
    threads: CONFIG.third_party_tools.mageck.test.threads if CONFIG.third_party_tools.mageck.test else 1
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.MAGECK
    log: OUTPUT_DIR + "/logs/readcounts/mageck/mageck_test.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/mageck/mageck_test.tsv",
    message: "Running MAGeCK test (RRA) analysis",
    shell: """
    mageck test \
    -k {input.counts} \
    -t {params.treatment_samples} \
    -c {params.control_samples} \
    -n {params.output_prefix} \
    {params.control_sgrna} \
    {params.control_gene} \
    {params.options} \
    > {log} 2>&1
    """


rule mageck_mle:
    input:
        counts=OUTPUT_DIR + "/readcounts/mageck/mageck_count.count.txt",
        design_matrix=lambda wildcards: (
            CONFIG.third_party_tools.mageck.design_matrix
            if CONFIG.third_party_tools.mageck.design_matrix
            else OUTPUT_DIR + "/readcounts/mageck/design_matrix.txt"
        ),
    output:
        mle_gene_summary=OUTPUT_DIR + "/readcounts/mageck/mageck_mle.gene_summary.txt",
        mle_sgrna_summary=OUTPUT_DIR + "/readcounts/mageck/mageck_mle.sgrna_summary.txt",
    params:
        options=str(CONFIG.third_party_tools.mageck.mle.command_line_arguments),
        output_prefix=OUTPUT_DIR + "/readcounts/mageck/mageck_mle",
        control_sgrna=CONFIG.third_party_tools.mageck.control_sgrna if CONFIG.third_party_tools.mageck.control_sgrna else "",
        control_gene=CONFIG.third_party_tools.mageck.control_gene if CONFIG.third_party_tools.mageck.control_gene else "",
    threads: CONFIG.third_party_tools.mageck.mle.threads
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: CONTAINERS.MAGECK
    log: OUTPUT_DIR + "/logs/readcounts/mageck/mageck_mle.log",
    benchmark: OUTPUT_DIR + "/.benchmark/readcounts/mageck/mageck_mle.tsv",
    message: "Running MAGeCK MLE analysis",
    shell: """
    mageck mle \
    -k {input.counts} \
    -d {input.design_matrix} \
    -n {params.output_prefix} \
    {params.control_sgrna} \
    {params.control_gene} \
    {params.options} \
    > {log} 2>&1
    """