################################
# Genome Index Building Rules
# Builds Bowtie2 and STAR indices for standard and composite genomes
################################


################################
# Standard Genome Index Rules
################################


rule samtools_faidx:
    input:
        fasta=OUTPUT_DIR + "/{genome}/sequence/{genome}.fa",
    output:
        fai=OUTPUT_DIR + "/{genome}/sequence/{genome}.fa.fai",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/index/{genome}_samtools_faidx.log",
    message:
        "Indexing FASTA with samtools faidx for {wildcards.genome}"
    shell:
        """
        samtools faidx {input.fasta} 2>&1 | tee {log}
        """


rule build_bowtie2:
    input:
        fasta=OUTPUT_DIR + "/{genome}/sequence/{genome}.fa",
    output:
        index=OUTPUT_DIR + "/{genome}/bt2_index/{genome}.1.bt2",
    params:
        prefix=OUTPUT_DIR + "/{genome}/bt2_index/{genome}",
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    threads: 16
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/index/{genome}_bowtie2_build.log",
    message:
        "Building Bowtie2 index for {wildcards.genome}"
    shell:
        """
        bowtie2-build --threads {threads} \
            {input.fasta} \
            {params.prefix} 2>&1 | tee {log}
        """


rule build_star:
    input:
        fasta=OUTPUT_DIR + "/{genome}/sequence/{genome}.fa",
        gtf=OUTPUT_DIR + "/{genome}/genes/{genome}.ncbiRefSeq.gtf",
    output:
        genome_dir=directory(OUTPUT_DIR + "/{genome}/STAR_2.7.10b"),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    threads: 4
    resources:
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    log:
        OUTPUT_DIR + "/logs/index/{genome}_star_build.log",
    message:
        "Building STAR index for {wildcards.genome}"
    shell:
        """
        mkdir -p {output.genome_dir}

        genome_length=$(awk '/^>/ {{next}} {{sum+=length($0)}} END {{print sum}}' {input.fasta})
        recommended_value=$(awk -v gl="$genome_length" 'BEGIN {{ v=int(log(gl)/log(2)/2)-1; print (v > 14 ? 14 : v) }}')

        echo "Genome length: $genome_length" | tee {log}
        echo "genomeSAindexNbases: $recommended_value" | tee -a {log}

        STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output.genome_dir} \
            --genomeFastaFiles {input.fasta} \
            --genomeSAindexNbases $recommended_value \
            --sjdbGTFfile {input.gtf} \
            --outTmpDir {output.genome_dir}/tmp 2>&1 | tee -a {log}

        rm -rf {output.genome_dir}/tmp
        """


################################
# Composite Genome Rules (only when SPIKEIN is specified)
################################

if SPIKEIN:
    rule concat_fasta:
        input:
            primary=lambda wc: OUTPUT_DIR + f"/{GENOME}/sequence/{GENOME}.fa",
            spikein=lambda wc: OUTPUT_DIR + f"/{SPIKEIN}/sequence/{SPIKEIN}.fa",
        output:
            fasta=OUTPUT_DIR + "/{primary}_{spikein}/sequence/{primary}_{spikein}.fa",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        log:
            OUTPUT_DIR + "/logs/composite/{primary}_{spikein}_concat_fasta.log",
        message:
            "Concatenating FASTA: {wildcards.primary} + {wildcards.spikein}"
        shell:
            """
            sed 's/^>chr/>{wildcards.spikein}_chr/' {input.spikein} > {output.fasta}.spikein_tmp 2>&1 | tee {log}
            cat {input.primary} {output.fasta}.spikein_tmp > {output.fasta} 2>&1 | tee -a {log}
            rm -f {output.fasta}.spikein_tmp
            """


    rule concat_gtf:
        input:
            primary=lambda wc: OUTPUT_DIR + f"/{GENOME}/genes/{GENOME}.ncbiRefSeq.gtf",
            spikein=lambda wc: OUTPUT_DIR + f"/{SPIKEIN}/genes/{SPIKEIN}.ncbiRefSeq.gtf",
        output:
            gtf=OUTPUT_DIR + "/{primary}_{spikein}/genes/{primary}_{spikein}.gtf",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        log:
            OUTPUT_DIR + "/logs/composite/{primary}_{spikein}_concat_gtf.log",
        message:
            "Concatenating GTF: {wildcards.primary} + {wildcards.spikein}"
        shell:
            """
            sed 's/^chr/{wildcards.spikein}_chr/' {input.spikein} > {output.gtf}.spikein_tmp 2>&1 | tee {log}
            cat {input.primary} {output.gtf}.spikein_tmp > {output.gtf} 2>&1 | tee -a {log}
            rm -f {output.gtf}.spikein_tmp
            """


    rule composite_samtools_faidx:
        input:
            fasta=OUTPUT_DIR + "/{primary}_{spikein}/sequence/{primary}_{spikein}.fa",
        output:
            fai=OUTPUT_DIR + "/{primary}_{spikein}/sequence/{primary}_{spikein}.fa.fai",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=1, attempts=attempt, scale=SCALE_RESOURCES),
        log:
            OUTPUT_DIR + "/logs/composite/{primary}_{spikein}_samtools_faidx.log",
        message:
            "Indexing composite FASTA with samtools faidx for {wildcards.primary}_{wildcards.spikein}"
        shell:
            """
            samtools faidx {input.fasta} 2>&1 | tee {log}
            """


    rule build_composite_bowtie2:
        input:
            fasta=OUTPUT_DIR + "/{primary}_{spikein}/sequence/{primary}_{spikein}.fa",
        output:
            index=OUTPUT_DIR + "/{primary}_{spikein}/bt2_index/{primary}_{spikein}.1.bt2",
        params:
            prefix=OUTPUT_DIR + "/{primary}_{spikein}/bt2_index/{primary}_{spikein}",
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        threads: 16
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=8, attempts=attempt, scale=SCALE_RESOURCES),
        log:
            OUTPUT_DIR + "/logs/composite/{primary}_{spikein}_bowtie2_build.log",
        message:
            "Building Bowtie2 index for composite {wildcards.primary}_{wildcards.spikein}"
        shell:
            """
            bowtie2-build --threads {threads} \
                {input.fasta} \
                {params.prefix} 2>&1 | tee {log}
            """


    rule build_composite_star:
        input:
            fasta=OUTPUT_DIR + "/{primary}_{spikein}/sequence/{primary}_{spikein}.fa",
            gtf=OUTPUT_DIR + "/{primary}_{spikein}/genes/{primary}_{spikein}.gtf",
        output:
            genome_dir=directory(OUTPUT_DIR + "/{primary}_{spikein}/STAR_2.7.10b"),
        container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
        threads: 4
        resources:
            mem=lambda wildcards, attempt: define_memory_requested(initial_value=32, attempts=attempt, scale=SCALE_RESOURCES),
            runtime=lambda wildcards, attempt: define_time_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
        log:
            OUTPUT_DIR + "/logs/composite/{primary}_{spikein}_star_build.log",
        message:
            "Building STAR index for composite {wildcards.primary}_{wildcards.spikein}"
        shell:
            """
            mkdir -p {output.genome_dir}

            genome_length=$(awk '/^>/ {{next}} {{sum+=length($0)}} END {{print sum}}' {input.fasta})
            recommended_value=$(awk -v gl="$genome_length" 'BEGIN {{ v=int(log(gl)/log(2)/2)-1; print (v > 14 ? 14 : v) }}')

            echo "Genome length: $genome_length" | tee {log}
            echo "genomeSAindexNbases: $recommended_value" | tee -a {log}

            STAR --runThreadN {threads} \
                --runMode genomeGenerate \
                --genomeDir {output.genome_dir} \
                --genomeFastaFiles {input.fasta} \
                --genomeSAindexNbases $recommended_value \
                --sjdbGTFfile {input.gtf} \
                --outTmpDir {output.genome_dir}/tmp 2>&1 | tee -a {log}

            rm -rf {output.genome_dir}/tmp
            """
