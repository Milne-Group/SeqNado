rule bcftools_call_snp:
    input:
        bam=OUTPUT_DIR + "/aligned/{sample}.bam",
        bai=OUTPUT_DIR + "/aligned/{sample}.bam.bai",
    output:
        vcf=temp(OUTPUT_DIR + "/variant/{sample}.raw.vcf.gz"),
        idx=temp(OUTPUT_DIR + "/variant/{sample}.raw.vcf.gz.tbi"),
        stats=OUTPUT_DIR + "/qc/variant/{sample}.stats.txt",
    params:
        fasta=CONFIG.genome.fasta,
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: CONFIG.third_party_tools.bcftools.call.threads
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/variant/{sample}.log",
    benchmark: OUTPUT_DIR + "/.benchmark/variant/{sample}.tsv",
    message: "Calling variants for sample {wildcards.sample} using bcftools"
    shell: """
    bcftools mpileup --threads {threads} -Ou -f {params.fasta} {input.bam} | bcftools call --threads {threads} -mv -Oz -o {output.vcf} > {log} 2>&1 &&
    tabix -f {output.vcf} > {output.idx} &&
    bcftools stats -F {params.fasta} -s - {output.vcf} > {output.stats}
    """


rule split_multiallelic:
    input:
        vcf=rules.bcftools_call_snp.output.vcf,
        idx=rules.bcftools_call_snp.output.idx,
    output:
        vcf=OUTPUT_DIR + "/variant/{sample}.vcf.gz",
        idx=OUTPUT_DIR + "/variant/{sample}.vcf.gz.tbi",
    params:
        outdir=OUTPUT_DIR + "/variant/",
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: 16
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/variant/{sample}_split.log",
    benchmark: OUTPUT_DIR + "/.benchmark/variant/{sample}_split.tsv",
    message: "Splitting multiallelic variants for sample {wildcards.sample}"
    shell: """
    bcftools norm --threads {threads} -m-any -Oz -o {output.vcf} {input.vcf} > {log} 2>&1 &&
    tabix -f {output.vcf} > {output.idx}
    """
