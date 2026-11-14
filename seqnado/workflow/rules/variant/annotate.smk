


rule bcftools_annotate:
    input:
        vcf=rules.split_multiallelic.output.vcf,
        idx=rules.split_multiallelic.output.idx,
    output:
        vcf=OUTPUT_DIR + "/variant/{sample}.anno.vcf.gz",
        idx=OUTPUT_DIR + "/variant/{sample}.anno.vcf.gz.tbi",
    params:
        db=CONFIG.assay_config.snp_calling.snp_database,
    resources:
        mem=lambda wildcards, attempt: f"{10 * 2 ** (attempt -1)}GB",
        runtime=lambda wildcards, attempt: f"{5 * 2 ** (attempt - 1)}h",
    threads: 16
    container: CONTAINERS.SEQNADO_PIPELINE
    log: OUTPUT_DIR + "/logs/variant/{sample}_anno.log",
    benchmark: OUTPUT_DIR + "/.benchmark/variant/{sample}_anno.tsv",
    message: "Annotating VCF for sample {wildcards.sample}"
    shell: """
    bcftools annotate --threads {threads} -c ID -a {params.db} -Oz -o {output.vcf} {input.vcf} 2> {log} &&
    tabix -f {output.vcf} > {output.idx}
    """


ruleorder: bcftools_annotate > split_multiallelic > bcftools_call_snp