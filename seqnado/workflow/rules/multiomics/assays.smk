
import yaml


# Load config for each assay upfront
LOADED_CONFIGS = {}
for assay in ASSAYS:
    config_path = ASSAY_CONFIGS[assay]
    with open(config_path, 'r') as f:
        assay_config = yaml.safe_load(f)
    assay_config["output_dir"] = f"{OUTPUT_DIR}/{assay.clean_name}"
    LOADED_CONFIGS[assay.clean_name] = assay_config


# Create individual module instances for each assay
if "atac" in LOADED_CONFIGS:
    module run_atac:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["atac"]

    use rule * from run_atac as atac_*

if "cat" in LOADED_CONFIGS:
    module run_cat:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["cat"]

    use rule * from run_cat as cat_*

if "chip" in LOADED_CONFIGS:
    module run_chip:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["chip"]

    use rule * from run_chip as chip_*

if "meth" in LOADED_CONFIGS:
    module run_meth:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["meth"]

    use rule * from run_meth as meth_*

if "rna" in LOADED_CONFIGS:
    module run_rna:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["rna"]

    use rule * from run_rna as rna_*

if "snp" in LOADED_CONFIGS:
    module run_snp:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["snp"]

    use rule * from run_snp as snp_*

if "mcc" in LOADED_CONFIGS:
    module run_mcc:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["mcc"]

    use rule * from run_mcc as mcc_*

if "crispr" in LOADED_CONFIGS:
    module run_crispr:
        snakefile: workflow.basedir + "/Snakefile"
        config: LOADED_CONFIGS["crispr"]

    use rule * from run_crispr as crispr_*