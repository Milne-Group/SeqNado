
from seqnado.workflow.helpers.multiomics import get_assay_all_inputs


rule multiomics_summary:
    """Generate a summary report of all assays."""
    input:
        get_assay_all_inputs(ASSAYS=ASSAYS, rules=rules)
    output:
        OUTPUT_DIR + "/multiomics_summary.txt"
    message:
        "Generating multiomics summary report"
    log:
        OUTPUT_DIR + "/logs/multiomics/multiomics_summary.log",
    benchmark:
        OUTPUT_DIR + "/.benchmark/multiomics/multiomics_summary.tsv",
    run:
        with open(output[0], 'w') as f:
            f.write("SeqNado Multiomics Project Summary\n")
            f.write("=" * 70 + "\n\n")
            f.write("Run Directory: " + str(Path.cwd()) + "\n")
            f.write(f"Total assays: {len(ASSAYS)}\n")
            f.write(f"Assays: {', '.join([assay.clean_name for assay in ASSAYS])}\n\n")

            for assay in ASSAYS:
                assay = assay.clean_name
                f.write(f"\n{assay.upper()}\n")
                f.write("-" * 50 + "\n")
                f.write(f"  Config:   config_{assay}.yaml\n")
                f.write(f"  Metadata: metadata_{assay}.csv\n")
                f.write(f"  Output:   seqnado_output/{assay}/\n")
                f.write(f"  Report:   seqnado_output/{assay}/seqnado_report.html\n")

                # Check for completeness
                complete_file = OUTPUT_DIR + f"/{assay}/seqnado_report.html"
                if os.path.exists(complete_file):
                    f.write(f"  STATUS:   COMPLETE\n")
                else:
                    f.write(f"  WARNING: SeqNado report not found! Check logs for errors.\n")
        with open(log[0], 'w') as log_handle:
            log_handle.write("Multiomics summary generated.\n")
        