suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rjson))
suppressPackageStartupMessages(library(edgeR))

# Setup logging
log_file <- snakemake@log[[1]]
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

tryCatch({
  cat("=== edgeR Spike-in Normalisation ===\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

  bam_dir <- normalizePath(snakemake@params$bam_dir)
  cat("BAM directory:", bam_dir, "\n")
  cat("Counts file:", snakemake@input[[1]], "\n")
  cat("Metadata file:", snakemake@input[[2]], "\n\n")

  # Load the data
  cat("Loading metadata...\n")
  metadata <- read_csv(snakemake@input[[2]], show_col_types = FALSE) %>%
    mutate(deseq2 = factor(deseq2))
  cat("  Samples:", nrow(metadata), "\n")
  cat("  Sample IDs:", paste(metadata$sample_id, collapse = ", "), "\n\n")

  cat("Loading count matrix...\n")
  counts <- read.delim(snakemake@input[[1]], comment.char = "#", check.names = FALSE) %>%
    rename_with(~ tools::file_path_sans_ext(basename(.)), everything()) %>%
    column_to_rownames(var = "Geneid")
  cat("  Initial dimensions:", nrow(counts), "genes x", ncol(counts), "columns\n")

  # Select only numeric columns and match them to sample names
  numeric_cols <- counts %>%
    dplyr::select(where(is.numeric)) %>%
    colnames()
  cat("  Numeric columns found:", length(numeric_cols), "\n")

  # Match columns to sample IDs
  cat("\nMatching sample columns...\n")
  sample_cols <- sapply(metadata$sample_id, function(sid) {
    matching_col <- numeric_cols[grep(sid, numeric_cols, fixed = TRUE)]
    if (length(matching_col) == 0) stop(paste("No matching column found for sample:", sid))
    return(matching_col[1])
  })
  cat("  Matched columns:\n")
  for (i in seq_along(sample_cols)) {
    cat("    ", names(sample_cols)[i], "->", sample_cols[i], "\n")
  }

  counts <- counts %>%
    dplyr::select(all_of(sample_cols)) %>%
    setNames(metadata$sample_id) %>%
    filter(rowSums(.) > 0)
  cat("\n  Final count matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

  # Create DGEList
  cat("\nCreating DGEList...\n")
  dge <- DGEList(counts = counts)

  # Handle spike-in genes
  spikein_gene_names <- snakemake@params$spikein_genes
  cat("\nSpike-in gene configuration:\n")
  if (length(spikein_gene_names) == 0 || all(spikein_gene_names == "")) {
    cat("  No spike-in genes provided - using TMM normalisation on all genes\n")
    dge <- calcNormFactors(dge, method = "TMM")
  } else {
    cat("  Spike-in genes requested:", paste(spikein_gene_names, collapse = ", "), "\n")
    spikein <- rownames(dge) %in% spikein_gene_names
    cat("  Spike-in genes found in data:", sum(spikein), "\n")
    if (sum(spikein) == 0) {
      cat("  WARNING: No spike-in genes found in count matrix - using TMM normalisation on all genes\n")
      dge <- calcNormFactors(dge, method = "TMM")
    } else {
      cat("  Using spike-in genes for TMM normalisation\n")
      dge <- calcNormFactors(dge, method = "TMM", subset = spikein)
      # Remove spike-ins
      dge <- dge[!spikein, , keep.lib.sizes = FALSE]
    }
  }

  # edgeR model fitting
  cat("\nFitting edgeR model...\n")
  design <- model.matrix(~ deseq2, data = metadata)
  dge <- estimateDisp(dge, design)
  fit <- glmQLFit(dge, design)

  # Extract normalization factors
  norm_factors <- dge$samples$norm.factors
  names(norm_factors) <- metadata$sample_id

  cat("\nNormalisation factors:\n")
  for (i in seq_along(norm_factors)) {
    cat("  ", names(norm_factors)[i], ":", round(norm_factors[i], 4), "\n")
  }

  # Write JSON
  sf_json <- toJSON(norm_factors)
  writeLines(sf_json, snakemake@output$normalisation_factors)
  cat("\nWrote JSON output:", snakemake@output$normalisation_factors, "\n")

  # Write TSV
  nf_df <- data.frame(
    sample = names(norm_factors),
    scale_factor = as.numeric(norm_factors)
  )
  write_tsv(nf_df, snakemake@output$normalisation_table)
  cat("Wrote TSV output:", snakemake@output$normalisation_table, "\n")

  cat("\n=== Completed successfully ===\n")

}, error = function(e) {
  cat("\n=== ERROR ===\n")
  cat("Error message:", conditionMessage(e), "\n")
  cat("\nTraceback:\n")
  traceback()
  stop(e)
})

sink(type = "message")
sink(type = "output")
close(log_con)
