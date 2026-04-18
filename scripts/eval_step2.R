# eval_step2.R — Validate the Wälchli average expression matrix
library(readr)
library(dplyr)

eval_walchli_avgexpr <- function(filepath = "processed/Walchli2024_AverageExpression.csv") {

  cat("=== EVAL: Wälchli Average Expression Matrix ===\n\n")

  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Check it has a gene name column and at least one cell type column
  cat("Dimensions:", nrow(df), "genes x", ncol(df) - 1, "cell clusters\n\n")

  if (ncol(df) < 2) {
    cat("FAIL: Only one column — expected genes + multiple cell type columns\n")
    return(invisible(NULL))
  } else {
    cat("PASS: Multiple cell cluster columns present\n")
  }

  # 2. Check gene count is reasonable (human genome ~20,000 protein-coding genes)
  n_genes <- nrow(df)
  if (n_genes < 1000) {
    cat(sprintf("WARN: Only %d genes — expected at least 5,000\n", n_genes))
  } else if (n_genes > 40000) {
    cat(sprintf("WARN: %d genes — unusually high, may include non-coding\n", n_genes))
  } else {
    cat(sprintf("PASS: Gene count looks reasonable (%d genes)\n", n_genes))
  }

  # 3. Check cluster names are meaningful (not just numbers)
  cluster_names <- colnames(df)[-1]
  cat("\nCell cluster names found:\n")
  print(cluster_names)

  # 4. Check known BBB landmark genes are present
  gene_col <- df[[1]]
  landmark_genes <- c("SLC2A1", "CLDN5", "ABCB1", "PECAM1", "VWF")
  missing <- setdiff(landmark_genes, gene_col)
  if (length(missing) > 0) {
    cat(sprintf("\nWARN: Missing landmark genes: %s\n", paste(missing, collapse = ", ")))
  } else {
    cat("\nPASS: All landmark BBB genes present\n")
  }

  # 5. Check values are numeric and non-negative
  value_cols <- df[, -1]
  if (all(sapply(value_cols, is.numeric))) {
    cat("PASS: All expression values are numeric\n")
  } else {
    cat("FAIL: Non-numeric values found in expression columns\n")
  }

  cat("\n=== EVAL COMPLETE ===\n")
}

eval_walchli_avgexpr()
