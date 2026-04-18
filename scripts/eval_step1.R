# eval_step1.R — Validate the master BBB gene list
library(readr)
library(dplyr)

eval_master_BBB_genelist <- function(filepath = "processed/master_BBB_genelist.csv") {

  cat("=== EVAL: Master BBB Gene List ===\n\n")
  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Check required columns exist
  required_cols <- c("GeneSymbol_Mouse", "GeneSymbol_Human", "Source")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    cat("FAIL: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  } else {
    cat("PASS: All required columns present\n")
  }

  # 2. Check gene count is in expected range
  n <- nrow(df)
  cat(sprintf("INFO: Total genes in master list: %d\n", n))
  if (n < 300) {
    cat("WARN: Fewer than 300 genes — something may have gone wrong with loading\n")
  } else if (n > 2500) {
    cat("WARN: More than 2500 genes — possible duplication issue\n")
  } else {
    cat("PASS: Gene count looks reasonable\n")
  }

  # 3. Check for duplicates
  n_dupes <- sum(duplicated(df$GeneSymbol_Human))
  if (n_dupes > 0) {
    cat(sprintf("FAIL: %d duplicate human gene symbols found — deduplication needed\n", n_dupes))
  } else {
    cat("PASS: No duplicate human gene symbols\n")
  }

  # 4. Check human gene symbols are uppercase (HGNC format)
  n_lowercase <- sum(df$GeneSymbol_Human != toupper(df$GeneSymbol_Human), na.rm = TRUE)
  if (n_lowercase > 0) {
    cat(sprintf("WARN: %d gene symbols are not uppercase — may not be in HGNC format\n", n_lowercase))
  } else {
    cat("PASS: All human gene symbols are uppercase (HGNC format)\n")
  }

  # 5. Check Source column has expected values
  valid_sources <- c("Daneman", "Munji", "Both")
  invalid_sources <- df %>% filter(!Source %in% valid_sources) %>% nrow()
  if (invalid_sources > 0) {
    cat(sprintf("WARN: %d rows have unexpected Source values\n", invalid_sources))
  } else {
    cat("PASS: Source column values are valid\n")
  }

  # 6. Check that known landmark BBB genes are present
  landmark_genes <- c("SLC2A1", "CLDN5", "ABCB1", "TJP1", "PECAM1", "VWF")
  missing_landmarks <- setdiff(landmark_genes, df$GeneSymbol_Human)
  if (length(missing_landmarks) > 0) {
    cat(sprintf("WARN: These known BBB landmark genes are missing: %s\n",
                paste(missing_landmarks, collapse = ", ")))
    cat("      (This may be okay — not all landmarks appear in every study)\n")
  } else {
    cat("PASS: All landmark BBB genes present (SLC2A1, CLDN5, ABCB1, TJP1, PECAM1, VWF)\n")
  }

  # 7. Summary table by source
  cat("\n--- Gene counts by source ---\n")
  print(df %>% count(Source))

  cat("\n=== EVAL COMPLETE ===\n")
}

# Run it
eval_master_BBB_genelist()
