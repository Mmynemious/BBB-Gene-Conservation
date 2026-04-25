# eval_step4c.R — Validate the liver control gene set
library(readr)
library(dplyr)

eval_liver_control <- function(filepath = "processed/liver_control_genelist.csv") {

  cat("=== EVAL: Liver Control Gene Set ===\n\n")
  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Required columns
  required <- c("Gene", "Liver_nTPM", "Brain_nTPM")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    cat("FAIL: Missing columns:", paste(missing, collapse = ", "), "\n")
  } else {
    cat("PASS: All required columns present\n")
  }

  # 2. Size check
  cat(sprintf("INFO: Control set size: %d genes\n", nrow(df)))
  if (nrow(df) < 3000) {
    cat("WARN: Control set smaller than 3x BBB list — consider lowering nTPM threshold\n")
  } else {
    cat("PASS: Control set is at least 3x the BBB list size\n")
  }

  # 3. No BBB gene overlap
  master <- read_csv("processed/master_BBB_genelist_complete.csv", show_col_types = FALSE)
  bbb_genes <- unique(na.omit(master$GeneSymbol_Human))
  overlap <- sum(df$Gene %in% bbb_genes)
  if (overlap > 0) {
    cat(sprintf("FAIL: %d genes overlap with BBB master list\n", overlap))
  } else {
    cat("PASS: No overlap with BBB master list\n")
  }

  # 4. Expression range
  cat(sprintf("\nLiver nTPM: median=%.1f  mean=%.1f  min=%.1f  max=%.1f\n",
      median(df$Liver_nTPM), mean(df$Liver_nTPM),
      min(df$Liver_nTPM), max(df$Liver_nTPM)))
  cat(sprintf("Brain nTPM: median=%.1f  mean=%.1f\n",
      median(df$Brain_nTPM), mean(df$Brain_nTPM)))

  # 5. Landmark liver genes present
  cat("\nLandmark liver gene check:\n")
  landmarks <- c("ALB", "APOB", "CYP3A4", "TTR")  # albumin, apoB, cytochrome P450, transthyretin
  for (g in landmarks) {
    row <- df[df$Gene == g, ]
    if (nrow(row) == 0) {
      cat(sprintf("  WARN: %s not found\n", g))
    } else {
      cat(sprintf("  PASS: %s — Liver nTPM=%.1f\n", g, row$Liver_nTPM[1]))
    }
  }

  cat("\n=== EVAL COMPLETE ===\n")
}

eval_liver_control()
