# eval_step3.R — Validate the cross-species Ensembl ID table
library(readr)
library(dplyr)

eval_crossspecies_ensembl <- function(filepath = "processed/BBB_genes_crossspecies_ensembl.csv") {

  cat("=== EVAL: Cross-Species Ensembl ID Table ===\n\n")
  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Required columns
  required <- c("GeneSymbol_Human", "HumanEnsembl", "MouseEnsembl", "MacaqueEnsembl")
  missing  <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    cat("FAIL: Missing columns:", paste(missing, collapse = ", "), "\n")
  } else {
    cat("PASS: All required columns present\n")
  }

  # 2. Coverage per species
  n_human   <- sum(!is.na(df$HumanEnsembl) & df$HumanEnsembl != "")
  n_mouse   <- sum(!is.na(df$MouseEnsembl) & df$MouseEnsembl != "")
  n_macaque <- sum(!is.na(df$MacaqueEnsembl) & df$MacaqueEnsembl != "")
  n_total   <- nrow(df)

  cat(sprintf("INFO: Total rows: %d\n", n_total))
  cat(sprintf("INFO: Human Ensembl IDs:   %d (%.1f%%)\n", n_human,   100 * n_human / n_total))
  cat(sprintf("INFO: Mouse orthologues:   %d (%.1f%%)\n", n_mouse,   100 * n_mouse / n_total))
  cat(sprintf("INFO: Macaque orthologues: %d (%.1f%%)\n", n_macaque, 100 * n_macaque / n_total))

  if (n_macaque / n_total < 0.5)
    cat("WARN: Fewer than 50% of genes have macaque orthologues — check Ensembl coverage\n")
  else
    cat("PASS: Macaque orthologue coverage looks reasonable\n")

  # 3. Check Ensembl ID format (should start with ENSG / ENSMUSG / ENSMMUG)
  bad_human   <- sum(!grepl("^ENSG", df$HumanEnsembl) & !is.na(df$HumanEnsembl))
  bad_mouse   <- sum(!grepl("^ENSMUSG", df$MouseEnsembl) & !is.na(df$MouseEnsembl))
  bad_macaque <- sum(!grepl("^ENSMMUG", df$MacaqueEnsembl) & !is.na(df$MacaqueEnsembl))

  if (bad_human + bad_mouse + bad_macaque > 0) {
    cat(sprintf("WARN: Unexpected Ensembl ID formats — human: %d, mouse: %d, macaque: %d\n",
                bad_human, bad_mouse, bad_macaque))
  } else {
    cat("PASS: Ensembl ID formats look correct\n")
  }

  # 4. Check landmark genes have all three IDs
  landmarks <- c("SLC2A1", "ABCB1", "ABCG2")
  for (g in landmarks) {
    row <- df[df$GeneSymbol_Human == g, ]
    if (nrow(row) == 0) {
      cat(sprintf("WARN: %s not found in table\n", g))
    } else {
      has_all <- !is.na(row$HumanEnsembl[1]) &
                 !is.na(row$MouseEnsembl[1]) &
                 !is.na(row$MacaqueEnsembl[1])
      cat(sprintf("%s %s: H=%s | Mo=%s | Ma=%s\n",
                  if (has_all) "PASS:" else "WARN:",
                  g,
                  ifelse(is.na(row$HumanEnsembl[1]),   "missing", row$HumanEnsembl[1]),
                  ifelse(is.na(row$MouseEnsembl[1]),    "missing", row$MouseEnsembl[1]),
                  ifelse(is.na(row$MacaqueEnsembl[1]),  "missing", row$MacaqueEnsembl[1])))
    }
  }

  cat("\n=== EVAL COMPLETE ===\n")
}

eval_crossspecies_ensembl()
