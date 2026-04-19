# eval_step4b.R — Validate the human-validated BBB gene list
library(readr)
library(dplyr)

eval_human_validation <- function(filepath = "processed/master_BBB_genelist_validated.csv") {

  cat("=== EVAL: Human Validation Cross-Reference ===\n\n")
  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Required columns
  required <- c("GeneSymbol_Human", "In_Walchli_EC", "In_Yang_EC",
                 "In_Winkler_EC", "Human_BBB_Datasets_N")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    cat("FAIL: Missing columns:", paste(missing, collapse = ", "), "\n")
  } else {
    cat("PASS: All required columns present\n")
  }

  # 2. Row count preserved
  if (nrow(df) == 1526) {
    cat("PASS: Row count preserved (1526 genes)\n")
  } else {
    cat(sprintf("FAIL: Expected 1526 rows, got %d\n", nrow(df)))
  }

  # 3. Coverage summary
  cat(sprintf("\nINFO: Wälchli EC:  %d / %d (%.1f%%)\n",
              sum(df$In_Walchli_EC), nrow(df), 100 * mean(df$In_Walchli_EC)))
  cat(sprintf("INFO: Yang EC:     %d / %d (%.1f%%)\n",
              sum(df$In_Yang_EC), nrow(df), 100 * mean(df$In_Yang_EC)))
  cat(sprintf("INFO: Winkler EC:  %d / %d (%.1f%%)\n",
              sum(df$In_Winkler_EC), nrow(df), 100 * mean(df$In_Winkler_EC)))

  # 4. Composite N distribution
  cat("\nINFO: Dataset confirmation distribution:\n")
  print(table(df$Human_BBB_Datasets_N))

  # 5. Wälchli coverage should be high (full expression matrix)
  if (mean(df$In_Walchli_EC) < 0.80) {
    cat("WARN: Wälchli coverage unexpectedly low — check expression matrix\n")
  } else {
    cat("\nPASS: Wälchli coverage looks reasonable for a full expression matrix\n")
  }

  # 6. Landmark gene check
  # Note: CLDN5 is intentionally absent — it is pan-endothelial and was not
  # enriched enough in either Daneman or Munji to make the source gene lists.
  cat("\nLandmark gene check:\n")
  landmarks <- c("SLC2A1", "ABCB1", "ABCG2", "SLC7A5")
  for (g in landmarks) {
    row <- df[df$GeneSymbol_Human == g, ]
    if (nrow(row) == 0) {
      cat(sprintf("WARN: %s not found\n", g))
    } else {
      cat(sprintf("  %s: Wälchli=%s | Yang=%s | Winkler=%s | N=%d\n",
                  g,
                  ifelse(row$In_Walchli_EC[1], "YES", "no"),
                  ifelse(row$In_Yang_EC[1],    "YES", "no"),
                  ifelse(row$In_Winkler_EC[1], "YES", "no"),
                  row$Human_BBB_Datasets_N[1]))
    }
  }

  cat("\n=== EVAL COMPLETE ===\n")
}

eval_human_validation()
