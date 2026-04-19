# eval_step3b.R — Validate the orthology confidence table
library(readr)
library(dplyr)

eval_orthology_confidence <- function(filepath = "processed/BBB_genes_orthology_confidence.csv") {

  cat("=== EVAL: Orthology Confidence Table ===\n\n")
  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Required columns
  required <- c("GeneSymbol_Human", "HumanEnsembl",
                 "MouseEnsembl", "Mouse_OrthoType", "Mouse_PercId",
                 "MacaqueEnsembl", "Macaque_OrthoType", "Macaque_PercId",
                 "Human_BBB_Datasets_N")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    cat("FAIL: Missing columns:", paste(missing, collapse = ", "), "\n")
  } else {
    cat("PASS: All required columns present\n")
  }

  # 2. Coverage
  n_human   <- sum(!is.na(df$HumanEnsembl))
  n_mouse   <- sum(!is.na(df$MouseEnsembl))
  n_macaque <- sum(!is.na(df$MacaqueEnsembl))
  cat(sprintf("\nINFO: Total rows:          %d\n", nrow(df)))
  cat(sprintf("INFO: Human Ensembl IDs:   %d (%.1f%%)\n", n_human,   100*n_human/nrow(df)))
  cat(sprintf("INFO: Mouse orthologues:   %d (%.1f%%)\n", n_mouse,   100*n_mouse/nrow(df)))
  cat(sprintf("INFO: Macaque orthologues: %d (%.1f%%)\n", n_macaque, 100*n_macaque/nrow(df)))

  # 3. Orthologue type distribution
  cat("\nMouse orthologue types:\n")
  print(table(df$Mouse_OrthoType, useNA="ifany"))
  cat("\nMacaque orthologue types:\n")
  print(table(df$Macaque_OrthoType, useNA="ifany"))

  # 4. % identity sanity checks
  cat(sprintf("\nMouse %%id:   median=%.1f  mean=%.1f\n",
      median(df$Mouse_PercId, na.rm=TRUE), mean(df$Mouse_PercId, na.rm=TRUE)))
  cat(sprintf("Macaque %%id: median=%.1f  mean=%.1f\n",
      median(df$Macaque_PercId, na.rm=TRUE), mean(df$Macaque_PercId, na.rm=TRUE)))

  if (median(df$Macaque_PercId, na.rm=TRUE) < median(df$Mouse_PercId, na.rm=TRUE)) {
    cat("WARN: Macaque % identity lower than mouse — check data\n")
  } else {
    cat("PASS: Macaque % identity higher than mouse (expected — macaque is closer to human)\n")
  }

  # 5. Landmark genes
  cat("\nLandmark gene check:\n")
  for (g in c("SLC2A1", "ABCB1", "ABCG2")) {
    rows <- df[!is.na(df$GeneSymbol_Human) & df$GeneSymbol_Human == g, ]
    if (nrow(rows) == 0) {
      cat(sprintf("  WARN: %s not found\n", g))
    } else {
      for (i in seq_len(nrow(rows))) {
        cat(sprintf("  %s: Mo=%s [%s, %.1f%%]  Ma=%s [%s, %.1f%%]\n",
            g,
            ifelse(is.na(rows$MouseEnsembl[i]),   "NA", substr(rows$MouseEnsembl[i], 1, 18)),
            ifelse(is.na(rows$Mouse_OrthoType[i]), "NA", rows$Mouse_OrthoType[i]),
            ifelse(is.na(rows$Mouse_PercId[i]),    NA,   rows$Mouse_PercId[i]),
            ifelse(is.na(rows$MacaqueEnsembl[i]),  "NA", substr(rows$MacaqueEnsembl[i], 1, 18)),
            ifelse(is.na(rows$Macaque_OrthoType[i]),"NA", rows$Macaque_OrthoType[i]),
            ifelse(is.na(rows$Macaque_PercId[i]),  NA,   rows$Macaque_PercId[i])))
      }
    }
  }

  cat("\n=== EVAL COMPLETE ===\n")
}

eval_orthology_confidence()
