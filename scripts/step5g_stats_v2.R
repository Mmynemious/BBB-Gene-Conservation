# ============================================================
# Step 5g: Statistical comparison using corrected (v2) dN/dS values
#
# Same analysis structure as step5e but on the codon-aligned
# dN/dS values from step5f. Adds Bonferroni multiple-testing
# correction across the primary comparisons.
#
# Output: processed/conservation_summary_v2.txt
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(readr)
library(dplyr)

df <- read_csv("processed/conservation_scores_combined_v2.csv",
               show_col_types = FALSE)

bbb   <- df[df$Gene_Set == "BBB", ]
liver <- df[df$Gene_Set == "Liver_Control", ]

sink("processed/conservation_summary_v2.txt")
cat("=======================================================\n")
cat("BBB Gene Conservation — v2 (codon-aligned dN/dS)\n")
cat("=======================================================\n\n")

p_values <- c()
labels   <- c()

compare_groups <- function(x, y, label_x, label_y, metric_name, store = TRUE) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 5 || length(y) < 5) {
    cat(sprintf("--- %s --- skipped (too few non-NA)\n\n", metric_name))
    return(invisible(NULL))
  }
  wt <- wilcox.test(x, y, exact = FALSE)
  r  <- 1 - (2 * wt$statistic) / (length(x) * length(y))
  cat(sprintf("--- %s ---\n", metric_name))
  cat(sprintf("  %s: n=%d  median=%.3f  IQR=[%.3f, %.3f]\n",
      label_x, length(x), median(x), quantile(x, 0.25), quantile(x, 0.75)))
  cat(sprintf("  %s: n=%d  median=%.3f  IQR=[%.3f, %.3f]\n",
      label_y, length(y), median(y), quantile(y, 0.25), quantile(y, 0.75)))
  cat(sprintf("  Wilcoxon p = %.4g  |  effect size r = %.3f\n\n",
      wt$p.value, r))
  if (store) {
    p_values <<- c(p_values, wt$p.value)
    labels   <<- c(labels,   metric_name)
  }
  invisible(list(p = wt$p.value, r = r))
}

# ── 1. Full BBB vs Liver (4 primary tests) ───────────────────────────────────
cat("== 1. Full BBB list vs Liver control ==\n\n")
compare_groups(bbb$PctId_Human_Mouse,    liver$PctId_Human_Mouse,
               "BBB", "Liver", "% identity — Human vs Mouse")
compare_groups(bbb$PctId_Human_Macaque,  liver$PctId_Human_Macaque,
               "BBB", "Liver", "% identity — Human vs Macaque")
compare_groups(bbb$dNdS_Human_Mouse_v2,    liver$dNdS_Human_Mouse_v2,
               "BBB", "Liver", "dN/dS v2 — Human vs Mouse")
compare_groups(bbb$dNdS_Human_Macaque_v2,  liver$dNdS_Human_Macaque_v2,
               "BBB", "Liver", "dN/dS v2 — Human vs Macaque")

# ── 2. S3-only BBB vs Liver ──────────────────────────────────────────────────
cat("== 2. S3-only BBB (strictest) vs Liver control ==\n\n")
bbb_s3 <- bbb[!is.na(bbb$Daneman_Filter) & bbb$Daneman_Filter == "liver_and_lung", ]
cat(sprintf("S3-only genes: %d\n\n", nrow(bbb_s3)))
compare_groups(bbb_s3$PctId_Human_Mouse,   liver$PctId_Human_Mouse,
               "BBB_S3", "Liver", "S3 % id — Human vs Mouse", store = FALSE)
compare_groups(bbb_s3$PctId_Human_Macaque, liver$PctId_Human_Macaque,
               "BBB_S3", "Liver", "S3 % id — Human vs Macaque", store = FALSE)
compare_groups(bbb_s3$dNdS_Human_Mouse_v2,    liver$dNdS_Human_Mouse_v2,
               "BBB_S3", "Liver", "S3 dN/dS v2 — Human vs Mouse", store = FALSE)
compare_groups(bbb_s3$dNdS_Human_Macaque_v2,  liver$dNdS_Human_Macaque_v2,
               "BBB_S3", "Liver", "S3 dN/dS v2 — Human vs Macaque", store = FALSE)

# ── 3. Multiple testing correction on the 4 primary comparisons ───────────────
cat("== 3. Multiple testing correction (Bonferroni & FDR) ==\n\n")
adj_bonf <- p.adjust(p_values, method = "bonferroni")
adj_fdr  <- p.adjust(p_values, method = "BH")
for (i in seq_along(p_values)) {
  cat(sprintf("  %-45s  p=%.4g  Bonf=%.4g  FDR=%.4g\n",
      labels[i], p_values[i], adj_bonf[i], adj_fdr[i]))
}
cat("\n")

# ── 4. By orthologue type ────────────────────────────────────────────────────
cat("== 4. BBB: one2one vs one2many mouse orthologues ==\n\n")
bbb_121 <- bbb[!is.na(bbb$Mouse_OrthoType) & bbb$Mouse_OrthoType == "ortholog_one2one", ]
bbb_12m <- bbb[!is.na(bbb$Mouse_OrthoType) & bbb$Mouse_OrthoType == "ortholog_one2many", ]
compare_groups(bbb_121$PctId_Human_Mouse, bbb_12m$PctId_Human_Mouse,
               "one2one", "one2many", "% id — Mouse", store = FALSE)
compare_groups(bbb_121$dNdS_Human_Mouse_v2,  bbb_12m$dNdS_Human_Mouse_v2,
               "one2one", "one2many", "dN/dS v2 — Mouse", store = FALSE)

cat("=======================================================\n")
cat("Notes:\n")
cat(" - Wilcoxon rank-sum (non-parametric)\n")
cat(" - Effect size: rank-biserial r (|r|>0.1 small, >0.3 medium, >0.5 large)\n")
cat(" - dN/dS v2: codon-aligned (translate→align→back-map)\n")
cat(" - Bonferroni and FDR applied to the 4 primary tests\n")
cat("=======================================================\n")
sink()

cat(readLines("processed/conservation_summary_v2.txt"), sep="\n")
cat("\nSaved: processed/conservation_summary_v2.txt\n")
