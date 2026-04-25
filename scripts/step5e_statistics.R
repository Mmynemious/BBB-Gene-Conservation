# ============================================================
# Step 5e: Statistical comparison — BBB vs liver conservation
#
# Tests whether BBB genes are significantly more or less
# conserved than liver control genes using Wilcoxon rank-sum
# tests (non-parametric, appropriate for skewed distributions).
#
# Also runs stratified analyses:
#   - S3-only (strictest Daneman) vs full BBB list
#   - By human validation score (Human_BBB_Datasets_N)
#   - By orthologue type (one2one vs one2many)
#
# Output: processed/conservation_stats.csv
#         processed/conservation_summary.txt
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(readr)
library(dplyr)

df <- read_csv("processed/conservation_scores_combined.csv",
               show_col_types = FALSE)

bbb   <- df[df$Gene_Set == "BBB", ]
liver <- df[df$Gene_Set == "Liver_Control", ]

sink("processed/conservation_summary.txt")
cat("=======================================================\n")
cat("BBB Gene Conservation Analysis — Statistical Summary\n")
cat("=======================================================\n\n")

# ── Helper: Wilcoxon test + effect size (rank-biserial r) ────────────────────
compare_groups <- function(x, y, label_x, label_y, metric_name) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  wt <- wilcox.test(x, y, exact = FALSE)
  # Rank-biserial correlation as effect size
  r  <- 1 - (2 * wt$statistic) / (length(x) * length(y))
  cat(sprintf("--- %s ---\n", metric_name))
  cat(sprintf("  %s: n=%d  median=%.2f  IQR=[%.2f, %.2f]\n",
      label_x, length(x), median(x), quantile(x, 0.25), quantile(x, 0.75)))
  cat(sprintf("  %s: n=%d  median=%.2f  IQR=[%.2f, %.2f]\n",
      label_y, length(y), median(y), quantile(y, 0.25), quantile(y, 0.75)))
  cat(sprintf("  Wilcoxon p = %.4g  |  effect size r = %.3f\n\n", wt$p.value, r))
  list(p = wt$p.value, r = as.numeric(r), n_x = length(x), n_y = length(y),
       median_x = median(x), median_y = median(y))
}

# ── 1. Full BBB vs Liver ──────────────────────────────────────────────────────
cat("== 1. Full BBB list vs Liver control ==\n\n")
r1a <- compare_groups(bbb$PctId_Human_Mouse,    liver$PctId_Human_Mouse,
                      "BBB", "Liver", "% identity — Human vs Mouse")
r1b <- compare_groups(bbb$PctId_Human_Macaque,  liver$PctId_Human_Macaque,
                      "BBB", "Liver", "% identity — Human vs Macaque")
r1c <- compare_groups(bbb$dNdS_Human_Mouse,     liver$dNdS_Human_Mouse,
                      "BBB", "Liver", "dN/dS — Human vs Mouse")
r1d <- compare_groups(bbb$dNdS_Human_Macaque,   liver$dNdS_Human_Macaque,
                      "BBB", "Liver", "dN/dS — Human vs Macaque")

# ── 2. S3-only BBB vs Liver ───────────────────────────────────────────────────
cat("== 2. S3-only BBB (strictest) vs Liver control ==\n\n")
bbb_s3 <- bbb[!is.na(bbb$Daneman_Filter) & bbb$Daneman_Filter == "liver_and_lung", ]
cat(sprintf("S3-only genes: %d\n\n", nrow(bbb_s3)))
compare_groups(bbb_s3$PctId_Human_Mouse,   liver$PctId_Human_Mouse,
               "BBB_S3", "Liver", "% identity — Human vs Mouse")
compare_groups(bbb_s3$PctId_Human_Macaque, liver$PctId_Human_Macaque,
               "BBB_S3", "Liver", "% identity — Human vs Macaque")
compare_groups(bbb_s3$dNdS_Human_Mouse,    liver$dNdS_Human_Mouse,
               "BBB_S3", "Liver", "dN/dS — Human vs Mouse")
compare_groups(bbb_s3$dNdS_Human_Macaque,  liver$dNdS_Human_Macaque,
               "BBB_S3", "Liver", "dN/dS — Human vs Macaque")

# ── 3. By human validation score ─────────────────────────────────────────────
cat("== 3. BBB genes by human validation score (N=0 to 3) ==\n\n")
for (n in 0:3) {
  sub <- bbb[!is.na(bbb$Human_BBB_Datasets_N) & bbb$Human_BBB_Datasets_N == n, ]
  if (nrow(sub) < 5) next
  cat(sprintf("N=%d (%d genes): Mouse %%id median=%.1f%%  Macaque %%id median=%.1f%%\n",
      n, nrow(sub),
      median(sub$PctId_Human_Mouse, na.rm=TRUE),
      median(sub$PctId_Human_Macaque, na.rm=TRUE)))
}
cat("\n")

# ── 4. One-to-one vs one-to-many orthologues ─────────────────────────────────
cat("== 4. BBB: one2one vs one2many mouse orthologues ==\n\n")
bbb_121 <- bbb[!is.na(bbb$Mouse_OrthoType) & bbb$Mouse_OrthoType == "ortholog_one2one", ]
bbb_12m <- bbb[!is.na(bbb$Mouse_OrthoType) & bbb$Mouse_OrthoType == "ortholog_one2many", ]
compare_groups(bbb_121$PctId_Human_Mouse, bbb_12m$PctId_Human_Mouse,
               "one2one", "one2many", "% identity — Human vs Mouse")
compare_groups(bbb_121$dNdS_Human_Mouse,  bbb_12m$dNdS_Human_Mouse,
               "one2one", "one2many", "dN/dS — Human vs Mouse")

# ── 5. Munji-only vs Daneman-only vs Both ────────────────────────────────────
cat("== 5. BBB: by source ==\n\n")
for (src in c("Daneman", "Munji", "Both")) {
  sub <- bbb[!is.na(bbb$Source) & bbb$Source == src, ]
  if (nrow(sub) < 5) next
  cat(sprintf("%s (%d genes): Mouse %%id=%.1f%%  Macaque %%id=%.1f%%  Mouse dN/dS=%.3f\n",
      src, nrow(sub),
      median(sub$PctId_Human_Mouse, na.rm=TRUE),
      median(sub$PctId_Human_Macaque, na.rm=TRUE),
      median(sub$dNdS_Human_Mouse, na.rm=TRUE)))
}

cat("\n=======================================================\n")
cat("Note: Wilcoxon rank-sum test used throughout (non-parametric).\n")
cat("Effect size: rank-biserial r (|r|>0.1 small, >0.3 medium, >0.5 large).\n")
cat("=======================================================\n")
sink()

# Print to console too
cat(readLines("processed/conservation_summary.txt"), sep="\n")
cat("\nSaved: processed/conservation_summary.txt\n")
