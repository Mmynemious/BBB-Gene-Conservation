# ============================================================
# Step 2: Process the Wälchli 2024 Seurat object
# Loads the .rds.gz file, computes average gene expression
# per cell cluster, and saves as a human-readable CSV.
# ============================================================

# ── Set working directory ────────────────────────────────────────────────────
setwd("~/Documents/Claude/Projects/BBB")

# ── Install Seurat if missing ─────────────────────────────────────────────────
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")

library(Seurat)
library(readr)

# ── STEP A: Load the Seurat object ───────────────────────────────────────────
# The file is a compressed R object (.rds.gz).
# gzcon() + file() lets R read it without fully decompressing to disk first.
# This is a large file (~600k cells) — expect this step to take 1-3 minutes.

cat("Loading Wälchli Seurat object (may take 1-3 min)...\n")

seurat_path <- "raw_data/walchli_2024/Walchli2024_AdultBrain_TemporalLobe_SortedEndothelial_SeuratObject.rds.gz"

walchli <- readRDS(pipe(paste("gunzip -c", shQuote(seurat_path), "| gunzip -c")))

cat("Loaded successfully.\n")
cat("Object class:", class(walchli), "\n")
cat("Number of cells:", ncol(walchli), "\n")
cat("Number of genes:", nrow(walchli), "\n")

# ── STEP B: Inspect cell cluster labels ──────────────────────────────────────
# Seurat objects group cells into clusters with identity labels.
# We need to know what those labels are called before averaging.

cat("\nCell identity labels (first 20):\n")
print(head(Idents(walchli), 20))

cat("\nAll unique cluster names:\n")
print(levels(Idents(walchli)))

cat("\nCells per cluster:\n")
print(table(Idents(walchli)))

# ── STEP C: Compute average expression per cluster ───────────────────────────
# AverageExpression() collapses all cells of the same type into one average.
# Output: genes as rows, cell clusters as columns, mean expression as values.
# We use the RNA assay (raw gene counts, normalised).

cat("\nComputing average expression per cluster...\n")

avg_expr <- AverageExpression(walchli, assays = "RNA", slot = "data",
                              group.by = "ECclusters")$RNA

cat("Average expression matrix dimensions:", dim(avg_expr), "\n")
cat("(rows = genes, columns = cell clusters)\n")

# ── STEP D: Save as CSV ───────────────────────────────────────────────────────
dir.create("processed", showWarnings = FALSE)
output_path <- "processed/Walchli2024_AverageExpression.csv"

write.csv(avg_expr, output_path)
cat("\nSaved:", output_path, "\n")

# ── STEP E: Quick preview ─────────────────────────────────────────────────────
cat("\nFirst 5 genes, first 3 clusters:\n")
print(avg_expr[1:5, 1:min(3, ncol(avg_expr))])

cat("\nStep 2 complete. Run scripts/eval_step2.R to validate.\n")
