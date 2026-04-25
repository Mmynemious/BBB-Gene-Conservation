# ============================================================
# Step 4c: Build liver control gene set from Human Protein Atlas
#
# Downloads HPA RNA tissue consensus expression (nTPM across tissues),
# selects genes highly expressed in liver but NOT in our BBB master list.
#
# Control set criteria:
#   - nTPM in liver >= 10 (meaningfully expressed)
#   - Not already in the BBB master list
#   - Not highly expressed in brain endothelium (nTPM brain < 10)
#     to avoid accidentally including BBB-adjacent genes
#
# Output: processed/liver_control_genelist.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(readr)
library(dplyr)

# ── Download HPA RNA tissue consensus file ───────────────────────────────────
hpa_zip  <- "raw_data/hpa_liver/rna_tissue_consensus.tsv.zip"
hpa_file <- "raw_data/hpa_liver/rna_tissue_consensus.tsv"

dir.create("raw_data/hpa_liver", showWarnings = FALSE)

if (!file.exists(hpa_file)) {
  cat("Downloading Human Protein Atlas tissue expression data...\n")
  download.file(
    url     = "https://www.proteinatlas.org/download/tsv/rna_tissue_consensus.tsv.zip",
    destfile = hpa_zip,
    mode    = "wb",
    method  = "curl"
  )
  unzip(hpa_zip, exdir = "raw_data/hpa_liver")
  cat("Downloaded and extracted.\n")
} else {
  cat("HPA file already present, skipping download.\n")
}

# ── Load and inspect ─────────────────────────────────────────────────────────
hpa <- read_tsv(hpa_file, show_col_types = FALSE)
cat("HPA data dimensions:", nrow(hpa), "rows x", ncol(hpa), "cols\n")
cat("Columns:", paste(colnames(hpa), collapse = ", "), "\n")
cat("Sample tissues:", paste(unique(hpa$Tissue)[1:10], collapse = ", "), "\n\n")

# ── Filter for liver-expressed genes ─────────────────────────────────────────
# Liver tissue label in HPA
liver_genes <- hpa |>
  filter(Tissue == "liver", nTPM >= 10) |>
  select(Gene = `Gene name`, Liver_nTPM = nTPM) |>
  arrange(desc(Liver_nTPM))

cat("Genes with liver nTPM >= 10:", nrow(liver_genes), "\n")

# ── Get brain nTPM for the same genes (to exclude brain-high genes) ───────────
brain_expr <- hpa |>
  filter(Tissue == "brain", `Gene name` %in% liver_genes$Gene) |>
  select(Gene = `Gene name`, Brain_nTPM = nTPM)

liver_genes <- liver_genes |>
  left_join(brain_expr, by = "Gene") |>
  mutate(Brain_nTPM = ifelse(is.na(Brain_nTPM), 0, Brain_nTPM))

# ── Exclude genes in the BBB master list ─────────────────────────────────────
master <- read_csv("processed/master_BBB_genelist_complete.csv",
                   show_col_types = FALSE)
bbb_genes <- unique(na.omit(master$GeneSymbol_Human))

liver_control <- liver_genes |>
  filter(
    !Gene %in% bbb_genes,      # not already a BBB gene
    Brain_nTPM < 10             # not highly expressed in brain
  )

cat("Liver control genes after exclusions:", nrow(liver_control), "\n")
cat(sprintf("  (removed %d BBB overlap genes)\n",
    sum(liver_genes$Gene %in% bbb_genes)))
cat(sprintf("  (removed %d brain-high genes)\n",
    sum(liver_genes$Brain_nTPM >= 10 & !liver_genes$Gene %in% bbb_genes)))

# ── Summary stats ─────────────────────────────────────────────────────────────
cat(sprintf("\nLiver control set: %d genes\n", nrow(liver_control)))
cat(sprintf("BBB gene list:     %d genes\n", length(bbb_genes)))
cat(sprintf("Ratio (control/BBB): %.1fx\n", nrow(liver_control) / length(bbb_genes)))
cat(sprintf("\nLiver nTPM: median=%.1f  mean=%.1f  min=%.1f  max=%.1f\n",
    median(liver_control$Liver_nTPM),
    mean(liver_control$Liver_nTPM),
    min(liver_control$Liver_nTPM),
    max(liver_control$Liver_nTPM)))

# ── Save ─────────────────────────────────────────────────────────────────────
write_csv(liver_control, "processed/liver_control_genelist.csv")
cat("\nSaved: processed/liver_control_genelist.csv\n")
cat("Run scripts/eval_step4c.R to validate.\n")
