# ============================================================
# Step 3b: Re-query biomaRt with orthology confidence scores
#
# Replaces the HomoloGene-based Step 1 mouse→human conversion
# for genes that were dropped (NA HumanEnsembl) and adds:
#   - orthology_type: one2one vs one2many vs many2many
#   - mouse_perc_id: % sequence identity to mouse orthologue
#   - macaque_perc_id: % sequence identity to macaque orthologue
#
# Input:  processed/master_BBB_genelist_validated.csv
#         processed/BBB_genes_crossspecies_ensembl.csv
# Output: processed/BBB_genes_orthology_confidence.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(biomaRt)
library(readr)
library(dplyr)

master    <- read_csv("processed/master_BBB_genelist_validated.csv", show_col_types = FALSE)
ensembl   <- read_csv("processed/BBB_genes_crossspecies_ensembl.csv",  show_col_types = FALSE)

cat("Master list:", nrow(master), "genes\n")
cat("Ensembl table:", nrow(ensembl), "rows\n")

# ── Connect to Ensembl release 113 ───────────────────────────────────────────
cat("\nConnecting to Ensembl BioMart release 113...\n")
human_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                          mirror = "useast")
cat("Connected.\n\n")

# ── Query 1: human gene symbol → Ensembl ID (catches HomoloGene misses) ──────
# Use all unique human gene symbols from master list
human_symbols <- unique(na.omit(master$GeneSymbol_Human))

cat("Querying human gene symbols → Ensembl IDs for", length(human_symbols), "genes...\n")
human_ids <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = human_symbols,
  mart       = human_mart
)
# Keep only one Ensembl ID per symbol (primary mapping)
human_ids <- human_ids |>
  filter(hgnc_symbol != "", ensembl_gene_id != "") |>
  distinct(hgnc_symbol, .keep_all = TRUE)

cat("Human IDs retrieved:", nrow(human_ids), "\n")

# ── Query 2: mouse orthologues + confidence ───────────────────────────────────
cat("Querying mouse orthologues + % identity...\n")
mouse_ortho <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmusculus_homolog_ensembl_gene",
    "mmusculus_homolog_orthology_type",
    "mmusculus_homolog_perc_id"           # % identity: human→mouse direction
  ),
  filters = "ensembl_gene_id",
  values  = human_ids$ensembl_gene_id,
  mart    = human_mart
)
mouse_ortho <- mouse_ortho |>
  filter(mmusculus_homolog_ensembl_gene != "") |>
  rename(
    HumanEnsembl      = ensembl_gene_id,
    MouseEnsembl      = mmusculus_homolog_ensembl_gene,
    Mouse_OrthoType   = mmusculus_homolog_orthology_type,
    Mouse_PercId      = mmusculus_homolog_perc_id
  )

cat("Mouse orthologue rows:", nrow(mouse_ortho), "\n")

# ── Query 3: macaque orthologues + confidence ─────────────────────────────────
cat("Querying macaque orthologues + % identity...\n")
macaque_ortho <- getBM(
  attributes = c(
    "ensembl_gene_id",
    "mmulatta_homolog_ensembl_gene",
    "mmulatta_homolog_orthology_type",
    "mmulatta_homolog_perc_id"
  ),
  filters = "ensembl_gene_id",
  values  = human_ids$ensembl_gene_id,
  mart    = human_mart
)
macaque_ortho <- macaque_ortho |>
  filter(mmulatta_homolog_ensembl_gene != "") |>
  rename(
    HumanEnsembl        = ensembl_gene_id,
    MacaqueEnsembl      = mmulatta_homolog_ensembl_gene,
    Macaque_OrthoType   = mmulatta_homolog_orthology_type,
    Macaque_PercId      = mmulatta_homolog_perc_id
  )

cat("Macaque orthologue rows:", nrow(macaque_ortho), "\n")

# ── Join everything ───────────────────────────────────────────────────────────
# Start from human IDs, join mouse then macaque
combined <- human_ids |>
  rename(GeneSymbol_Human = hgnc_symbol, HumanEnsembl = ensembl_gene_id) |>
  left_join(mouse_ortho,   by = "HumanEnsembl") |>
  left_join(macaque_ortho, by = "HumanEnsembl")

cat("\nCombined table:", nrow(combined), "rows\n")

# ── Join back to master list ──────────────────────────────────────────────────
# Attach validation columns and source metadata from master list
master_slim <- master |>
  select(GeneSymbol_Mouse, GeneSymbol_Human, Source, Daneman_Filter,
         In_Walchli_EC, In_Yang_EC, In_Winkler_EC, Human_BBB_Datasets_N)

final <- master_slim |>
  left_join(combined, by = "GeneSymbol_Human")

cat("Final table:", nrow(final), "rows\n")

# ── Coverage summary ─────────────────────────────────────────────────────────
n_human   <- sum(!is.na(final$HumanEnsembl) & final$HumanEnsembl != "")
n_mouse   <- sum(!is.na(final$MouseEnsembl) & final$MouseEnsembl != "")
n_macaque <- sum(!is.na(final$MacaqueEnsembl) & final$MacaqueEnsembl != "")

cat(sprintf("\n--- Coverage ---\n"))
cat(sprintf("Human Ensembl IDs:   %d / %d (%.1f%%)\n", n_human,   nrow(final), 100 * n_human   / nrow(final)))
cat(sprintf("Mouse orthologues:   %d / %d (%.1f%%)\n", n_mouse,   nrow(final), 100 * n_mouse   / nrow(final)))
cat(sprintf("Macaque orthologues: %d / %d (%.1f%%)\n", n_macaque, nrow(final), 100 * n_macaque / nrow(final)))

# ── Orthologue type distribution ─────────────────────────────────────────────
cat("\n--- Mouse orthologue type breakdown ---\n")
print(table(final$Mouse_OrthoType, useNA = "ifany"))

cat("\n--- Macaque orthologue type breakdown ---\n")
print(table(final$Macaque_OrthoType, useNA = "ifany"))

# ── % identity summary ───────────────────────────────────────────────────────
cat(sprintf("\nMouse %%id:   median=%.1f  mean=%.1f  min=%.1f\n",
    median(final$Mouse_PercId, na.rm=TRUE),
    mean(final$Mouse_PercId, na.rm=TRUE),
    min(final$Mouse_PercId, na.rm=TRUE)))
cat(sprintf("Macaque %%id: median=%.1f  mean=%.1f  min=%.1f\n",
    median(final$Macaque_PercId, na.rm=TRUE),
    mean(final$Macaque_PercId, na.rm=TRUE),
    min(final$Macaque_PercId, na.rm=TRUE)))

# ── Save ─────────────────────────────────────────────────────────────────────
write_csv(final, "processed/master_BBB_genelist_complete.csv")
cat("\nSaved: processed/master_BBB_genelist_complete.csv\n")
cat("Run scripts/eval_step3b.R to validate.\n")
